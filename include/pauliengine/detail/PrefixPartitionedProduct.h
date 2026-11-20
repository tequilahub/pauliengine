#pragma once
/// Internal: the prefix-partitioned product behind products, commutators and squares of
/// QubitHamiltonians with complex coefficients, and of unweighted sets.

#include <algorithm>
#include <bit>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "pauliengine/detail/FlatTerms.h"
#include "pauliengine/detail/Parallel.h"

namespace pauliengine_detail {

/// Products and commutators of two flat sums, formed in pieces that fit the cache.
///
/// The operator of a product is the XOR of the operators of its factors, so the
/// leading bits of a product follow from the leading bits of the factors alone.
/// With both sums ordered by operator, the terms sharing some leading bits are a
/// range of rows. A block is a range of a times a range of b; a task is a list of
/// blocks whose products all share their leading bits. Splitting a task by the
/// next d bits gives 2^d tasks, in the order of the operators they produce, so
/// the results of the tasks only have to be written one after the other.
///
/// A task is split until it is small. Its pairs are then formed, sorted by the
/// remaining bits and merged, all within the cache; where the remaining bits
/// allow few distinct operators the coefficients are added up in a table instead,
/// without sorting. No array of all n1 * n2 products ever exists.
///
/// Only bit positions that differ somewhere within a or within b can differ
/// between products. These are packed, most significant first, into 64-bit keys
/// (one per level of 64 positions); the key of a product is the XOR of the keys
/// of its factors, and its words are formed only for the terms that remain.
/// The result does not depend on the number of threads.
///
/// A commutator keeps the anti-commuting pairs only, each counted twice. A square
/// a * a keeps the pairs i <= j that commute: P_i P_j and P_j P_i are equal if the
/// strings commute and cancel if they anti-commute, so
///
///     a * a = sum_i a_i^2 I + 2 sum_{i < j, P_i P_j = P_j P_i} a_i a_j P_i P_j .
///
/// About a quarter of the pairs remain, and the cancelling ones leave no rounding
/// residue. Blocks of rows i > j are dropped when a task is split.
///
/// first_pairs() serves sums without coefficients, which keep the first pair of
/// every operator instead of a sum: the coefficients of the rows of a and b are then
/// (index, phase exponent), and every operator gets the (rank, phase exponent) of
/// its pair of lowest rank i * size of b + j, where i and j are those indices.
class PrefixPartitionedProduct {
    public:
        /// a * b, or [a, b] if `commutator`: compacted and ordered by operator.
        static FlatTerms run(const FlatTerms& a, const FlatTerms& b, bool commutator) {
                if (a.size() == 0 || b.size() == 0) return FlatTerms(0, std::max(a.width(), b.width()));
                if (a.size() > UINT32_MAX || b.size() > UINT32_MAX) {
                        throw std::length_error("QubitHamiltonian: more than 2^32 - 1 terms in a product");
                }
                return PrefixPartitionedProduct(a, b, commutator ? Mode::Commutator : Mode::Product).result();
        }

        /// For every operator among the products (or the commutators, if
        /// `anticommuting`) of the rows of a and b, the pair of lowest rank, as
        /// described above. Ordered by operator.
        static FlatTerms first_pairs(const FlatTerms& a, const FlatTerms& b, bool anticommuting) {
                if (a.size() == 0 || b.size() == 0) return FlatTerms(0, std::max(a.width(), b.width()));
                if (a.size() > UINT32_MAX || b.size() > UINT32_MAX) {
                        throw std::length_error("more than 2^32 - 1 terms in a product");
                }
                const Mode mode = anticommuting ? Mode::FirstAnticommuting : Mode::First;
                return PrefixPartitionedProduct(a, b, mode).result();
        }

        /// a * a: compacted and ordered by operator.
        static FlatTerms square(const FlatTerms& a) {
                if (a.size() == 0) return FlatTerms(0, a.width());
                if (a.size() > UINT32_MAX) {
                        throw std::length_error("QubitHamiltonian: more than 2^32 - 1 terms in a product");
                }
                return PrefixPartitionedProduct(a, a, Mode::Square).result();
        }

        /// Whether at least an eighth of the pairs anti-commute, judged by a grid of
        /// at most 32 x 32 pairs. run() visits every pair of a commutator twice, which
        /// pays off only if a fair share of the pairs contribute.
        static bool many_anticommute(const FlatTerms& a, const FlatTerms& b) {
                // a smaller grid for small sums, where it is a visible share of the work
                const size_t side = a.size() * b.size() >= 65536 ? 32 : 16;
                const size_t rows_a = std::min<size_t>(a.size(), side);
                const size_t rows_b = std::min<size_t>(b.size(), side);
                size_t found = 0;
                for (size_t s = 0; s < rows_a; ++s) {
                        const TermView left = a.view(s * a.size() / rows_a);
                        for (size_t t = 0; t < rows_b; ++t) {
                                found += anticommute(left, b.view(t * b.size() / rows_b));
                        }
                }
                return 8 * found >= rows_a * rows_b;
        }

    private:
        enum class Mode { Product, Commutator, Square, First, FirstAnticommuting };

        /// A task of at most this many pairs is not split further. Its records
        /// (32 bytes each, 128 KB in all) fit even a 256 KB L2 cache; leaves of
        /// 2048 to 8192 pairs ran equally fast, larger ones were slower with many threads.
        static constexpr size_t kLeafPairs = 4096;
        /// Size aimed at when choosing how many bits to split by.
        static constexpr size_t kLeafTarget = 2048;
        /// Blocks are kept at this many pairs on average; smaller ones cost more to visit than to compute.
        static constexpr size_t kMinBlockPairs = 32;
        /// Tasks of at most this many pairs are handed to the threads as a whole.
        static constexpr size_t kMaxUnitPairs = 32768;
        static constexpr unsigned kMaxDigitBits = 8;
        /// Largest table of coefficient sums: 2^13 entries of 24 bytes (192 KB),
        /// within a 256 KB L2 cache as well.
        static constexpr unsigned kMaxTableBits = 13;
        static constexpr unsigned kMaxScatterBits = 13;
        static constexpr size_t kInsertionLimit = 32;
        /// Results of fewer terms have their gaps closed in place by one thread.
        static constexpr size_t kMoveInPlaceRows = size_t{1} << 17;
        static constexpr uint32_t kNone = UINT32_MAX;

        /// Rows [a0, a1) of a times rows [b0, b1) of b.
        /// In a square the two ranges are either the same, and stand for the pairs
        /// i <= j, or the rows of the first come before those of the second.
        struct Block {
                uint32_t a0, a1, b0, b1;
        };

        /// One product: its key at the current level, coefficient and factors.
        struct Record {
                uint64_t key;
                double re, im;
                uint32_t i, j;
        };

        /// Rows [first, last) sharing the digit `digit`.
        struct Run {
                uint32_t first, last, digit;
        };

        /// How the bits [src_shift, src_shift + len) of a key word enter a packed key.
        struct Piece {
                uint32_t key_word;
                uint32_t level;
                uint8_t src_shift, dst_shift, len;
        };

        /// Where the terms of a unit are written: rows of 2 * nw words and
        /// (real, imaginary) pairs, with room for all terms the unit can give.
        struct Output {
                uint64_t* words;
                double* coeffs;
                size_t n;
        };

        /// A task processed by one thread; its results come before those of the next unit.
        struct Unit {
                std::vector<Block> blocks;
                size_t pos;
                size_t pairs;
                size_t max_rows = 0;  // bound on the number of its terms
                size_t offset = 0;    // sum of the bounds of the units before it
                size_t rows = 0;      // number of its terms
                size_t first = 0;     // row of its first term in the result
        };

        /// Working memory of one thread.
        struct Scratch {
                explicit Scratch(size_t positions) : children(positions + 1) {}

                std::unique_ptr<Record[]> records;
                size_t record_cap = 0;
                std::vector<uint32_t> ends;
                std::vector<double> sums;
                std::vector<uint32_t> factors;
                std::vector<uint64_t> kept;  // (i << 32) | j of the pairs that count
                std::vector<Run> runs_a, runs_b;
                /// children[pos]: the tasks a task at bit position pos was split into.
                std::vector<std::vector<std::vector<Block>>> children;
        };

        PrefixPartitionedProduct(const FlatTerms& a, const FlatTerms& b, Mode mode)
                : mode_(mode), nw_(std::max(a.width(), b.width())) {
                a_ = ordered(a, a_sorted_);
                b_ = mode == Mode::Square ? a_ : ordered(b, b_sorted_);
                n1_ = a_->size();
                n2_ = b_->size();
                nwa_ = a_->width();
                nwb_ = b_->width();
                one_word_ = nwa_ == 1 && nwb_ == 1;
                aw_ = a_->x(0);
                bw_ = b_->x(0);

                // Positions that differ within a or within b, as masks of the key words
                // (x words, then y words, both of the common width).
                std::vector<uint64_t> varying(2 * nw_, 0);
                mark_varying(*a_, varying);
                if (mode_ != Mode::Square) mark_varying(*b_, varying);

                // Runs of such positions, most significant first, cut at the level bounds.
                std::vector<Piece> pieces;
                size_t pos = 0;
                for (size_t q = 0; q < 2 * nw_; ++q) {
                        uint64_t m = varying[q];
                        while (m != 0) {
                                const unsigned high = 63 - static_cast<unsigned>(std::countl_zero(m));
                                const uint64_t upto_high = high == 63 ? ~uint64_t{0} : ((uint64_t{1} << (high + 1)) - 1);
                                const uint64_t zeros = ~m & upto_high;
                                const unsigned low = zeros == 0 ? 0 : 64 - static_cast<unsigned>(std::countl_zero(zeros));
                                m &= low == 0 ? 0 : ((uint64_t{1} << low) - 1);
                                unsigned len = high - low + 1;
                                unsigned top = high + 1;  // one past the highest bit not yet placed
                                while (len > 0) {
                                        const unsigned p = static_cast<unsigned>(pos % 64);
                                        const unsigned take = std::min(len, 64 - p);
                                        pieces.push_back({static_cast<uint32_t>(q), static_cast<uint32_t>(pos / 64),
                                                          static_cast<uint8_t>(top - take), static_cast<uint8_t>(64 - p - take),
                                                          static_cast<uint8_t>(take)});
                                        pos += take;
                                        len -= take;
                                        top -= take;
                                }
                        }
                }
                nvar_ = pos;
                levels_ = std::max<size_t>(1, (nvar_ + 63) / 64);
                ka_ = pack_keys(*a_, pieces);
                kb_ = mode_ == Mode::Square ? ka_ : pack_keys(*b_, pieces);
        }

        static bool row_less(const FlatTerms& t, size_t i, size_t j) {
                const size_t row_words = 2 * t.width();
                return std::lexicographical_compare(t.x(i), t.x(i) + row_words, t.x(j), t.x(j) + row_words);
        }

        /// `t` if its rows are ordered by operator, as the results of all operations
        /// are; otherwise an ordered copy, kept in `storage`.
        static const FlatTerms* ordered(const FlatTerms& t, FlatTerms& storage) {
                bool sorted = true;
                for (size_t i = 1; i < t.size() && sorted; ++i) sorted = !row_less(t, i, i - 1);
                if (sorted) return &t;
                std::vector<uint32_t> order(t.size());
                std::iota(order.begin(), order.end(), uint32_t{0});
                std::stable_sort(order.begin(), order.end(), [&](uint32_t i, uint32_t j) { return row_less(t, i, j); });
                storage = FlatTerms(t.size(), t.width());
                for (size_t r = 0; r < order.size(); ++r) {
                        std::copy(t.x(order[r]), t.x(order[r]) + 2 * t.width(), storage.x(r));
                        storage.set_coeff(r, t.coeff(order[r]));
                }
                return &storage;
        }

        void mark_varying(const FlatTerms& t, std::vector<uint64_t>& varying) const {
                const size_t w = t.width();
                if (w == 0) return;
                std::vector<uint64_t> any(2 * w, 0);
                std::vector<uint64_t> all(2 * w, ~uint64_t{0});
                for (size_t i = 0; i < t.size(); ++i) {
                        const uint64_t* row = t.x(i);
                        for (size_t k = 0; k < 2 * w; ++k) {
                                any[k] |= row[k];
                                all[k] &= row[k];
                        }
                }
                for (size_t k = 0; k < w; ++k) {
                        varying[k] |= any[k] ^ all[k];
                        varying[nw_ + k] |= any[w + k] ^ all[w + k];
                }
        }

        /// Packed keys of all rows: level * size() + row.
        std::vector<uint64_t> pack_keys(const FlatTerms& t, const std::vector<Piece>& pieces) const {
                const size_t n = t.size();
                const size_t w = t.width();
                std::vector<uint64_t> keys(levels_ * n, 0);
                for (const Piece& piece : pieces) {
                        const size_t k = piece.key_word < nw_ ? piece.key_word : piece.key_word - nw_;
                        if (k >= w) continue;  // a word this sum does not have: zero
                        const size_t offset = piece.key_word < nw_ ? k : w + k;
                        const uint64_t mask = piece.len == 64 ? ~uint64_t{0} : ((uint64_t{1} << piece.len) - 1);
                        uint64_t* out = keys.data() + piece.level * n;
                        const uint64_t* word = t.x(0) + offset;
                        for (size_t i = 0; i < n; ++i, word += 2 * w) {
                                out[i] |= ((*word >> piece.src_shift) & mask) << piece.dst_shift;
                        }
                }
                return keys;
        }

        unsigned bits_in_level(size_t level) const noexcept {
                return static_cast<unsigned>(std::min<size_t>(64, nvar_ - 64 * level));
        }

        // ------------------------------------------------------------------
        // Splitting
        // ------------------------------------------------------------------

        /// Whether splitting leaves blocks that are worth visiting.
        static bool can_split(size_t nblocks, size_t pairs) noexcept {
                return pairs / (4 * kMinBlockPairs) >= nblocks;
        }

        /// Number of bits to split a task by, given that `remaining` are left in its level.
        static unsigned digit_bits(size_t nblocks, size_t pairs, unsigned remaining) noexcept {
                const unsigned by_blocks = (static_cast<unsigned>(std::bit_width(pairs / (kMinBlockPairs * nblocks))) - 1) / 2;
                const unsigned by_size = static_cast<unsigned>(std::bit_width((pairs - 1) / kLeafTarget));
                return std::clamp(std::min(by_blocks, by_size), 1u, std::min(remaining, kMaxDigitBits));
        }

        static void find_runs(const uint64_t* keys, uint32_t first, uint32_t last, unsigned shift, uint64_t mask,
                              std::vector<Run>& runs) {
                runs.clear();
                for (uint32_t i = first; i < last;) {
                        const uint32_t digit = static_cast<uint32_t>((keys[i] >> shift) & mask);
                        const uint32_t start = i;
                        do { ++i; } while (i < last && ((keys[i] >> shift) & mask) == digit);
                        runs.push_back({start, i, digit});
                }
        }

        /// Splits the blocks by the d bits at position pos; children[c] gets the blocks
        /// whose products have digit c.
        void split(const Block* blocks, size_t nblocks, size_t pos, unsigned d,
                   std::vector<std::vector<Block>>& children, std::vector<Run>& runs_a, std::vector<Run>& runs_b) const {
                const size_t fanout = size_t{1} << d;
                if (children.size() < fanout) children.resize(fanout);
                for (size_t c = 0; c < fanout; ++c) children[c].clear();
                const size_t level = pos / 64;
                const unsigned shift = 64 - static_cast<unsigned>(pos % 64) - d;
                const uint64_t mask = fanout - 1;
                const uint64_t* ka = ka_.data() + level * n1_;
                const uint64_t* kb = kb_.data() + level * n2_;
                const bool square = mode_ == Mode::Square;
                for (size_t r = 0; r < nblocks; ++r) {
                        find_runs(ka, blocks[r].a0, blocks[r].a1, shift, mask, runs_a);
                        find_runs(kb, blocks[r].b0, blocks[r].b1, shift, mask, runs_b);
                        for (const Run& ra : runs_a) {
                                for (const Run& rb : runs_b) {
                                        if (square && ra.first > rb.first) continue;  // the pairs i > j
                                        children[ra.digit ^ rb.digit].push_back({ra.first, ra.last, rb.first, rb.last});
                                }
                        }
                }
        }

        size_t pairs_of(const Block& block) const noexcept {
                const size_t rows = block.a1 - block.a0;
                if (mode_ == Mode::Square && block.a0 == block.b0) return rows * (rows + 1) / 2;
                return rows * (block.b1 - block.b0);
        }

        size_t count_pairs(const std::vector<Block>& blocks) const noexcept {
                size_t pairs = 0;
                for (const Block& block : blocks) pairs += pairs_of(block);
                return pairs;
        }

        /// Splits a task into units, appended in the order of their operators.
        void expand(std::vector<Block>&& blocks, size_t pos, size_t pairs, size_t unit_pairs,
                    std::vector<Unit>& units, std::vector<Run>& runs_a, std::vector<Run>& runs_b) const {
                if (pairs <= unit_pairs || pos == nvar_ || !can_split(blocks.size(), pairs)) {
                        units.push_back(Unit{std::move(blocks), pos, pairs});
                        return;
                }
                const unsigned d = digit_bits(blocks.size(), pairs, bits_in_level(pos / 64) - static_cast<unsigned>(pos % 64));
                std::vector<std::vector<Block>> children;
                split(blocks.data(), blocks.size(), pos, d, children, runs_a, runs_b);
                std::vector<Block>().swap(blocks);
                for (auto& child : children) {
                        if (!child.empty()) expand(std::move(child), pos + d, count_pairs(child), unit_pairs, units, runs_a, runs_b);
                }
        }

        // ------------------------------------------------------------------
        // Pairs
        // ------------------------------------------------------------------

        /// Calls f(i, j) for every pair of the blocks, the longer range innermost.
        template<typename F>
        void for_pairs(const Block* blocks, size_t nblocks, F&& f) const {
                const bool square = mode_ == Mode::Square;
                for (size_t r = 0; r < nblocks; ++r) {
                        const Block block = blocks[r];
                        if (square && block.a0 == block.b0) {
                                for (uint32_t i = block.a0; i < block.a1; ++i) {
                                        for (uint32_t j = i; j < block.a1; ++j) f(i, j);
                                }
                        } else if (block.b1 - block.b0 >= block.a1 - block.a0) {
                                for (uint32_t i = block.a0; i < block.a1; ++i) {
                                        for (uint32_t j = block.b0; j < block.b1; ++j) f(i, j);
                                }
                        } else {
                                for (uint32_t j = block.b0; j < block.b1; ++j) {
                                        for (uint32_t i = block.a0; i < block.a1; ++i) f(i, j);
                                }
                        }
                }
        }

        /// Exponent k of the phase i^k of row i of a times row j of b.
        template<bool OneWord>
        int phase(uint32_t i, uint32_t j) const noexcept {
                if constexpr (OneWord) {
                        const uint64_t* ra = aw_ + 2 * static_cast<size_t>(i);
                        const uint64_t* rb = bw_ + 2 * static_cast<size_t>(j);
                        const uint64_t x1y2 = ra[0] & rb[1];
                        const uint64_t anti = (rb[0] & ra[1]) ^ x1y2;
                        const uint64_t c2 = (ra[0] ^ rb[0] ^ ra[1] ^ rb[1] ^ x1y2) & anti;
                        return 2 * std::popcount(c2) - std::popcount(anti);
                } else {
                        return symplectic_product_phase(a_->x(i), a_->y(i), nwa_, b_->x(j), b_->y(j), nwb_);
                }
        }

        template<bool OneWord>
        bool anticommutes(uint32_t i, uint32_t j) const noexcept {
                if constexpr (OneWord) {
                        const uint64_t* ra = aw_ + 2 * static_cast<size_t>(i);
                        const uint64_t* rb = bw_ + 2 * static_cast<size_t>(j);
                        return std::popcount((ra[0] & rb[1]) ^ (rb[0] & ra[1])) & 1;
                } else {
                        return anticommute(a_->view(i), b_->view(j));
                }
        }

        /// Whether a pair contributes, given whether it anti-commutes.
        template<Mode M>
        static constexpr bool keeps(bool anticommuting) noexcept {
                if constexpr (M == Mode::Product || M == Mode::First) return true;
                else if constexpr (M == Mode::Square) return !anticommuting;
                else return anticommuting;
        }

        /// Whether equal operators keep their first pair instead of adding up.
        template<Mode M>
        static constexpr bool keeps_first() noexcept {
                return M == Mode::First || M == Mode::FirstAnticommuting;
        }

        /// Value of the pair (i, j) with phase exponent k: its coefficient, or
        /// (rank, phase exponent) when equal operators keep their first pair.
        template<Mode M>
        std::complex<double> pair_value(uint32_t i, uint32_t j, int k) const noexcept {
                if constexpr (keeps_first<M>()) {
                        const std::complex<double> left = a_->coeff(i);
                        const std::complex<double> right = b_->coeff(j);
                        const double rank = left.real() * static_cast<double>(n2_) + right.real();
                        const int phase = (static_cast<int>(left.imag()) + static_cast<int>(right.imag()) + k) & 3;
                        return {rank, static_cast<double>(phase)};
                } else {
                        return product_coeff(a_->coeff(i), b_->coeff(j), k, scale_of<M>(i, j));
                }
        }

        /// Adds the value (re2, im2) of a pair with the same operator to (re, im).
        template<Mode M>
        static void combine(double& re, double& im, double re2, double im2) noexcept {
                if constexpr (keeps_first<M>()) {
                        if (re2 < re) {
                                re = re2;
                                im = im2;
                        }
                } else {
                        re += re2;
                        im += im2;
                }
        }

        /// Factor of the pair: both orders of a commutator or of a square give the same term.
        template<Mode M>
        static constexpr double scale_of(uint32_t i, uint32_t j) noexcept {
                return (M == Mode::Commutator || (M == Mode::Square && i != j)) ? 2.0 : 1.0;
        }

        /// Appends the term coeff * (row i of a)(row j of b), given without its phase.
        void emit(uint32_t i, uint32_t j, double re, double im, Output& out) const noexcept {
                if (re == 0.0 && im == 0.0 && mode_ != Mode::First && mode_ != Mode::FirstAnticommuting) return;
                uint64_t* row = out.words + out.n * 2 * nw_;
                const uint64_t* ra = a_->x(i);
                const uint64_t* rb = b_->x(j);
                if (one_word_) {
                        row[0] = ra[0] ^ rb[0];
                        row[1] = ra[1] ^ rb[1];
                } else {
                        for (size_t k = 0; k < nw_; ++k) {
                                row[k] = (k < nwa_ ? ra[k] : 0) ^ (k < nwb_ ? rb[k] : 0);
                                row[nw_ + k] = (k < nwa_ ? ra[nwa_ + k] : 0) ^ (k < nwb_ ? rb[nwb_ + k] : 0);
                        }
                }
                out.coeffs[2 * out.n] = re;
                out.coeffs[2 * out.n + 1] = im;
                ++out.n;
        }

        // ------------------------------------------------------------------
        // Leaves
        // ------------------------------------------------------------------

        /// A task whose products differ in the lowest `bits` positions of the last
        /// level only: adds every coefficient to the entry of its operator.
        template<bool OneWord, Mode M>
        void sum_in_table(const Block* blocks, size_t nblocks, unsigned bits, Scratch& s, Output& out) const {
                const size_t level = levels_ - 1;
                const size_t entries = size_t{1} << bits;
                const unsigned shift = bits == 0 ? 0 : 64 - bits_in_level(level);
                const uint64_t mask = entries - 1;
                const uint64_t* ka = ka_.data() + level * n1_;
                const uint64_t* kb = kb_.data() + level * n2_;
                s.sums.assign(2 * entries, 0.0);
                s.factors.assign(2 * entries, kNone);
                double* const sums = s.sums.data();
                uint32_t* const factors = s.factors.data();
                for_pairs(blocks, nblocks, [&](uint32_t i, uint32_t j) {
                        const int k = phase<OneWord>(i, j);  // odd if the pair anti-commutes
                        if (!keeps<M>(k & 1)) return;
                        const size_t entry = static_cast<size_t>(((ka[i] ^ kb[j]) >> shift) & mask);
                        const std::complex<double> c = pair_value<M>(i, j, k);
                        if (keeps_first<M>() && factors[2 * entry] == kNone) {
                                sums[2 * entry] = c.real();
                                sums[2 * entry + 1] = c.imag();
                        } else {
                                combine<M>(sums[2 * entry], sums[2 * entry + 1], c.real(), c.imag());
                        }
                        factors[2 * entry] = i;
                        factors[2 * entry + 1] = j;
                });
                for (size_t entry = 0; entry < entries; ++entry) {
                        if (factors[2 * entry] != kNone) {
                                emit(factors[2 * entry], factors[2 * entry + 1], sums[2 * entry], sums[2 * entry + 1], out);
                        }
                }
        }

        /// Sorts by key. `nearly_sorted`: no record is far from its place.
        static void sort_records(Record* first, Record* last, bool nearly_sorted) {
                if (nearly_sorted || static_cast<size_t>(last - first) <= kInsertionLimit) {
                        for (Record* p = first + 1; p < last; ++p) {
                                if (p[-1].key <= p->key) continue;
                                const Record moved = *p;
                                Record* q = p;
                                do {
                                        *q = q[-1];
                                        --q;
                                } while (q > first && q[-1].key > moved.key);
                                *q = moved;
                        }
                } else {
                        std::sort(first, last, [](const Record& x, const Record& y) { return x.key < y.key; });
                }
        }

        /// Records of equal keys up to `level`; appends their terms.
        template<Mode M>
        void merge_equal(Record* first, Record* last, size_t level, Output& out) const {
                if (level == levels_) {  // the same operator
                        double re = first->re;
                        double im = first->im;
                        for (const Record* p = first + 1; p != last; ++p) combine<M>(re, im, p->re, p->im);
                        emit(first->i, first->j, re, im, out);
                        return;
                }
                const uint64_t* ka = ka_.data() + level * n1_;
                const uint64_t* kb = kb_.data() + level * n2_;
                for (Record* p = first; p != last; ++p) p->key = ka[p->i] ^ kb[p->j];
                sort_records(first, last, false);
                for (Record* p = first; p != last;) {
                        Record* next = p + 1;
                        while (next != last && next->key == p->key) ++next;
                        if (next - p == 1) {
                                emit(p->i, p->j, p->re, p->im, out);
                        } else {
                                merge_equal<M>(p, next, level + 1, out);
                        }
                        p = next;
                }
        }

        /// A task whose products differ in the lowest `bits` positions of `level` and
        /// in all later levels: forms the pairs grouped by their leading bits, sorts
        /// the groups and merges equal operators.
        template<bool OneWord, Mode M>
        void sort_and_merge(const Block* blocks, size_t nblocks, size_t pairs, size_t level, unsigned bits,
                            Scratch& s, Output& out) const {
                if (pairs > UINT32_MAX) throw std::length_error("QubitHamiltonian: product too large");
                const unsigned scatter = std::min(bits, std::clamp(static_cast<unsigned>(std::bit_width(pairs)), 3u, kMaxScatterBits + 2) - 2);
                const size_t groups = size_t{1} << scatter;
                const unsigned shift = 64 - bits_in_level(level) + bits - scatter;
                const uint64_t mask = groups - 1;
                const uint64_t* ka = ka_.data() + level * n1_;
                const uint64_t* kb = kb_.data() + level * n2_;

                // First the group sizes, then every record straight to its place.
                s.ends.assign(groups, 0);
                uint32_t* const ends = s.ends.data();
                // Where only some pairs count (commutators, squares), they are picked out
                // here, so the second pass visits them alone.
                constexpr bool filtered = !keeps<M>(false) || !keeps<M>(true);
                if constexpr (filtered) s.kept.clear();
                for_pairs(blocks, nblocks, [&](uint32_t i, uint32_t j) {
                        if constexpr (filtered) {
                                if (!keeps<M>(anticommutes<OneWord>(i, j))) return;
                                s.kept.push_back((static_cast<uint64_t>(i) << 32) | j);
                        }
                        ++ends[((ka[i] ^ kb[j]) >> shift) & mask];
                });
                uint32_t count = 0;
                uint32_t largest = 0;
                for (size_t g = 0; g < groups; ++g) {
                        const uint32_t size = ends[g];
                        ends[g] = count;
                        count += size;
                        largest = std::max(largest, size);
                }
                if (count == 0) return;
                if (s.record_cap < count) {
                        s.records.reset(new Record[count]);
                        s.record_cap = count;
                }
                Record* const records = s.records.get();
                const auto place = [&](uint32_t i, uint32_t j) {
                        const uint64_t key = ka[i] ^ kb[j];
                        const std::complex<double> c = pair_value<M>(i, j, phase<OneWord>(i, j));
                        // Field by field: a whole record built first would be read back
                        // from the stack in pieces other than it was written in.
                        Record& record = records[ends[(key >> shift) & mask]++];
                        record.key = key;
                        record.re = c.real();
                        record.im = c.imag();
                        record.i = i;
                        record.j = j;
                };
                if constexpr (filtered) {
                        for (const uint64_t pair : s.kept) place(static_cast<uint32_t>(pair >> 32), static_cast<uint32_t>(pair));
                } else {
                        for_pairs(blocks, nblocks, place);
                }

                // The groups are in order; what is left is the order within each.
                if (scatter != bits) {
                        if (largest <= kInsertionLimit) {
                                sort_records(records, records + count, true);
                        } else {
                                uint32_t begin = 0;
                                for (size_t g = 0; g < groups; ++g) {
                                        if (ends[g] - begin > 1) sort_records(records + begin, records + ends[g], false);
                                        begin = ends[g];
                                }
                        }
                }
                const bool last_level = level + 1 == levels_;
                for (uint32_t r = 0; r < count;) {
                        uint32_t next = r + 1;
                        while (next < count && records[next].key == records[r].key) ++next;
                        if (next - r == 1) {
                                emit(records[r].i, records[r].j, records[r].re, records[r].im, out);
                        } else if (last_level) {
                                double re = records[r].re;
                                double im = records[r].im;
                                for (uint32_t q = r + 1; q < next; ++q) combine<M>(re, im, records[q].re, records[q].im);
                                emit(records[r].i, records[r].j, re, im, out);
                        } else {
                                merge_equal<M>(records + r, records + next, level + 1, out);
                        }
                        r = next;
                }
        }

        template<bool OneWord, Mode M>
        void process(const Block* blocks, size_t nblocks, size_t pos, size_t pairs, Scratch& s, Output& out) const {
                if (pos == nvar_) {
                        sum_in_table<OneWord, M>(blocks, nblocks, 0, s, out);
                        return;
                }
                const size_t level = pos / 64;
                const unsigned bits = bits_in_level(level) - static_cast<unsigned>(pos % 64);
                if (level + 1 == levels_ && bits <= kMaxTableBits && (size_t{1} << bits) <= 2 * pairs) {
                        sum_in_table<OneWord, M>(blocks, nblocks, bits, s, out);
                        return;
                }
                if (pairs <= kLeafPairs || !can_split(nblocks, pairs)) {
                        sort_and_merge<OneWord, M>(blocks, nblocks, pairs, level, bits, s, out);
                        return;
                }
                const unsigned d = digit_bits(nblocks, pairs, bits);
                std::vector<std::vector<Block>>& children = s.children[pos];
                split(blocks, nblocks, pos, d, children, s.runs_a, s.runs_b);
                for (size_t c = 0; c < (size_t{1} << d); ++c) {
                        if (!children[c].empty()) {
                                process<OneWord, M>(children[c].data(), children[c].size(), pos + d, count_pairs(children[c]), s, out);
                        }
                }
        }

        template<bool OneWord>
        void process_as(const Block* blocks, size_t nblocks, const Unit& unit, Scratch& s, Output& out) const {
                switch (mode_) {
                        case Mode::Product: process<OneWord, Mode::Product>(blocks, nblocks, unit.pos, unit.pairs, s, out); break;
                        case Mode::Commutator: process<OneWord, Mode::Commutator>(blocks, nblocks, unit.pos, unit.pairs, s, out); break;
                        case Mode::Square: process<OneWord, Mode::Square>(blocks, nblocks, unit.pos, unit.pairs, s, out); break;
                        case Mode::First: process<OneWord, Mode::First>(blocks, nblocks, unit.pos, unit.pairs, s, out); break;
                        case Mode::FirstAnticommuting:
                                process<OneWord, Mode::FirstAnticommuting>(blocks, nblocks, unit.pos, unit.pairs, s, out);
                                break;
                }
        }

        void process_unit(Unit& unit, Scratch& s, Output& out) const {
                const Block* blocks = unit.blocks.data();
                const size_t nblocks = unit.blocks.size();
                if (one_word_) process_as<true>(blocks, nblocks, unit, s, out);
                else process_as<false>(blocks, nblocks, unit, s, out);
                std::vector<Block>().swap(unit.blocks);
                unit.rows = out.n;
        }

        /// Bound on the number of terms of a task: one per pair and one per operator.
        size_t max_rows(size_t pos, size_t pairs) const noexcept {
                const size_t remaining = nvar_ - pos;
                return remaining < 63 ? std::min(pairs, size_t{1} << remaining) : pairs;
        }

        FlatTerms result() const {
                const size_t total = mode_ == Mode::Square ? n1_ * (n1_ + 1) / 2 : n1_ * n2_;
                // Independent of the number of threads, so that the result is too.
                const size_t unit_pairs = std::clamp(total / 64, kLeafPairs, kMaxUnitPairs);
                std::vector<Unit> units;
                {
                        std::vector<Run> runs_a, runs_b;
                        std::vector<Block> root{{0, static_cast<uint32_t>(n1_), 0, static_cast<uint32_t>(n2_)}};
                        expand(std::move(root), 0, total, unit_pairs, units, runs_a, runs_b);
                }
                size_t capacity = 0;
                for (Unit& unit : units) {
                        unit.max_rows = max_rows(unit.pos, unit.pairs);
                        unit.offset = capacity;
                        capacity += unit.max_rows;
                }
                // Room for every term there can be. Only the rows written are ever
                // touched, so the rest costs address space, not memory.
                FlatTerms out(capacity, nw_);
                const std::ptrdiff_t count = static_cast<std::ptrdiff_t>(units.size());  // signed for OpenMP 2.0
#ifdef PAULIENGINE_HAS_OPENMP
                // Two or three units are done sooner in turn than by a team of threads.
                if (count >= 4 && use_threads(total * std::max<size_t>(nw_, 1))) {
                        // Units finish in any order, so each writes to the range its bound
                        // reserves. The threads share nothing but the list of units.
                        #pragma omp parallel
                        {
                                Scratch scratch(nvar_);
                                #pragma omp for schedule(dynamic)
                                for (std::ptrdiff_t u = 0; u < count; ++u) {
                                        Unit& unit = units[u];
                                        Output unit_out{out.x(unit.offset), out.coeff_data(unit.offset), 0};
                                        process_unit(unit, scratch, unit_out);
                                }
                        }
                        size_t rows = 0;
                        for (Unit& unit : units) {
                                unit.first = rows;
                                rows += unit.rows;
                        }
                        // Every unit reached its bound (no two products were the same
                        // operator, none vanished): the ranges are the result as they are.
                        if (rows == capacity) return out;
                        // Otherwise the gaps are closed. Small gaps, or few terms: every
                        // unit moves down a little in turn, within memory that is in use
                        // already. Large gaps: the units are copied to a new block by all
                        // threads, which is faster than one thread moving them far.
                        if (capacity - rows <= rows / 4 || rows < kMoveInPlaceRows) {
                                for (const Unit& unit : units) {
                                        if (unit.first != unit.offset) out.move_rows(unit.first, unit.offset, unit.rows);
                                }
                                out.shrink_to(rows);
                                return out;
                        }
                        FlatTerms packed(rows, nw_);
                        for_each_index(units.size(), true, [&](size_t u) {
                                const Unit& unit = units[u];
                                packed.assign_rows(unit.first, out.x(unit.offset), out.coeff_data(unit.offset), unit.rows);
                        });
                        return packed;
                }
#endif
                // In order, every unit knows its place: written where it stays.
                size_t rows = 0;
                Scratch scratch(nvar_);
                for (std::ptrdiff_t u = 0; u < count; ++u) {
                        Output unit_out{out.x(rows), out.coeff_data(rows), 0};
                        process_unit(units[u], scratch, unit_out);
                        rows += unit_out.n;
                }
                out.shrink_to(rows);
                return out;
        }

        Mode mode_;
        size_t nw_;
        const FlatTerms* a_ = nullptr;
        const FlatTerms* b_ = nullptr;
        FlatTerms a_sorted_, b_sorted_;
        size_t n1_ = 0, n2_ = 0;
        size_t nwa_ = 0, nwb_ = 0;
        bool one_word_ = false;
        const uint64_t* aw_ = nullptr;
        const uint64_t* bw_ = nullptr;
        size_t nvar_ = 0;
        size_t levels_ = 1;
        std::vector<uint64_t> ka_, kb_;
};

}  // namespace pauliengine_detail
