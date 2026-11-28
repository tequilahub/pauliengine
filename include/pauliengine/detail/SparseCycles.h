#pragma once
/// Internal: products and commutators of PauliCycleSums whose bases act on few
/// sites, computed on sorted site lists instead of words.

#include <algorithm>
#include <bit>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <vector>

#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/detail/Parallel.h"
#include "pauliengine/PauliString.h"

namespace pauliengine_detail {

/// The sites of a string in ascending order with their Paulis (x | y << 1).
struct SparseString {
        std::vector<uint32_t> pos;
        std::vector<uint8_t> pauli;

        static SparseString of(const WordVec& x, const WordVec& y) {
                SparseString s;
                const size_t words = std::min(x.size(), y.size());
                for (size_t w = 0; w < words; ++w) {
                        uint64_t bits = x[w] | y[w];
                        while (bits) {
                                const int b = std::countr_zero(bits);
                                s.pos.push_back(static_cast<uint32_t>(w * 64 + static_cast<size_t>(b)));
                                s.pauli.push_back(static_cast<uint8_t>(((x[w] >> b) & 1ULL) | (((y[w] >> b) & 1ULL) << 1)));
                                bits &= bits - 1;
                        }
                }
                return s;
        }
};

/// sum_ij a_i b_j or sum_ij [a_i, b_j] for cycles a_i, b_j with sparse bases.
///
/// The cycle product a b = sum_s a R^s(b) and the commutator are formed shift by
/// shift as in PauliCycle, but on site lists: R^s(b) is b's list cut at one place,
/// the product a merge of two lists, O(weight) instead of O(n / 64) per shift. The
/// rotation class of a product is the least rotation of its (gap to the next site,
/// Pauli) items, which is what PauliCycle's unified form is built from, so it is
/// found without forming words. The summands are merged in flat tables: first those
/// of one pair (i, j), then the pairs in order, with the same additions in the same
/// order as merging the PauliCycleSums of the pairs, so the summands, their order and
/// their coefficients are the same. No words are formed: a summand is its class.
///
/// With several threads the pairs are formed in parallel and the classes are split
/// over partitions by their hash; every partition takes its records in pair order,
/// so sums and first occurrences, and with them the result, do not depend on the
/// number of threads.
template <typename Coeff>
class SparseCycleProducts {
    public:
        enum class Op { Product, Commutator };

        /// A summand of the result: its class (canonical items at items[off, off + len))
        /// and its merged coefficient.
        struct Summand {
                uint32_t len;
                size_t off;
                Coeff coeff;
        };
        std::vector<Summand> summands;  ///< in the order of the merged sum
        std::vector<uint64_t> items;    ///< canonical items of the summands

        /// Bases act on few sites for their length: the site lists are cheaper than
        /// the words. The inverse of PauliCycle::prefers_signature().
        static bool suitable(size_t weight, size_t n_words) { return weight <= 4 * n_words + 16; }

        /// Site visits per unit of work for use_threads().
        static constexpr size_t kVisitsPerUnit = 64;

        void run(const std::vector<SparseString>& a, const std::vector<Coeff>& ca,
                 const std::vector<SparseString>& b, const std::vector<Coeff>& cb, uint32_t n, Op op) {
                summands.clear();
                items.clear();
                if (a.empty() || b.empty()) {
                        return;
                }
                // Site visits, a bound on the work: n shifts per pair for the product,
                // at most |supp a| |supp b| for the commutator; kVisitsPerUnit visits
                // count as one unit of use_threads().
                size_t visits = 0;
                for (const auto& x : a) {
                        for (const auto& y : b) {
                                const size_t shifts = op == Op::Product ? n : std::min<size_t>(n, x.pos.size() * y.pos.size());
                                visits += shifts * (x.pos.size() + y.pos.size() + 1);
                        }
                }
                const size_t pairs = a.size() * b.size();
                if (pairs < 2 || !use_threads(visits / kVisitsPerUnit)) {
                        run_serial(a, ca, b, cb, n, op);
                } else {
                        run_parallel(a, ca, b, cb, n, op);
                }
        }

        /// Unified words (see PauliCycle::unify_words) from canonical items.
        static void unified_words(const uint64_t* it, size_t len, uint64_t* ux, uint64_t* uy, size_t nw) {
                std::fill(ux, ux + nw, 0ULL);
                std::fill(uy, uy + nw, 0ULL);
                uint64_t q = 0;
                for (size_t t = 0; t < len; ++t) {
                        if (it[t] & 1ULL) ux[q >> 6] |= 1ULL << (q & 63);
                        if (it[t] & 2ULL) uy[q >> 6] |= 1ULL << (q & 63);
                        q += it[t] >> 2;
                }
        }

    private:
        struct Phases {
                int power[4][4];
                bool anti[4][4];
                Phases() {
                        for (int u = 0; u < 4; ++u) {
                                for (int v = 0; v < 4; ++v) {
                                        const uint64_t xa = u & 1, ya = u >> 1, xb = v & 1, yb = v >> 1;
                                        uint64_t ox = 0, oy = 0;
                                        power[u][v] = multiply_symplectic_words(&xa, &ya, 1, &xb, &yb, 1, &ox, &oy);
                                        anti[u][v] = ((xa & yb) ^ (xb & ya)) != 0;
                                }
                        }
                }
        };
        static const Phases& phases() {
                static const Phases p;
                return p;
        }

        /// One product a R^s(b): its sites, the phase exponent, whether a and R^s(b)
        /// anticommute, and its canonical items.
        struct Product {
                std::vector<uint32_t> pos;
                std::vector<uint8_t> pauli;
                std::vector<uint64_t> raw;
                std::vector<uint64_t> items;
                int power = 0;
                bool anti = false;
        };

        static void multiply(const SparseString& a, const SparseString& b, uint32_t n, uint32_t s, Product& out) {
                const Phases& ph = phases();
                out.pos.clear();
                out.pauli.clear();
                int power = 0;
                bool anti = false;
                const size_t wa = a.pos.size();
                const size_t wb = b.pos.size();
                // R^s(b) in ascending order: the sites that wrap around come first
                const size_t split = static_cast<size_t>(std::lower_bound(b.pos.begin(), b.pos.end(), n - s) - b.pos.begin());
                const size_t wrapped = wb - split;
                size_t ia = 0, kb = 0;
                while (ia < wa || kb < wb) {
                        const size_t idx = kb < wrapped ? split + kb : kb - wrapped;
                        const uint32_t pa = ia < wa ? a.pos[ia] : UINT32_MAX;
                        const uint32_t pb = kb < wb ? (kb < wrapped ? b.pos[idx] + s - n : b.pos[idx] + s) : UINT32_MAX;
                        if (pa < pb) {
                                out.pos.push_back(pa);
                                out.pauli.push_back(a.pauli[ia++]);
                        } else if (pb < pa) {
                                out.pos.push_back(pb);
                                out.pauli.push_back(b.pauli[idx]);
                                ++kb;
                        } else {
                                const uint8_t u = a.pauli[ia++], v = b.pauli[idx];
                                ++kb;
                                power += ph.power[u][v];
                                anti ^= ph.anti[u][v];
                                if (u ^ v) {
                                        out.pos.push_back(pa);
                                        out.pauli.push_back(static_cast<uint8_t>(u ^ v));
                                }
                        }
                }
                out.power = power;
                out.anti = anti;
        }

        /// Least rotation of the product's items, as in PauliCycle::least_rotation.
        static void canonical(uint32_t n, Product& p) {
                const size_t w = p.pos.size();
                p.raw.resize(w);
                for (size_t k = 0; k + 1 < w; ++k) {
                        p.raw[k] = (static_cast<uint64_t>(p.pos[k + 1] - p.pos[k]) << 2) | p.pauli[k];
                }
                if (w) {
                        p.raw[w - 1] = (static_cast<uint64_t>(p.pos[0] + n - p.pos[w - 1]) << 2) | p.pauli[w - 1];
                }
                size_t i = 0, j = 1, k = 0;
                while (i < w && j < w && k < w) {
                        const size_t ia = i + k < w ? i + k : i + k - w;
                        const size_t jb = j + k < w ? j + k : j + k - w;
                        const uint64_t x = p.raw[ia], y = p.raw[jb];
                        if (x == y) {
                                ++k;
                                continue;
                        }
                        if (x > y) i += k + 1; else j += k + 1;
                        if (i == j) ++j;
                        k = 0;
                }
                const size_t start = w ? std::min(i, j) : 0;
                p.items.resize(w);
                for (size_t t = 0; t < w; ++t) {
                        p.items[t] = p.raw[start + t < w ? start + t : start + t - w];
                }
        }

        static uint64_t hash_items(const uint64_t* v, size_t m) {
                uint64_t h = 0x9e3779b97f4a7c15ULL ^ m;
                for (size_t k = 0; k < m; ++k) {
                        h ^= v[k] + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                }
                h ^= h >> 33;
                h *= 0xff51afd7ed558ccdULL;
                h ^= h >> 33;
                return h;
        }

        /// Classes with their merged coefficient, in order of first occurrence.
        struct Table {
                struct Entry {
                        uint64_t hash;
                        size_t off;
                        uint32_t len;
                        uint64_t seq;  // position of the first occurrence in the stream of pair summands
                        Coeff coeff;
                };
                std::vector<Entry> entries;
                std::vector<uint64_t> arena;
                std::vector<int64_t> slots;
                size_t mask = 0;

                void reset(size_t expected) {
                        entries.clear();
                        arena.clear();
                        size_t cap = 16;
                        while (cap < 2 * expected + 2) cap <<= 1;
                        if (slots.size() != cap) {
                                slots.assign(cap, -1);
                        } else {
                                std::fill(slots.begin(), slots.end(), int64_t{-1});
                        }
                        mask = cap - 1;
                }
                void add(const uint64_t* it, uint32_t len, uint64_t h, const Coeff& c, uint64_t seq) {
                        size_t p = h & mask;
                        while (slots[p] >= 0) {
                                Entry& e = entries[static_cast<size_t>(slots[p])];
                                if (e.hash == h && e.len == len && std::equal(it, it + len, arena.data() + e.off)) {
                                        e.coeff = e.coeff + c;
                                        return;
                                }
                                p = (p + 1) & mask;
                        }
                        slots[p] = static_cast<int64_t>(entries.size());
                        entries.push_back({h, arena.size(), len, seq, c});
                        arena.insert(arena.end(), it, it + len);
                        if (2 * entries.size() > mask + 1) grow();
                }
                void grow() {
                        const size_t cap = 2 * (mask + 1);
                        slots.assign(cap, -1);
                        mask = cap - 1;
                        for (size_t e = 0; e < entries.size(); ++e) {
                                size_t p = entries[e].hash & mask;
                                while (slots[p] >= 0) p = (p + 1) & mask;
                                slots[p] = static_cast<int64_t>(e);
                        }
                }
        };

        /// The shifts in the order of PauliCycle::overlap_shifts (commutator) or 0..n-1.
        static void shifts_of(const SparseString& a, const SparseString& b, uint32_t n, Op op,
                              std::vector<uint32_t>& shifts, std::vector<uint8_t>& seen) {
                shifts.clear();
                if (op == Op::Product) {
                        for (uint32_t s = 0; s < n; ++s) shifts.push_back(s);
                        return;
                }
                seen.assign(n, 0);
                for (uint32_t qp : a.pos) {
                        for (uint32_t qq : b.pos) {
                                const uint32_t k = qp >= qq ? qp - qq : qp + n - qq;
                                if (!seen[k]) {
                                        seen[k] = 1;
                                        shifts.push_back(k);
                                        if (shifts.size() == n) return;
                                }
                        }
                }
        }

        /// Scratch of one thread for forming pairs.
        struct Scratch {
                Product p;
                Table pair;
                std::vector<uint32_t> shifts;
                std::vector<uint8_t> seen;
        };

        /// The summands of pair (i, j), merged within the pair, into s.pair (zeros kept;
        /// the caller skips them as PauliCycleSum's constructor drops them).
        static void form_pair(const SparseString& a, const Coeff& ca, const SparseString& b, const Coeff& cb,
                              uint32_t n, Op op, Scratch& s) {
                shifts_of(a, b, n, op, s.shifts, s.seen);
                s.pair.reset(s.shifts.size());
                const Coeff ab = ca * cb;  // as in PauliCycle::operator*
                for (uint32_t k : s.shifts) {
                        multiply(a, b, n, k, s.p);
                        Coeff c;
                        if (op == Op::Product) {
                                c = ab;
                                if (s.p.power & 1) c *= std::complex<double>(0.0, 1.0);
                                if (s.p.power & 2) c = -c;
                        } else {
                                if (!s.p.anti) continue;  // commuting: no commutator term
                                c = PauliString<Coeff>::commutator_coeff(ca, cb, s.p.power);
                                if (coeff_is_zero(c)) continue;
                        }
                        canonical(n, s.p);
                        const uint32_t len = static_cast<uint32_t>(s.p.items.size());
                        s.pair.add(s.p.items.data(), len, hash_items(s.p.items.data(), len), c, 0);
                }
        }

        void finish(const Table& t) {
                for (const auto& e : t.entries) {
                        if (coeff_is_zero(e.coeff)) continue;
                        summands.push_back({e.len, items.size(), e.coeff});
                        items.insert(items.end(), t.arena.begin() + static_cast<std::ptrdiff_t>(e.off),
                                     t.arena.begin() + static_cast<std::ptrdiff_t>(e.off + e.len));
                }
        }

        void run_serial(const std::vector<SparseString>& a, const std::vector<Coeff>& ca,
                        const std::vector<SparseString>& b, const std::vector<Coeff>& cb, uint32_t n, Op op) {
                Scratch s;
                Table global;
                global.reset(1024);
                for (uint32_t i = 0; i < a.size(); ++i) {
                        for (uint32_t j = 0; j < b.size(); ++j) {
                                form_pair(a[i], ca[i], b[j], cb[j], n, op, s);
                                for (const auto& e : s.pair.entries) {
                                        if (coeff_is_zero(e.coeff)) continue;  // dropped within the pair
                                        global.add(s.pair.arena.data() + e.off, e.len, e.hash, e.coeff, 0);
                                }
                        }
                }
                finish(global);
        }

        /// A summand of one pair, handed from the pair stage to the partitions.
        struct Record {
                uint64_t hash;
                uint32_t len;
                size_t off;  // into the arena of the thread that formed it
                Coeff coeff;
        };

        void run_parallel(const std::vector<SparseString>& a, const std::vector<Coeff>& ca,
                          const std::vector<SparseString>& b, const std::vector<Coeff>& cb, uint32_t n, Op op) {
#ifdef PAULIENGINE_HAS_OPENMP
                const size_t nb = b.size(), pairs = a.size() * nb;
                const size_t threads = static_cast<size_t>(omp_get_max_threads());
                // partitions: a few per thread, a power of two, chosen by the top hash bits
                unsigned bits = 0;
                while ((size_t{1} << bits) < 4 * threads && bits < 10) ++bits;
                const size_t parts = size_t{1} << bits;
                std::vector<Table> tables(parts);
                for (auto& t : tables) t.reset(1024);
                // Pairs are taken in chunks, so the records of one chunk stay in memory.
                const size_t per_pair = op == Op::Product ? n : 8;
                const size_t chunk = std::max(threads, (size_t{1} << 20) / std::max<size_t>(per_pair, 1));
                std::vector<std::vector<Record>> recs(threads);
                std::vector<std::vector<uint64_t>> arenas(threads);
                std::vector<Scratch> scratch(threads);
                std::vector<std::vector<uint32_t>> bucket(threads * parts);  // [t * parts + q]: records of thread t
                uint64_t seq = 0;  // records so far
                for (size_t first = 0; first < pairs; first += chunk) {
                        const size_t last = std::min(pairs, first + chunk);
                        const std::ptrdiff_t count = static_cast<std::ptrdiff_t>(last - first);
                        #pragma omp parallel
                        {
                                const size_t t = static_cast<size_t>(omp_get_thread_num());
                                auto& mine = recs[t];
                                auto& arena = arenas[t];
                                mine.clear();
                                arena.clear();
                                Scratch& s = scratch[t];
                                // static schedule: thread t takes ascending pairs, so the
                                // records of all threads in thread order are in pair order
                                #pragma omp for schedule(static)
                                for (std::ptrdiff_t q = 0; q < count; ++q) {
                                        const size_t pair = first + static_cast<size_t>(q);
                                        form_pair(a[pair / nb], ca[pair / nb], b[pair % nb], cb[pair % nb], n, op, s);
                                        for (const auto& e : s.pair.entries) {
                                                if (coeff_is_zero(e.coeff)) continue;  // dropped within the pair
                                                mine.push_back({e.hash, e.len, arena.size(), e.coeff});
                                                arena.insert(arena.end(), s.pair.arena.begin() + static_cast<std::ptrdiff_t>(e.off),
                                                             s.pair.arena.begin() + static_cast<std::ptrdiff_t>(e.off + e.len));
                                        }
                                }
                                // the thread's own records to the partitions, in order
                                for (size_t q = 0; q < parts; ++q) bucket[t * parts + q].clear();
                                for (size_t r = 0; r < mine.size(); ++r) {
                                        bucket[t * parts + (mine[r].hash >> (64 - bits))].push_back(static_cast<uint32_t>(r));
                                }
                        }
                        // The records of all threads in thread order are the stream of pair
                        // summands; offset[t] is where those of thread t start in it.
                        std::vector<uint64_t> offset(threads + 1, seq);
                        for (size_t t = 0; t < threads; ++t) offset[t + 1] = offset[t] + recs[t].size();
                        seq = offset[threads];
                        // every partition takes its records in stream order
                        for_each_index(parts, true, [&](size_t q) {
                                Table& tab = tables[q];
                                for (size_t t = 0; t < threads; ++t) {
                                        for (const uint32_t r : bucket[t * parts + q]) {
                                                const Record& rec = recs[t][r];
                                                tab.add(arenas[t].data() + rec.off, rec.len, rec.hash, rec.coeff, offset[t] + r);
                                        }
                                }
                        });
                }
                // The classes in order of first occurrence: every range of the stream is
                // collected from all partitions (each in stream order) and sorted.
                const size_t ranges = 4 * threads;
                std::vector<std::vector<Summand>> part_summands(ranges);
                std::vector<std::vector<uint64_t>> part_items(ranges);
                for_each_index(ranges, true, [&](size_t g) {
                        const uint64_t lo = seq * g / ranges, hi = seq * (g + 1) / ranges;
                        std::vector<std::pair<uint64_t, std::pair<size_t, size_t>>> mine;  // (seq, (partition, entry))
                        for (size_t q = 0; q < parts; ++q) {
                                const auto& entries = tables[q].entries;
                                auto it = std::lower_bound(entries.begin(), entries.end(), lo,
                                                           [](const auto& e, uint64_t v) { return e.seq < v; });
                                for (; it != entries.end() && it->seq < hi; ++it) {
                                        if (!coeff_is_zero(it->coeff)) {
                                                mine.push_back({it->seq, {q, static_cast<size_t>(it - entries.begin())}});
                                        }
                                }
                        }
                        std::sort(mine.begin(), mine.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
                        auto& out = part_summands[g];
                        auto& its = part_items[g];
                        out.reserve(mine.size());
                        for (const auto& [sq, at] : mine) {
                                const auto& en = tables[at.first].entries[at.second];
                                out.push_back({en.len, its.size(), en.coeff});
                                its.insert(its.end(), tables[at.first].arena.begin() + static_cast<std::ptrdiff_t>(en.off),
                                           tables[at.first].arena.begin() + static_cast<std::ptrdiff_t>(en.off + en.len));
                        }
                });
                size_t total = 0, total_items = 0;
                for (size_t g = 0; g < ranges; ++g) {
                        total += part_summands[g].size();
                        total_items += part_items[g].size();
                }
                summands.reserve(total);
                items.reserve(total_items);
                for (size_t g = 0; g < ranges; ++g) {
                        const size_t base = items.size();
                        for (auto sm : part_summands[g]) {
                                sm.off += base;
                                summands.push_back(sm);
                        }
                        items.insert(items.end(), part_items[g].begin(), part_items[g].end());
                }
#else
                run_serial(a, ca, b, cb, n, op);
#endif
        }
};

}  // namespace pauliengine_detail
