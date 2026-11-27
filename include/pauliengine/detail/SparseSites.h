#pragma once
/// Internal: Pauli strings as sorted lists of sites, and the accumulation of many
/// such terms into a compacted sum. Used by SparsePauliString and
/// SparsePauliStringSum, whose cost depends on the weight of the strings
/// instead of on the highest qubit they act on.

#include <algorithm>
#include <bit>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "pauliengine/PauliString.h"
#include "pauliengine/detail/Coefficients.h"
#include "pauliengine/detail/Parallel.h"

namespace pauliengine_detail {

/// A site is (qubit << 2) | pauli with pauli = x | y << 1, i.e. X = 1, Y = 2, Z = 3,
/// the bits of the symplectic form of PauliString. A string lists its sites by
/// ascending qubit; identities are not listed.
using Site = uint32_t;

/// Largest qubit index a site can hold.
inline constexpr uint32_t kMaxSparseQubit = (uint32_t{1} << 30) - 1;

inline uint32_t site_qubit(Site s) { return s >> 2; }
inline uint32_t site_pauli(Site s) { return s & 3u; }
inline Site make_site(uint32_t qubit, uint32_t pauli) { return (qubit << 2) | pauli; }

/// Exponent of the phase i^k of the single-qubit product p * q (0 for equal or identity factors).
inline constexpr int kSitePhase[4][4] = {
        {0, 0, 0, 0},
        {0, 0, 1, 3},   // X*Y = iZ, X*Z = -iY
        {0, 3, 0, 1},   // Y*X = -iZ, Y*Z = iX
        {0, 1, 3, 0},   // Z*X = iY, Z*Y = -iX
};

/// Pauli letter of a site's pauli bits.
inline char pauli_letter(uint32_t p) {
        static constexpr char kLetters[4] = {'I', 'X', 'Y', 'Z'};
        return kLetters[p & 3u];
}

/// Pauli bits of a letter; throws for anything else than I, X, Y, Z.
inline uint32_t pauli_bits(char c) {
        switch (c) {
                case 'I': case 'i': return 0;
                case 'X': case 'x': return 1;
                case 'Y': case 'y': return 2;
                case 'Z': case 'z': return 3;
                default: throw std::invalid_argument(std::string("invalid Pauli operator '") + c + "'");
        }
}

/// Product of two strings: writes the sites of the operator part to `out` (room for
/// na + nb sites) and returns the exponent k of the phase i^k. Odd k means the
/// strings anti-commute.
inline int multiply_sites(const Site* a, size_t na, const Site* b, size_t nb, Site* out, size_t& nout) {
        // Branch-free merge: for qubits in random order the three-way branch of a
        // plain merge is mispredicted on most steps.
        size_t i = 0, j = 0, k = 0;
        unsigned phase = 0;
        while (i < na && j < nb) {
                const Site x = a[i];
                const Site y = b[j];
                const uint32_t qa = site_qubit(x);
                const uint32_t qb = site_qubit(y);
                const bool lt = qa < qb;
                const bool gt = qb < qa;
                const bool eq = !(lt | gt);
                const uint32_t p = (x ^ y) & 3u;                          // product Pauli on a shared qubit
                phase += eq ? static_cast<unsigned>(kSitePhase[x & 3u][y & 3u]) : 0u;
                out[k] = lt ? x : (gt ? y : (x & ~3u) | p);
                k += static_cast<size_t>(!(eq & (p == 0)));
                i += static_cast<size_t>(!gt);
                j += static_cast<size_t>(!lt);
        }
        while (i < na) out[k++] = a[i++];
        while (j < nb) out[k++] = b[j++];
        nout = k;
        return static_cast<int>(phase & 3u);
}

/// Exponent of the phase of the product, without forming it.
inline int product_phase_sites(const Site* a, size_t na, const Site* b, size_t nb) {
        size_t i = 0, j = 0;
        unsigned phase = 0;
        while (i < na && j < nb) {
                const uint32_t qa = site_qubit(a[i]);
                const uint32_t qb = site_qubit(b[j]);
                phase += qa == qb ? static_cast<unsigned>(kSitePhase[a[i] & 3u][b[j] & 3u]) : 0u;
                i += static_cast<size_t>(qa <= qb);
                j += static_cast<size_t>(qb <= qa);
        }
        return static_cast<int>(phase & 3u);
}

/// Whether two strings anti-commute: an odd number of shared qubits carry different Paulis.
inline bool anticommute_sites(const Site* a, size_t na, const Site* b, size_t nb) {
        size_t i = 0, j = 0;
        unsigned odd = 0;
        while (i < na && j < nb) {
                const uint32_t qa = site_qubit(a[i]);
                const uint32_t qb = site_qubit(b[j]);
                odd ^= static_cast<unsigned>((qa == qb) & (a[i] != b[j]));
                i += static_cast<size_t>(qa <= qb);
                j += static_cast<size_t>(qb <= qa);
        }
        return odd != 0;
}

/// Hash of an operator given by its sites; equal operators hash alike.
inline uint64_t hash_sites(const Site* s, size_t n) {
        uint64_t h = 0x9e3779b97f4a7c15ULL ^ (static_cast<uint64_t>(n) * 0xff51afd7ed558ccdULL);
        for (size_t k = 0; k < n; ++k) {
                h = (h ^ s[k]) * 0x100000001b3ULL;
                h ^= h >> 29;
        }
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdULL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53ULL;
        h ^= h >> 33;
        return h;
}

/// The sites of the symplectic words x, y (n words each).
inline void sites_of_words(const uint64_t* x, const uint64_t* y, size_t n, std::vector<Site>& out) {
        out.clear();
        for (size_t w = 0; w < n; ++w) {
                uint64_t bits = x[w] | y[w];
                while (bits) {
                        const unsigned b = static_cast<unsigned>(std::countr_zero(bits));
                        const uint64_t q = 64 * w + b;
                        if (q > kMaxSparseQubit) throw std::length_error("sparse Pauli string: qubit index beyond 2^30 - 1");
                        out.push_back(make_site(static_cast<uint32_t>(q),
                                static_cast<uint32_t>(((x[w] >> b) & 1u) | (((y[w] >> b) & 1u) << 1))));
                        bits &= bits - 1;
                }
        }
}

/// The symplectic words of a string given by its sites, trimmed after the last site.
inline void words_of_sites(const Site* s, size_t n, std::vector<uint64_t>& x, std::vector<uint64_t>& y) {
        const size_t nw = n == 0 ? 0 : site_qubit(s[n - 1]) / 64 + 1;
        x.assign(nw, 0);
        y.assign(nw, 0);
        for (size_t k = 0; k < n; ++k) {
                const uint32_t q = site_qubit(s[k]);
                const uint32_t p = site_pauli(s[k]);
                if (p & 1u) x[q / 64] |= uint64_t{1} << (q % 64);
                if (p & 2u) y[q / 64] |= uint64_t{1} << (q % 64);
        }
}

/// Whether the sites form a valid string: ascending qubits, no identities.
inline void check_sites(const Site* s, size_t n) {
        for (size_t k = 0; k < n; ++k) {
                if (site_pauli(s[k]) == 0) throw std::invalid_argument("sparse Pauli string: identity site");
                if (k > 0 && site_qubit(s[k]) <= site_qubit(s[k - 1])) {
                        throw std::invalid_argument("sparse Pauli string: sites must have ascending qubits");
                }
        }
}

/// When operators are stored sparse: an operator spanning `span` qubits (highest
/// qubit index + 1) whose heaviest string has `max_weight` sites is sparse if
/// span >= min_qubits and every string is light, max_weight <= max_density * span
/// (a single heavy string makes the whole operator dense). The
/// Python factories build sparse operators by this rule, and a sparse sum is turned
/// back into bit vectors once it misses the rule by more than kLeaveFactor (the
/// Python bindings return products and commutators with bit vectors anyway).
struct SparseThreshold {
        double min_qubits = 256.0;
        // 1/256: for weight-4 strings the sparse product only beats the bit vectors
        // from about 1024 qubits on (measured with 100 to 1000 terms, one thread).
        double max_density = 1.0 / 256.0;
        static constexpr double kLeaveFactor = 1.0;
};

/// The threshold in use (process-wide; set from Python with set_sparse_threshold).
inline SparseThreshold& sparse_threshold() {
        static SparseThreshold threshold;
        return threshold;
}

/// Whether an operator should be built sparse.
inline bool sparse_enters(double span, double max_weight) {
        const SparseThreshold& t = sparse_threshold();
        return span > 0 && span >= t.min_qubits && max_weight <= t.max_density * span;
}

/// Whether a sparse result should stay sparse (the rule, relaxed by kLeaveFactor;
/// the zero operator and the identity always stay).
inline bool sparse_stays(double span, double max_weight) {
        const SparseThreshold& t = sparse_threshold();
        if (max_weight == 0) return true;
        return span >= t.min_qubits / SparseThreshold::kLeaveFactor
            && max_weight <= SparseThreshold::kLeaveFactor * t.max_density * span;
}

/// scale * i^k * a * b for either coefficient type.
template <typename Coeff>
Coeff phase_product(const Coeff& a, const Coeff& b, int k, double scale) {
        if constexpr (std::is_same_v<Coeff, std::complex<double>>) {
                static constexpr double kRealSign[4] = {1.0, -1.0, -1.0, 1.0};
                static constexpr double kImagSign[4] = {1.0, 1.0, -1.0, -1.0};
                const double parts[2] = {a.real() * b.real() - a.imag() * b.imag(),
                                         a.real() * b.imag() + a.imag() * b.real()};
                k &= 3;
                return {scale * kRealSign[k] * parts[k & 1], scale * kImagSign[k] * parts[(k & 1) ^ 1]};
        } else {
                Coeff c = a * b;
                if (k & 1) c = c * Coeff(std::complex<double>(0.0, 1.0));
                if (k & 2) c = coeff_zero<Coeff>() - c;
                if (scale != 1.0) c = c * Coeff(scale);
                return c;
        }
}

/// Complex conjugate of a coefficient.
inline std::complex<double> conj_of(const std::complex<double>& c) { return std::conj(c); }
inline pauliengine::symbolic::Expression conj_of(const pauliengine::symbolic::Expression& c) {
        return c.conjugate();
}

/// Whether two coefficients are equal (polynomials are canonical).
inline bool coeff_equal_to(const std::complex<double>& a, const std::complex<double>& b) { return a == b; }
inline bool coeff_equal_to(const pauliengine::symbolic::Expression& a, const pauliengine::symbolic::Expression& b) {
        return a == b;
}

/// std::allocator whose value-less construction leaves trivial types uninitialised,
/// so that resizing a vector that is about to be overwritten does not write it twice.
template <typename T>
struct DefaultInitAllocator : std::allocator<T> {
        template <typename U>
        struct rebind {
                using other = DefaultInitAllocator<U>;
        };
        DefaultInitAllocator() noexcept = default;
        template <typename U>
        DefaultInitAllocator(const DefaultInitAllocator<U>&) noexcept {}
        template <typename U, typename... Args>
        void construct(U* p, Args&&... args) {
                if constexpr (sizeof...(Args) == 0) {
                        ::new (static_cast<void*>(p)) U;
                } else {
                        ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
                }
        }
};

/// The terms of a sum: term i has the sites [start[i], start[i + 1]) and coeffs[i].
template <typename Coeff>
struct SparseTerms {
        std::vector<uint64_t, DefaultInitAllocator<uint64_t>> start{0};
        std::vector<Site, DefaultInitAllocator<Site>> sites;
        std::vector<Coeff> coeffs;

        size_t size() const { return coeffs.size(); }
        const Site* sites_of(size_t i) const { return sites.data() + start[i]; }
        size_t weight_of(size_t i) const { return static_cast<size_t>(start[i + 1] - start[i]); }

        void push_back(const Site* s, size_t n, const Coeff& c) {
                sites.insert(sites.end(), s, s + n);
                start.push_back(sites.size());
                coeffs.push_back(c);
        }
};

/// The operator of an emitted term: up to Inline sites are kept inline, longer
/// ones in a pool (offset `spill`). 16 + 4 * Inline bytes.
template <size_t Inline>
struct SiteKey {
        static constexpr size_t kInlineSites = Inline;
        uint64_t hash;
        uint32_t length;
        uint32_t spill;
        Site inline_sites[kInlineSites];

        template <typename Pool>
        const Site* sites(const Pool& pool) const {
                return length <= kInlineSites ? inline_sites : pool.data() + spill;
        }
};

/// A key for the sites s (n of them) with hash h; long keys are appended to pool.
template <typename Key, typename Pool>
inline Key make_key(uint64_t h, const Site* s, size_t n, Pool& pool) {
        Key k;
        k.hash = h;
        k.length = static_cast<uint32_t>(n);
        k.spill = 0;
        if (n <= Key::kInlineSites) {
                std::copy_n(s, n, k.inline_sites);
                std::fill(k.inline_sites + n, k.inline_sites + Key::kInlineSites, Site{0});
        } else {
                if (pool.size() > UINT32_MAX - n) throw std::length_error("sparse sum: site pool beyond 2^32");
                k.spill = static_cast<uint32_t>(pool.size());
                pool.insert(pool.end(), s, s + n);
        }
        return k;
}

/// Adds up terms with equal operators, in the order they are emitted, and drops
/// the ones whose coefficient sums to zero.
///
/// run() calls gen(row, emit) for the rows in blocks. Within a block the rows are
/// formed on all threads, each thread appending its terms to its own buffer. The
/// terms are then routed by the top bits of the hash of their operator into
/// buckets, each with its own segment of one large open-addressing table: every
/// thread sorts its terms by bucket (a stable counting sort), and the buckets are
/// folded in parallel, each taking the terms of the threads in thread order. The
/// number of buckets grows with the expected number of terms so that a segment
/// stays within the cache while it is filled.
///
/// A static schedule gives every thread an ascending range of rows, so every bucket
/// sees its terms in row order and the sums, and the result, do not depend on the
/// number of threads. The result lists the buckets in turn, each in the order of
/// first occurrence. Only one block of emitted terms is held at a time, and all
/// storage is a few large arrays, reused from block to block.
template <typename Coeff, size_t Inline = 4>
class SparseAccumulator {
        using SiteKey = pauliengine_detail::SiteKey<Inline>;
    public:
        static constexpr bool kParallelCoeff = std::is_same_v<Coeff, std::complex<double>>;

        /// Emitted terms per block, summed over all threads.
        static constexpr size_t kBlockTerms = size_t{1} << 19;

        /// Emits the terms of one row; `scratch` is per thread.
        class Emitter {
            public:
                void operator()(const Site* s, size_t n, const Coeff& c) {
                        keys_.push_back(make_key<SiteKey>(hash_sites(s, n), s, n, pool_));
                        coeffs_.push_back(c);
                }
                std::vector<Site>& scratch() { return scratch_; }

            private:
                friend class SparseAccumulator;
                std::vector<SiteKey, DefaultInitAllocator<SiteKey>> keys_;
                std::vector<Coeff, DefaultInitAllocator<Coeff>> coeffs_;
                std::vector<Site> pool_;
                std::vector<Site> scratch_;
                // the terms again, sorted by bucket (stable), and where each bucket starts
                std::vector<SiteKey, DefaultInitAllocator<SiteKey>> sorted_keys_;
                std::vector<Coeff, DefaultInitAllocator<Coeff>> sorted_coeffs_;
                std::vector<size_t> bucket_start_;

                void clear() {
                        keys_.clear();
                        coeffs_.clear();
                        pool_.clear();
                }

                size_t bytes() const {
                        return (keys_.capacity() + sorted_keys_.capacity()) * sizeof(SiteKey)
                             + (coeffs_.capacity() + sorted_coeffs_.capacity()) * sizeof(Coeff)
                             + (pool_.capacity() + scratch_.capacity()) * sizeof(Site);
                }
        };

        /// Calls gen(row, emitter) for row in [0, rows); terms_per_row estimates
        /// the emitted terms of a row (for the block size and the number of buckets),
        /// work_per_row the work (for starting threads).
        template <typename Gen>
        static SparseTerms<Coeff> run(size_t rows, size_t terms_per_row, size_t work_per_row, Gen&& gen) {
                const double expected = static_cast<double>(rows) * static_cast<double>(std::max<size_t>(terms_per_row, 1));
                unsigned bits = 0;
                while (bits < 12 && expected > static_cast<double>(size_t{4096} << bits)) ++bits;
                const bool parallel = kParallelCoeff && use_threads(rows * std::max<size_t>(work_per_row, 1));
                int threads = 1;
#ifdef PAULIENGINE_HAS_OPENMP
                if (parallel) threads = omp_get_max_threads();
#endif
                // The buffers of the previous call on this thread are reused: fresh
                // memory costs a page fault per 4 KB, often more than the work itself.
                Workspace& ws = workspace();
                SparseAccumulator acc(bits);
                acc.slots_.swap(ws.slots);
                acc.keys_.swap(ws.keys);
                acc.coeffs_.swap(ws.coeffs);
                std::vector<Emitter>& emitters = ws.emitters;
                if (emitters.size() < static_cast<size_t>(threads)) emitters.resize(static_cast<size_t>(threads));
                emitters.resize(static_cast<size_t>(threads));
                for (Emitter& em : emitters) em.clear();
                // room for the expected terms up front: growing copies every bucket
                const double per_bucket = expected / static_cast<double>(acc.buckets_);
                acc.ensure_capacity(static_cast<size_t>(std::min(per_bucket * 1.125 + 64.0, 1048576.0)), parallel);
                const size_t block = std::max<size_t>(1, kBlockTerms / std::max<size_t>(terms_per_row, 1));
                for (size_t first = 0; first < rows; first += block) {
                        const size_t last = std::min(rows, first + block);
                        acc.form_block(first, last, parallel, emitters, gen);
                        acc.fold(parallel, emitters);
                }
                SparseTerms<Coeff> out = acc.collect(parallel);
                acc.slots_.swap(ws.slots);
                acc.keys_.swap(ws.keys);
                acc.coeffs_.swap(ws.coeffs);
                ws.release_if_large();
                return out;
        }

    private:
        static constexpr uint32_t kEmpty = UINT32_MAX;

        /// Buffers kept between calls on the same thread, up to kKeepBytes.
        struct Workspace {
                static constexpr size_t kKeepBytes = size_t{1} << 28;
                std::vector<uint32_t, DefaultInitAllocator<uint32_t>> slots;
                std::vector<SiteKey, DefaultInitAllocator<SiteKey>> keys;
                std::vector<Coeff, DefaultInitAllocator<Coeff>> coeffs;
                std::vector<Emitter> emitters;

                size_t bytes() const {
                        size_t b = slots.capacity() * sizeof(uint32_t) + keys.capacity() * sizeof(SiteKey)
                                 + coeffs.capacity() * sizeof(Coeff);
                        for (const Emitter& em : emitters) b += em.bytes();
                        return b;
                }

                void release_if_large() {
                        if (bytes() <= kKeepBytes) return;
                        *this = Workspace();
                }
        };

        static Workspace& workspace() {
                thread_local Workspace ws;
                return ws;
        }

        const unsigned bits_;
        const size_t buckets_;
        size_t slots_per_bucket_ = 0;   // power of two
        size_t capacity_ = 0;           // entries per bucket: half the slots
        std::vector<uint32_t, DefaultInitAllocator<uint32_t>> slots_;      // buckets_ * slots_per_bucket_
        std::vector<SiteKey, DefaultInitAllocator<SiteKey>> keys_;        // buckets_ * capacity_
        std::vector<Coeff, DefaultInitAllocator<Coeff>> coeffs_;          // buckets_ * capacity_
        std::vector<uint32_t> used_;                                       // entries per bucket
        std::vector<std::vector<Site>> pools_;                             // long keys, per bucket

        explicit SparseAccumulator(unsigned bits)
                : bits_(bits), buckets_(size_t{1} << bits), used_(size_t{1} << bits, 0), pools_(size_t{1} << bits) {}

        size_t bucket_of(uint64_t h) const { return bits_ == 0 ? 0 : static_cast<size_t>(h >> (64 - bits_)); }

        template <typename Gen>
        void form_block(size_t first, size_t last, bool parallel, std::vector<Emitter>& emitters, Gen& gen) {
                const std::ptrdiff_t lo = static_cast<std::ptrdiff_t>(first);
                const std::ptrdiff_t hi = static_cast<std::ptrdiff_t>(last);
#ifdef PAULIENGINE_HAS_OPENMP
                if (parallel) {
                        #pragma omp parallel
                        {
                                Emitter& mine = emitters[static_cast<size_t>(omp_get_thread_num())];
                                #pragma omp for schedule(static)
                                for (std::ptrdiff_t r = lo; r < hi; ++r) gen(static_cast<size_t>(r), mine);
                        }
                        return;
                }
#else
                (void)parallel;
#endif
                for (std::ptrdiff_t r = lo; r < hi; ++r) gen(static_cast<size_t>(r), emitters[0]);
        }

        /// Room for `needed` entries in every bucket: doubles the segments and moves
        /// the entries over, bucket by bucket.
        void ensure_capacity(size_t needed, bool parallel) {
                if (needed <= capacity_) return;
                size_t cap = std::max<size_t>(capacity_, 64);
                while (cap < needed) cap *= 2;
                if (capacity_ == 0) cap = needed;   // first allocation: exactly as asked
                size_t spb = 64;
                while (spb < 2 * cap) spb *= 2;
                if (capacity_ == 0) {
                        // nothing to move: size the (reused) arrays in place
                        slots_.resize(buckets_ * spb);
                        keys_.resize(buckets_ * cap);
                        coeffs_.resize(buckets_ * cap);
                        for_each_index(buckets_, parallel, [&](size_t b) {
                                std::fill(slots_.data() + b * spb, slots_.data() + (b + 1) * spb, kEmpty);
                        });
                        capacity_ = cap;
                        slots_per_bucket_ = spb;
                        return;
                }
                decltype(slots_) slots(buckets_ * spb);
                decltype(keys_) keys(buckets_ * cap);
                decltype(coeffs_) coeffs(buckets_ * cap);
                for_each_index(buckets_, parallel, [&](size_t b) {
                        uint32_t* seg = slots.data() + b * spb;
                        std::fill(seg, seg + spb, kEmpty);
                        for (size_t e = 0; e < used_[b]; ++e) {
                                keys[b * cap + e] = keys_[b * capacity_ + e];
                                coeffs[b * cap + e] = std::move(coeffs_[b * capacity_ + e]);
                                size_t pos = static_cast<size_t>(keys[b * cap + e].hash) & (spb - 1);
                                while (seg[pos] != kEmpty) pos = (pos + 1) & (spb - 1);
                                seg[pos] = static_cast<uint32_t>(e);
                        }
                });
                slots_.swap(slots);
                keys_.swap(keys);
                coeffs_.swap(coeffs);
                capacity_ = cap;
                slots_per_bucket_ = spb;
        }

        void add(size_t b, const SiteKey& key, const std::vector<Site>& key_pool, const Coeff& c) {
                uint32_t* seg = slots_.data() + b * slots_per_bucket_;
                SiteKey* keys = keys_.data() + b * capacity_;
                const size_t mask = slots_per_bucket_ - 1;
                const std::vector<Site>& pool = pools_[b];
                const Site* s = key.sites(key_pool);
                size_t pos = static_cast<size_t>(key.hash) & mask;
                while (seg[pos] != kEmpty) {
                        const uint32_t e = seg[pos];
                        const SiteKey& other = keys[e];
                        if (other.hash == key.hash && other.length == key.length
                            && std::equal(s, s + key.length, other.sites(pool))) {
                                Coeff& target = coeffs_[b * capacity_ + e];
                                target = target + c;
                                return;
                        }
                        pos = (pos + 1) & mask;
                }
                const uint32_t e = used_[b]++;
                seg[pos] = e;
                keys[e] = make_key<SiteKey>(key.hash, s, key.length, pools_[b]);
                coeffs_[b * capacity_ + e] = c;
        }

        void fold(bool parallel, std::vector<Emitter>& emitters) {
                // every thread moves its terms into bucket order (a stable counting
                // sort), so that each bucket later reads them sequentially
                for_each_index(emitters.size(), parallel && emitters.size() > 1, [&](size_t t) {
                        Emitter& em = emitters[t];
                        const size_t n = em.keys_.size();
                        em.bucket_start_.assign(buckets_ + 1, 0);
                        for (const SiteKey& k : em.keys_) ++em.bucket_start_[bucket_of(k.hash) + 1];
                        for (size_t b = 0; b < buckets_; ++b) em.bucket_start_[b + 1] += em.bucket_start_[b];
                        em.sorted_keys_.resize(n);
                        em.sorted_coeffs_.resize(n);
                        std::vector<size_t> next(em.bucket_start_.begin(), em.bucket_start_.end() - 1);
                        for (size_t r = 0; r < n; ++r) {
                                const size_t to = next[bucket_of(em.keys_[r].hash)]++;
                                em.sorted_keys_[to] = em.keys_[r];
                                em.sorted_coeffs_[to] = std::move(em.coeffs_[r]);
                        }
                });
                // room for all terms of this block in every bucket, at worst all new
                size_t needed = 0;
                for (size_t b = 0; b < buckets_; ++b) {
                        size_t incoming = used_[b];
                        for (const Emitter& em : emitters) incoming += em.bucket_start_[b + 1] - em.bucket_start_[b];
                        needed = std::max(needed, incoming);
                }
                if (needed > UINT32_MAX / 2) throw std::length_error("sparse sum: too many terms in one bucket");
                ensure_capacity(needed, parallel);
                for_each_index(buckets_, parallel, [&](size_t b) {
                        for (const Emitter& em : emitters) {
                                for (size_t p = em.bucket_start_[b]; p < em.bucket_start_[b + 1]; ++p) {
                                        add(b, em.sorted_keys_[p], em.pool_, em.sorted_coeffs_[p]);
                                }
                        }
                });
                for (Emitter& em : emitters) em.clear();
        }

        SparseTerms<Coeff> collect(bool parallel) {
                // kept terms per bucket, then their offsets in the result
                std::vector<size_t> count(buckets_ + 1, 0), sites(buckets_ + 1, 0);
                for_each_index(buckets_, parallel, [&](size_t b) {
                        size_t kept = 0, kept_sites = 0;  // local: neighbouring counters share cache lines
                        for (size_t e = 0; e < used_[b]; ++e) {
                                if (!coeff_is_zero(coeffs_[b * capacity_ + e])) {
                                        ++kept;
                                        kept_sites += keys_[b * capacity_ + e].length;
                                }
                        }
                        count[b + 1] = kept;
                        sites[b + 1] = kept_sites;
                });
                for (size_t b = 0; b < buckets_; ++b) {
                        count[b + 1] += count[b];
                        sites[b + 1] += sites[b];
                }
                SparseTerms<Coeff> out;
                out.start.resize(count[buckets_] + 1);
                out.start[0] = 0;
                out.sites.resize(sites[buckets_]);
                out.coeffs.resize(count[buckets_]);
                for_each_index(buckets_, parallel, [&](size_t b) {
                        size_t term = count[b];
                        size_t site = sites[b];
                        for (size_t e = 0; e < used_[b]; ++e) {
                                Coeff& c = coeffs_[b * capacity_ + e];
                                if (coeff_is_zero(c)) continue;
                                const SiteKey& k = keys_[b * capacity_ + e];
                                std::copy_n(k.sites(pools_[b]), k.length, out.sites.data() + site);
                                site += k.length;
                                out.start[term + 1] = site;
                                out.coeffs[term] = std::move(c);
                                ++term;
                        }
                });
                return out;
        }
};

/// SparseAccumulator::run with keys that hold max_weight sites inline (4, 8 or 16;
/// heavier terms spill into a pool).
template <typename Coeff, typename Gen>
SparseTerms<Coeff> accumulate(size_t max_weight, size_t rows, size_t terms_per_row, size_t work_per_row, Gen&& gen) {
        if (max_weight <= 4) return SparseAccumulator<Coeff, 4>::run(rows, terms_per_row, work_per_row, gen);
        if (max_weight <= 8) return SparseAccumulator<Coeff, 8>::run(rows, terms_per_row, work_per_row, gen);
        return SparseAccumulator<Coeff, 16>::run(rows, terms_per_row, work_per_row, gen);
}

}  // namespace pauliengine_detail
