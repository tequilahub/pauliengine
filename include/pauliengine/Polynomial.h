#pragma once

/// @file
/// Symbolic coefficients of PauliEngine: sparse multivariate polynomials in named
/// symbols with std::complex<double> coefficients.
///
/// This replaces SymEngine. PauliEngine only ever forms sums, products and integer
/// powers of its coefficients, takes derivatives, substitutes numbers and
/// conjugates, so polynomials are all that is needed. The representation is
/// canonical (terms sorted by monomial, equal monomials merged, zero terms
/// dropped), hence a zero test or a comparison is a plain look at the data, with no
/// expand step.
///
/// An Expression is a value type that owns one flat buffer and shares nothing with
/// other expressions, so expressions may be created, copied and combined on any
/// number of threads at once. Only the symbol table (name <-> id) is shared; it is
/// guarded by a spinlock and is touched when a symbol is created by name, printed,
/// substituted or conjugated, not during arithmetic.
///
/// Conjugation does not assume anything about a symbol: conjugate(a) is the partner
/// symbol "conjugate(a)", a flag that a was conjugated. Substituting a value v (or a
/// polynomial p) for "a" puts conj(v) (conjugate(p)) in for "conjugate(a)", so
/// dagger() is right whatever is substituted later. A symbol has one of three kinds:
///
/// - unknown (the default): conjugated as above; diff("a") differentiates
///   conjugate(a) like a (d/da conjugate(a) = 1), right for real parameters.
/// - real (set_symbol_real("a")): conjugate(a) == a, so a Hermitian sum with real
///   symbols stays equal to its dagger before anything is substituted.
/// - complex (set_symbol_complex("a")): conjugated as above; diff treats a and
///   conjugate(a) as independent variables (Wirtinger calculus).

#include <algorithm>
#include <bit>
#include <cctype>
#include <charconv>
#include <cmath>
#include <complex>
#include <cstdint>
#include <initializer_list>
#include <atomic>
#include <thread>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace pauliengine::symbolic {

// ---------------------------------------------------------------------------
// Symbol table
// ---------------------------------------------------------------------------

/// Process-wide table of symbol names. Ids are dense and never reused.
class SymbolTable {
    public:
        static SymbolTable& instance() {
                static SymbolTable table;
                return table;
        }

        /// Id of the symbol `name`, created (real) on first use.
        uint32_t id(const std::string& name) {
                Lock lock(flag_);
                return id_locked(name);
        }

        /// Id of `name` if the symbol exists.
        bool find(const std::string& name, uint32_t& out) {
                Lock lock(flag_);
                const auto it = ids_.find(name);
                if (it == ids_.end()) return false;
                out = it->second;
                return true;
        }

        std::string name(uint32_t id) {
                Lock lock(flag_);
                return names_.at(id);
        }

        /// Names of all symbols, indexed by id (one lock for many lookups).
        std::vector<std::string> names() {
                Lock lock(flag_);
                return names_;
        }

        enum class Kind : uint8_t { Unknown, Real, Complex };

        /// The conjugation map, indexed by id: a symbol's partner, or the symbol
        /// itself if it is real. Partners are created here for every symbol that
        /// has none yet, so the map covers all symbols that exist.
        std::vector<uint32_t> conjugates() {
                Lock lock(flag_);
                for (uint32_t a = 0, n = static_cast<uint32_t>(names_.size()); a < n; ++a) {
                        if (partner_[a] == kNone && kind_[a] != Kind::Real) make_partner_locked(a);
                }
                return conj_;
        }

        /// Partners by id (kNone if a symbol has none yet), for subs() and diff().
        std::vector<uint32_t> partners() {
                Lock lock(flag_);
                return partner_;
        }

        /// The kind of the symbol `id`; a partner has the kind of its symbol.
        std::vector<Kind> kinds() {
                Lock lock(flag_);
                return kind_;
        }

        void set_kind(const std::string& name, Kind kind) {
                Lock lock(flag_);
                const uint32_t a = id_locked(name);
                const uint32_t b = partner_[a];
                kind_[a] = kind;
                if (b != kNone) kind_[b] = kind;
                if (kind == Kind::Real) {
                        conj_[a] = a;                   // conjugate(a) == a from now on
                        if (b != kNone) conj_[b] = a;   // and an existing conjugate(a) goes back to a
                } else {
                        if (b == kNone) {
                                make_partner_locked(a);
                        } else {
                                conj_[a] = b;
                                conj_[b] = a;
                        }
                }
        }

        Kind kind(const std::string& name) {
                Lock lock(flag_);
                const auto it = ids_.find(name);
                return it == ids_.end() ? Kind::Unknown : kind_[it->second];
        }

        static constexpr uint32_t kNone = UINT32_MAX;

    private:
        SymbolTable() = default;

        /// Spinlock instead of std::mutex: the table is only touched outside of
        /// arithmetic, and std::mutex crashes (access violation in lock()) when
        /// an older msvcp140.dll than the one we were built against is already
        /// loaded into the process, which happens next to other extensions.
        struct Lock {
                std::atomic_flag& f;
                explicit Lock(std::atomic_flag& flag) : f(flag) {
                        while (f.test_and_set(std::memory_order_acquire)) std::this_thread::yield();
                }
                ~Lock() { f.clear(std::memory_order_release); }
        };

        uint32_t id_locked(const std::string& name) {
                const auto it = ids_.find(name);
                if (it != ids_.end()) return it->second;
                const auto id = static_cast<uint32_t>(names_.size());
                names_.push_back(name);
                conj_.push_back(id);
                partner_.push_back(kNone);
                kind_.push_back(Kind::Unknown);
                ids_.emplace(name, id);
                return id;
        }

        /// Creates the partner "conjugate(name)" of symbol a and links the two.
        void make_partner_locked(uint32_t a) {
                const uint32_t b = id_locked("conjugate(" + names_[a] + ")");
                partner_[a] = b;
                partner_[b] = a;
                kind_[b] = kind_[a];
                conj_[a] = b;
                conj_[b] = a;
        }

        std::atomic_flag flag_ = ATOMIC_FLAG_INIT;
        std::vector<std::string> names_;
        std::vector<uint32_t> conj_;
        std::vector<uint32_t> partner_;
        std::vector<Kind> kind_;
        std::unordered_map<std::string, uint32_t> ids_;
};

/// Declare `name` complex (Wirtinger derivatives), or with false real.
inline void set_symbol_complex(const std::string& name, bool is_complex = true) {
        SymbolTable::instance().set_kind(name, is_complex ? SymbolTable::Kind::Complex : SymbolTable::Kind::Real);
}

/// Declare `name` real (conjugate(name) == name), or with false back to unknown.
inline void set_symbol_real(const std::string& name, bool is_real = true) {
        SymbolTable::instance().set_kind(name, is_real ? SymbolTable::Kind::Real : SymbolTable::Kind::Unknown);
}

inline bool is_symbol_complex(const std::string& name) {
        return SymbolTable::instance().kind(name) == SymbolTable::Kind::Complex;
}

inline bool is_symbol_real(const std::string& name) {
        return SymbolTable::instance().kind(name) == SymbolTable::Kind::Real;
}

// ---------------------------------------------------------------------------
// WordBuf: the storage of an Expression
// ---------------------------------------------------------------------------

/// A vector of 64-bit words with room for kInline words inside the object, so a
/// single term with up to three factors (the common case of a product of Pauli
/// strings) needs no heap allocation. Only the operations Expression uses.
class WordBuf {
    public:
        static constexpr uint32_t kInline = 6;

        WordBuf() noexcept : ptr_(inline_) {}
        WordBuf(std::initializer_list<uint64_t> init) : WordBuf() {
                reserve(init.size());
                std::copy(init.begin(), init.end(), ptr_);
                size_ = static_cast<uint32_t>(init.size());
        }
        WordBuf(const WordBuf& o) : WordBuf() {
                reserve(o.size_);
                std::copy(o.ptr_, o.ptr_ + o.size_, ptr_);
                size_ = o.size_;
        }
        WordBuf(WordBuf&& o) noexcept : WordBuf() { take(std::move(o)); }
        WordBuf& operator=(const WordBuf& o) {
                if (this != &o) {
                        size_ = 0;
                        reserve(o.size_);
                        std::copy(o.ptr_, o.ptr_ + o.size_, ptr_);
                        size_ = o.size_;
                }
                return *this;
        }
        WordBuf& operator=(WordBuf&& o) noexcept {
                if (this != &o) {
                        release();
                        take(std::move(o));
                }
                return *this;
        }
        WordBuf& operator=(std::initializer_list<uint64_t> init) { return *this = WordBuf(init); }
        ~WordBuf() { release(); }

        uint64_t* data() noexcept { return ptr_; }
        const uint64_t* data() const noexcept { return ptr_; }
        size_t size() const noexcept { return size_; }
        bool empty() const noexcept { return size_ == 0; }
        uint64_t* begin() noexcept { return ptr_; }
        uint64_t* end() noexcept { return ptr_ + size_; }
        const uint64_t* begin() const noexcept { return ptr_; }
        const uint64_t* end() const noexcept { return ptr_ + size_; }
        uint64_t& operator[](size_t i) noexcept { return ptr_[i]; }
        uint64_t operator[](size_t i) const noexcept { return ptr_[i]; }
        void clear() noexcept { size_ = 0; }

        void reserve(size_t n) {
                if (n <= cap_) return;
                const size_t cap = std::max<size_t>(n, 2 * static_cast<size_t>(cap_));
                uint64_t* fresh = new uint64_t[cap];
                std::copy(ptr_, ptr_ + size_, fresh);
                if (ptr_ != inline_) delete[] ptr_;
                ptr_ = fresh;
                cap_ = static_cast<uint32_t>(cap);
        }
        void resize(size_t n) {
                reserve(n);
                if (n > size_) std::fill(ptr_ + size_, ptr_ + n, uint64_t{0});
                size_ = static_cast<uint32_t>(n);
        }
        void push_back(uint64_t w) {
                if (size_ == cap_) reserve(size_ + 1);
                ptr_[size_++] = w;
        }
        /// Appends [first, last); `pos` must be end().
        void insert(const uint64_t* /*pos*/, const uint64_t* first, const uint64_t* last) {
                const size_t n = static_cast<size_t>(last - first);
                reserve(size_ + n);
                std::copy(first, last, ptr_ + size_);
                size_ += static_cast<uint32_t>(n);
        }

        friend bool operator==(const WordBuf& a, const WordBuf& b) noexcept {
                return a.size_ == b.size_ && std::equal(a.ptr_, a.ptr_ + a.size_, b.ptr_);
        }

    private:
        void release() noexcept {
                if (ptr_ != inline_) delete[] ptr_;
                ptr_ = inline_;
                cap_ = kInline;
                size_ = 0;
        }
        void take(WordBuf&& o) noexcept {
                if (o.ptr_ == o.inline_) {
                        std::copy(o.inline_, o.inline_ + o.size_, inline_);
                        ptr_ = inline_;
                        cap_ = kInline;
                } else {
                        ptr_ = o.ptr_;
                        cap_ = o.cap_;
                        o.ptr_ = o.inline_;
                        o.cap_ = kInline;
                }
                size_ = o.size_;
                o.size_ = 0;
        }

        uint64_t* ptr_;
        uint32_t size_ = 0;
        uint32_t cap_ = kInline;
        uint64_t inline_[kInline];
};

// ---------------------------------------------------------------------------
// Expression
// ---------------------------------------------------------------------------

/// A polynomial sum_k c_k * m_k with complex c_k and monomials m_k.
///
/// Storage: one buffer of 64-bit words. Term k is [nf, re, im, f_1 .. f_nf], with
/// re and im the bits of the coefficient and f_i = (symbol id << 32) | exponent,
/// sorted by id. Terms are sorted by monomial (lexicographic over the factor
/// words, a prefix first), so the constant term comes first.
class Expression {
    public:
        using complex = std::complex<double>;

        Expression() = default;
        Expression(int v) : Expression(complex(static_cast<double>(v), 0.0)) {}
        Expression(double v) : Expression(complex(v, 0.0)) {}
        Expression(complex v) {
                if (v != complex(0.0, 0.0)) {
                        buf_ = {0, bits(v.real()), bits(v.imag())};
                        n_ = 1;
                }
        }
        /// Parse "2*a*b + 3", "a**2 - b/4", "(1 + 2*I)*c", "2a", "conjugate(a)".
        explicit Expression(const std::string& text);
        explicit Expression(const char* text) : Expression(std::string(text)) {}

        static Expression symbol(const std::string& name) {
                return symbol(SymbolTable::instance().id(name));
        }
        static Expression symbol(uint32_t id) {
                Expression e;
                e.buf_ = {1, bits(1.0), bits(0.0), factor(id, 1)};
                e.n_ = 1;
                return e;
        }

        /// Number of terms; 0 for the zero polynomial.
        size_t n_terms() const { return n_; }
        bool is_zero() const { return n_ == 0; }
        /// True if the expression has no symbols.
        bool is_number() const { return n_ == 0 || (n_ == 1 && buf_[0] == 0); }

        /// The value of a numeric expression; throws if symbols are left.
        complex to_complex() const {
                if (n_ == 0) return {0.0, 0.0};
                if (!is_number()) {
                        throw std::runtime_error("Expression cannot be converted to std::complex<double>: "
                                                 "it still contains symbols (" + str() + ").");
                }
                return coeff_at(buf_.data());
        }

        bool operator==(const Expression& o) const { return n_ == o.n_ && buf_ == o.buf_; }
        bool operator!=(const Expression& o) const { return !(*this == o); }

        // -- arithmetic -------------------------------------------------------

        friend Expression operator+(const Expression& a, const Expression& b) { return combine(a, b, false); }
        friend Expression operator-(const Expression& a, const Expression& b) { return combine(a, b, true); }
        friend Expression operator*(const Expression& a, const Expression& b) { return multiply(a, b); }
        Expression operator-() const { return scaled(complex(-1.0, 0.0)); }

        friend Expression operator*(const Expression& a, complex s) { return a.scaled(s); }
        friend Expression operator*(complex s, const Expression& a) { return a.scaled(s); }
        friend Expression operator*(const Expression& a, double s) { return a.scaled(complex(s, 0.0)); }
        friend Expression operator*(double s, const Expression& a) { return a.scaled(complex(s, 0.0)); }
        friend Expression operator*(const Expression& a, int s) { return a.scaled(complex(s, 0.0)); }
        friend Expression operator*(int s, const Expression& a) { return a.scaled(complex(s, 0.0)); }
        friend Expression operator+(const Expression& a, complex s) { return a + Expression(s); }
        friend Expression operator+(complex s, const Expression& a) { return Expression(s) + a; }
        friend Expression operator-(const Expression& a, complex s) { return a - Expression(s); }
        friend Expression operator-(complex s, const Expression& a) { return Expression(s) - a; }
        friend Expression operator+(const Expression& a, double s) { return a + Expression(s); }
        friend Expression operator+(double s, const Expression& a) { return Expression(s) + a; }
        friend Expression operator-(const Expression& a, double s) { return a - Expression(s); }
        friend Expression operator-(double s, const Expression& a) { return Expression(s) - a; }

        Expression& operator+=(const Expression& o) { return *this = *this + o; }
        Expression& operator-=(const Expression& o) { return *this = *this - o; }
        Expression& operator*=(const Expression& o) { return *this = *this * o; }
        Expression& operator*=(complex s) { scale_in_place(s); return *this; }
        Expression& operator*=(double s) { scale_in_place(complex(s, 0.0)); return *this; }
        Expression& operator*=(int s) { scale_in_place(complex(s, 0.0)); return *this; }

        /// get(0) + ... + get(n - 1) in one merge; get(k) returns a const Expression&.
        /// Adds equal monomials in index order, like repeated + would.
        template<typename Get>
        static Expression sum_of(size_t n, Get&& get) {
                if (n == 0) return {};
                if (n == 1) return get(0);
                if (n == 2) return get(0) + get(1);
                Builder b;
                for (size_t k = 0; k < n; ++k) {
                        const Expression& e = get(k);
                        for (const uint64_t* t = e.buf_.data(), *end = t + e.buf_.size(); t < end; t = next(t)) {
                                b.begin(coeff_at(t));
                                for (uint64_t i = 0; i < t[0]; ++i) b.push(t[3 + i]);
                        }
                }
                return b.finish();
        }

        /// this**n for n >= 0, by repeated squaring.
        Expression pow(unsigned n) const {
                Expression result(1);
                Expression base = *this;
                while (n) {
                        if (n & 1) result = result * base;
                        n >>= 1;
                        if (n) base = base * base;
                }
                return result;
        }

        // -- calculus and evaluation -------------------------------------------

        /// The symbols a derivative by one name acts on: the symbol, and its partner
        /// conjugate(a) unless the symbol is complex (then the two are independent).
        struct DiffTarget {
                uint32_t id = SymbolTable::kNone;
                uint32_t partner = SymbolTable::kNone;
        };

        /// Resolves `name` for diff(); false if the symbol does not exist.
        static bool prepare_diff(const std::string& name, DiffTarget& target) {
                SymbolTable& table = SymbolTable::instance();
                uint32_t id;
                if (!table.find(name, id)) return false;
                target.id = id;
                const std::vector<uint32_t> partners = table.partners();
                const std::vector<SymbolTable::Kind> kinds = table.kinds();
                if (partners[id] != SymbolTable::kNone && kinds[id] != SymbolTable::Kind::Complex) {
                        target.partner = partners[id];
                }
                return true;
        }

        /// Derivative with respect to the symbol `name` (0 if it does not occur).
        /// conjugate(name) counts as the same variable unless the symbol is complex.
        Expression diff(const std::string& name) const {
                DiffTarget target;
                if (!prepare_diff(name, target)) return {};
                return diff(target);
        }

        Expression diff(const DiffTarget& target) const {
                if (target.partner == SymbolTable::kNone) return diff(target.id);
                return diff(target.id) + diff(target.partner);
        }

        Expression diff(uint32_t id) const {
                Builder out;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        const uint64_t* f = factors(t);
                        const uint64_t nf = t[0];
                        for (uint64_t i = 0; i < nf; ++i) {
                                if (factor_id(f[i]) != id) continue;
                                const uint32_t e = factor_exp(f[i]);
                                out.begin(coeff_at(t) * static_cast<double>(e));
                                for (uint64_t j = 0; j < nf; ++j) {
                                        if (j != i) out.push(f[j]);
                                        else if (e > 1) out.push(factor(id, e - 1));
                                }
                                break;
                        }
                }
                return out.finish();
        }

        /// Values for symbols, resolved to ids once for many subs() calls.
        struct Substitution {
                std::vector<std::pair<uint32_t, complex>> by_id;  // sorted by id
                bool empty() const { return by_id.empty(); }
                const complex* find(uint32_t id) const {
                        const auto it = std::lower_bound(by_id.begin(), by_id.end(), id,
                                [](const std::pair<uint32_t, complex>& p, uint32_t key) { return p.first < key; });
                        return (it != by_id.end() && it->first == id) ? &it->second : nullptr;
                }
        };

        /// Resolves {name: value}. Unknown names are ignored. For a complex symbol a,
        /// conjugate(a) gets conj(value) unless it is given too.
        static Substitution prepare_subs(const std::unordered_map<std::string, complex>& values) {
                Substitution s;
                if (values.empty()) return s;
                SymbolTable& table = SymbolTable::instance();
                std::unordered_map<uint32_t, complex> by_id;
                for (const auto& [name, v] : values) {
                        uint32_t id;
                        if (table.find(name, id)) by_id[id] = v;
                }
                const std::vector<uint32_t> partners = table.partners();
                for (const auto& [name, v] : values) {
                        uint32_t id;
                        if (table.find(name, id) && partners[id] != SymbolTable::kNone) {
                                by_id.try_emplace(partners[id], std::conj(v));
                        }
                }
                s.by_id.assign(by_id.begin(), by_id.end());
                std::sort(s.by_id.begin(), s.by_id.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });
                return s;
        }

        /// Substitute numbers for symbols, given by name (see prepare_subs).
        Expression subs(const std::unordered_map<std::string, complex>& values) const {
                if (n_ == 0 || values.empty()) return *this;
                return subs(prepare_subs(values));
        }

        Expression subs(const Substitution& s) const {
                if (n_ == 0 || s.empty()) return *this;
                bool touched = false;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end && !touched; t = next(t)) {
                        for (uint64_t i = 0; i < t[0] && !touched; ++i) touched = s.find(factor_id(t[3 + i])) != nullptr;
                }
                if (!touched) return *this;
                Builder out;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        complex c = coeff_at(t);
                        out.begin({});
                        const uint64_t* f = factors(t);
                        for (uint64_t i = 0; i < t[0]; ++i) {
                                const complex* v = s.find(factor_id(f[i]));
                                if (v == nullptr) {
                                        out.push(f[i]);
                                } else {
                                        c *= ipow(*v, factor_exp(f[i]));
                                }
                        }
                        out.set_coeff(c);
                }
                return out.finish();
        }

        /// Polynomials for symbols, resolved to ids once for many subs() calls.
        struct PolySubstitution {
                std::vector<std::pair<uint32_t, Expression>> by_id;  // sorted by id
                bool empty() const { return by_id.empty(); }
                const Expression* find(uint32_t id) const {
                        const auto it = std::lower_bound(by_id.begin(), by_id.end(), id,
                                [](const std::pair<uint32_t, Expression>& p, uint32_t key) { return p.first < key; });
                        return (it != by_id.end() && it->first == id) ? &it->second : nullptr;
                }
        };

        /// Resolves {name: polynomial}. Unknown names are ignored. For a complex
        /// symbol a, conjugate(a) gets the conjugated polynomial unless it is given too.
        static PolySubstitution prepare_subs(const std::unordered_map<std::string, Expression>& values) {
                PolySubstitution s;
                if (values.empty()) return s;
                SymbolTable& table = SymbolTable::instance();
                std::unordered_map<uint32_t, Expression> by_id;
                for (const auto& [name, v] : values) {
                        uint32_t id;
                        if (table.find(name, id)) by_id[id] = v;
                }
                const std::vector<uint32_t> partners = table.partners();
                std::vector<uint32_t> conj_map;
                for (const auto& [name, v] : values) {
                        uint32_t id;
                        if (table.find(name, id) && partners[id] != SymbolTable::kNone) {
                                if (conj_map.empty()) conj_map = conjugation_map();
                                by_id.try_emplace(partners[id], v.conjugate(conj_map));
                        }
                }
                s.by_id.assign(by_id.begin(), by_id.end());
                std::sort(s.by_id.begin(), s.by_id.end(),
                          [](const auto& a, const auto& b) { return a.first < b.first; });
                return s;
        }

        /// Substitute polynomials for symbols, given by name, all at once: {"a": b,
        /// "b": a} swaps a and b. Numbers are constant polynomials, so this also
        /// covers subs with values.
        Expression subs(const std::unordered_map<std::string, Expression>& values) const {
                if (n_ == 0 || values.empty()) return *this;
                return subs(prepare_subs(values));
        }

        Expression subs(const PolySubstitution& s) const {
                if (n_ == 0 || s.empty()) return *this;
                bool touched = false;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end && !touched; t = next(t)) {
                        for (uint64_t i = 0; i < t[0] && !touched; ++i) touched = s.find(factor_id(t[3 + i])) != nullptr;
                }
                if (!touched) return *this;
                // Powers of the substituted polynomials, each formed once.
                std::vector<uint64_t> power_of;   // factor words (symbol, exponent) ...
                std::vector<Expression> powers;   // ... and polynomial**exponent for them
                const auto power = [&](uint64_t f, const Expression& p) -> size_t {
                        for (size_t k = 0; k < power_of.size(); ++k) {
                                if (power_of[k] == f) return k;
                        }
                        power_of.push_back(f);
                        powers.push_back(p.pow(factor_exp(f)));
                        return powers.size() - 1;
                };
                std::vector<Expression> parts;
                parts.reserve(n_);
                std::vector<uint64_t> kept;
                std::vector<size_t> inserted;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        kept.clear();
                        inserted.clear();
                        const uint64_t* f = factors(t);
                        for (uint64_t i = 0; i < t[0]; ++i) {
                                const Expression* p = s.find(factor_id(f[i]));
                                if (p == nullptr) kept.push_back(f[i]);
                                else inserted.push_back(power(f[i], *p));
                        }
                        Expression term;
                        term.append_term_unchecked(coeff_at(t), kept.data(), kept.size());
                        for (size_t k : inserted) term = term * powers[k];
                        parts.push_back(std::move(term));
                }
                return sum_of(parts.size(), [&](size_t k) -> const Expression& { return parts[k]; });
        }

        /// The conjugation map for conjugate(map); fetch it once for many
        /// conjugations. It covers every symbol that exists when it is fetched.
        static std::vector<uint32_t> conjugation_map() {
                return SymbolTable::instance().conjugates();
        }

        /// Complex conjugate: conjugated coefficients, every symbol swapped with its
        /// partner conjugate(a), except real symbols, which stay.
        Expression conjugate() const {
                if (is_number()) return conjugate(std::vector<uint32_t>{});
                return conjugate(conjugation_map());
        }

        /// conjugate() with the map from conjugation_map(); empty means all real.
        Expression conjugate(const std::vector<uint32_t>& conj) const {
                Expression out = *this;
                bool renamed = false;
                for (uint64_t* t = out.buf_.data(), *end = t + out.buf_.size(); t < end; t = next_mut(t)) {
                        t[2] = bits(-std::bit_cast<double>(t[2]) + 0.0);
                        if (t[0] == 0 || conj.empty()) continue;
                        for (uint64_t i = 0; i < t[0]; ++i) {
                                const uint32_t id = factor_id(t[3 + i]);
                                if (id < conj.size() && conj[id] != id) {
                                        t[3 + i] = factor(conj[id], factor_exp(t[3 + i]));
                                        renamed = true;
                                }
                        }
                }
                if (!renamed) return out;
                // Renaming changes the order of factors and terms: rebuild.
                Builder b;
                for (const uint64_t* t = out.buf_.data(), *end = t + out.buf_.size(); t < end; t = next(t)) {
                        b.begin(coeff_at(t));
                        for (uint64_t i = 0; i < t[0]; ++i) b.push(t[3 + i]);
                        b.sort_last_factors();
                }
                return b.finish();
        }

        /// Names of the symbols that occur, sorted by name.
        std::vector<std::string> free_symbols() const {
                std::vector<uint32_t> ids;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        for (uint64_t i = 0; i < t[0]; ++i) ids.push_back(factor_id(t[3 + i]));
                }
                std::sort(ids.begin(), ids.end());
                ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                std::vector<std::string> out;
                if (ids.empty()) return out;
                const std::vector<std::string> names = SymbolTable::instance().names();
                for (uint32_t id : ids) out.push_back(names[id]);
                std::sort(out.begin(), out.end());
                return out;
        }

        /// Readable form that the parser reads back exactly, e.g. "2 + 3*a*b**2 - I*c".
        std::string str() const {
                if (n_ == 0) return "0";
                std::vector<std::string> names;
                std::string out;
                bool first = true;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        std::string mono;
                        if (t[0] > 0 && names.empty()) names = SymbolTable::instance().names();
                        for (uint64_t i = 0; i < t[0]; ++i) {
                                if (i) mono += '*';
                                mono += names[factor_id(t[3 + i])];
                                const uint32_t e = factor_exp(t[3 + i]);
                                if (e != 1) mono += "**" + std::to_string(e);
                        }
                        bool negative = false;
                        const std::string term = format_term(coeff_at(t), mono, negative);
                        if (first) {
                                out += negative ? "-" + term : term;
                        } else {
                                out += negative ? " - " + term : " + " + term;
                        }
                        first = false;
                }
                return out;
        }

        // -- low-level access for kernels that work on whole Hamiltonians ---------

        static uint32_t symbol_of(uint64_t f) { return factor_id(f); }
        static uint32_t exponent_of(uint64_t f) { return factor_exp(f); }
        static uint64_t make_factor(uint32_t id, uint32_t exp) { return factor(id, exp); }

        /// Calls f(coefficient, factors, n_factors) for every term in canonical order;
        /// a factor is a word read with symbol_of() and exponent_of().
        template<typename F>
        void for_each_term(F&& f) const {
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        f(coeff_at(t), t + 3, static_cast<size_t>(t[0]));
                }
        }

        /// Appends a term without any check: the coefficient must be non-zero, the
        /// factors sorted by symbol, and the monomial greater than all present ones.
        void append_term_unchecked(complex c, const uint64_t* f, size_t nf) {
                buf_.reserve(buf_.size() + 3 + nf);
                buf_.push_back(nf);
                buf_.push_back(bits(c.real()));
                buf_.push_back(bits(c.imag()));
                buf_.insert(buf_.end(), f, f + nf);
                ++n_;
        }

        size_t hash() const {
                size_t h = 0xcbf29ce484222325ULL ^ n_;
                for (uint64_t w : buf_) h ^= std::hash<uint64_t>{}(w) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                return h;
        }

    private:
        WordBuf buf_;
        uint32_t n_ = 0;

        // -- word helpers -----------------------------------------------------

        static uint64_t bits(double d) { return std::bit_cast<uint64_t>(d + 0.0); }  // + 0.0 turns -0.0 into 0.0
        static uint64_t factor(uint32_t id, uint32_t exp) { return (static_cast<uint64_t>(id) << 32) | exp; }
        static uint32_t factor_id(uint64_t f) { return static_cast<uint32_t>(f >> 32); }
        static uint32_t factor_exp(uint64_t f) { return static_cast<uint32_t>(f); }
        static complex coeff_at(const uint64_t* t) { return {std::bit_cast<double>(t[1]), std::bit_cast<double>(t[2])}; }
        static const uint64_t* factors(const uint64_t* t) { return t + 3; }
        static const uint64_t* next(const uint64_t* t) { return t + 3 + t[0]; }
        static uint64_t* next_mut(uint64_t* t) { return t + 3 + t[0]; }

        /// Monomial order: lexicographic in the exponents, the symbol with the highest
        /// id being the most significant. The factor lists are compared from the back:
        /// a higher id, or the same id with a higher exponent, is the greater monomial.
        /// This is a true monomial order (m1 < m2 implies m1*m < m2*m), so multiplying
        /// a sorted polynomial by a monomial keeps it sorted, and it is the integer
        /// order of packed monomials (MonomialPacking). The constant term comes first.
        static int compare_mono(const uint64_t* a, size_t na, const uint64_t* b, size_t nb) {
                while (na > 0 && nb > 0) {
                        --na;
                        --nb;
                        if (a[na] != b[nb]) return a[na] < b[nb] ? -1 : 1;
                }
                return na == nb ? 0 : (na < nb ? -1 : 1);
        }

        /// Fixed array with the append interface multiply_mono() writes to.
        struct WordSpan {
                uint64_t* p;
                size_t n = 0;
                void push_back(uint64_t w) { p[n++] = w; }
                const uint64_t* end() const { return p + n; }
                void insert(const uint64_t*, const uint64_t* first, const uint64_t* last) {
                        std::copy(first, last, p + n);
                        n += static_cast<size_t>(last - first);
                }
        };

        static complex ipow(complex v, uint32_t e) {
                complex r(1.0, 0.0);
                while (e) {
                        if (e & 1) r *= v;
                        e >>= 1;
                        if (e) v *= v;
                }
                return r;
        }

        /// Appends the product of two sorted factor lists to `out`.
        template<typename Words>
        static void multiply_mono(const uint64_t* a, size_t na, const uint64_t* b, size_t nb, Words& out) {
                size_t i = 0, j = 0;
                while (i < na && j < nb) {
                        const uint32_t ia = factor_id(a[i]), ib = factor_id(b[j]);
                        if (ia < ib) out.push_back(a[i++]);
                        else if (ib < ia) out.push_back(b[j++]);
                        else {
                                const uint64_t e = static_cast<uint64_t>(factor_exp(a[i])) + factor_exp(b[j]);
                                if (e > 0xffffffffULL) throw std::overflow_error("Expression: exponent overflow");
                                out.push_back(factor(ia, static_cast<uint32_t>(e)));
                                ++i;
                                ++j;
                        }
                }
                out.insert(out.end(), a + i, a + na);
                out.insert(out.end(), b + j, b + nb);
        }

        /// Collects terms in any order and builds the canonical expression.
        class Builder {
            public:
                void begin(complex c) {
                        starts_.push_back(static_cast<uint32_t>(factors_.size()));
                        coeffs_.push_back(c);
                }
                void push(uint64_t f) { factors_.push_back(f); }
                void set_coeff(complex c) { coeffs_.back() = c; }
                std::vector<uint64_t>& factor_buffer() { return factors_; }
                /// Sorts the factors of the last term by symbol id (after renaming).
                void sort_last_factors() {
                        std::sort(factors_.begin() + starts_.back(), factors_.end());
                        // Renaming a and conjugate(a) can make two factors share an id only
                        // if both occurred, which renaming keeps distinct; nothing to merge.
                }

                Expression finish() {
                        const size_t n = coeffs_.size();
                        starts_.push_back(static_cast<uint32_t>(factors_.size()));
                        Expression e;
                        if (n == 0) return e;
                        const auto mono = [&](size_t k) { return factors_.data() + starts_[k]; };
                        const auto len = [&](size_t k) { return static_cast<size_t>(starts_[k + 1] - starts_[k]); };

                        std::vector<uint32_t> order(n);
                        for (size_t k = 0; k < n; ++k) order[k] = static_cast<uint32_t>(k);
                        bool sorted = true;
                        for (size_t k = 1; k < n && sorted; ++k) {
                                sorted = compare_mono(mono(k - 1), len(k - 1), mono(k), len(k)) < 0;
                        }
                        if (!sorted) {
                                std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
                                        const int c = compare_mono(mono(a), len(a), mono(b), len(b));
                                        return c != 0 ? c < 0 : a < b;
                                });
                        }

                        e.buf_.reserve(3 * n + factors_.size());
                        for (size_t k = 0; k < n;) {
                                const uint32_t head = order[k];
                                complex sum = coeffs_[head];
                                size_t j = k + 1;
                                while (j < n && compare_mono(mono(head), len(head), mono(order[j]), len(order[j])) == 0) {
                                        sum += coeffs_[order[j]];
                                        ++j;
                                }
                                if (sum != complex(0.0, 0.0)) {
                                        e.buf_.push_back(len(head));
                                        e.buf_.push_back(bits(sum.real()));
                                        e.buf_.push_back(bits(sum.imag()));
                                        e.buf_.insert(e.buf_.end(), mono(head), mono(head) + len(head));
                                        ++e.n_;
                                }
                                k = j;
                        }
                        return e;
                }

            private:
                std::vector<complex> coeffs_;
                std::vector<uint32_t> starts_;
                std::vector<uint64_t> factors_;
        };

        // -- arithmetic kernels -------------------------------------------------

        static Expression combine(const Expression& a, const Expression& b, bool subtract) {
                if (b.n_ == 0) return a;
                if (a.n_ == 0) return subtract ? -b : b;
                Expression out;
                out.buf_.reserve(a.buf_.size() + b.buf_.size());
                const uint64_t* p = a.buf_.data();
                const uint64_t* pe = p + a.buf_.size();
                const uint64_t* q = b.buf_.data();
                const uint64_t* qe = q + b.buf_.size();
                const double sign = subtract ? -1.0 : 1.0;
                const auto emit = [&out](const uint64_t* t, complex c) {
                        if (c == complex(0.0, 0.0)) return;
                        out.buf_.push_back(t[0]);
                        out.buf_.push_back(bits(c.real()));
                        out.buf_.push_back(bits(c.imag()));
                        out.buf_.insert(out.buf_.end(), t + 3, t + 3 + t[0]);
                        ++out.n_;
                };
                while (p < pe && q < qe) {
                        const int c = compare_mono(p + 3, p[0], q + 3, q[0]);
                        if (c < 0) {
                                emit(p, coeff_at(p));
                                p = next(p);
                        } else if (c > 0) {
                                emit(q, sign * coeff_at(q));
                                q = next(q);
                        } else {
                                emit(p, coeff_at(p) + sign * coeff_at(q));
                                p = next(p);
                                q = next(q);
                        }
                }
                for (; p < pe; p = next(p)) emit(p, coeff_at(p));
                for (; q < qe; q = next(q)) emit(q, sign * coeff_at(q));
                return out;
        }

        Expression scaled(complex s) const {
                Expression out = *this;
                out.scale_in_place(s);
                return out;
        }

        void scale_in_place(complex s) {
                if (n_ == 0) return;
                if (s == complex(0.0, 0.0)) {
                        buf_.clear();
                        n_ = 0;
                        return;
                }
                if (s == complex(1.0, 0.0)) return;
                bool any_zero = false;
                for (uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next_mut(t)) {
                        const complex c = coeff_at(t) * s;
                        t[1] = bits(c.real());
                        t[2] = bits(c.imag());
                        any_zero |= (c == complex(0.0, 0.0));
                }
                if (any_zero) drop_zero_terms();
        }

        void drop_zero_terms() {
                WordBuf kept;
                kept.reserve(buf_.size());
                uint32_t n = 0;
                for (const uint64_t* t = buf_.data(), *end = t + buf_.size(); t < end; t = next(t)) {
                        if (coeff_at(t) == complex(0.0, 0.0)) continue;
                        kept.insert(kept.end(), t, next(t));
                        ++n;
                }
                buf_ = std::move(kept);
                n_ = n;
        }

        static Expression multiply(const Expression& a, const Expression& b) {
                if (a.n_ == 0 || b.n_ == 0) return {};
                if (a.n_ == 1 && a.buf_[0] == 0) return b.scaled(coeff_at(a.buf_.data()));
                if (b.n_ == 1 && b.buf_[0] == 0) return a.scaled(coeff_at(b.buf_.data()));
                if (a.n_ == 1 && b.n_ == 1) {
                        const uint64_t* p = a.buf_.data();
                        const uint64_t* q = b.buf_.data();
                        const complex c = coeff_at(p) * coeff_at(q);
                        Expression out;
                        if (c == complex(0.0, 0.0)) return out;
                        out.buf_.reserve(3 + p[0] + q[0]);
                        out.buf_.resize(3);
                        multiply_mono(p + 3, p[0], q + 3, q[0], out.buf_);
                        out.buf_[0] = out.buf_.size() - 3;
                        out.buf_[1] = bits(c.real());
                        out.buf_[2] = bits(c.imag());
                        out.n_ = 1;
                        return out;
                }
                // Words of all pair products together (an upper bound for the result).
                const size_t count = static_cast<size_t>(a.n_) * b.n_;
                const size_t words = a.buf_.size() * b.n_ + b.buf_.size() * a.n_ - 3 * count;
                if (a.n_ == 1 || b.n_ == 1) {
                        // Monomial times polynomial: the order is kept and no two results
                        // are equal, so the terms are written straight to the result.
                        const Expression& mono = a.n_ == 1 ? a : b;
                        const Expression& poly = a.n_ == 1 ? b : a;
                        const uint64_t* m = mono.buf_.data();
                        const complex cm = coeff_at(m);
                        Expression out;
                        out.buf_.reserve(words);
                        for (const uint64_t* q = poly.buf_.data(), *qe = q + poly.buf_.size(); q < qe; q = next(q)) {
                                const complex c = cm * coeff_at(q);
                                if (c == complex(0.0, 0.0)) continue;  // underflow
                                const size_t head = out.buf_.size();
                                out.buf_.resize(head + 3);
                                multiply_mono(m + 3, m[0], q + 3, q[0], out.buf_);
                                out.buf_[head] = out.buf_.size() - head - 3;
                                out.buf_[head + 1] = bits(c.real());
                                out.buf_[head + 2] = bits(c.imag());
                                ++out.n_;
                        }
                        return out;
                }
                constexpr size_t kSmallTerms = 32;
                constexpr size_t kSmallWords = 320;
                if (count <= kSmallTerms && words <= kSmallWords) {
                        // Small product: all pair products on the stack, an insertion sort of
                        // their indices, then one pass that merges equal monomials.
                        uint64_t scratch[kSmallWords];
                        uint16_t start[kSmallTerms + 1];
                        uint8_t order[kSmallTerms];
                        WordSpan span{scratch};
                        size_t n = 0;
                        for (const uint64_t* p = a.buf_.data(), *pe = p + a.buf_.size(); p < pe; p = next(p)) {
                                for (const uint64_t* q = b.buf_.data(), *qe = q + b.buf_.size(); q < qe; q = next(q)) {
                                        const complex c = coeff_at(p) * coeff_at(q);
                                        start[n] = static_cast<uint16_t>(span.n);
                                        span.push_back(bits(c.real()));
                                        span.push_back(bits(c.imag()));
                                        multiply_mono(p + 3, p[0], q + 3, q[0], span);
                                        ++n;
                                }
                        }
                        start[n] = static_cast<uint16_t>(span.n);
                        const auto mono = [&](size_t k) { return scratch + start[k] + 2; };
                        const auto len = [&](size_t k) { return static_cast<size_t>(start[k + 1] - start[k] - 2); };
                        for (size_t k = 0; k < n; ++k) {
                                size_t j = k;
                                for (; j > 0 && compare_mono(mono(k), len(k), mono(order[j - 1]), len(order[j - 1])) < 0; --j) {
                                        order[j] = order[j - 1];
                                }
                                order[j] = static_cast<uint8_t>(k);
                        }
                        Expression out;
                        out.buf_.reserve(span.n + n);
                        for (size_t k = 0; k < n;) {
                                const size_t head = order[k];
                                complex sum(std::bit_cast<double>(scratch[start[head]]), std::bit_cast<double>(scratch[start[head] + 1]));
                                size_t j = k + 1;
                                for (; j < n && compare_mono(mono(head), len(head), mono(order[j]), len(order[j])) == 0; ++j) {
                                        const size_t o = start[order[j]];
                                        sum += complex(std::bit_cast<double>(scratch[o]), std::bit_cast<double>(scratch[o + 1]));
                                }
                                if (sum != complex(0.0, 0.0)) {
                                        out.buf_.push_back(len(head));
                                        out.buf_.push_back(bits(sum.real()));
                                        out.buf_.push_back(bits(sum.imag()));
                                        out.buf_.insert(out.buf_.end(), mono(head), mono(head) + len(head));
                                        ++out.n_;
                                }
                                k = j;
                        }
                        return out;
                }
                Builder builder;
                for (const uint64_t* p = a.buf_.data(), *pe = p + a.buf_.size(); p < pe; p = next(p)) {
                        for (const uint64_t* q = b.buf_.data(), *qe = q + b.buf_.size(); q < qe; q = next(q)) {
                                builder.begin(coeff_at(p) * coeff_at(q));
                                multiply_mono(p + 3, p[0], q + 3, q[0], builder.factor_buffer());
                        }
                }
                return builder.finish();
        }

        // -- printing -----------------------------------------------------------

        static std::string format_number(double x) {
                if (x == std::floor(x) && std::abs(x) < 1e15) {
                        return std::to_string(static_cast<long long>(x));
                }
                char buf[64];
                const auto res = std::to_chars(buf, buf + sizeof(buf), x);
                return std::string(buf, res.ptr);
        }

        /// The term without its leading sign; `negative` tells whether it has one.
        static std::string format_term(complex c, const std::string& mono, bool& negative) {
                std::string coeff;
                if (c.imag() == 0.0) {
                        negative = c.real() < 0.0;
                        const double m = std::abs(c.real());
                        if (mono.empty()) return format_number(m);
                        if (m == 1.0) return mono;
                        coeff = format_number(m);
                } else if (c.real() == 0.0) {
                        negative = c.imag() < 0.0;
                        const double m = std::abs(c.imag());
                        coeff = m == 1.0 ? "I" : format_number(m) + "*I";
                } else {
                        negative = false;
                        coeff = "(" + format_number(c.real()) + (c.imag() < 0.0 ? " - " : " + ")
                                + format_number(std::abs(c.imag())) + "*I)";
                }
                return mono.empty() ? coeff : coeff + "*" + mono;
        }

        friend class Parser;
};

// ---------------------------------------------------------------------------
// MonomialPacking
// ---------------------------------------------------------------------------

/// Packs the monomials of a product of two operands into one 64-bit word each:
/// every symbol gets a bit field wide enough for the highest exponent the product
/// can reach, the symbol with the highest id in the most significant bits. Then
/// the integer order of packed monomials is the canonical monomial order, and the
/// product of two monomials is the sum of their words. This lets a symbolic
/// Hamiltonian be handled as flat records (operator, monomial) -> complex number.
class MonomialPacking {
    public:
        /// Raises max_exp[id] to the highest exponent of each symbol in `e`.
        static void note_exponents(const Expression& e, std::vector<uint32_t>& max_exp) {
                e.for_each_term([&](std::complex<double>, const uint64_t* f, size_t nf) {
                        for (size_t i = 0; i < nf; ++i) {
                                const uint32_t id = Expression::symbol_of(f[i]);
                                if (id >= max_exp.size()) max_exp.resize(id + 1, 0);
                                max_exp[id] = std::max(max_exp[id], Expression::exponent_of(f[i]));
                        }
                });
        }

        /// Lays out the bit fields for products of monomials bounded by max_a and
        /// max_b (from note_exponents). Returns false if 64 bits are not enough.
        bool init(const std::vector<uint32_t>& max_a, const std::vector<uint32_t>& max_b) {
                const size_t n = std::max(max_a.size(), max_b.size());
                shift_.assign(n, kAbsent);
                slots_.clear();
                unsigned position = 0;
                for (size_t id = 0; id < n; ++id) {
                        const uint64_t bound = static_cast<uint64_t>(id < max_a.size() ? max_a[id] : 0)
                                + (id < max_b.size() ? max_b[id] : 0);
                        if (bound == 0) continue;
                        if (bound > 0xffffffffULL) return false;  // exponents are 32-bit
                        const unsigned width = static_cast<unsigned>(std::bit_width(bound));
                        if (position + width > 64) return false;
                        shift_[id] = static_cast<uint8_t>(position);
                        slots_.push_back({static_cast<uint32_t>(id), static_cast<uint8_t>(position),
                                          width == 64 ? ~uint64_t{0} : (uint64_t{1} << width) - 1});
                        position += width;
                }
                return true;
        }

        uint64_t pack(const uint64_t* f, size_t nf) const {
                uint64_t word = 0;
                for (size_t i = 0; i < nf; ++i) {
                        word += static_cast<uint64_t>(Expression::exponent_of(f[i])) << shift_[Expression::symbol_of(f[i])];
                }
                return word;
        }

        /// Appends coefficient * monomial to `out`; call with ascending words.
        void append(uint64_t word, std::complex<double> coefficient, Expression& out) const {
                uint64_t factors[64];
                size_t nf = 0;
                for (const Slot& s : slots_) {
                        const uint64_t e = (word >> s.shift) & s.mask;
                        if (e != 0) factors[nf++] = Expression::make_factor(s.id, static_cast<uint32_t>(e));
                }
                out.append_term_unchecked(coefficient, factors, nf);
        }

    private:
        static constexpr uint8_t kAbsent = 255;
        struct Slot {
                uint32_t id;
                uint8_t shift;
                uint64_t mask;
        };
        std::vector<uint8_t> shift_;  // by symbol id
        std::vector<Slot> slots_;     // ascending id
};

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Recursive-descent parser for polynomial expressions:
///   sum     := product (('+' | '-') product)*
///   product := unary (('*' | '/') unary | implicit unary)*     "2a" = 2*a
///   unary   := ('+' | '-') unary | power
///   power   := atom (('**' | '^') unary)?
///   atom    := number | name | name '(' sum ')' | '(' sum ')'
/// Names are symbols, except I (imaginary unit), pi and E; conjugate(...) and
/// conj(...) conjugate their argument. Division and non-integer powers are only
/// allowed for numbers, since the result must stay a polynomial.
class Parser {
    public:
        explicit Parser(const std::string& text) : s_(text) {}

        Expression parse() {
                Expression e = sum();
                ws();
                if (p_ != s_.size()) fail("unexpected '" + std::string(1, s_[p_]) + "'");
                return e;
        }

    private:
        const std::string& s_;
        size_t p_ = 0;

        [[noreturn]] void fail(const std::string& what) const {
                throw std::invalid_argument("cannot parse expression '" + s_ + "' at position "
                                            + std::to_string(p_) + ": " + what);
        }

        void ws() {
                while (p_ < s_.size() && (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n')) ++p_;
        }
        char peek() {
                ws();
                return p_ < s_.size() ? s_[p_] : '\0';
        }
        static bool name_start(char c) { return std::isalpha(static_cast<unsigned char>(c)) || c == '_'; }
        static bool name_char(char c) { return std::isalnum(static_cast<unsigned char>(c)) || c == '_'; }
        static bool digit(char c) { return c >= '0' && c <= '9'; }

        Expression sum() {
                Expression acc = product();
                for (;;) {
                        const char c = peek();
                        if (c == '+') {
                                ++p_;
                                acc = acc + product();
                        } else if (c == '-') {
                                ++p_;
                                acc = acc - product();
                        } else {
                                return acc;
                        }
                }
        }

        Expression product() {
                Expression acc = unary();
                for (;;) {
                        const char c = peek();
                        if (c == '*' && !(p_ + 1 < s_.size() && s_[p_ + 1] == '*')) {
                                ++p_;
                                acc = acc * unary();
                        } else if (c == '/') {
                                ++p_;
                                const Expression d = unary();
                                if (!d.is_number() || d.is_zero()) fail("division is only possible by a non-zero number");
                                acc = acc * (std::complex<double>(1.0, 0.0) / d.to_complex());
                        } else if (name_start(c) || c == '(') {
                                acc = acc * unary();  // implicit multiplication, as in "2a"
                        } else {
                                return acc;
                        }
                }
        }

        Expression unary() {
                const char c = peek();
                if (c == '-') {
                        ++p_;
                        return -unary();
                }
                if (c == '+') {
                        ++p_;
                        return unary();
                }
                return power();
        }

        Expression power() {
                Expression base = atom();
                const char c = peek();
                bool is_pow = false;
                if (c == '^') {
                        ++p_;
                        is_pow = true;
                } else if (c == '*' && p_ + 1 < s_.size() && s_[p_ + 1] == '*') {
                        p_ += 2;
                        is_pow = true;
                }
                if (!is_pow) return base;
                const Expression e = unary();
                if (!e.is_number()) fail("exponents must be numbers");
                const std::complex<double> v = e.to_complex();
                if (v.imag() == 0.0 && v.real() >= 0.0 && v.real() == std::floor(v.real()) && v.real() < 4294967296.0) {
                        return base.pow(static_cast<unsigned>(v.real()));
                }
                if (base.is_number()) return Expression(std::pow(base.to_complex(), v));
                fail("only non-negative integer powers of symbols are polynomials");
        }

        Expression atom() {
                const char c = peek();
                if (c == '(') {
                        ++p_;
                        Expression e = sum();
                        if (peek() != ')') fail("missing ')'");
                        ++p_;
                        return e;
                }
                if (digit(c) || c == '.') return number();
                if (name_start(c)) {
                        const size_t start = p_;
                        while (p_ < s_.size() && name_char(s_[p_])) ++p_;
                        const std::string name = s_.substr(start, p_ - start);
                        if ((name == "conjugate" || name == "conj") && peek() == '(') {
                                ++p_;
                                Expression e = sum();
                                if (peek() != ')') fail("missing ')'");
                                ++p_;
                                return e.conjugate();
                        }
                        if (name == "I") return Expression(std::complex<double>(0.0, 1.0));
                        if (name == "pi") return Expression(3.141592653589793238462643383279502884);
                        if (name == "E") return Expression(2.718281828459045235360287471352662498);
                        return Expression::symbol(name);
                }
                if (c == '\0') fail("unexpected end");
                fail("unexpected '" + std::string(1, c) + "'");
        }

        Expression number() {
                const size_t start = p_;
                while (p_ < s_.size() && (digit(s_[p_]) || s_[p_] == '.')) ++p_;
                if (p_ < s_.size() && (s_[p_] == 'e' || s_[p_] == 'E')) {
                        // exponent only if digits follow, so "2E" stays 2*E
                        size_t q = p_ + 1;
                        if (q < s_.size() && (s_[q] == '+' || s_[q] == '-')) ++q;
                        if (q < s_.size() && digit(s_[q])) {
                                p_ = q;
                                while (p_ < s_.size() && digit(s_[p_])) ++p_;
                        }
                }
                double v = 0.0;
                const auto res = std::from_chars(s_.data() + start, s_.data() + p_, v);
                if (res.ec != std::errc() || res.ptr != s_.data() + p_) fail("invalid number");
                return Expression(v);
        }
};

inline Expression::Expression(const std::string& text) {
        *this = Parser(text).parse();
}

}  // namespace pauliengine::symbolic
