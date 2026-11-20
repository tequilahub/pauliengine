#pragma once
/// @file
/// Storage of the x and y words of a Pauli string.

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

/// Small-buffer word storage for the symplectic bit representation: up to
/// kInlineWords 64-bit words (= 128 qubits) live inline without touching the
/// heap; longer strings fall back to a single heap allocation. Pauli string
/// multiplication creates one short-lived string per pair, so avoiding
/// malloc/free on that path matters. Size and capacity are 32-bit, which keeps
/// the object at 32 bytes; a sum of millions of strings is written faster.
class WordVec {
    public:
        static constexpr size_t kInlineWords = 2;

        WordVec() noexcept : size_(0), cap_(kInlineWords), ptr_(buf_) {}

        explicit WordVec(size_t n) : WordVec() { resize(n); }

        WordVec(const WordVec& o) : WordVec() { assign(o.ptr_, o.ptr_ + o.size_); }

        WordVec(WordVec&& o) noexcept : size_(o.size_), cap_(kInlineWords), ptr_(buf_) {
                if (o.ptr_ != o.buf_) {
                        ptr_ = o.ptr_;
                        cap_ = o.cap_;
                } else {
                        for (size_t i = 0; i < size_; ++i) buf_[i] = o.buf_[i];
                }
                o.ptr_ = o.buf_;
                o.size_ = 0;
                o.cap_ = kInlineWords;
        }

        WordVec(const std::vector<uint64_t>& v) : WordVec() { assign(v.data(), v.data() + v.size()); }

        WordVec& operator=(const WordVec& o) {
                if (this != &o) assign(o.ptr_, o.ptr_ + o.size_);
                return *this;
        }

        WordVec& operator=(WordVec&& o) noexcept {
                if (this == &o) return *this;
                if (ptr_ != buf_) delete[] ptr_;
                if (o.ptr_ != o.buf_) {
                        ptr_ = o.ptr_;
                        cap_ = o.cap_;
                } else {
                        ptr_ = buf_;
                        cap_ = kInlineWords;
                        for (size_t i = 0; i < o.size_; ++i) buf_[i] = o.buf_[i];
                }
                size_ = o.size_;
                o.ptr_ = o.buf_;
                o.size_ = 0;
                o.cap_ = kInlineWords;
                return *this;
        }

        ~WordVec() {
                if (ptr_ != buf_) delete[] ptr_;
        }

        size_t size() const noexcept { return size_; }
        bool empty() const noexcept { return size_ == 0; }
        uint64_t* data() noexcept { return ptr_; }
        const uint64_t* data() const noexcept { return ptr_; }
        uint64_t& operator[](size_t i) noexcept { return ptr_[i]; }
        const uint64_t& operator[](size_t i) const noexcept { return ptr_[i]; }
        uint64_t* begin() noexcept { return ptr_; }
        uint64_t* end() noexcept { return ptr_ + size_; }
        const uint64_t* begin() const noexcept { return ptr_; }
        const uint64_t* end() const noexcept { return ptr_ + size_; }

        /// Grows zero-filled; shrinking never reallocates.
        void resize(size_t n) {
                const uint32_t m = checked_count(n);
                if (m > cap_) {
                        uint64_t* p = new uint64_t[m];
                        for (size_t i = 0; i < size_; ++i) p[i] = ptr_[i];
                        if (ptr_ != buf_) delete[] ptr_;
                        ptr_ = p;
                        cap_ = m;
                }
                for (size_t i = size_; i < m; ++i) ptr_[i] = 0;
                size_ = m;
        }

        void assign(const uint64_t* first, const uint64_t* last) {
                const uint32_t n = checked_count(static_cast<size_t>(last - first));
                if (n > cap_) {
                        uint64_t* p = new uint64_t[n];
                        if (ptr_ != buf_) delete[] ptr_;
                        ptr_ = p;
                        cap_ = n;
                }
                for (size_t i = 0; i < n; ++i) ptr_[i] = first[i];
                size_ = n;
        }

        std::vector<uint64_t> to_vector() const { return std::vector<uint64_t>(ptr_, ptr_ + size_); }

        friend bool operator==(const WordVec& a, const WordVec& b) {
                if (a.size_ != b.size_) return false;
                for (size_t i = 0; i < a.size_; ++i) {
                        if (a.ptr_[i] != b.ptr_[i]) return false;
                }
                return true;
        }
        friend bool operator!=(const WordVec& a, const WordVec& b) { return !(a == b); }

    private:
        /// A word count as stored; more than 2^32 - 1 words (about 2.7e11 qubits) is refused.
        static uint32_t checked_count(size_t n) {
                if (n > UINT32_MAX) {
                        throw std::length_error("WordVec: more than 2^32 - 1 words");
                }
                return static_cast<uint32_t>(n);
        }

        uint32_t size_;
        uint32_t cap_;
        uint64_t* ptr_;
        uint64_t buf_[kInlineWords];
};

static_assert(sizeof(WordVec) == 32, "WordVec should stay at 32 bytes");
