#pragma once

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "pauliengine/PauliOrbitSum.h"

/// Permutation-invariant Pauli orbit basis (arXiv:2604.16701, Section V).
///
/// B_{p,q,r} on n qubits is the normalized sum of all Pauli strings with p X-,
/// q Y- and r Z-sites. This header provides the orbit data structure and its
/// algebra; the PauliOrbitSum data structure lives in PauliOrbitSum.h.
class PauliOrbit {
        public:
                using Coeff = std::complex<double>;

                /// The structure constants are built from n! and its divisors in
                /// float64, which overflows at 171! (> 1.8e308). Beyond this the
                /// coefficients would silently become inf/nan, so construction is
                /// refused instead.
                static constexpr int kMaxQubits = 170;

                /// Fixed at construction (exposed read-only to Python): they define
                /// equality and the hash, so changing them would corrupt sets/dicts.
                int n_qubits;
                int p;  ///< number of X sites
                int q;  ///< number of Y sites
                int r;  ///< number of Z sites

                /// Orbit B_{p,q,r} on n qubits. Throws std::invalid_argument unless
                /// 0 <= n <= kMaxQubits, p, q, r >= 0 and p + q + r <= n.
                PauliOrbit(int n, int p_, int q_, int r_)
                        : n_qubits(n), p(p_), q(q_), r(r_) {
                        if (n < 0 || n > kMaxQubits) {
                                throw std::invalid_argument("PauliOrbit: n_qubits must be in [0, "
                                        + std::to_string(kMaxQubits) + "], got " + std::to_string(n));
                        }
                        if (!is_valid()) {
                                throw std::invalid_argument("PauliOrbit: invalid label (p=" + std::to_string(p)
                                        + ", q=" + std::to_string(q) + ", r=" + std::to_string(r)
                                        + "); need p, q, r >= 0 and p + q + r <= n_qubits = "
                                        + std::to_string(n));
                        }
                }

                /// Number of identity sites, n - p - q - r.
                int identities() const { return n_qubits - p - q - r; }

                /// Always true for a constructed orbit (the constructor rejects
                /// anything else); kept as a public query.
                bool is_valid() const {
                        return p >= 0 && q >= 0 && r >= 0 && identities() >= 0;
                }

                // Number of distinct Pauli strings in this orbit (Eq. 52).
                /// Number of distinct Pauli strings in the orbit, n! / (p! q! r! (n-p-q-r)!) (Eq. 52).
                double n_terms() const {
                        const int s = identities();
                        if (p < 0 || q < 0 || r < 0 || s < 0) {
                                return 0.0;
                        }
                        const std::vector<double>& fact = factorial_table(n_qubits);
                        return fact[n_qubits] / (fact[p] * fact[q] * fact[r] * fact[s]);
                }

                /// Equal iff n_qubits and the label (p, q, r) agree.
                bool operator==(const PauliOrbit& o) const {
                        return n_qubits == o.n_qubits && p == o.p && q == o.q && r == o.r;
                }

                /// Negation of operator==.
                bool operator!=(const PauliOrbit& o) const { return !(*this == o); }

                /// Fixed-size hash, consistent with operator==.
                size_t hash() const {
                        size_t h = std::hash<int>{}(n_qubits);
                        for (int v : {p, q, r}) {
                                h ^= std::hash<int>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                        }
                        return h;
                }


                /// Invertible key of the orbit label, consistent with operator==:
                /// n * 2^24 + p * 2^16 + q * 2^8 + r (all four are <= kMaxQubits < 256).
                uint64_t get_hash() const {
                        return (static_cast<uint64_t>(n_qubits) << 24) | (static_cast<uint64_t>(p) << 16)
                                | (static_cast<uint64_t>(q) << 8) | static_cast<uint64_t>(r);
                }

                /// Inverse of get_hash(); the constructor validates the decoded label.
                static PauliOrbit from_hash(uint64_t key) {
                        if (key >> 32) {
                                throw std::invalid_argument("PauliOrbit.from_hash: invalid hash");
                        }
                        return PauliOrbit(static_cast<int>((key >> 24) & 0xFF), static_cast<int>((key >> 16) & 0xFF),
                                static_cast<int>((key >> 8) & 0xFF), static_cast<int>(key & 0xFF));
                }

                /// Commutator [this, other] in the orbit basis (Theorem 1 / Prop. 7),
                /// via the contingency tables of the two labels; never expands the strings.
                PauliOrbitSum commutator(const PauliOrbit& other) const {
                        require_same_n_qubits(other, "commutator");
                        PauliOrbitSum out(n_qubits);
                        const int n = n_qubits;
                        const std::vector<double>& fact = factorial_table(n);
                        const double fn = fact[n];
                        const auto n_terms_of = [&](int P, int Q, int R) {
                                const int S = n - P - Q - R;
                                return fn / (fact[P] * fact[Q] * fact[R] * fact[S]);
                        };
                        const double pref = 2.0 / (n_terms_of(p, q, r)
                                * n_terms_of(other.p, other.q, other.r));

                        std::vector<std::array<int, 3>> keys;
                        std::vector<double> vals;
                        std::unordered_map<int, size_t> index;
                        // Only anticommuting tables (D odd) contribute; each carries the
                        // real sign (-1)^{E + (D-1)/2 + 1} (verified against the expanded
                        // operator algebra).
                        for_each_overlap_table(other,
                                [&](int pt, int qt, int rt, int E, int D, double W) {
                                        if ((D & 1) == 0) {
                                                return;
                                        }
                                        const int expo = E + (D - 1) / 2 + 1;
                                        const double sgn = (expo & 1) ? -1.0 : 1.0;
                                        accumulate_term(keys, vals, index, n, pt, qt, rt, sgn * W);
                                });

                        for (size_t k = 0; k < keys.size(); ++k) {
                                if (vals[k] != 0.0) {
                                        out.add(PauliOrbit(n, keys[k][0], keys[k][1], keys[k][2]), pref * vals[k]);
                                }
                        }
                        return out;
                }


                /// Product this * other in the orbit basis (complex coefficients).
                PauliOrbitSum multiply(const PauliOrbit& other) const {
                        require_same_n_qubits(other, "multiply");
                        PauliOrbitSum out(n_qubits);
                        const int n = n_qubits;
                        const std::vector<double>& fact = factorial_table(n);
                        const double fn = fact[n];
                        const auto n_terms_of = [&](int P, int Q, int R) {
                                const int S = n - P - Q - R;
                                return fn / (fact[P] * fact[Q] * fact[R] * fact[S]);
                        };
                        const double pref = 1.0 / (n_terms_of(p, q, r)
                                * n_terms_of(other.p, other.q, other.r));

                        std::vector<std::array<int, 3>> keys;
                        std::vector<std::complex<double>> vals;
                        std::unordered_map<int, size_t> index;
                        for_each_overlap_table(other,
                                [&](int pt, int qt, int rt, int E, int D, double W) {
                                        // phase = (-1)^E * i^{D+1}
                                        std::complex<double> phase;
                                        switch ((D + 1) & 3) {
                                                case 0:  phase = {1.0, 0.0};  break;
                                                case 1:  phase = {0.0, 1.0};  break;
                                                case 2:  phase = {-1.0, 0.0}; break;
                                                default: phase = {0.0, -1.0}; break;
                                        }
                                        if (E & 1) {
                                                phase = -phase;
                                        }
                                        accumulate_term(keys, vals, index, n, pt, qt, rt, phase * W);
                                });

                        for (size_t k = 0; k < keys.size(); ++k) {
                                if (vals[k] != std::complex<double>(0.0, 0.0)) {
                                        out.add(PauliOrbit(n, keys[k][0], keys[k][1], keys[k][2]), pref * vals[k]);
                                }
                        }
                        return out;
                }


                /// The structure constant is real; it is returned as Coeff (complex)
                /// for a uniform interface with PauliCycle and the sum classes.
                Coeff targeted_commutator(const PauliOrbit& other, const PauliOrbit& target) const {
                        require_same_n_qubits(other, "targeted_commutator");
                        require_same_n_qubits(target, "targeted_commutator");
                        const int n = n_qubits;

                        const int p1 = p, q1 = q, r1 = r, s1 = n - p1 - q1 - r1;
                        const int p2 = other.p, q2 = other.q, r2 = other.r, s2 = n - p2 - q2 - r2;
                        const int tp = target.p, tq = target.q, tr = target.r;

                        const std::vector<double>& fact = factorial_table(n);
                        const double fn = fact[n];
                        const auto n_terms_of = [&](int P, int Q, int R) {
                                const int S = n - P - Q - R;
                                return fn / (fact[P] * fact[Q] * fact[R] * fact[S]);
                        };
                        const double pref = 2.0 / (n_terms_of(p1, q1, r1) * n_terms_of(p2, q2, r2));

                        double acc = 0.0;

                        // Two exact routes to the same constant: split the target's
                        // letters over the table cells (cost ~ (tp tq tr)^3, Corollary 5;
                        // best for a light target), or enumerate all overlap tables of the
                        // operands and keep those producing the target (cost ~ (p q r)^3;
                        // best when the operands are lighter than the target). Pick the
                        // cheaper one by its loop-count bound.
                        const double split_cost = compositions4(tp) * compositions4(tq) * compositions4(tr);
                        const double table_cost = compositions4(p1) * compositions4(q1) * compositions4(r1);
                        if (table_cost <= split_cost) {
                                for_each_overlap_table(other,
                                        [&](int pt, int qt, int rt, int E, int D, double W) {
                                                if ((D & 1) == 0 || pt != tp || qt != tq || rt != tr) {
                                                        return;
                                                }
                                                const int expo = E + (D - 1) / 2 + 1;
                                                acc += (expo & 1) ? -W : W;
                                        });
                                return Coeff(pref * acc, 0.0);
                        }

                        // Split target X-count over (XI, IX, YZ, ZY).
                        for (int nXI = 0; nXI <= tp; ++nXI)
                          for (int nIX = 0; nIX <= tp - nXI; ++nIX)
                            for (int nYZ = 0; nYZ <= tp - nXI - nIX; ++nYZ) {
                              const int nZY = tp - nXI - nIX - nYZ;
                              // Split target Y-count over (YI, IY, ZX, XZ).
                              for (int nYI = 0; nYI <= tq; ++nYI)
                                for (int nIY = 0; nIY <= tq - nYI; ++nIY)
                                  for (int nZX = 0; nZX <= tq - nYI - nIY; ++nZX) {
                                    const int nXZ = tq - nYI - nIY - nZX;
                                    // Split target Z-count over (ZI, IZ, XY, YX).
                                    for (int nZI = 0; nZI <= tr; ++nZI)
                                      for (int nIZ = 0; nIZ <= tr - nZI; ++nIZ)
                                        for (int nXY = 0; nXY <= tr - nZI - nIZ; ++nXY) {
                                          const int nYX = tr - nZI - nIZ - nXY;

                                          // Diagonals from the row sums.
                                          const int nXX = p1 - nXY - nXZ - nXI;
                                          const int nYY = q1 - nYX - nYZ - nYI;
                                          const int nZZ = r1 - nZX - nZY - nZI;
                                          const int nII = s1 - nIX - nIY - nIZ;
                                          if (nXX < 0 || nYY < 0 || nZZ < 0 || nII < 0) continue;

                                          // Column sums must match `other`.
                                          if (nXX + nYX + nZX + nIX != p2) continue;
                                          if (nXY + nYY + nZY + nIY != q2) continue;
                                          if (nXZ + nYZ + nZZ + nIZ != r2) continue;
                                          if (nXI + nYI + nZI + nII != s2) continue;

                                          // Parity: only anticommuting patterns (D odd).
                                          const int D = nXY + nXZ + nYX + nYZ + nZX + nZY;
                                          if ((D & 1) == 0) continue;

                                          const double W = fn / (
                                                  fact[nXX] * fact[nXY] * fact[nXZ] * fact[nXI] *
                                                  fact[nYX] * fact[nYY] * fact[nYZ] * fact[nYI] *
                                                  fact[nZX] * fact[nZY] * fact[nZZ] * fact[nZI] *
                                                  fact[nIX] * fact[nIY] * fact[nIZ] * fact[nII]);

                                          const int E = nYX + nZY + nXZ;
                                          const int expo = E + (D - 1) / 2 + 1;
                                          const double sgn = (expo & 1) ? -1.0 : 1.0;
                                          acc += sgn * W;
                                        }
                                  }
                            }

                        return Coeff(pref * acc, 0.0);
                }

        private:

                void require_same_n_qubits(const PauliOrbit& other, const char* op) const {
                        if (other.n_qubits != n_qubits) {
                                throw std::invalid_argument(std::string("PauliOrbit.") + op
                                        + ": n_qubits mismatch (" + std::to_string(n_qubits)
                                        + " vs " + std::to_string(other.n_qubits) + ")");
                        }
                }

                /// Adds w to the output label (pt, qt, rt), O(1) via a label index.
                template <typename T>
                static void accumulate_term(std::vector<std::array<int, 3>>& keys,
                                            std::vector<T>& vals,
                                            std::unordered_map<int, size_t>& index,
                                            int n, int pt, int qt, int rt, T w) {
                        const int key = (pt * (n + 1) + qt) * (n + 1) + rt;
                        const auto [it, inserted] = index.try_emplace(key, keys.size());
                        if (!inserted) {
                                vals[it->second] += w;
                                return;
                        }
                        keys.push_back({pt, qt, rt});
                        vals.push_back(w);
                }


                template <typename F>
                void for_each_overlap_table(const PauliOrbit& other, F&& fn) const {
                        const int n = n_qubits;
                        const int p1 = p, q1 = q, r1 = r;
                        const int p2 = other.p, q2 = other.q, r2 = other.r;
                        const int s2 = n - p2 - q2 - r2;
                        const double fn_n = factorial_table(n)[n];
                        const std::vector<double>& inv = inv_factorial_table(n);

                        // W = n! / prod(cell!) is built from inverse factorials, with
                        // the X- and Y-row factors hoisted out of the inner loops. All
                        // partial products are >= 1/n! >= 1/170!, so nothing underflows.
                        for (int nXX = std::max(0, p1 - (n - p2)); nXX <= std::min(p1, p2); ++nXX) {
                          const int rX1 = p1 - nXX;
                          for (int nXY = std::max(0, rX1 - (n - q2)); nXY <= std::min(rX1, q2); ++nXY) {
                            const int rX2 = rX1 - nXY;
                            for (int nXZ = std::max(0, rX2 - (n - r2)); nXZ <= std::min(rX2, r2); ++nXZ) {
                              const int nXI = rX2 - nXZ;
                              if (nXI > s2) continue;
                              const double rowX = inv[nXX] * inv[nXY] * inv[nXZ] * inv[nXI];
                              const int cp2 = p2 - nXX, cq2 = q2 - nXY, cr2 = r2 - nXZ, cs2 = s2 - nXI;

                              for (int nYX = std::max(0, q1 - (cq2 + cr2 + cs2)); nYX <= std::min(q1, cp2); ++nYX) {
                                const int rY1 = q1 - nYX;
                                for (int nYY = std::max(0, rY1 - (cr2 + cs2)); nYY <= std::min(rY1, cq2); ++nYY) {
                                  const int rY2 = rY1 - nYY;
                                  for (int nYZ = std::max(0, rY2 - cs2); nYZ <= std::min(rY2, cr2); ++nYZ) {
                                    const int nYI = rY2 - nYZ;
                                    if (nYI > cs2) continue;
                                    const double rowXY = rowX * inv[nYX] * inv[nYY] * inv[nYZ] * inv[nYI];
                                    const int dp2 = cp2 - nYX, dq2 = cq2 - nYY, dr2 = cr2 - nYZ, ds2 = cs2 - nYI;

                                    for (int nZX = std::max(0, r1 - (dq2 + dr2 + ds2)); nZX <= std::min(r1, dp2); ++nZX) {
                                      const int rZ1 = r1 - nZX;
                                      for (int nZY = std::max(0, rZ1 - (dr2 + ds2)); nZY <= std::min(rZ1, dq2); ++nZY) {
                                        const int rZ2 = rZ1 - nZY;
                                        for (int nZZ = std::max(0, rZ2 - ds2); nZZ <= std::min(rZ2, dr2); ++nZZ) {
                                          const int nZI = rZ2 - nZZ;
                                          const int nIX = dp2 - nZX;
                                          const int nIY = dq2 - nZY;
                                          const int nIZ = dr2 - nZZ;
                                          const int nII = ds2 - nZI;

                                          const int D = nXY + nXZ + nYX + nYZ + nZX + nZY;
                                          const int pt = nXI + nIX + nYZ + nZY;
                                          const int qt = nYI + nIY + nZX + nXZ;
                                          const int rt = nZI + nIZ + nXY + nYX;
                                          const double W = fn_n * (rowXY
                                                  * inv[nZX] * inv[nZY] * inv[nZZ] * inv[nZI]
                                                  * inv[nIX] * inv[nIY] * inv[nIZ] * inv[nII]);
                                          const int E = nYX + nZY + nXZ;
                                          fn(pt, qt, rt, E, D, W);
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                }

                /// Factorials [0..n] as float64, grown on demand and cached across
                /// calls
                static const std::vector<double>& factorial_table(int n) {
                        thread_local std::vector<double> fact = {1.0};
                        if (static_cast<int>(fact.size()) <= n) {
                                const size_t old = fact.size();
                                fact.resize(static_cast<size_t>(n) + 1);
                                for (size_t i = old; i < fact.size(); ++i) {
                                        fact[i] = fact[i - 1] * static_cast<double>(i);
                                }
                        }
                        return fact;
                }

                /// 1/k! for k in [0..n], grown on demand and cached like factorial_table.
                static const std::vector<double>& inv_factorial_table(int n) {
                        thread_local std::vector<double> inv = {1.0};
                        if (static_cast<int>(inv.size()) <= n) {
                                const std::vector<double>& fact = factorial_table(n);
                                const size_t old = inv.size();
                                inv.resize(static_cast<size_t>(n) + 1);
                                for (size_t i = old; i < inv.size(); ++i) {
                                        inv[i] = 1.0 / fact[i];
                                }
                        }
                        return inv;
                }

                /// Number of ways to split k sites over 4 table cells, C(k + 3, 3):
                /// the loop count of one row (or one target letter) in the enumerations.
                static double compositions4(int k) {
                        return (k + 1.0) * (k + 2.0) * (k + 3.0) / 6.0;
                }
};


inline PauliOrbitSum::PauliOrbitSum(int n) : n_qubits(n) {
        if (n < 0 || n > PauliOrbit::kMaxQubits) {
                throw std::invalid_argument("PauliOrbitSum: n_qubits must be in [0, "
                        + std::to_string(PauliOrbit::kMaxQubits) + "], got " + std::to_string(n));
        }
}

inline int PauliOrbitSum::label_key(int p, int q, int r) {
        // p, q, r <= kMaxQubits, so this packing is collision-free and fits an int.
        constexpr int base = PauliOrbit::kMaxQubits + 1;
        return (p * base + q) * base + r;
}

inline void PauliOrbitSum::add(const PauliOrbit& orbit, PauliOrbitSum::Coeff coeff) {
        if (orbit.n_qubits != n_qubits) {
                throw std::invalid_argument("PauliOrbitSum.add: n_qubits mismatch ("
                        + std::to_string(n_qubits) + " vs " + std::to_string(orbit.n_qubits) + ")");
        }
        const auto [it, inserted] = index_.try_emplace(label_key(orbit.p, orbit.q, orbit.r), data.size());
        if (!inserted) {
                coeffs[it->second] += coeff;
                return;
        }
        data.push_back(orbit);
        coeffs.push_back(coeff);
}

inline PauliOrbitSum::Coeff PauliOrbitSum::contains(const PauliOrbit& orbit) const {
        if (orbit.n_qubits != n_qubits) {
                return Coeff(0.0, 0.0);
        }
        const auto it = index_.find(label_key(orbit.p, orbit.q, orbit.r));
        return it == index_.end() ? Coeff(0.0, 0.0) : coeffs[it->second];
}

inline PauliOrbitSum PauliOrbitSum::commutator(const PauliOrbitSum& other) const {
        // [ sum_i a_i B_i , sum_j b_j B_j ] = sum_{i,j} a_i b_j [B_i, B_j].
        PauliOrbitSum result(n_qubits);
        for (size_t i = 0; i < data.size(); ++i) {
                for (size_t j = 0; j < other.data.size(); ++j) {
                        const PauliOrbitSum part = data[i].commutator(other.data[j]);
                        const Coeff scale = coeffs[i] * other.coeffs[j];
                        for (size_t k = 0; k < part.data.size(); ++k) {
                                result.add(part.data[k], scale * part.coeffs[k]);
                        }
                }
        }
        result.prune_zeros();
        return result;
}

inline PauliOrbitSum PauliOrbitSum::multiply(const PauliOrbitSum& other) const {
        // ( sum_i a_i B_i )( sum_j b_j B_j ) = sum_{i,j} a_i b_j B_i B_j.
        PauliOrbitSum result(n_qubits);
        for (size_t i = 0; i < data.size(); ++i) {
                for (size_t j = 0; j < other.data.size(); ++j) {
                        const PauliOrbitSum part = data[i].multiply(other.data[j]);
                        const Coeff scale = coeffs[i] * other.coeffs[j];
                        for (size_t k = 0; k < part.data.size(); ++k) {
                                result.add(part.data[k], scale * part.coeffs[k]);
                        }
                }
        }
        result.prune_zeros();
        return result;
}

inline PauliOrbitSum::Coeff PauliOrbitSum::targeted_commutator(const PauliOrbitSum& other,
                                                               const PauliOrbit& target) const {
        Coeff total(0.0, 0.0);
        for (size_t i = 0; i < data.size(); ++i) {
                for (size_t j = 0; j < other.data.size(); ++j) {
                        total += coeffs[i] * other.coeffs[j]
                                * data[i].targeted_commutator(other.data[j], target);
                }
        }
        return total;
}

inline void PauliOrbitSum::prune_zeros() {
        std::vector<PauliOrbit> kept_data;
        std::vector<Coeff> kept_coeffs;
        for (size_t i = 0; i < data.size(); ++i) {
                if (coeffs[i] != Coeff(0.0, 0.0)) {
                        kept_data.push_back(data[i]);
                        kept_coeffs.push_back(coeffs[i]);
                }
        }
        data = std::move(kept_data);
        coeffs = std::move(kept_coeffs);

        index_.clear();
        for (size_t i = 0; i < data.size(); ++i) {
                index_.emplace(label_key(data[i].p, data[i].q, data[i].r), i);
        }
}
