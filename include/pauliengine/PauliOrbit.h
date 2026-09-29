#pragma once

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

#include "pauliengine/PauliOrbitSum.h"

// Permutation-invariant Pauli orbit basis (arXiv:2604.16701, Section V).
//
// This header provides only the orbit data structure and the commutator; the
// PauliOrbitSum data structure lives in PauliOrbitSum.h.

class PauliOrbit {
        public:
                int n_qubits;
                int p;  // number of X sites
                int q;  // number of Y sites
                int r;  // number of Z sites

                PauliOrbit(int n, int p_, int q_, int r_)
                        : n_qubits(n), p(p_), q(q_), r(r_) {}

                // Number of identity sites, n - p - q - r.
                int identities() const { return n_qubits - p - q - r; }

                bool is_valid() const {
                        return p >= 0 && q >= 0 && r >= 0 && identities() >= 0;
                }

                // Number of distinct Pauli strings in this orbit (Eq. 52).
                double n_terms() const {
                        const int s = identities();
                        if (p < 0 || q < 0 || r < 0 || s < 0) {
                                return 0.0;
                        }
                        const std::vector<double>& fact = factorial_table(n_qubits);
                        return fact[n_qubits] / (fact[p] * fact[q] * fact[r] * fact[s]);
                }

                bool operator==(const PauliOrbit& o) const {
                        return n_qubits == o.n_qubits && p == o.p && q == o.q && r == o.r;
                }

                bool operator!=(const PauliOrbit& o) const { return !(*this == o); }

                size_t hash() const {
                        size_t h = std::hash<int>{}(n_qubits);
                        for (int v : {p, q, r}) {
                                h ^= std::hash<int>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                        }
                        return h;
                }


                PauliOrbitSum commutator(const PauliOrbit& other) const {
                        PauliOrbitSum out(n_qubits);
                        const int n = n_qubits;
                        if (other.n_qubits != n || !is_valid() || !other.is_valid()) {
                                return out;
                        }
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
                                        accumulate_term(keys, vals, pt, qt, rt, sgn * W);
                                });

                        for (size_t k = 0; k < keys.size(); ++k) {
                                if (vals[k] != 0.0) {
                                        out.add(PauliOrbit(n, keys[k][0], keys[k][1], keys[k][2]), pref * vals[k]);
                                }
                        }
                        return out;
                }


                PauliOrbitSum multiply(const PauliOrbit& other) const {
                        PauliOrbitSum out(n_qubits);
                        const int n = n_qubits;
                        if (other.n_qubits != n || !is_valid() || !other.is_valid()) {
                                return out;
                        }
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
                                        accumulate_term(keys, vals, pt, qt, rt, phase * W);
                                });

                        for (size_t k = 0; k < keys.size(); ++k) {
                                if (vals[k] != std::complex<double>(0.0, 0.0)) {
                                        out.add(PauliOrbit(n, keys[k][0], keys[k][1], keys[k][2]), pref * vals[k]);
                                }
                        }
                        return out;
                }


                double targeted_commutator(const PauliOrbit& other, const PauliOrbit& target) const {
                        const int n = n_qubits;
                        if (other.n_qubits != n || target.n_qubits != n
                                || !is_valid() || !other.is_valid() || !target.is_valid()) {
                                return 0.0;
                        }

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

                        return pref * acc;
                }

        private:

                template <typename T>
                static void accumulate_term(std::vector<std::array<int, 3>>& keys,
                                            std::vector<T>& vals,
                                            int pt, int qt, int rt, T w) {
                        for (size_t k = 0; k < keys.size(); ++k) {
                                if (keys[k][0] == pt && keys[k][1] == qt && keys[k][2] == rt) {
                                        vals[k] += w;
                                        return;
                                }
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
                        const std::vector<double>& fact = factorial_table(n);
                        const double fn_n = fact[n];

                        for (int nXX = std::max(0, p1 - (n - p2)); nXX <= std::min(p1, p2); ++nXX) {
                          const int rX1 = p1 - nXX;
                          for (int nXY = std::max(0, rX1 - (n - q2)); nXY <= std::min(rX1, q2); ++nXY) {
                            const int rX2 = rX1 - nXY;
                            for (int nXZ = std::max(0, rX2 - (n - r2)); nXZ <= std::min(rX2, r2); ++nXZ) {
                              const int nXI = rX2 - nXZ;
                              if (nXI > s2) continue;
                              const int cp2 = p2 - nXX, cq2 = q2 - nXY, cr2 = r2 - nXZ, cs2 = s2 - nXI;

                              for (int nYX = std::max(0, q1 - (cq2 + cr2 + cs2)); nYX <= std::min(q1, cp2); ++nYX) {
                                const int rY1 = q1 - nYX;
                                for (int nYY = std::max(0, rY1 - (cr2 + cs2)); nYY <= std::min(rY1, cq2); ++nYY) {
                                  const int rY2 = rY1 - nYY;
                                  for (int nYZ = std::max(0, rY2 - cs2); nYZ <= std::min(rY2, cr2); ++nYZ) {
                                    const int nYI = rY2 - nYZ;
                                    if (nYI > cs2) continue;
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
                                          const double W = fn_n / (
                                                  fact[nXX] * fact[nXY] * fact[nXZ] * fact[nXI] *
                                                  fact[nYX] * fact[nYY] * fact[nYZ] * fact[nYI] *
                                                  fact[nZX] * fact[nZY] * fact[nZZ] * fact[nZI] *
                                                  fact[nIX] * fact[nIY] * fact[nIZ] * fact[nII]);
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

                // Factorials [0..n] as float64, grown on demand and cached across
                // calls
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
};


inline void PauliOrbitSum::add(const PauliOrbit& orbit, PauliOrbitSum::Coeff coeff) {
        for (size_t i = 0; i < data.size(); ++i) {
                if (data[i] == orbit) {
                        coeffs[i] += coeff;
                        return;
                }
        }
        data.push_back(orbit);
        coeffs.push_back(coeff);
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
}
