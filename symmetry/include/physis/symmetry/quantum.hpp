// The recoupling theory of the quantum group U_q(sl₂) at the root of unity q = e^{iπ/r}: the
// quantum 6j symbols of the Turaev–Viro state sum.
//
// At q = e^{iπ/r} (r ≥ 3) the quantum integers [n] = (qⁿ − q⁻ⁿ)/(q − q⁻¹) = sin(nπ/r) / sin(π/r)
// are real, positive for 1 ≤ n ≤ r − 1, and [r] = 0. The representations that survive (the
// Verlinde truncation) are the spins j = 0, ½, …, (r − 2)/2, with quantum dimension [2j + 1];
// a triad (j₁, j₂, j₃) is **admissible** when it satisfies the classical triangle conditions
// and j₁ + j₂ + j₃ ≤ r − 2. The total quantum dimension is
//   w² = Σ_j [2j + 1]² = r / (2 sin²(π/r)).
// The quantum 6j symbol (Racah–Wigner normalisation, tetrahedrally symmetric, real) is
// Racah's formula with every factorial n! replaced by the quantum factorial [n]! = [1][2]…[n]
// (Kirillov and Reshetikhin 1989, in Infinite-Dimensional Lie Algebras and Groups, 285–339,
// eq. 4.12; Kauffman and Lins 1994, Temperley–Lieb Recoupling Theory, ch. 10). For an
// admissible symbol the denominators never contain [r]; terms whose numerator [z + 1]! does
// are zero. As r → ∞ it tends to the classical 6j symbol (recoupling.hpp).
//
// Sign conventions of state sums. The Turaev–Viro invariant of a triangulated closed
// 3-manifold with V vertices is, in this normalisation,
//   TV(M) = w^{−2V} Σ_colourings Π_edges (−1)^{2j_e} [2j_e + 1] Π_tetrahedra (−1)^{Σ 2j} {6j}_q
// up to conventions for the signs, which the state sum (physis quantum_gravity) fixes and
// checks by the invariance under change of triangulation (Turaev and Viro 1992, Topology 31,
// 865–902). This header gives the symbols; quantum_dimension_signed is (−1)^{2j}[2j + 1].
//
// Accuracy: double precision; the alternating sum loses about log₁₀ of its largest term over
// the result in digits, a few digits for spins below twenty.

#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include "recoupling.hpp"

namespace physis::symmetry {

// [n] at q = e^{iπ/r}.
inline double quantum_integer(int n, int r) {
  return std::sin(n * std::numbers::pi / r) / std::sin(std::numbers::pi / r);
}
// [n]! = [1]…[n]; exactly 0 when n ≥ r (the product then contains [r] = 0).
inline double quantum_factorial(int n, int r) {
  if (n < 0) throw std::domain_error("physis::symmetry::quantum_factorial: negative argument");
  if (n >= r) return 0;
  double f = 1;
  for (int k = 2; k <= n; k++) f *= quantum_integer(k, r);
  return f;
}
// [2j + 1].
inline double quantum_dimension(int two_j, int r) { return quantum_integer(two_j + 1, r); }
// (−1)^{2j} [2j + 1].
inline double quantum_dimension_signed(int two_j, int r) { return (two_j % 2 ? -1 : 1) * quantum_dimension(two_j, r); }
// w² = Σ_{j ≤ (r−2)/2} [2j + 1]² = r / (2 sin²(π/r)).
inline double total_quantum_dimension_squared(int r) {
  const double s = std::sin(std::numbers::pi / r);
  return r / (2 * s * s);
}

// Whether the triad (a, b, c), as two_j, is admissible at level r.
inline bool admissible(int a, int b, int c, int r) { return triangle(a, b, c) && (a + b + c) / 2 <= r - 2; }

// Δ_q(abc) = √([a+b−c]! [a−b+c]! [−a+b+c]! / [a+b+c+1]!) for an admissible triad.
inline double quantum_delta(int a, int b, int c, int r) {
  return std::sqrt(quantum_factorial((a + b - c) / 2, r) * quantum_factorial((a - b + c) / 2, r) * quantum_factorial((-a + b + c) / 2, r) /
                   quantum_factorial((a + b + c) / 2 + 1, r));
}

// The quantum 6j symbol {a b c; d e f}_q at q = e^{iπ/r}, arguments as two_j; zero unless its
// four triads (abc), (aef), (dbf), (dec) are admissible.
inline double quantum_6j(int r, int a, int b, int c, int d, int e, int f) {
  if (r < 3) throw std::domain_error("physis::symmetry::quantum_6j: r must be at least 3");
  if (!admissible(a, b, c, r) || !admissible(a, e, f, r) || !admissible(d, b, f, r) || !admissible(d, e, c, r)) return 0;
  const int al[4] = {(a + b + c) / 2, (a + e + f) / 2, (d + b + f) / 2, (d + e + c) / 2};
  const int be[3] = {(a + b + d + e) / 2, (b + c + e + f) / 2, (c + a + f + d) / 2};
  const int lo = *std::max_element(al, al + 4), hi = *std::min_element(be, be + 3);
  double s = 0;
  for (int z = lo; z <= hi; z++) {
    const double num = quantum_factorial(z + 1, r);
    if (num == 0) continue;
    double den = 1;
    for (int i = 0; i < 4; i++) den *= quantum_factorial(z - al[i], r);
    for (int k = 0; k < 3; k++) den *= quantum_factorial(be[k] - z, r);
    s += (z % 2 ? -1 : 1) * num / den;
  }
  return quantum_delta(a, b, c, r) * quantum_delta(a, e, f, r) * quantum_delta(d, b, f, r) * quantum_delta(d, e, c, r) * s;
}

}  // namespace physis::symmetry
