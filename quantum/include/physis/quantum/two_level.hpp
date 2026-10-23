// The open two-level system: a driven qubit that decays by emitting photons.
//
// A two-level system (an atom, a quantum dot, a superconducting qubit) has a ground state |g>
// and an excited state |e>. Driven at the Rabi frequency Omega, detuned by Delta from its
// transition, its Hamiltonian in the frame rotating with the drive is (hbar = 1, basis (g, e))
//   H = -Delta |e><e| + (Omega / 2)(|e><g| + |g><e|).
// It decays by spontaneous emission, L_1 = sqrt(gamma) sigma_-, and loses phase,
// L_2 = sqrt(gamma_phi / 2) sigma_z, and its density matrix obeys the GKSL master equation
// (Gorini, Kossakowski and Sudarshan 1976; Lindblad 1976)
//   rho' = -i [H, rho] + sum_k (L_k rho L_k^+ - {L_k^+ L_k, rho} / 2).
// On the Bloch vector r = (u, v, w) = (2 Re rho_eg, 2 Im rho_eg, rho_ee - rho_gg) the equation
// is affine (the optical Bloch equations), with the coherence decaying at Gamma_2 = gamma / 2 +
// gamma_phi:
//   u' = -Gamma_2 u - Delta v,   v' = Delta u - Gamma_2 v + Omega w,   w' = -Omega v - gamma (w + 1).
// Its stationary state has w = -1 / (1 + s), s = Omega^2 Gamma_2 / (gamma (Gamma_2^2 + Delta^2)),
// so that without dephasing rho_ee = (Omega^2 / 4) / (Delta^2 + gamma^2 / 4 + Omega^2 / 2).
//
// The same evolution unravels into quantum trajectories (Dalibard, Castin and Mølmer 1992;
// Carmichael 1993; Plenio and Knight 1998): with every emitted photon detected, the state
// between detections evolves by the effective Hamiltonian H_eff = H - (i gamma / 2) |e><e|,
// unnormalised, the square of its norm being the probability that no photon has yet been
// detected; at a detection it jumps to |g>. Averaged over trajectories (with the dephasing
// unravelled likewise), the state is the density matrix of the master equation. From |g>, the
// waiting time to the first photon has the density gamma |c_e(t)|^2, and its mean is the
// inverse of the stationary emission rate gamma rho_ee: the photons of one emitter form a
// renewal process, antibunched, since a photon is never followed at once by another
// (Cohen-Tannoudji and Dalibard 1986).

#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include <physis/dynamics/hybrid.hpp>

namespace physis::quantum {

using Complex = std::complex<double>;
// A 2x2 complex matrix, row-major: m[i][j].
using Matrix2 = std::array<std::array<Complex, 2>, 2>;
// A state vector (c_g, c_e).
using Vector2 = std::array<Complex, 2>;

struct TwoLevel {
  double Omega = 0;      // Rabi frequency, rad/s
  double Delta = 0;      // detuning of the drive from the transition, rad/s
  double gamma = 0;      // rate of spontaneous emission, 1/s
  double gamma_phi = 0;  // rate of pure dephasing, 1/s
};

// m c.
inline Vector2 apply(const Matrix2& m, const Vector2& c) {
  return {m[0][0] * c[0] + m[0][1] * c[1], m[1][0] * c[0] + m[1][1] * c[1]};
}

// exp(M) of a 2x2 complex matrix. With m = tr M / 2 and q^2 = m^2 - det M, the eigenvalues of
// M are m +- q, and by the theorem of Cayley and Hamilton exp M = e^m (cosh q I +
// (sinh q / q)(M - m I)), sinh q / q taken as its series where q is small.
inline Matrix2 exp2(const Matrix2& M) {
  const Complex m = (M[0][0] + M[1][1]) / 2.0;
  const Complex det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
  const Complex q = std::sqrt(m * m - det);
  const Complex c = std::cosh(q);
  const Complex s = std::abs(q) < 1e-4 ? 1.0 + q * q / 6.0 + q * q * q * q / 120.0 : std::sinh(q) / q;
  const Complex e = std::exp(m);
  return {{{e * (c + s * (M[0][0] - m)), e * s * M[0][1]}, {e * s * M[1][0], e * (c + s * (M[1][1] - m))}}};
}

// H_eff = H - (i gamma / 2) |e><e|, in the basis (g, e).
inline Matrix2 effective_hamiltonian(const TwoLevel& s) {
  return {{{Complex(0, 0), Complex(s.Omega / 2, 0)}, {Complex(s.Omega / 2, 0), Complex(-s.Delta, -s.gamma / 2)}}};
}

// The propagator between detections, exp(-i H_eff t).
inline Matrix2 no_jump_propagator(const TwoLevel& s, double t) {
  Matrix2 M = effective_hamiltonian(s);
  for (auto& row : M)
    for (auto& x : row) x *= Complex(0, -t);
  return exp2(M);
}

// The optical Bloch equations as a linear system of physis/dynamics, r' = A r + c, with no input:
// linear_propagator<Bloch>(s, h) integrates the master equation exactly over h.
struct Bloch {
  static constexpr std::size_t n = 3;
  using Params = TwoLevel;
  static void system(const TwoLevel& s, dynamics::Matrix& A, dynamics::State<3>& B, dynamics::State<3>& c) {
    const double G = s.gamma / 2 + s.gamma_phi;
    A(0, 0) = -G; A(0, 1) = -s.Delta;
    A(1, 0) = s.Delta; A(1, 1) = -G; A(1, 2) = s.Omega;
    A(2, 1) = -s.Omega; A(2, 2) = -s.gamma;
    B = {0, 0, 0};
    c = {0, 0, -s.gamma};
  }
};

// The stationary Bloch vector of the master equation (gamma > 0).
inline std::array<double, 3> stationary(const TwoLevel& s) {
  const double G = s.gamma / 2 + s.gamma_phi;
  const double sat = s.Omega * s.Omega * G / (s.gamma * (G * G + s.Delta * s.Delta));
  const double w = -1 / (1 + sat);
  const double v = s.Omega * G * w / (G * G + s.Delta * s.Delta);
  const double u = -s.Delta * v / G;
  return {u, v, w};
}

// The excited population of a Bloch vector, rho_ee = (1 + w) / 2.
inline double excited(const std::array<double, 3>& r) { return (1 + r[2]) / 2; }

// The stationary rate of emission, gamma rho_ee.
inline double emission_rate(const TwoLevel& s) { return s.gamma * excited(stationary(s)); }

}  // namespace physis::quantum
