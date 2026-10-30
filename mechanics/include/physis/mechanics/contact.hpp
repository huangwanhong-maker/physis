// Contact with the ground, compliant: a point pressed into a ground line feels a normal force of
// Hunt–Crossley's law and a regularised Coulomb friction (decision 0014; docs/theory/discrete-mechanics.md
// §7).
//
// The ground. A line in the plane (x, z): z = z₀ + s x, its outward normal n = (−s, 1)/√(1 + s²) and its
// tangent t = (1, s)/√(1 + s²). A point p is pressed into it by δ = −n·(p − (0, z₀)) when that is positive.
//
// The normal force (Hunt and Crossley 1975, J. Appl. Mech. 42, 440–445): f_n = k δⁿ (1 + c δ̇), never
// pulling (f_n ≥ 0). The exponent n = 3/2 is Hertz's for a sphere on a half-space; c (s/m) sets the loss
// of an impact: for small c v₀ the coefficient of restitution is e ≈ 1 − (2/3) c v₀ (Marhefka and Orin
// 1999, IEEE Trans. Syst. Man Cybern. A 29, 566–572), so a landing's loss grows with its speed, and the
// force has no jump at first touch, unlike a spring and a dashpot.
//
// Friction. Coulomb's law f_t = −μ f_n sign(v_t), regularised as −μ f_n tanh(v_t / v_s): below the slip
// speed v_s the ground holds like a stiff viscous damper (f_t ≈ −μ f_n v_t / v_s), above it the full
// Coulomb force acts. A point held by it creeps at a speed of the order of v_s times the fraction of the
// friction cone used; v_s is the model's error and is stated with it (Stewart 2000, SIAM Rev. 42, 3–39, on
// compliant and rigid models of contact).
//
// Both are smooth enough to be differentiated (calx jets), so a variational integrator can take them as
// forces (Lagrange–d'Alembert), and both are evaluated on doubles by the explicit integrators.
//
// What is not modelled: the ground's own motion and deformation beyond δ (a soil's sinkage and bulldozing,
// terramechanics), rolling resistance, and the stiction of a body at rest (the regularisation creeps).

#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

namespace physis::mechanics {

struct GroundLine {
  double z0 = 0;     // the ground's height at x = 0 (m)
  double slope = 0;  // dz/dx

  std::array<double, 2> normal() const {
    const double r = std::sqrt(1 + slope * slope);
    return {-slope / r, 1 / r};
  }
  std::array<double, 2> tangent() const {
    const double r = std::sqrt(1 + slope * slope);
    return {1 / r, slope / r};
  }
  double height(double x) const { return z0 + slope * x; }
  // The signed distance of p above the ground (negative when pressed in).
  template <class S>
  S distance(const std::array<S, 2>& p) const {
    const auto n = normal();
    return p[0] * n[0] + (p[1] - z0) * n[1];
  }
};

struct CompliantContact {
  double stiffness = 1e6;   // k, in N/mⁿ
  double exponent = 1.5;    // n
  double damping = 1.0;     // c, in s/m
  double friction = 0.7;    // μ
  double slip_speed = 0.01; // v_s, in m/s

  void validate() const {
    if (!(stiffness > 0) || !(exponent >= 1) || !(damping >= 0) || !(friction >= 0) || !(slip_speed > 0))
      throw std::invalid_argument("physis::mechanics::CompliantContact: a positive stiffness and slip speed, n ≥ 1, c and μ not negative");
  }
};

// The force on a point at p moving at pdot, and its normal part (0 when it does not touch).
template <class S>
struct ContactForce {
  std::array<S, 2> force{};
  S normal{};
  S penetration{};
  bool touching = false;
};

template <class S>
ContactForce<S> contact_force(const GroundLine& ground, const CompliantContact& law, const std::array<S, 2>& p, const std::array<S, 2>& pdot) {
  using std::pow, std::tanh;
  ContactForce<S> out;
  const S d = ground.distance(p);
  if (!(d < 0.0)) {
    out.force = {S(0.0), S(0.0)};
    out.normal = S(0.0);
    out.penetration = S(0.0);
    return out;
  }
  const auto n = ground.normal();
  const auto t = ground.tangent();
  const S delta = -d;
  const S delta_rate = -(pdot[0] * n[0] + pdot[1] * n[1]);
  S fn = pow(delta, law.exponent) * law.stiffness * (1.0 + delta_rate * law.damping);
  if (fn < 0.0) fn = S(0.0);
  const S vt = pdot[0] * t[0] + pdot[1] * t[1];
  const S ft = -(fn * law.friction) * tanh(vt / law.slip_speed);
  out.force = {fn * n[0] + ft * t[0], fn * n[1] + ft * t[1]};
  out.normal = fn;
  out.penetration = delta;
  out.touching = true;
  return out;
}

}  // namespace physis::mechanics
