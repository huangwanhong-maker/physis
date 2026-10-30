// Checks of a tree in a fluid (mechanics/fluid.hpp) and of the Euler–Lagrange equations solved for the
// acceleration (mechanics/lagrangian.hpp).
//
// - lagrangian_acceleration on a tree alone is the tree's forward dynamics (Newton–Euler), to rounding.
// - A plate (one link, 0.5 m long, 0.2 m high, 1 kg) in water:
//   - pushed across itself from rest it accelerates as F/(m + ρπh²ℓ/4): the added mass;
//   - moving across itself it slows as v₀/(1 + k v₀ t/(m + m_a)), k = ½ρC_n hℓ; along itself as
//     v₀/(1 + k_t v₀ t/m), k_t = ½ρC_f 2hℓ, with no added mass: the closed forms of quadratic drag;
//   - at rest in a current across it, it is carried toward the current's speed, as quadratic drag allows;
//   - without drag, turning and moving, it keeps the energy of the body and the fluid together;
//   - the variational integrator of the tree and the fluid, its drag as a force, agrees with Runge–Kutta.

#include <physis/mechanics/fluid.hpp>
#include <physis/mechanics/lagrangian.hpp>
#include <physis/mechanics/variational.hpp>

#include <calx/verification/check.hpp>

#include <array>
#include <cmath>
#include <numbers>
#include <vector>

namespace pm = physis::mechanics;

namespace {

struct Drag {
  const pm::PlanarTree<3>* tree;
  const std::vector<pm::SlenderPart>* parts;
  const pm::Fluid* fluid;
  template <class S>
  std::array<S, 3> operator()(double, const std::array<S, 3>& q, const std::array<S, 3>& v) const {
    return pm::fluid_drag<3>(*tree, tree->kinematics(q, v), *parts, *fluid);
  }
};

}  // namespace

int main() {
  calx::verification::Suite s;
  {
    // A floating tree of three links, at a general state: the two ways to its acceleration.
    std::vector<pm::PlanarLink> links = {{"a", -1, {0, 0}, 0.0, 1.2, {0.05, 0.01}, 0.02, 0},
                                         {"b", 0, {0.2, 0.0}, 0.1, 0.4, {0.1, 0.0}, 0.003, 1e-3},
                                         {"c", 1, {0.2, 0.0}, -0.2, 0.3, {0.08, 0.0}, 0.002, 5e-4}};
    const pm::PlanarTree<5> tree(pm::Base::floating, links, 9.81);
    const std::array<double, 5> q{0.3, -0.2, 0.4, 0.7, -0.9}, v{0.5, -0.3, 1.1, -2.0, 0.6}, Q{0.4, -1.0, 0.2, 0.3, -0.1};
    const auto a = tree.forward_dynamics(q, v, Q);
    const auto b = pm::lagrangian_acceleration<5>(pm::TreeLagrangian<5, pm::NoPotential>{&tree, {}}, 0.0, q, v, Q);
    double worst = 0, size = 0;
    for (int k = 0; k < 5; k++) worst = std::max(worst, std::fabs(a[k] - b[k])), size = std::max(size, std::fabs(a[k]));
    s.check("lagrangian: the Euler–Lagrange equations by jets give the tree's Newton–Euler acceleration (relative)", worst / size, 1e-12);
  }
  // The plate: one link along x, its centre at its joint, 1 kg, I = 0.03 kg m²; in water.
  const std::vector<pm::PlanarLink> plate = {{"plate", -1, {0, 0}, 0.0, 1.0, {0.0, 0.0}, 0.03, 0}};
  const pm::PlanarTree<3> tree(pm::Base::floating, plate, 0.0);
  const std::vector<pm::SlenderPart> parts = {{0, {-0.25, 0.0}, 0.5, 0.2}};
  pm::Fluid water;
  water.density = 1000, water.normal_drag = 1.2, water.friction_drag = 0.02;
  const double m = 1.0, ma = 1000 * std::numbers::pi * 0.2 * 0.2 / 4 * 0.5, kn = 0.5 * 1000 * 1.2 * 0.2 * 0.5, kt = 0.5 * 1000 * 0.02 * 2 * 0.2 * 0.5;
  const auto L = [&](const pm::Fluid& f) { return pm::TreeInFluid<3>{&tree, &parts, &f}; };
  {
    const pm::Fluid still = [&] { auto f = water; f.normal_drag = 0, f.friction_drag = 0; return f; }();
    const auto a = pm::lagrangian_acceleration<3>(L(still), 0.0, {0, 0, 0}, {0, 0, 0}, {0, 2.0, 0});
    s.check("added mass: pushed across itself from rest, it accelerates as F/(m + ρπh²ℓ/4) (relative)", a[1] / (2.0 / (m + ma)) - 1, 1e-12);
    const auto b = pm::lagrangian_acceleration<3>(L(still), 0.0, {0, 0, 0}, {0, 0, 0}, {2.0, 0, 0});
    s.check("added mass: along itself, as F/m: no added mass (relative)", b[0] / 2.0 - 1, 1e-12);
  }
  const Drag drag{&tree, &parts, &water};
  const auto Q = [&](double t, const std::array<double, 3>& q, const std::array<double, 3>& v) { return drag(t, q, v); };
  {
    double t = 0;
    std::array<double, 3> q{0, 0, 0}, v{0, 0.5, 0};
    for (int i = 0; i < 1000; i++) pm::runge_kutta4_lagrangian<3>(L(water), Q, t, q, v, 1e-3);
    const double exact = 0.5 / (1 + kn * 0.5 * 1.0 / (m + ma));
    s.check("drag: moving across itself it slows as v₀/(1 + k v₀ t/(m + m_a)) (after 1 s, relative)", v[1] / exact - 1, 1e-6);
  }
  {
    double t = 0;
    std::array<double, 3> q{0, 0, 0}, v{0.5, 0, 0};
    for (int i = 0; i < 1000; i++) pm::runge_kutta4_lagrangian<3>(L(water), Q, t, q, v, 1e-3);
    const double exact = 0.5 / (1 + kt * 0.5 * 1.0 / m);
    s.check("drag: along itself, as v₀/(1 + k_t v₀ t/m), its skin friction alone (after 1 s, relative)", v[0] / exact - 1, 1e-6);
  }
  {
    pm::Fluid current = water;
    current.current = {0.0, 0.3};
    const Drag dc{&tree, &parts, &current};
    const auto Qc = [&](double t, const std::array<double, 3>& q, const std::array<double, 3>& v) { return dc(t, q, v); };
    double t = 0;
    std::array<double, 3> q{0, 0, 0}, v{0, 0, 0};
    const auto a0 = pm::lagrangian_acceleration<3>(L(current), 0.0, q, v, Qc(0, q, v));
    s.check("current: at rest in a current across it, it is pushed by k u²/(m + m_a) (relative)", a0[1] / (kn * 0.09 / (m + ma)) - 1, 1e-6);
    for (int i = 0; i < 20000; i++) pm::runge_kutta4_lagrangian<3>(L(current), Qc, t, q, v, 1e-3);
    // (m + m_a) dw/dt = −k|w|w for w = v − u: w(t) = w₀/(1 + k|w₀|t/(m + m_a)), hyperbolically to the current's speed.
    const double exact = 0.3 - 0.3 / (1 + kn * 0.3 * 20.0 / (m + ma));
    s.check("current: and carried toward the current's speed as quadratic drag allows, u − u/(1 + k u t/(m + m_a)) (after 20 s, relative)", v[1] / exact - 1, 1e-6);
  }
  {
    pm::Fluid still = water;
    still.normal_drag = 0, still.friction_drag = 0;
    const auto none = [](double, const std::array<double, 3>&, const std::array<double, 3>&) { return std::array<double, 3>{0, 0, 0}; };
    const auto energy = [&](const std::array<double, 3>& q, const std::array<double, 3>& v) { return L(still)(0.0, q, v); };  // no potential: L = T + T_f
    double t = 0;
    std::array<double, 3> q{0, 0, 0}, v{0.3, 0.2, 2.0};
    const double e0 = energy(q, v);
    for (int i = 0; i < 2000; i++) pm::runge_kutta4_lagrangian<3>(L(still), none, t, q, v, 1e-3);
    s.check("energy: turning and moving without drag, the body and the fluid keep their energy together (relative, 2 s)", energy(q, v) / e0 - 1, 1e-9);
  }
  {
    // The variational integrator: the tree and the fluid's kinetic energy in its Lagrangian, the drag a force.
    using Lag = pm::TreeInFluid<3>;
    pm::VariationalIntegrator<3, Lag, Drag> vi(pm::DiscreteLagrangian<3, Lag>(L(water), 1e-3), drag, 1e-12, 50);
    auto vs = vi.start(0.0, {0, 0, 0}, {0.2, 0.5, 0.5});
    for (int i = 0; i < 1000; i++) vi.step(vs);
    double t = 0;
    std::array<double, 3> q{0, 0, 0}, v{0.2, 0.5, 0.5};
    for (int i = 0; i < 1000; i++) pm::runge_kutta4_lagrangian<3>(L(water), Q, t, q, v, 1e-3);
    double worst = 0;
    for (int k = 0; k < 3; k++) worst = std::max(worst, std::fabs(vs.q[k] - q[k]));
    s.check("variational: the discrete action of the plate and the fluid, its drag by Lagrange–d'Alembert, follows Runge–Kutta within its order-2 error (1 s at 1 ms, m and rad)", worst, 1e-5);
  }
  return s.done();
}
