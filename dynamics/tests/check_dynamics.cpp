// Checks of physis/dynamics: the exact propagator of a linear system against its closed
// form (a leaky membrane driven by an exponentially decaying current, Rotter and Diesmann
// 1999, eq. 15 ff.), at equal and unequal time constants; held components; events of a
// linear hybrid system against the closed-form first passage on the grid; and a nonlinear
// gated system integrated by Rush and Larsen's method, exact for frozen coefficients and of
// order one otherwise.

#include <physis/dynamics/hybrid.hpp>
#include <physis/dynamics/linear.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

namespace pd = physis::dynamics;

namespace {

// tau v' = -v + R I + R u,  I' = -I / tau_s: a membrane and its synaptic current.
struct Membrane {
  static constexpr std::size_t n = 2;
  struct Params { double tau = 10e-3, tau_s = 2e-3, R = 40e6; };
  static void system(const Params& p, pd::Matrix& A, pd::State<2>& B, pd::State<2>& c) {
    A(0, 0) = -1 / p.tau; A(0, 1) = p.R / p.tau;
    A(1, 1) = -1 / p.tau_s;
    B = {p.R / p.tau, 0};
    c = {0, 0};
  }
};

// The integrate-and-fire neuron in its simplest form: tau x' = -x + u, x >= 1 fires, x <- 0,
// held for 2 ms.
struct Fire {
  static constexpr std::size_t n = 1;
  struct Params { double tau = 20e-3, t_ref = 2e-3; };
  static void system(const Params& p, pd::Matrix& A, pd::State<1>& B, pd::State<1>& c) { A(0, 0) = -1 / p.tau; B = {1 / p.tau}; c = {0}; }
  static bool guard(const Params&, const pd::State<1>& x) { return x[0] >= 1; }
  static void jump(const Params&, pd::State<1>& x) { x[0] = 0; }
  static double dwell(const Params& p) { return p.t_ref; }
  static constexpr std::array<bool, 1> held{true};
};

// A gate m' = (m_inf(v) - m) / tau_m(v) driven by v' = -v / tau_v: Rush and Larsen's form.
struct Gate {
  static constexpr std::size_t n = 2;
  struct Params { double tau_v = 5e-3; };
  static double inf(double v) { return 1 / (1 + std::exp(-v)); }
  static double tau(double v) { return 1e-3 * (1 + 0.5 * std::tanh(v)); }
  static void derivative(const Params& p, const pd::State<2>& x, double, pd::State<2>& dx) {
    dx[0] = -x[0] / p.tau_v;
    dx[1] = (inf(x[0]) - x[1]) / tau(x[0]);
  }
  static void gates(const Params&, const pd::State<2>& x, pd::State<2>& m_inf, pd::State<2>& m_tau) {
    m_inf = {0, inf(x[0])};
    m_tau = {1, tau(x[0])};
  }
  static constexpr std::array<bool, 2> gate{false, true};
};

}  // namespace

int main() {
  calx::verification::Suite s;
  // ---- the propagator against its closed form -------------------------------------------------------
  for (const double tau_s : {2e-3, 10e-3}) {
    Membrane::Params p;
    p.tau_s = tau_s;
    const double h = 0.1e-3;
    pd::Matrix A(2);
    pd::State<2> B{}, c{};
    Membrane::system(p, A, B, c);
    // The potential in volts, the current in amperes: their typical sizes keep the exponential accurate.
    const auto P = pd::propagator<2, 1>(A, B, c, h, {}, {1e-3, 1e-10}, {1e-10});
    const double e_m = std::exp(-h / p.tau), e_s = std::exp(-h / p.tau_s);
    // v(h) from I(0) = 1 is R tau_s (e_m - e_s) / (tau - tau_s), and R h e_m / tau when tau_s = tau;
    // v(h) from u = 1 is R (1 - e_m).
    const double p21 = tau_s == p.tau ? p.R * h / p.tau * e_m : p.R * p.tau_s / (p.tau - p.tau_s) * (e_m - e_s);
    const double err = std::max({std::fabs(P.P[0] - e_m), std::fabs(P.P[3] - e_s), std::fabs(P.P[2]), calx::verification::relative(P.P[1], p21),
                                 calx::verification::relative(P.Q[0], p.R * (1 - e_m))});
    s.check(std::string("dynamics: the exact propagator of a membrane and its current, tau_s ") + (tau_s == p.tau ? "= tau" : "!= tau"), err, 1e-13);
    auto held = pd::propagator<2, 1>(A, B, c, h, {true, false}, {1e-3, 1e-10}, {1e-10});
    // Without the sizes the norm of the augmented matrix is 4e5 and seventeen squarings cost accuracy.
    const auto raw = pd::propagator<2, 1>(A, B, c, h);
    s.require(std::string("dynamics: scaling by the components' sizes is more accurate than none") + (tau_s == p.tau ? " (equal)" : ""),
              calx::verification::relative(raw.P[1], p21) > calx::verification::relative(P.P[1], p21));
    pd::State<2> x{0.3, 1};
    held.apply(x, 5.0);
    s.require(std::string("dynamics: a held component stays still while the others flow") + (tau_s == p.tau ? " (equal)" : ""),
              x[0] == 0.3 && std::fabs(x[1] - e_s) < 1e-15);
  }
  // ---- events of a linear hybrid system ---------------------------------------------------------------
  {
    Fire::Params p;
    const double h = 0.1e-3, u = 1.5, T = 0.5;
    const auto tr = pd::simulate<Fire>(p, {0}, [&](double) { return u; }, T, h, {}, 0);
    // From 0, x(t) = u (1 - e^{-t/tau}) first reaches 1 at t* = tau ln(u / (u - 1)); on the grid the
    // event is at the first k with x(k h) >= 1, then 2 ms held: the period is that many steps plus 20.
    const double tstar = p.tau * std::log(u / (u - 1));
    long long k = static_cast<long long>(std::ceil(tstar / h - 1e-9));
    while (u * (1 - std::exp(-(k - 1) * h / p.tau)) >= 1) k--;
    while (u * (1 - std::exp(-k * h / p.tau)) < 1) k++;
    const long long period = k + 20;
    bool regular = !tr.events.empty();
    for (std::size_t i = 0; i < tr.events.size(); i++)
      regular &= std::llround(tr.events[i] / h) == static_cast<long long>(k + i * period);
    s.require("dynamics: a hybrid system's events fall on the grid where the closed form first crosses the guard",
              regular && tr.events.size() == static_cast<std::size_t>((std::llround(T / h) - k) / period + 1),
              std::to_string(tr.events.size()) + " events, period " + std::to_string(period) + " steps");
  }
  // ---- Rush and Larsen ---------------------------------------------------------------------------------------
  {
    Gate::Params p;
    // With v held still, m relaxes with constant coefficients and the method is exact.
    pd::State<2> x{0.7, 0.1};
    double sub = 0;
    const pd::Stepping rl{pd::Method::rush_larsen, 1, {}};
    pd::step<Gate>(p, x, 0, 0.5e-3, rl, {true, false}, sub, {1, 1});
    const double want = Gate::inf(0.7) + (0.1 - Gate::inf(0.7)) * std::exp(-0.5e-3 / Gate::tau(0.7));
    s.check("dynamics: Rush and Larsen's step is exact for a gate with frozen coefficients", x[1] - want, 1e-15);
    // With v moving, the method is of order one: errors at h and h/2 against a fine RK4 solution.
    auto run = [&](pd::Method m, int steps) {
      pd::State<2> y{2.0, 0.0};
      double sb = 0;
      const pd::Stepping how{m, 1, {}};
      for (int i = 0; i < steps; i++) pd::step<Gate>(p, y, 0, 10e-3 / steps, how, {false, false}, sb, {1, 1});
      return y[1];
    };
    const double ref = run(pd::Method::rk4, 20000);
    const double e1 = std::fabs(run(pd::Method::rush_larsen, 400) - ref), e2 = std::fabs(run(pd::Method::rush_larsen, 800) - ref);
    s.check("dynamics: Rush and Larsen's method converges with order one, observed order minus 1", calx::verification::order(e1, e2) - 1, 0.1);
  }
  return s.done();
}
