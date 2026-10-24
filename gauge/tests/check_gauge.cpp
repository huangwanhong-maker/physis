// Checks of physis/gauge: gauge invariance of the action and of Wilson loops, on lattices and
// on one-vertex complexes, for ℤ_3, S_3, U(1), SU(2) and SU(3); staples against the local
// action; Monte Carlo (heat bath with overrelaxation, Metropolis where the action is not
// linear, Cabibbo–Marinari for SU(3)) against Migdal and Witten's exact partition function
// on closed surfaces and against the area law on a disc; exact enumeration of ℤ_N theories
// against the character formula; the integer U(1) topological charge; and flat connections
// counted through the fundamental group and by enumeration.

#include <physis/gauge/exact.hpp>
#include <physis/gauge/field.hpp>
#include <physis/gauge/observables.hpp>
#include <physis/gauge/run.hpp>
#include <physis/gauge/updates.hpp>

#include <calx/verification/check.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/dual.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/surfaces.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <string>
#include <vector>

namespace pg = physis::gauge;
namespace mt = morphe::topology;
namespace ma = morphe::algebra;
namespace ft = physis::field_theory;

namespace {

// The closed path through the given vertices (each consecutive pair joined by an edge).
mt::Path path_through(const mt::CWComplex& X, const std::vector<std::size_t>& vertices) {
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge;
  for (std::size_t e = 0; e < X.edges.size(); e++) edge[{X.edges[e][0], X.edges[e][1]}] = e;
  mt::Path p;
  for (std::size_t i = 0; i < vertices.size(); i++) {
    const std::size_t u = vertices[i], v = vertices[(i + 1) % vertices.size()];
    if (edge.count({u, v})) p.push_back({edge.at({u, v}), 1});
    else p.push_back({edge.at({v, u}), -1});
  }
  return p;
}

// The boundary of the w × h rectangle at the origin of a product of intervals with `side` + 1
// vertices along y (vertex (x, y) is x (side + 1) + y).
mt::Path rectangle(const mt::CWComplex& X, std::size_t side, std::size_t w, std::size_t h) {
  std::vector<std::size_t> v;
  auto at = [&](std::size_t x, std::size_t y) { return x * (side + 1) + y; };
  for (std::size_t x = 0; x < w; x++) v.push_back(at(x, 0));
  for (std::size_t y = 0; y < h; y++) v.push_back(at(w, y));
  for (std::size_t x = w; x > 0; x--) v.push_back(at(x, h));
  for (std::size_t y = h; y > 0; y--) v.push_back(at(0, y));
  return path_through(X, v);
}

template <pg::GaugeGroup G>
double invariance_defect(const mt::CWComplex& X, G group, calx::random::Stream& z) {
  pg::GaugeField<G> U(group, X);
  U.randomise(z);
  const auto loops = mt::fundamental_group(X).loops;
  std::vector<double> before{U.action(1.7)};
  for (std::size_t f = 0; f < X.faces.size(); f++) before.push_back(U.plaquette(f));
  for (const auto& l : loops) before.push_back(U.wilson_loop(l));
  std::vector<typename G::Element> g(X.vertices);
  for (auto& x : g) x = U.group().random(z);
  U.transform(g);
  std::vector<double> after{U.action(1.7)};
  for (std::size_t f = 0; f < X.faces.size(); f++) after.push_back(U.plaquette(f));
  for (const auto& l : loops) after.push_back(U.wilson_loop(l));
  double d = 0;
  for (std::size_t i = 0; i < before.size(); i++) d = std::max(d, std::fabs(before[i] - after[i]));
  return d;
}

template <pg::StapleGroup G>
double staple_defect(const mt::CWComplex& X, G group, calx::random::Stream& z) {
  pg::GaugeField<G> U(group, X);
  U.randomise(z);
  double d = 0;
  for (std::size_t e = 0; e < X.edges.size(); e++)
    if (U.linear(e)) d = std::max(d, std::fabs(U.group().re_trace_with(U.link(e), U.staple(e)) - U.local_plaquettes(e)));
  return d;
}

// The mean plaquette by Monte Carlo, in errors from the exact value.
double deviation(const pg::RunResult& r, double exact) { return std::fabs(r.plaquette.mean - exact) / r.plaquette.error; }

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(21, 1, 0, 0);
  const auto lattice = mt::torus_lattice({3, 4});
  const auto genus2 = mt::orientable_surface(2);
  const auto dual = mt::dual_two_complex(mt::simplex_boundary(4)).complex;
  const auto S3 = ma::symmetric(3);

  // ---- gauge invariance and staples ----------------------------------------------------------------------------------
  {
    double worst = 0;
    for (const auto* X : {&lattice, &genus2, &dual}) {
      worst = std::max(worst, invariance_defect(*X, pg::ZN{3}, z));
      worst = std::max(worst, invariance_defect(*X, pg::FiniteGroupAdapter{&S3, {}}, z));
      worst = std::max(worst, invariance_defect(*X, pg::U1Group{}, z));
      worst = std::max(worst, invariance_defect(*X, pg::SU2Group{}, z));
      worst = std::max(worst, invariance_defect(*X, pg::SUNGroup{3}, z));
    }
    s.check("gauge: action, plaquettes and Wilson loops of π1's generators are gauge invariant (Z3, S3, U(1), SU(2), SU(3); lattice, genus 2, dual of ∂Δ⁴)", worst, 1e-12);
    double st = 0;
    for (const auto* X : {&lattice, &dual}) {
      st = std::max(st, staple_defect(*X, pg::U1Group{}, z));
      st = std::max(st, staple_defect(*X, pg::SU2Group{}, z));
      st = std::max(st, staple_defect(*X, pg::SUNGroup{3}, z));
    }
    s.check("gauge: Re tr(U_e Σ_e)/N equals the local plaquette sum on every linear edge", st, 1e-12);
    pg::GaugeField<pg::U1Group> T(pg::U1Group{}, mt::polygon_complex("abAB"));
    s.require("gauge: the edges of the one-vertex torus a b a⁻¹ b⁻¹ are not linear", !T.linear(0) && !T.linear(1));
  }
  // ---- the heat-bath samplers ------------------------------------------------------------------------------------------------
  {
    // ⟨x0⟩ under √(1 - x0²) e^{a x0} is I_2(a)/I_1(a); ⟨cos φ⟩ under e^{κ cos φ} is I_1(κ)/I_0(κ);
    // both ratios here by quadrature, on both of the SU(2) sampler's branches (a < 2, a >= 2).
    double worst = 0;
    for (double a : {0.3, 1.5, 2.5, 8.0}) {
      auto w = [&](double x) { return std::sqrt(1 - x * x) * std::exp(a * x); };
      const double want = calx::quadrature::integrate([&](double x) { return x * w(x); }, -1, 1, 1e-15, 1e-13).value /
                          calx::quadrature::integrate(w, -1, 1, 1e-15, 1e-13).value;
      std::vector<double> x(200000);
      for (auto& v : x) v = pg::detail::su2_x0(a, z);
      const double m = ft::mean(x), e = std::sqrt(ft::variance(x) / static_cast<double>(x.size()));
      worst = std::max(worst, std::fabs(m - want) / e);
    }
    for (double kappa : {0.2, 2.0, 10.0}) {
      const auto X = pg::expansion_u1(kappa, 1);
      const double want = X.reps[2].c / X.reps[1].c;  // I_1 / I_0 (reps: n = -1, 0, 1)
      std::vector<double> c(200000);
      for (auto& v : c) v = std::cos(pg::detail::von_mises(kappa, z));
      const double m = ft::mean(c), e = std::sqrt(ft::variance(c) / static_cast<double>(c.size()));
      worst = std::max(worst, std::fabs(m - want) / e);
    }
    s.check("gauge: the SU(2) heat-bath x0 (Creutz below a = 2, Kennedy–Pendleton above) and von Mises samples have their exact means (errors)", worst, 4.0, "σ");
  }
  // ---- Monte Carlo against Migdal–Witten ----------------------------------------------------------------------------------
  {
    const auto T44 = mt::torus_lattice({4, 4});
    pg::RunOptions o;
    o.thermalise = 300, o.sweeps = 6000, o.seed = 3;
    o.beta = 1.5;
    auto r = pg::run(T44, {"U1", 0}, o);
    s.check("gauge: U(1) on the 4 × 4 torus at β = 1.5, heat bath and overrelaxation, against Migdal–Witten (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    s.require("gauge: ... its topological charge is an integer", r.charge.has_value());
    o.beta = 2.5;
    r = pg::run(T44, {"SU2", 0}, o);
    s.check("gauge: SU(2) on the 4 × 4 torus at β = 2.5, heat bath and overrelaxation, against Migdal–Witten (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    o.beta = 0.7;
    r = pg::run(T44, {"SU2", 0}, o);
    s.check("gauge: SU(2) on the 4 × 4 torus at β = 0.7 (strong coupling), against Migdal–Witten (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    o.update = "metropolis", o.step = 1.2, o.sweeps = 20000;
    o.beta = 2.0;
    r = pg::run(mt::polygon_complex("abAB"), {"SU2", 0}, o);
    s.check("gauge: SU(2) on the one-vertex torus (non-linear links) at β = 2, Metropolis, against Migdal–Witten (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    r = pg::run(genus2, {"SU2", 0}, o);
    s.check("gauge: SU(2) on the one-vertex surface of genus 2 (χ = -2) at β = 2, Metropolis, against Migdal–Witten (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    o.sweeps = 6000;
    o.beta = 1.0;
    r = pg::run(mt::torus_lattice({3, 3}), {"Z", 3}, o);
    s.check("gauge: Z3 on the 3 × 3 torus at β = 1, Metropolis, against the character formula (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    s.require("gauge: ... with an acceptance between 0.1 and 0.9", r.acceptance > 0.1 && r.acceptance < 0.9, std::to_string(r.acceptance));
    o.update = "heatbath";
    r = pg::run(mt::torus_lattice({3, 3}), {"Z", 3}, o);
    s.check("gauge: Z3 on the 3 × 3 torus at β = 1, the heat bath of a finite group, against the character formula (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
    o.beta = 0.8;
    r = pg::run(mt::torus_lattice({4, 4}), {"Z", 2}, o);
    s.check("gauge: Z2 on the 4 × 4 torus at β = 0.8, heat bath, against Migdal–Witten (errors)", deviation(r, *r.exact_plaquette), 4.0, "σ");
  }
  // ---- discs: plaquettes and the area law ------------------------------------------------------------------------------------
  {
    const auto D = mt::product(mt::interval(4), mt::interval(4));
    pg::RunOptions o;
    o.thermalise = 200, o.sweeps = 8000, o.seed = 5, o.beta = 2.0;
    o.loops = {rectangle(D, 4, 1, 1), rectangle(D, 4, 2, 2), rectangle(D, 4, 1, 3)};
    const auto u1 = pg::expansion_u1(2.0);
    auto r = pg::run(D, {"U1", 0}, o);
    double worst = deviation(r, pg::disc(u1, 16).mean_plaquette);
    const std::size_t areas[] = {1, 4, 3};
    for (std::size_t l = 0; l < 3; l++) worst = std::max(worst, std::fabs(r.loops[l].mean - pg::disc_wilson_loop(u1, areas[l])) / r.loops[l].error);
    s.check("gauge: U(1) on a 4 × 4 disc at β = 2: the plaquette and Wilson loops of areas 1, 4, 3 against (I1/I0)^A (errors)", worst, 4.0, "σ");
    o.beta = 3.0;
    const auto su2 = pg::expansion_su2(3.0);
    r = pg::run(D, {"SU2", 0}, o);
    worst = deviation(r, pg::disc(su2, 16).mean_plaquette);
    for (std::size_t l = 0; l < 3; l++) worst = std::max(worst, std::fabs(r.loops[l].mean - pg::disc_wilson_loop(su2, areas[l])) / r.loops[l].error);
    s.check("gauge: SU(2) on the disc at β = 3: the plaquette and Wilson loops against (I2/I1)^A (errors)", worst, 4.0, "σ");
    o.loops.clear();
    o.sweeps = 3000;
    o.beta = 4.0;
    r = pg::run(mt::product(mt::interval(3), mt::interval(3)), {"SU", 3}, o);
    s.check("gauge: SU(3) on a 3 × 3 disc at β = 4, Cabibbo–Marinari, against the Weyl integral of one plaquette (errors)",
            deviation(r, pg::free_plaquette_sun(3, 4.0)), 4.0, "σ");
    s.check("gauge: ... and run() found that exact value itself", *r.exact_plaquette - pg::free_plaquette_sun(3, 4.0), 1e-12);
  }
  // ---- exact enumeration of ℤ_N against the character formula ----------------------------------------------------------------
  {
    double worst = 0;
    for (const mt::CWComplex* X : std::initializer_list<const mt::CWComplex*>{&genus2, nullptr}) {
      const mt::CWComplex Y = X ? *X : mt::torus_lattice({2, 2});
      for (std::uint32_t N : {2u, 3u}) {
        const double beta = 0.8;
        pg::GaugeField<pg::ZN> U(pg::ZN{N}, Y);
        const std::size_t E = Y.edges.size();
        std::vector<std::uint32_t> k(E, 0);
        double Z = 0, P = 0;
        for (;;) {
          for (std::size_t e = 0; e < E; e++) U.set(e, k[e]);
          const double w = std::exp(-U.action(beta) + beta * static_cast<double>(Y.faces.size()));
          Z += w;
          P += w * U.mean_plaquette();
          std::size_t e = 0;
          while (e < E && ++k[e] == N) k[e++] = 0;
          if (e == E) break;
        }
        const auto mw = pg::migdal_witten(pg::expansion_zn(N, beta), Y.euler_characteristic(), Y.faces.size());
        worst = std::max(worst, std::fabs(std::log(Z / std::pow(N, E)) - mw.log_z));
        worst = std::max(worst, std::fabs(P / Z - mw.mean_plaquette));
      }
    }
    s.check("gauge: Z2 and Z3 theories on the genus-2 surface and the 2 × 2 torus, enumerated, against Σ_n c_n^F (log Z and plaquette)", worst, 1e-12);
  }
  // ---- the topological charge ------------------------------------------------------------------------------------------------
  {
    const std::size_t L = 5;
    const auto T = mt::torus_lattice({L, L});
    const auto S = mt::classify_surface(T);
    pg::GaugeField<pg::U1Group> U(pg::U1Group{}, T);
    double worst = 0;
    for (int t = 0; t < 50; t++) {
      U.randomise(z);
      const double q = pg::topological_charge(U, S.orientation);
      worst = std::max(worst, std::fabs(q - std::round(q)));
    }
    s.check("gauge: the U(1) topological charge of 50 random fields on the 5 × 5 torus is an integer", worst, 1e-9);
    // A field of uniform flux 2π / L² per plaquette has |Q| = 1 (vertex (x, y) is x L + y).
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> edge;
    for (std::size_t e = 0; e < T.edges.size(); e++) edge[{T.edges[e][0], T.edges[e][1]}] = e;
    const double a = 2 * std::numbers::pi / static_cast<double>(L * L);
    for (std::size_t x = 0; x < L; x++)
      for (std::size_t y = 0; y < L; y++) {
        // The y-link at (x, y) carries a x; the x-link from (L-1, y) to (0, y) carries -a L y.
        U.set(edge.at({x * L + y, x * L + (y + 1) % L}), physis::symmetry::U1::from_angle(a * static_cast<double>(x)));
        U.set(edge.at({x * L + y, ((x + 1) % L) * L + y}), physis::symmetry::U1::from_angle(x == L - 1 ? -a * static_cast<double>(L * y) : 0.0));
      }
    const double q = pg::topological_charge(U, S.orientation);
    s.check("gauge: a uniform flux of 2π/L² per plaquette on the 5 × 5 torus has |Q| = 1", std::fabs(q) - 1, 1e-12);
  }
  // ---- flat connections ------------------------------------------------------------------------------------------------------
  {
    const auto T22 = mt::torus_lattice({2, 2});
    const auto RP2 = mt::polygon_complex("aa");
    const ma::FiniteGroup Z2 = ma::cyclic(2);
    struct Case { const mt::CWComplex* X; const ma::FiniteGroup* G; long long want; };
    const Case cases[] = {{&genus2, &S3, 486}, {&T22, &S3, 3888}, {&dual, &Z2, 16}, {&RP2, &S3, 4}};
    bool ok = true;
    std::string detail;
    for (const auto& c : cases) {
      const auto by_pi1 = pg::count_flat_connections(*c.X, *c.G);
      const auto by_enum = pg::count_flat_connections_by_enumeration(*c.X, *c.G);
      detail += by_pi1.str() + " ";
      ok = ok && by_pi1 == morphe::arithmetic::Integer(c.want) && by_enum == static_cast<std::uint64_t>(c.want);
    }
    s.require("gauge: flat S3, Z2 connections on the genus-2 surface (486), the 2 × 2 torus (3888), the dual of ∂Δ⁴ (16) and ℝP² (4), through π1 and by enumeration",
              ok, detail);
  }
  return s.done();
}
