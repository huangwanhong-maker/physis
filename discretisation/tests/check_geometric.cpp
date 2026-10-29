// Checks of cellulations in space (geometric.hpp) and of cochains on boxes as fields (cochain_field.hpp).
//
// - A box of uneven cuts: every cell's measure is the product of its extents, its centroid their
//   midpoints; its counts are the products of the cuts'.
// - d_{k+1} d_k = 0 as integer matrices, in every degree, on boxes of three, two and one dimensions.
// - Geometry against topology: every incidence of d_{n−1} is o_i σ(i, f), the orientation of the top cell
//   times the side of the facet's conormal it lies on.
// - A box's facets: two cells on each interior facet, one on each of the 2(n₀n₁ + n₁n₂ + n₂n₀) boundary
//   facets.
// - Points located in their cells.
// - Cochain fields: an affine 0-cochain is read exactly between its vertices, a top cochain as a
//   density; other degrees refuse a point; pairing with chains and d.
// - Transfer to a refinement (uniform, and by uneven new cuts): T d = d T in every degree, the total of
//   a top cochain kept, an affine 0-cochain reproduced (formalism Proposition 6.2).

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/discretisation/cochain_field.hpp>
#include <physis/discretisation/geometric.hpp>

using calx::verification::Suite;
namespace dz = physis::discretisation;
using dz::Vec3;

namespace {

// The largest entry of d_{k+1} d_k, as integers.
double dd(const dz::CWComplex& X, std::size_t k) {
  const auto a = dz::coboundary(X, k), b = dz::coboundary(X, k + 1);
  double worst = 0;
  for (std::size_t i = 0; i < b.pattern.rows; i++) {
    std::vector<double> row(a.pattern.cols, 0.0);
    for (std::uint64_t e = b.pattern.begin(i); e < b.pattern.end(i); e++) {
      const std::size_t j = b.pattern.col[e];
      for (std::uint64_t g = a.pattern.begin(j); g < a.pattern.end(j); g++) row[a.pattern.col[g]] += b.value[e] * a.value[g];
    }
    for (double x : row) worst = std::max(worst, std::fabs(x));
  }
  return worst;
}

// Every incidence of d_{n−1} against the geometry: d(i, f) = o_i σ(i, f).
bool orientations_agree(const dz::GeometricComplex& K) {
  const auto M = dz::measure(K);
  const auto F = dz::facets(K, M);
  const auto d = dz::coboundary(K.topology, M.dimension - 1);
  for (std::size_t i = 0; i < d.pattern.rows; i++)
    for (std::uint64_t e = d.pattern.begin(i); e < d.pattern.end(i); e++) {
      const std::size_t f = d.pattern.col[e];
      const int sigma = F.cells[f][0] == i ? 1 : -1;
      if (d.value[e] != M.orientation[i] * sigma) return false;
    }
  return true;
}

std::vector<double> random_cochain(calx::random::Stream& rng, std::size_t n) {
  std::vector<double> c(n);
  for (double& x : c) x = rng.uniform(-1, 1);
  return c;
}

}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261029, 3, 0, 0);
  const std::array<std::vector<double>, 3> cuts{{{0, 0.3, 1.0, 1.2, 2.0}, {0, 0.5, 1.5}, {-1, 0.25, 1}}};
  const dz::BoxComplex box = dz::box_complex(cuts);
  const auto M = dz::measure(box.complex);

  // ---- measures --------------------------------------------------------------------------------------------------
  {
    const auto& X = box.complex.topology;
    s.require("box: 45 vertices, 4·3·3 + 5·2·3 + 5·3·2 edges, and so on to 16 cubes; dimension 3",
              X.vertices == 45 && X.edges.size() == 36 + 30 + 30 && X.faces.size() == 4 * 2 * 3 + 4 * 3 * 2 + 5 * 2 * 2 && X.higher[0].size() == 16 && box.dimension() == 3);
    double worst = 0;
    for (std::size_t k = 1; k <= 3; k++) {
      const auto m = M.of(k);
      for (std::size_t i = 0; i < box.cells[k].size(); i++) {
        const auto& c = box.cells[k][i];
        double product = 1;
        Vec3 mid{};
        for (int a = 0; a < 3; a++) {
          const auto& x = cuts[static_cast<std::size_t>(a)];
          const int j = c.corner[static_cast<std::size_t>(a)];
          const bool spans = c.axes & (1u << a);
          if (spans) product *= x[static_cast<std::size_t>(j) + 1] - x[static_cast<std::size_t>(j)];
          mid[static_cast<std::size_t>(a)] = spans ? 0.5 * (x[static_cast<std::size_t>(j)] + x[static_cast<std::size_t>(j) + 1]) : x[static_cast<std::size_t>(j)];
        }
        worst = std::max({worst, std::fabs(m[i] - product), physis::geometry::norm({M.centroid[k][i][0] - mid[0], M.centroid[k][i][1] - mid[1], M.centroid[k][i][2] - mid[2]})});
      }
    }
    s.check("box: every cell's measure is the product of its extents and its centroid their midpoints (edges, faces, cubes)", worst, 1e-15);
  }

  // ---- d d = 0, and geometry against topology ------------------------------------------------------------------------
  {
    const dz::BoxComplex flat = dz::box_complex({{{0, 1, 2.5, 3}, {0, 2, 3}, {0.5}}});
    const dz::BoxComplex line = dz::box_complex({{{0, 1, 3, 4, 7}, {2}, {1}}});
    s.check("d: d₁d₀ = 0 and d₂d₁ = 0 on the box of three dimensions, as integer matrices", std::max(dd(box.complex.topology, 0), dd(box.complex.topology, 1)), 0);
    s.check("d: d₁d₀ = 0 on a box of two dimensions", dd(flat.complex.topology, 0), 0);
    s.require("d: a box of one dimension has d₀ only", line.dimension() == 1 && flat.dimension() == 2);
    s.require("orientation: every incidence of d is the cell's orientation times the side of the facet's conormal (three, two and one dimensions)",
              orientations_agree(box.complex) && orientations_agree(flat.complex) && orientations_agree(line.complex));
    const auto F = dz::facets(box.complex, M);
    std::size_t boundary = 0, interior = 0;
    for (std::size_t f = 0; f < F.cells.size(); f++) (F.boundary(f) ? boundary : interior)++;
    s.require("facets: 2(4·2 + 2·2 + 2·4) = 40 on the boundary, each with one cell; the other 3·2·2 + 4·1·2 + 4·2·1 = 28 with two",
              boundary == 40 && interior == 28 && F.cells.size() == 68);
  }

  // ---- location ------------------------------------------------------------------------------------------------------
  {
    bool inside = true;
    for (int i = 0; i < 200; i++) {
      const Vec3 x{rng.uniform(0, 2), rng.uniform(0, 1.5), rng.uniform(-1, 1)};
      const auto L = box.locate(x);
      if (!L || L->cell == dz::no_cell) { inside = false; break; }
      const auto& c = box.cells[3][L->cell];
      for (int a = 0; a < 3; a++) {
        const auto& cut = cuts[static_cast<std::size_t>(a)];
        const int j = c.corner[static_cast<std::size_t>(a)];
        if (!(x[static_cast<std::size_t>(a)] >= cut[static_cast<std::size_t>(j)] && x[static_cast<std::size_t>(a)] <= cut[static_cast<std::size_t>(j) + 1]) ||
            L->local[static_cast<std::size_t>(a)] < 0 || L->local[static_cast<std::size_t>(a)] > 1)
          inside = false;
      }
    }
    s.require("locate: 200 points each in the cube that holds them, at coordinates in [0, 1]; a point outside none",
              inside && !box.locate({2.5, 0, 0}) && box.locate({2.0, 1.5, 1.0}).has_value());
  }

  // ---- cochain fields ------------------------------------------------------------------------------------------------
  const auto shared = std::make_shared<const dz::BoxComplex>(box);
  {
    auto affine = [](const Vec3& x) { return 1 + 2 * x[0] - x[1] + 0.5 * x[2]; };
    std::vector<double> phi(box.cells[0].size());
    for (std::size_t v = 0; v < phi.size(); v++) phi[v] = affine(box.complex.position[v]);
    physis::fields::Descriptor d;
    d.name = "phi";
    const dz::CochainField f0(d, shared, 0, phi);
    double worst = 0;
    for (int i = 0; i < 100; i++) {
      const Vec3 x{rng.uniform(0, 2), rng.uniform(0, 1.5), rng.uniform(-1, 1)};
      worst = std::max(worst, std::fabs(f0.at({0, x}).value - affine(x)));
    }
    s.check("cochain field: an affine 0-cochain read between its vertices, multilinearly, exactly (100 points)", worst, 1e-14);
    std::vector<double> mass(box.cells[3].size());
    const auto vol = M.of(3);
    for (std::size_t i = 0; i < mass.size(); i++) mass[i] = 2.5 * vol[i];
    const dz::CochainField f3(d, shared, 3, mass);
    const auto rho = f3.at({0, {1.1, 0.7, 0.1}});
    const dz::CochainField f1(d, shared, 1, random_cochain(rng, box.cells[1].size()));
    s.require("cochain field: a top cochain read as its density, with its note; an edge cochain refuses a point; outside the box, outside",
              std::fabs(rho.value - 2.5) < 1e-14 && !rho.note.empty() && f1.at({0, {1, 1, 0}}).quality == physis::fields::Quality::invalid &&
                  f0.at({0, {3, 0, 0}}).quality == physis::fields::Quality::outside);
    const auto d0 = f0.coboundary();
    const double along = d0.pair({{0, 1.0}});
    const auto& e0 = box.complex.topology.edges[0];
    s.check("cochain field: d of the 0-cochain paired with an edge is φ(head) − φ(tail)", along - (phi[e0[1]] - phi[e0[0]]), 1e-15);
    s.require("cochain field: pairing, d and transfer among its capabilities; it is a 1-form after d",
              physis::fields::described(f0).capabilities.has(physis::fields::Capability::chain) && physis::fields::described(f0).capabilities.has(physis::fields::Capability::coboundary) &&
                  d0.degree() == 1 && d0.descriptor().fibre.degree == 1);
  }

  // ---- transfer ------------------------------------------------------------------------------------------------------
  {
    const dz::BoxComplex fine2 = dz::refine(box, 2);
    std::array<std::vector<double>, 3> uneven = cuts;
    uneven[0].insert(uneven[0].begin() + 1, 0.1);
    uneven[0].insert(uneven[0].begin() + 4, 1.13);
    uneven[2].insert(uneven[2].begin() + 1, -0.6);
    const dz::BoxComplex fine3 = dz::box_complex(uneven);
    for (const dz::BoxComplex* fine : {&fine2, &fine3}) {
      const std::string which = fine == &fine2 ? "halved" : "cut unevenly";
      double commute = 0;
      for (std::size_t k = 0; k < 3; k++) {
        const auto c = random_cochain(rng, box.cells[k].size());
        const auto Tdc = dz::transfer(box, *fine, k + 1, dz::coboundary(box.complex.topology, k).apply(c));
        const auto dTc = dz::coboundary(fine->complex.topology, k).apply(dz::transfer(box, *fine, k, c));
        for (std::size_t i = 0; i < Tdc.size(); i++) commute = std::max(commute, std::fabs(Tdc[i] - dTc[i]));
      }
      s.check("transfer (" + which + "): T d = d T in degrees 0, 1 and 2", commute, 1e-14);
      const auto u = random_cochain(rng, box.cells[3].size());
      const auto Tu = dz::transfer(box, *fine, 3, u);
      double a = 0, b = 0;
      for (double x : u) a += x;
      for (double x : Tu) b += x;
      s.check("transfer (" + which + "): a top cochain keeps its total", b - a, 1e-13);
      std::vector<double> phi(box.cells[0].size());
      for (std::size_t v = 0; v < phi.size(); v++) phi[v] = 3 - box.complex.position[v][0] + 2 * box.complex.position[v][2];
      const auto Tphi = dz::transfer(box, *fine, 0, phi);
      double worst = 0;
      for (std::size_t v = 0; v < Tphi.size(); v++) worst = std::max(worst, std::fabs(Tphi[v] - (3 - fine->complex.position[v][0] + 2 * fine->complex.position[v][2])));
      s.check("transfer (" + which + "): an affine 0-cochain reproduced at the fine vertices", worst, 1e-14);
    }
  }
  return s.done();
}
