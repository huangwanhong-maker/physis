// Checks of the H3 grid (h3.hpp) against the reference library's own answers (tests/fixtures/h3.txt, made by
// tests/make_h3_fixtures.py from the `h3` package on the reference C library): for 800 points uniform on the sphere at random
// resolutions, seven places at every resolution (the poles, the antimeridian, three cities) and the twelve pentagons of five
// resolutions, each point's cell, the cell's centre, its boundary and its parent; the areas of 40 cells.
// And what holds of the grid itself: every point lies in the cell found for it (a round trip through a cell's centre finds that
// cell again); a cell has six or five vertices, more only where a class III edge crosses a face; twelve pentagons at every
// resolution, and a cell's index round-trips through its text.

#include <physis/computational_geometry/h3.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace H = physis::computational_geometry::h3;

int main() {
  calx::verification::Suite s;
  const double DEG = std::numbers::pi / 180;
  std::ifstream file(std::string(PHYSIS_GEOMETRY_FIXTURES) + "/h3.txt");
  s.require("the fixture is there", static_cast<bool>(file));
  std::string line;
  int n = 0, wrong_cell = 0, wrong_parent = 0, wrong_pent = 0, wrong_base = 0, wrong_count = 0, wrong_text = 0, lost = 0, pentagons = 0, hexes_many = 0;
  double worst_centre = 0, worst_boundary = 0, worst_area = 0;
  std::string first_wrong;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream in(line);
    std::string kind;
    in >> kind;
    if (kind == "A") {
      std::string cell;
      double area;
      in >> cell >> area;
      worst_area = std::max(worst_area, std::fabs(H::cell_area_m2(H::from_string(cell)) - area) / area);
      continue;
    }
    double lat, lon, clat, clon;
    int res, pent, base, m;
    std::string cell, parent;
    in >> lat >> lon >> res >> cell >> parent >> pent >> base >> clat >> clon >> m;
    n++;
    const std::uint64_t c = H::cell_at({lat * DEG, lon * DEG}, res);
    if (H::to_string(c) != cell) { wrong_cell++; if (first_wrong.empty()) first_wrong = line.substr(0, 90) + " found " + H::to_string(c); continue; }
    if (H::from_string(cell) != c) wrong_text++;
    if (H::to_string(H::parent(c, std::max(0, res - 1))) != parent) wrong_parent++;
    if (static_cast<int>(H::is_pentagon(c)) != pent) wrong_pent++;
    pentagons += pent;
    if (H::base_cell(c) != base) wrong_base++;
    const auto ctr = H::centre(c);
    const double dlon = std::remainder(ctr.lng / DEG - clon, 360.0);
    worst_centre = std::max({worst_centre, std::fabs(ctr.lat / DEG - clat), std::fabs(dlon) * std::cos(clat * DEG)});
    // The cell holding its own centre is the cell.
    if (H::cell_at(ctr, res) != c) lost++;
    const auto b = H::boundary(c);
    if (static_cast<int>(b.size()) != m) { wrong_count++; for (int i = 0; i < 2 * m; i++) { double x; in >> x; } continue; }
    if (b.size() > 6) hexes_many++;
    for (int i = 0; i < m; i++) {
      double la, lo;
      in >> la >> lo;
      worst_boundary = std::max({worst_boundary, std::fabs(b[i].lat / DEG - la), std::fabs(std::remainder(b[i].lng / DEG - lo, 360.0)) * std::cos(la * DEG)});
    }
  }
  s.require("every point's cell is the reference's (" + std::to_string(n) + " points at all resolutions)", wrong_cell == 0, first_wrong);
  s.require("each index round-trips through its text", wrong_text == 0);
  s.require("every parent is the reference's", wrong_parent == 0, std::to_string(wrong_parent));
  s.require("pentagons and base cells are the reference's", wrong_pent == 0 && wrong_base == 0, std::to_string(wrong_pent) + " and " + std::to_string(wrong_base));
  s.require("the pentagons of five resolutions are all there, twelve each", pentagons >= 60, std::to_string(pentagons));
  // The reference since 4.5 projects through unit vectors and this port through spherical trigonometry, as it did before: they part
  // by rounding, under 2 mm on the ground.
  s.check("a cell's centre against the reference (degrees)", worst_centre, 2e-8);
  s.require("a cell's centre lies in it, at every resolution", lost == 0, std::to_string(lost));
  s.require("a cell's boundary has the reference's vertices", wrong_count == 0, std::to_string(wrong_count));
  s.require("class III cells crossing a face have the extra vertices", hexes_many > 0, std::to_string(hexes_many));
  s.check("a cell's boundary against the reference (degrees)", worst_boundary, 2e-8);
  s.check("a cell's area against the reference (relative)", worst_area, 1e-6);
  return s.done();
}
