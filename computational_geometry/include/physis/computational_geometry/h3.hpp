// H3: a hierarchical grid of hexagons on the sphere (Brodsky et al., H3, Uber Technologies; https://h3geo.org). Twenty
// faces of an icosahedron, each projected gnomonically; on each an aperture-7 hierarchy of hexagonal grids, alternately
// aligned with the face (class II) and turned by asin(sqrt(3/28)) ≈ 19.1° (class III); twelve pentagons at the
// icosahedron's vertices, at every resolution. A cell's index is 64 bits: a mode, a resolution, one of 122 base cells, and
// a digit in 0..6 for each finer resolution. It serves a global grid on which any quantity can be binned and compared.
//
// This is a port of the reference implementation's core (uber/h3, Apache-2.0; the tables in h3_tables.hpp are copied from
// it by tools/make_h3_tables.py). It holds what a value at a place needs: a point to its cell, a cell to its centre and
// its boundary, its parent, and the walk through its index's digits. Functions keep the reference's order of steps, so
// that they can be read beside it. It is held to the reference library's own answers (tests/check_h3.cpp, tests/fixtures/
// h3.json, made by tests/make_h3_fixtures.py).
//
// Geometry. The grid lives on a sphere; a geodetic latitude is taken as a latitude on it, as the reference does, so a cell's
// place on WGS 84 follows from reading its latitude and longitude as geodetic. Latitudes and longitudes are in radians here,
// as everywhere in physis.
//
// Caveat. Hexagons of one resolution are not of one area: from the middle of a face to its corner they shrink by about a
// factor of two. Densities must divide by each cell's own area (cell_area_m2), never by a resolution's mean.

#pragma once

#include <physis/computational_geometry/h3_tables.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace physis::computational_geometry::h3 {

struct LatLng { double lat = 0, lng = 0; };   // radians

namespace detail {

using tables::Ijk;

inline constexpr double PI = 3.14159265358979323846, TWO_PI = 6.28318530717958647692528676655900576839433, HALF_PI = 1.5707963267948966;
inline constexpr double EPS = 0.0000000000000001;
inline constexpr double SQRT3_2 = 0.8660254037844386467637231707529361834714, RSIN60 = 1.1547005383792515290182975610039149112953;
inline constexpr double ONE_THIRD = 0.333333333333333333333333333333333333333, ONE_SEVENTH = 0.14285714285714285714285714285714285;
inline constexpr double AP7_ROT = 0.333473172251832115336090755351601070065900389;
inline constexpr double SQRT7 = 2.6457513110645905905016157536392604257102, RSQRT7 = 0.37796447300922722721451653623418006081576;
inline constexpr double RES0_U_GNOMONIC = 0.38196601125010500003, INV_RES0_U_GNOMONIC = 2.61803398874989588842;
inline constexpr int MAX_FACE_COORD = 2;
inline constexpr std::uint64_t INIT = 35184372088831ULL;
inline constexpr int MODE_OFFSET = 59, RES_OFFSET = 52, BC_OFFSET = 45, DIGIT_BITS = 3, MAX_RES = 15;

// The digits, named as the reference names them.
enum Digit { CENTER = 0, K_AXES = 1, J_AXES = 2, JK_AXES = 3, I_AXES = 4, IK_AXES = 5, IJ_AXES = 6, INVALID = 7 };

struct FaceIjk { int face = 0; Ijk coord{0, 0, 0}; };
struct Vec2 { double x = 0, y = 0; };

inline int res_of(std::uint64_t h) { return static_cast<int>((h >> RES_OFFSET) & 15); }
inline int base_of(std::uint64_t h) { return static_cast<int>((h >> BC_OFFSET) & 127); }
inline int digit_shift(int res, int r) { return (res - r) * DIGIT_BITS; }
inline int digit_of(std::uint64_t h, int r) { return static_cast<int>((h >> ((MAX_RES - r) * DIGIT_BITS)) & 7); }
inline void set_digit(std::uint64_t& h, int r, int d) {
  const int s = (MAX_RES - r) * DIGIT_BITS;
  h = (h & ~(std::uint64_t{7} << s)) | (static_cast<std::uint64_t>(d) << s);
}
inline bool class3(int res) { return res % 2 == 1; }
inline bool base_is_pentagon(int b) { return b >= 0 && b < 122 && tables::base_cell_data[b].is_pentagon == 1; }

inline void normalize(Ijk& c) {
  if (c.i < 0) { c.j -= c.i; c.k -= c.i; c.i = 0; }
  if (c.j < 0) { c.i -= c.j; c.k -= c.j; c.j = 0; }
  if (c.k < 0) { c.i -= c.k; c.j -= c.k; c.k = 0; }
  int m = c.i;
  if (c.j < m) m = c.j;
  if (c.k < m) m = c.k;
  if (m > 0) { c.i -= m; c.j -= m; c.k -= m; }
}
inline Ijk add(Ijk a, Ijk b) { return {a.i + b.i, a.j + b.j, a.k + b.k}; }
inline Ijk sub(Ijk a, Ijk b) { return {a.i - b.i, a.j - b.j, a.k - b.k}; }
inline Ijk scale(Ijk a, int f) { return {a.i * f, a.j * f, a.k * f}; }
inline Ijk combine(Ijk c, Ijk iv, Ijk jv, Ijk kv) {
  Ijk r = add(add(scale(iv, c.i), scale(jv, c.j)), scale(kv, c.k));
  normalize(r);
  return r;
}
inline void up_ap7(Ijk& c) {
  const int i = c.i - c.k, j = c.j - c.k;
  c = {static_cast<int>(std::lround((3 * i - j) * ONE_SEVENTH)), static_cast<int>(std::lround((i + 2 * j) * ONE_SEVENTH)), 0};
  normalize(c);
}
inline void up_ap7r(Ijk& c) {
  const int i = c.i - c.k, j = c.j - c.k;
  c = {static_cast<int>(std::lround((2 * i + j) * ONE_SEVENTH)), static_cast<int>(std::lround((3 * j - i) * ONE_SEVENTH)), 0};
  normalize(c);
}
inline void down_ap7(Ijk& c) { c = combine(c, {3, 0, 1}, {1, 3, 0}, {0, 1, 3}); }
inline void down_ap7r(Ijk& c) { c = combine(c, {3, 1, 0}, {0, 3, 1}, {1, 0, 3}); }
inline void down_ap3(Ijk& c) { c = combine(c, {2, 0, 1}, {1, 2, 0}, {0, 1, 2}); }
inline void down_ap3r(Ijk& c) { c = combine(c, {2, 1, 0}, {0, 2, 1}, {1, 0, 2}); }
inline void rotate60ccw(Ijk& c) { c = combine(c, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}); }
inline void rotate60cw(Ijk& c) { c = combine(c, {1, 0, 1}, {1, 1, 0}, {0, 1, 1}); }

inline constexpr Ijk UNIT[7] = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}};
inline int unit_to_digit(Ijk c) {
  normalize(c);
  for (int d = 0; d < 7; d++) if (c.i == UNIT[d].i && c.j == UNIT[d].j && c.k == UNIT[d].k) return d;
  return INVALID;
}
inline void neighbor(Ijk& c, int digit) {
  if (digit > CENTER && digit < INVALID) { c = add(c, UNIT[digit]); normalize(c); }
}
inline int digit_ccw(int d) {
  switch (d) { case K_AXES: return IK_AXES; case IK_AXES: return I_AXES; case I_AXES: return IJ_AXES; case IJ_AXES: return J_AXES; case J_AXES: return JK_AXES; case JK_AXES: return K_AXES; default: return d; }
}
inline int digit_cw(int d) {
  switch (d) { case K_AXES: return JK_AXES; case JK_AXES: return J_AXES; case J_AXES: return IJ_AXES; case IJ_AXES: return I_AXES; case I_AXES: return IK_AXES; case IK_AXES: return K_AXES; default: return d; }
}

inline int leading_digit(std::uint64_t h) {
  for (int r = 1, res = res_of(h); r <= res; r++) if (digit_of(h, r)) return digit_of(h, r);
  return CENTER;
}
inline std::uint64_t rotate_cell_ccw(std::uint64_t h) {
  for (int r = 1, res = res_of(h); r <= res; r++) set_digit(h, r, digit_ccw(digit_of(h, r)));
  return h;
}
inline std::uint64_t rotate_cell_cw(std::uint64_t h) {
  for (int r = 1, res = res_of(h); r <= res; r++) set_digit(h, r, digit_cw(digit_of(h, r)));
  return h;
}
inline std::uint64_t rotate_pent_ccw(std::uint64_t h) {
  bool found = false;
  for (int r = 1, res = res_of(h); r <= res; r++) {
    set_digit(h, r, digit_ccw(digit_of(h, r)));
    if (!found && digit_of(h, r) != 0) {
      found = true;
      if (leading_digit(h) == K_AXES) h = rotate_cell_ccw(h);
    }
  }
  return h;
}

inline double pos_angle(double a) {
  double t = a < 0.0 ? a + TWO_PI : a;
  if (a >= TWO_PI) t -= TWO_PI;
  return t;
}
inline double constrain_lng(double l) {
  while (l > PI) l -= TWO_PI;
  while (l < -PI) l += TWO_PI;
  return l;
}
inline double azimuth(LatLng a, LatLng b) {
  return std::atan2(std::cos(b.lat) * std::sin(b.lng - a.lng), std::cos(a.lat) * std::sin(b.lat) - std::sin(a.lat) * std::cos(b.lat) * std::cos(b.lng - a.lng));
}
inline LatLng az_distance(LatLng p, double az, double distance) {
  if (distance < EPS) return p;
  LatLng q;
  az = pos_angle(az);
  if (az < EPS || std::fabs(az - PI) < EPS) {
    q.lat = az < EPS ? p.lat + distance : p.lat - distance;
    if (std::fabs(q.lat - HALF_PI) < EPS) { q.lat = HALF_PI; q.lng = 0.0; }
    else if (std::fabs(q.lat + HALF_PI) < EPS) { q.lat = -HALF_PI; q.lng = 0.0; }
    else q.lng = constrain_lng(p.lng);
    return q;
  }
  double sinlat = std::sin(p.lat) * std::cos(distance) + std::cos(p.lat) * std::sin(distance) * std::cos(az);
  sinlat = std::max(-1.0, std::min(1.0, sinlat));
  q.lat = std::asin(sinlat);
  if (std::fabs(q.lat - HALF_PI) < EPS) { q.lat = HALF_PI; q.lng = 0.0; return q; }
  if (std::fabs(q.lat + HALF_PI) < EPS) { q.lat = -HALF_PI; q.lng = 0.0; return q; }
  const double inv = 1.0 / std::cos(q.lat);
  double sinlng = std::sin(az) * std::sin(distance) * inv, coslng = (std::cos(distance) - std::sin(p.lat) * std::sin(q.lat)) / std::cos(p.lat) * inv;
  sinlng = std::max(-1.0, std::min(1.0, sinlng));
  coslng = std::max(-1.0, std::min(1.0, coslng));
  q.lng = constrain_lng(p.lng + std::atan2(sinlng, coslng));
  return q;
}

// ---- the projection -----------------------------------------------------------------------------------------------

inline void closest_face(LatLng g, int& face, double& sqd) {
  const double r = std::cos(g.lat), v[3] = {std::cos(g.lng) * r, std::sin(g.lng) * r, std::sin(g.lat)};
  face = 0;
  sqd = 5.0;
  for (int f = 0; f < 20; f++) {
    const auto& c = tables::face_center_point[f];
    const double d = (c[0] - v[0]) * (c[0] - v[0]) + (c[1] - v[1]) * (c[1] - v[1]) + (c[2] - v[2]) * (c[2] - v[2]);
    if (d < sqd) { face = f; sqd = d; }
  }
}
inline Vec2 geo_to_hex2d(LatLng g, int res, int& face) {
  double sqd;
  closest_face(g, face, sqd);
  double r = std::acos(1 - sqd * 0.5);
  if (r < EPS) return {0, 0};
  const auto& c = tables::face_center_geo[face];
  double theta = pos_angle(tables::face_axes_az[face][0] - pos_angle(azimuth({c.lat, c.lng}, g)));
  if (class3(res)) theta = pos_angle(theta - AP7_ROT);
  r = std::tan(r) * INV_RES0_U_GNOMONIC;
  for (int i = 0; i < res; i++) r *= SQRT7;
  return {r * std::cos(theta), r * std::sin(theta)};
}
inline LatLng hex2d_to_geo(Vec2 v, int face, int res, bool substrate) {
  double r = std::hypot(v.x, v.y);
  const auto& c = tables::face_center_geo[face];
  if (r < EPS) return {c.lat, c.lng};
  double theta = std::atan2(v.y, v.x);
  for (int i = 0; i < res; i++) r *= RSQRT7;
  if (substrate) {
    r *= ONE_THIRD;
    if (class3(res)) r *= RSQRT7;
  }
  r = std::atan(r * RES0_U_GNOMONIC);
  if (!substrate && class3(res)) theta = pos_angle(theta + AP7_ROT);
  theta = pos_angle(tables::face_axes_az[face][0] - theta);
  return az_distance({c.lat, c.lng}, theta, r);
}
inline Ijk hex2d_to_ijk(Vec2 v) {
  Ijk h{0, 0, 0};
  const double a1 = std::fabs(v.x), a2 = std::fabs(v.y), x2 = a2 * RSIN60, x1 = a1 + x2 / 2.0;
  const int m1 = static_cast<int>(x1), m2 = static_cast<int>(x2);
  const double r1 = x1 - m1, r2 = x2 - m2;
  if (r1 < 0.5) {
    if (r1 < 1.0 / 3.0) {
      h.i = m1;
      h.j = r2 < (1.0 + r1) / 2.0 ? m2 : m2 + 1;
    } else {
      h.j = r2 < (1.0 - r1) ? m2 : m2 + 1;
      h.i = ((1.0 - r1) <= r2 && r2 < (2.0 * r1)) ? m1 + 1 : m1;
    }
  } else {
    if (r1 < 2.0 / 3.0) {
      h.j = r2 < (1.0 - r1) ? m2 : m2 + 1;
      h.i = ((2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1)) ? m1 : m1 + 1;
    } else {
      h.i = m1 + 1;
      h.j = r2 < (r1 / 2.0) ? m2 : m2 + 1;
    }
  }
  if (v.x < 0.0) {
    if (h.j % 2 == 0) {
      const long long axis = h.j / 2, diff = h.i - axis;
      h.i = static_cast<int>(h.i - 2.0 * diff);
    } else {
      const long long axis = (h.j + 1) / 2, diff = h.i - axis;
      h.i = static_cast<int>(h.i - (2.0 * diff + 1));
    }
  }
  if (v.y < 0.0) {
    h.i = h.i - (2 * h.j + 1) / 2;
    h.j = -h.j;
  }
  normalize(h);
  return h;
}
inline Vec2 ijk_to_hex2d(Ijk h) {
  const int i = h.i - h.k, j = h.j - h.k;
  return {i - 0.5 * j, j * SQRT3_2};
}

// ---- crossing the faces ----------------------------------------------------------------------------------------------

enum Overage { NO_OVERAGE = 0, FACE_EDGE = 1, NEW_FACE = 2 };

inline Overage adjust_overage(FaceIjk& f, int res, bool pent_leading4, bool substrate) {
  Overage overage = NO_OVERAGE;
  Ijk& c = f.coord;
  int max_dim = tables::max_dim[res];
  if (substrate) max_dim *= 3;
  if (substrate && c.i + c.j + c.k == max_dim) overage = FACE_EDGE;
  else if (c.i + c.j + c.k > max_dim) {
    overage = NEW_FACE;
    const tables::FaceOrient* orient;
    if (c.k > 0) {
      if (c.j > 0) orient = &tables::face_neighbors[f.face][tables::JK];
      else {
        orient = &tables::face_neighbors[f.face][tables::KI];
        if (pent_leading4) {
          const Ijk origin{max_dim, 0, 0};
          Ijk t = sub(c, origin);
          rotate60cw(t);
          c = add(t, origin);
        }
      }
    } else orient = &tables::face_neighbors[f.face][tables::IJ];
    f.face = orient->face;
    for (int i = 0; i < orient->ccw_rot60; i++) rotate60ccw(c);
    int unit = tables::unit_scale[res];
    if (substrate) unit *= 3;
    c = add(c, scale(orient->translate, unit));
    normalize(c);
    if (substrate && c.i + c.j + c.k == max_dim) overage = FACE_EDGE;
  }
  return overage;
}
inline Overage adjust_pent_overage(FaceIjk& f, int res) {
  Overage o;
  do { o = adjust_overage(f, res, false, true); } while (o == NEW_FACE);
  return o;
}
inline Vec2 intersect(Vec2 p0, Vec2 p1, Vec2 p2, Vec2 p3) {
  const Vec2 s1{p1.x - p0.x, p1.y - p0.y}, s2{p3.x - p2.x, p3.y - p2.y};
  const double t = (s2.x * (p0.y - p2.y) - s2.y * (p0.x - p2.x)) / (-s2.x * s1.y + s1.x * s2.y);
  return {p0.x + t * s1.x, p0.y + t * s1.y};
}
inline bool almost_equal(Vec2 a, Vec2 b) {
  constexpr double FLOAT_EPS = 1.19209290e-07;
  return std::fabs(a.x - b.x) < FLOAT_EPS && std::fabs(a.y - b.y) < FLOAT_EPS;
}

// A cell's digits, followed down from its base cell's position on its home face.
inline FaceIjk cell_to_face_ijk(std::uint64_t h) {
  const int bc = base_of(h);
  if (bc < 0 || bc >= 122) throw std::invalid_argument("h3: no such base cell");
  if (base_is_pentagon(bc) && leading_digit(h) == IK_AXES) h = rotate_cell_cw(h);
  FaceIjk f{tables::base_cell_data[bc].home.face, tables::base_cell_data[bc].home.coord};
  const int cell_res = res_of(h);
  bool possible = true;
  if (!base_is_pentagon(bc) && (cell_res == 0 || (f.coord.i == 0 && f.coord.j == 0 && f.coord.k == 0))) possible = false;
  for (int r = 1; r <= cell_res; r++) {
    if (class3(r)) down_ap7(f.coord); else down_ap7r(f.coord);
    neighbor(f.coord, digit_of(h, r));
  }
  if (!possible) return f;
  const Ijk orig = f.coord;
  int res = cell_res;
  if (class3(res)) { down_ap7r(f.coord); res++; }
  const bool leading4 = base_is_pentagon(bc) && leading_digit(h) == I_AXES;
  if (adjust_overage(f, res, leading4, false) != NO_OVERAGE) {
    if (base_is_pentagon(bc)) while (adjust_overage(f, res, false, false) != NO_OVERAGE) continue;
    if (res != cell_res) up_ap7r(f.coord);
  } else if (res != cell_res) f.coord = orig;
  return f;
}

}  // namespace detail

// The resolution of a cell (0..15), its base cell (0..121), and whether it is one of the twelve pentagons of its resolution.
inline int resolution(std::uint64_t h) { return detail::res_of(h); }
inline int base_cell(std::uint64_t h) { return detail::base_of(h); }
inline bool is_pentagon(std::uint64_t h) { return detail::base_is_pentagon(detail::base_of(h)) && detail::leading_digit(h) == detail::CENTER; }

// A cell's index as the reference writes it (lower-case hexadecimal), and back.
inline std::string to_string(std::uint64_t h) {
  char buf[17];
  std::snprintf(buf, sizeof buf, "%llx", static_cast<unsigned long long>(h));
  return buf;
}
inline std::uint64_t from_string(std::string_view s) {
  std::uint64_t h = 0;
  for (const char c : s) {
    int d = c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1;
    if (d < 0) throw std::invalid_argument("h3: not a hexadecimal index");
    h = (h << 4) | static_cast<std::uint64_t>(d);
  }
  return h;
}

// The cell of a resolution holding a point (radians).
inline std::uint64_t cell_at(LatLng g, int res) {
  using namespace detail;
  if (res < 0 || res > MAX_RES) throw std::invalid_argument("h3: a resolution is 0..15");
  if (!std::isfinite(g.lat) || !std::isfinite(g.lng)) throw std::invalid_argument("h3: a point has finite coordinates");
  FaceIjk f;
  f.coord = hex2d_to_ijk(geo_to_hex2d(g, res, f.face));
  std::uint64_t h = INIT;
  h = (h & ~(std::uint64_t{15} << MODE_OFFSET)) | (std::uint64_t{1} << MODE_OFFSET);
  h = (h & ~(std::uint64_t{15} << RES_OFFSET)) | (static_cast<std::uint64_t>(res) << RES_OFFSET);
  const auto base_rotation = [](const FaceIjk& q) { return tables::face_ijk_base_cells[q.face][q.coord.i][q.coord.j][q.coord.k]; };
  const auto set_base = [&](int b) { h = (h & ~(std::uint64_t{127} << BC_OFFSET)) | (static_cast<std::uint64_t>(b) << BC_OFFSET); };
  if (res == 0) {
    if (f.coord.i > MAX_FACE_COORD || f.coord.j > MAX_FACE_COORD || f.coord.k > MAX_FACE_COORD) throw std::invalid_argument("h3: beyond the face");
    set_base(base_rotation(f).base_cell);
    return h;
  }
  FaceIjk bc = f;
  Ijk& c = bc.coord;
  for (int r = res - 1; r >= 0; r--) {
    const Ijk last = c;
    Ijk centre;
    if (class3(r + 1)) { up_ap7(c); centre = c; down_ap7(centre); }
    else { up_ap7r(c); centre = c; down_ap7r(centre); }
    Ijk diff = sub(last, centre);
    normalize(diff);
    set_digit(h, r + 1, unit_to_digit(diff));
  }
  if (bc.coord.i > MAX_FACE_COORD || bc.coord.j > MAX_FACE_COORD || bc.coord.k > MAX_FACE_COORD) throw std::invalid_argument("h3: beyond the face");
  const auto rot = base_rotation(bc);
  set_base(rot.base_cell);
  if (base_is_pentagon(rot.base_cell)) {
    if (leading_digit(h) == K_AXES) {
      const auto& off = tables::base_cell_data[rot.base_cell].cw_offset_pent;
      h = (off[0] == bc.face || off[1] == bc.face) ? rotate_cell_cw(h) : rotate_cell_ccw(h);
    }
    for (int i = 0; i < rot.ccw_rot60; i++) h = rotate_pent_ccw(h);
  } else {
    for (int i = 0; i < rot.ccw_rot60; i++) h = rotate_cell_ccw(h);
  }
  return h;
}

// A cell's centre (radians).
inline LatLng centre(std::uint64_t h) {
  using namespace detail;
  const FaceIjk f = cell_to_face_ijk(h);
  return hex2d_to_geo(ijk_to_hex2d(f.coord), f.face, res_of(h), false);
}

// A cell's boundary: its vertices in the reference's order (counter-clockwise seen from outside the sphere), six for a
// hexagon, five for a pentagon, more where a class III cell's edge crosses a face of the icosahedron.
inline std::vector<LatLng> boundary(std::uint64_t h) {
  using namespace detail;
  const int res = res_of(h);
  FaceIjk centre_f = cell_to_face_ijk(h);
  // The cell's vertices lie on a finer grid, of aperture 3 (twice, a further 7 for class III).
  const bool pent = is_pentagon(h);
  const int n = pent ? 5 : 6;
  static constexpr Ijk V2[6] = {{2, 1, 0}, {1, 2, 0}, {0, 2, 1}, {0, 1, 2}, {1, 0, 2}, {2, 0, 1}};
  static constexpr Ijk V3[6] = {{5, 4, 0}, {1, 5, 0}, {0, 5, 4}, {0, 1, 5}, {4, 0, 5}, {5, 0, 1}};
  down_ap3(centre_f.coord);
  down_ap3r(centre_f.coord);
  int adj = res;
  if (class3(res)) { down_ap7r(centre_f.coord); adj += 1; }
  std::array<FaceIjk, 6> verts;
  for (int v = 0; v < n; v++) {
    verts[v].face = centre_f.face;
    verts[v].coord = add(centre_f.coord, class3(res) ? V3[v] : V2[v]);
    normalize(verts[v].coord);
  }
  const int max_dim = tables::max_dim[adj];
  const Vec2 v0{3.0 * max_dim, 0.0}, v1{-1.5 * max_dim, 3.0 * SQRT3_2 * max_dim}, v2{-1.5 * max_dim, -3.0 * SQRT3_2 * max_dim};
  const auto edge_of = [&](int dir, Vec2& e0, Vec2& e1) {
    if (dir == tables::IJ) { e0 = v0; e1 = v1; }
    else if (dir == tables::JK) { e0 = v1; e1 = v2; }
    else { e0 = v2; e1 = v0; }
  };
  std::vector<LatLng> out;
  if (pent) {
    FaceIjk last;
    for (int vert = 0; vert < n + 1; vert++) {
      const int v = vert % n;
      FaceIjk f = verts[v];
      adjust_pent_overage(f, adj);
      if (class3(res) && vert > 0) {
        FaceIjk t = f;
        const Vec2 a = ijk_to_hex2d(last.coord);
        const int back = tables::adjacent_face_dir[t.face][last.face];
        const auto& orient = tables::face_neighbors[t.face][back];
        t.face = orient.face;
        for (int i = 0; i < orient.ccw_rot60; i++) rotate60ccw(t.coord);
        t.coord = add(t.coord, scale(orient.translate, tables::unit_scale[adj] * 3));
        normalize(t.coord);
        const Vec2 b = ijk_to_hex2d(t.coord);
        Vec2 e0, e1;
        edge_of(tables::adjacent_face_dir[t.face][f.face], e0, e1);
        out.push_back(hex2d_to_geo(intersect(a, b, e0, e1), t.face, adj, true));
      }
      if (vert < n) out.push_back(hex2d_to_geo(ijk_to_hex2d(f.coord), f.face, adj, true));
      last = f;
    }
    return out;
  }
  int last_face = -1;
  Overage last_overage = NO_OVERAGE;
  for (int vert = 0; vert < n + 1; vert++) {
    const int v = vert % n;
    FaceIjk f = verts[v];
    const Overage overage = adjust_overage(f, adj, false, true);
    if (class3(res) && vert > 0 && f.face != last_face && last_overage != FACE_EDGE) {
      const int lastv = (v + 5) % n;
      const Vec2 a = ijk_to_hex2d(verts[lastv].coord), b = ijk_to_hex2d(verts[v].coord);
      const int face2 = last_face == centre_f.face ? f.face : last_face;
      Vec2 e0, e1;
      edge_of(tables::adjacent_face_dir[centre_f.face][face2], e0, e1);
      const Vec2 inter = intersect(a, b, e0, e1);
      if (!almost_equal(a, inter) && !almost_equal(b, inter)) out.push_back(hex2d_to_geo(inter, centre_f.face, adj, true));
    }
    if (vert < n) out.push_back(hex2d_to_geo(ijk_to_hex2d(f.coord), f.face, adj, true));
    last_face = f.face;
    last_overage = overage;
  }
  return out;
}

// The cell of a coarser resolution that holds a cell.
inline std::uint64_t parent(std::uint64_t h, int res) {
  using namespace detail;
  const int own = res_of(h);
  if (res < 0 || res > own) throw std::invalid_argument("h3: a parent's resolution is no finer than its cell's");
  h = (h & ~(std::uint64_t{15} << RES_OFFSET)) | (static_cast<std::uint64_t>(res) << RES_OFFSET);
  for (int r = res + 1; r <= own; r++) set_digit(h, r, 7);
  return h;
}

// A cell's area on the sphere of radius R (m), from the triangles it makes with its centre (each by l'Huilier's theorem).
inline double cell_area_m2(std::uint64_t h, double radius_m = 6371007.180918475) {
  const LatLng c = centre(h);
  const auto bnd = boundary(h);
  const auto arc = [](LatLng a, LatLng b) {
    const double sa = std::sin((b.lat - a.lat) / 2), so = std::sin((b.lng - a.lng) / 2);
    return 2 * std::atan2(std::sqrt(sa * sa + std::cos(a.lat) * std::cos(b.lat) * so * so), std::sqrt(std::max(0.0, 1 - sa * sa - std::cos(a.lat) * std::cos(b.lat) * so * so)));
  };
  double total = 0;
  for (std::size_t i = 0; i < bnd.size(); i++) {
    const LatLng a = bnd[i], b = bnd[(i + 1) % bnd.size()];
    const double la = arc(c, a), lb = arc(c, b), lc = arc(a, b), s = (la + lb + lc) / 2;
    total += 4 * std::atan(std::sqrt(std::max(0.0, std::tan(s / 2) * std::tan((s - la) / 2) * std::tan((s - lb) / 2) * std::tan((s - lc) / 2))));
  }
  return total * radius_m * radius_m;
}

}  // namespace physis::computational_geometry::h3
