// Point groups: the isometries of a finite set of labelled points in ℝ³, their classification,
// and the responses they allow.
//
// The symmetry group of a finite labelled set A ⊂ ℝ³ is the group of isometries of ℝ³ that map A
// onto itself, labels kept. Every such isometry fixes the centroid c (it permutes the points, so
// it fixes their mean), hence is x ↦ c + Q(x − c) with Q ∈ O(3): the group is a finite subgroup of
// O(3), a point group, as long as the points do not lie on one line (then the group contains every
// rotation about that line, and is C∞v or D∞h). The finite subgroups of O(3) are the seven infinite
// families Cn, Cnv, Cnh, S2n, Dn, Dnh, Dnd and the seven groups T, Td, Th, O, Oh, I, Ih (Hessel
// 1830; Cotton 1990, Chemical Applications of Group Theory, 3rd ed., ch. 3; Conway and Smith 2003,
// On Quaternions and Octonions, ch. 3). point_group() finds the group and names it in Schoenflies's
// notation.
//
// Method. With the points centred, an orthogonal Q is fixed by the images of three linearly
// independent points a, b, c: Q = [Qa Qb Qc][a b c]⁻¹. So the group is found by trying, for a chosen
// independent triple, every triple of images of the same labels, norms and inner products, and
// keeping the Q that are orthogonal and permute the whole set. When the points span only a plane,
// a and b fix Q on the plane and Q(a × b) = det Q · (Qa × Qb) fixes it on the normal, for each sign
// of det Q; a planar set is thus always fixed by the reflection in its plane. congruences() does
// the same between two sets, which decides whether two labelled configurations are images of each
// other by a rotation, by a reflection, or not at all (the two configurations of a tetrahedron
// with four distinct labels are images only by improper isometries: they are enantiomorphic).
//
// Classification. With n the highest order of a proper rotation: two or more 3-fold axes make the
// group cubic or icosahedral (I if a 5-fold rotation, O if a 4-fold, T otherwise; then h with the
// inversion, d with mirrors and no inversion). Otherwise n C₂ axes perpendicular to the principal
// axis make it dihedral (Dnh with the horizontal mirror, Dnd with n vertical mirrors, else Dn), and
// their absence cyclic (Cnh, Cnv, S2n, or Cn; for n = 1: Cs, Ci or C1). The order is checked against
// the symbol's (|Cn| = n, |Cnv| = |Cnh| = |Dn| = |S2n| = 2n, |Dnh| = |Dnd| = 4n, |T| = 12,
// |Td| = |Th| = |O| = 24, |Oh| = 48, |I| = 60, |Ih| = 120).
//
// Responses. A response of a state fixed by H, transforming under H by a representation ρ, must
// lie in the fixed subspace of ρ(H) (Neumann's principle; Nye 1957, Physical Properties of
// Crystals, ch. 1; Curie 1894). Its dimension is the mean of the character,
//
//   dim V^H = (1/|H|) Σ_{h ∈ H} χ_ρ(h),
//
// for the vector (polar) representation χ(Q) = tr Q and the axial (pseudovector) one
// χ(Q) = det Q · tr Q. A static polar vector, such as an electric dipole, is allowed exactly in the
// polar groups C1, Cs, Cn, Cnv; an axial one in C1, Ci, Cs, Cn, Cnh, S2n. A set is chiral (not
// superposable on its mirror image) exactly when its group holds only proper rotations: C1, Cn,
// Dn, T, O, I.
//
// The group is returned with its multiplication table as a morphe::algebra::FiniteGroup, the
// identity first, so that the exact methods of morphe (conjugacy classes, character tables)
// apply to it. The table is built from the matrices, not the permutations of the points, which
// need not be faithful (the reflection in the plane of a planar set moves no point).
//
// Tolerances are relative: lengths to tol times the largest distance from the centroid.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/algebra/finite_group.hpp>

namespace physis::symmetry {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;  // rows

struct LabelledPoints {
  std::vector<std::string> label;
  std::vector<Vec3> x;
  std::size_t size() const { return x.size(); }
  void add(std::string l, Vec3 p) {
    label.push_back(std::move(l));
    x.push_back(p);
  }
};

namespace pg_detail {

inline Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec3 cross(const Vec3& a, const Vec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 scale(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
inline Vec3 apply(const Mat3& Q, const Vec3& v) { return {dot(Q[0], v), dot(Q[1], v), dot(Q[2], v)}; }
inline Mat3 mul(const Mat3& A, const Mat3& B) {
  Mat3 C{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++) C[i][j] += A[i][k] * B[k][j];
  return C;
}
inline Mat3 transpose(const Mat3& A) {
  Mat3 T{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) T[i][j] = A[j][i];
  return T;
}
inline double det(const Mat3& A) { return dot(A[0], cross(A[1], A[2])); }
inline Mat3 columns(const Vec3& a, const Vec3& b, const Vec3& c) { return {{{a[0], b[0], c[0]}, {a[1], b[1], c[1]}, {a[2], b[2], c[2]}}}; }
inline Mat3 inverse(const Mat3& A) {
  const double d = det(A);
  Mat3 I{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      const int a = (j + 1) % 3, b = (j + 2) % 3, c = (i + 1) % 3, e = (i + 2) % 3;
      I[i][j] = (A[a][c] * A[b][e] - A[a][e] * A[b][c]) / d;
    }
  return I;
}
inline double max_abs_diff(const Mat3& A, const Mat3& B) {
  double m = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m = std::max(m, std::fabs(A[i][j] - B[i][j]));
  return m;
}
inline Mat3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }

}  // namespace pg_detail

// x ↦ c + Q(x − c), and the permutation of the points it induces: point i goes to image[i].
struct Isometry {
  Mat3 Q{};
  std::vector<std::size_t> image;
  int det() const { return pg_detail::det(Q) > 0 ? 1 : -1; }
  double trace() const { return Q[0][0] + Q[1][1] + Q[2][2]; }
  bool proper() const { return det() > 0; }
};

struct Centred {
  std::vector<Vec3> y;
  Vec3 centroid{};
  double radius = 0;   // the largest distance from the centroid
  std::size_t rank = 0;  // the dimension of the span of the centred points
  std::array<std::size_t, 3> basis{};  // independent points, as many as the rank
};

inline Centred centre(const LabelledPoints& A, double tol) {
  using namespace pg_detail;
  Centred C;
  const std::size_t n = A.size();
  if (n == 0) throw std::invalid_argument("physis::symmetry::centre: no points");
  for (const auto& p : A.x)
    for (int k = 0; k < 3; k++) C.centroid[k] += p[k] / static_cast<double>(n);
  for (const auto& p : A.x) {
    C.y.push_back(sub(p, C.centroid));
    C.radius = std::max(C.radius, norm(C.y.back()));
  }
  const double r = C.radius;
  if (r == 0) return C;
  // The rarest labels first, so that few images are tried.
  std::vector<std::size_t> idx(n);
  for (std::size_t i = 0; i < n; i++) idx[i] = i;
  const auto count = [&](std::size_t i) { return std::count(A.label.begin(), A.label.end(), A.label[i]); };
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t i, std::size_t j) { return count(i) < count(j); });
  for (auto i : idx)
    if (norm(C.y[i]) > tol * r) { C.basis[0] = i; C.rank = 1; break; }
  if (C.rank < 1) return C;
  for (auto i : idx)
    if (norm(cross(C.y[C.basis[0]], C.y[i])) > tol * r * r) { C.basis[1] = i; C.rank = 2; break; }
  if (C.rank < 2) return C;
  const Vec3 n01 = cross(C.y[C.basis[0]], C.y[C.basis[1]]);
  for (auto i : idx)
    if (std::fabs(dot(n01, C.y[i])) > tol * r * r * r) { C.basis[2] = i; C.rank = 3; break; }
  return C;
}

// The isometries x ↦ c_B + Q(x − c_A) that carry A onto B, labels kept. Throws if the points lie
// on a line (infinitely many) or in one point.
inline std::vector<Isometry> congruences(const LabelledPoints& A, const LabelledPoints& B, double tol = 1e-6) {
  using namespace pg_detail;
  std::vector<Isometry> out;
  if (A.size() != B.size()) return out;
  const Centred a = centre(A, tol), b = centre(B, tol);
  if (std::fabs(a.radius - b.radius) > tol * std::max(1.0, a.radius)) return out;
  if (a.rank < 2) throw std::invalid_argument("physis::symmetry::congruences: the points lie on a line; the group is infinite");
  const double r = a.radius, eps = tol * r;
  const std::size_t n = A.size();
  const auto candidates = [&](std::size_t i, const std::vector<std::size_t>& chosen, const std::vector<std::size_t>& images) {
    std::vector<std::size_t> c;
    for (std::size_t j = 0; j < n; j++) {
      if (B.label[j] != A.label[i] || std::fabs(norm(b.y[j]) - norm(a.y[i])) > eps) continue;
      bool ok = true;
      for (std::size_t k = 0; k < chosen.size() && ok; k++)
        ok = images[k] != j && std::fabs(dot(b.y[j], b.y[images[k]]) - dot(a.y[i], a.y[chosen[k]])) <= 4 * eps * r;
      if (ok) c.push_back(j);
    }
    return c;
  };
  const auto accept = [&](const Mat3& Q) {
    if (max_abs_diff(mul(transpose(Q), Q), identity()) > 1e3 * tol) return;
    for (const auto& iso : out)
      if (max_abs_diff(iso.Q, Q) < 1e3 * tol) return;
    Isometry iso;
    iso.Q = Q;
    iso.image.assign(n, n);
    std::vector<char> used(n, 0);
    for (std::size_t i = 0; i < n; i++) {
      const Vec3 y = apply(Q, a.y[i]);
      std::size_t found = n;
      for (std::size_t j = 0; j < n; j++)
        if (!used[j] && B.label[j] == A.label[i] && norm(sub(b.y[j], y)) <= 10 * eps) { found = j; break; }
      if (found == n) return;
      used[found] = 1;
      iso.image[i] = found;
    }
    out.push_back(std::move(iso));
  };
  const std::size_t i0 = a.basis[0], i1 = a.basis[1];
  for (auto j0 : candidates(i0, {}, {}))
    for (auto j1 : candidates(i1, {i0}, {j0})) {
      if (a.rank == 3) {
        const std::size_t i2 = a.basis[2];
        const Mat3 Ainv = inverse(columns(a.y[i0], a.y[i1], a.y[i2]));
        for (auto j2 : candidates(i2, {i0, i1}, {j0, j1})) accept(mul(columns(b.y[j0], b.y[j1], b.y[j2]), Ainv));
      } else {
        const Mat3 Ainv = inverse(columns(a.y[i0], a.y[i1], cross(a.y[i0], a.y[i1])));
        const Vec3 m = cross(b.y[j0], b.y[j1]);
        for (double s : {1.0, -1.0}) accept(mul(columns(b.y[j0], b.y[j1], scale(s, m)), Ainv));
      }
    }
  return out;
}

// A symmetry element, for display: a rotation axis with its highest order, a rotoreflection axis,
// a mirror plane (by its normal), or the centre of inversion.
struct SymmetryElement {
  enum class Kind { rotation, rotoreflection, mirror, inversion } kind = Kind::rotation;
  int order = 1;
  Vec3 axis{};  // the axis, or the mirror's normal; unit length
};

struct PointGroup {
  std::string schoenflies;            // "C1", "C2v", "D6h", "Td", …; "C∞v" or "D∞h" for a line
  bool linear = false;                // the points lie on a line (the group is infinite)
  std::vector<Isometry> elements;     // the identity first (empty when linear)
  morphe::algebra::FiniteGroup group; // its multiplication table, elements in the same order
  int principal_order = 1;
  bool inversion = false;
  std::size_t mirrors = 0;
  std::vector<SymmetryElement> symmetry_elements;
  Vec3 centroid{};

  std::size_t order() const { return elements.size(); }
  // A set on a line is fixed by the mirrors through the line, so it is never chiral.
  bool chiral() const {
    if (linear) return false;
    for (const auto& g : elements)
      if (!g.proper()) return false;
    return true;
  }
  // dim of the vectors fixed by the whole group, polar (tr Q) or axial (det Q · tr Q).
  std::size_t fixed_vectors(bool axial = false) const {
    if (elements.empty()) return 0;
    double s = 0;
    for (const auto& g : elements) s += (axial ? g.det() : 1) * g.trace();
    return static_cast<std::size_t>(std::lround(s / static_cast<double>(elements.size())));
  }
  // The mean of an arbitrary character χ(Q): the multiplicity of the trivial representation.
  double average(const std::function<double(const Mat3&)>& chi) const {
    double s = 0;
    for (const auto& g : elements) s += chi(g.Q);
    return s / static_cast<double>(elements.size());
  }
};

namespace pg_detail {

// The rotation axis of a proper Q (unit; sign arbitrary for a half turn) and its angle.
inline std::pair<Vec3, double> axis_angle(const Mat3& Q) {
  const double c = std::clamp((Q[0][0] + Q[1][1] + Q[2][2] - 1) / 2, -1.0, 1.0);
  const double theta = std::acos(c);
  Vec3 a{Q[2][1] - Q[1][2], Q[0][2] - Q[2][0], Q[1][0] - Q[0][1]};
  if (norm(a) > 1e-6) return {scale(1 / norm(a), a), theta};
  if (c > 0) return {{0, 0, 1}, 0};
  // A half turn: Q + I = 2 a aᵀ; take its largest column.
  Vec3 best{};
  for (int j = 0; j < 3; j++) {
    const Vec3 col{Q[0][j] + (j == 0), Q[1][j] + (j == 1), Q[2][j] + (j == 2)};
    if (norm(col) > norm(best)) best = col;
  }
  return {scale(1 / norm(best), best), theta};
}
inline int element_order(const Mat3& Q, double tol) {
  Mat3 P = Q;
  for (int k = 1; k <= 240; k++) {
    if (max_abs_diff(P, identity()) < tol) return k;
    P = mul(P, Q);
  }
  return 0;
}
inline bool parallel(const Vec3& a, const Vec3& b, double tol) { return norm(cross(a, b)) < tol; }

}  // namespace pg_detail

inline PointGroup point_group(const LabelledPoints& A, double tol = 1e-6) {
  using namespace pg_detail;
  PointGroup G;
  const Centred c = centre(A, tol);
  G.centroid = c.centroid;
  if (c.rank < 2) {
    // On a line: C∞v, or D∞h when the inversion maps the set onto itself.
    G.linear = true;
    bool inv = true;
    for (std::size_t i = 0; i < A.size() && inv; i++) {
      bool found = false;
      for (std::size_t j = 0; j < A.size() && !found; j++)
        found = A.label[j] == A.label[i] && norm(sub(c.y[j], scale(-1, c.y[i]))) <= 10 * tol * std::max(1.0, c.radius);
      inv = found;
    }
    G.schoenflies = c.rank == 0 ? "K_h" : inv ? "D∞h" : "C∞v";
    G.inversion = inv;
    return G;
  }
  G.elements = congruences(A, A, tol);
  // The identity first.
  const double mtol = 1e3 * tol;
  for (std::size_t i = 0; i < G.elements.size(); i++)
    if (max_abs_diff(G.elements[i].Q, identity()) < mtol) { std::swap(G.elements[0], G.elements[i]); break; }
  const std::size_t n = G.elements.size();
  // The multiplication table, from the matrices.
  std::vector<std::uint32_t> table(n * n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      const Mat3 P = mul(G.elements[i].Q, G.elements[j].Q);
      std::size_t k = n;
      for (std::size_t m = 0; m < n && k == n; m++)
        if (max_abs_diff(P, G.elements[m].Q) < mtol) k = m;
      if (k == n) throw std::runtime_error("physis::symmetry::point_group: the isometries found are not closed under composition");
      table[i * n + j] = static_cast<std::uint32_t>(k);
    }
  // Proper rotations: orders and axes.
  int maxn = 1, c3 = 0, c4 = 0, c5 = 0;
  Vec3 principal{0, 0, 1};
  std::vector<std::pair<Vec3, int>> rot_axes;  // each axis line with its highest order
  for (const auto& g : G.elements) {
    if (!g.proper()) {
      if (max_abs_diff(g.Q, Mat3{{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}}) < mtol) G.inversion = true;
      else if (std::fabs(g.trace() - 1) < 1e-6) G.mirrors++;
      continue;
    }
    const int k = element_order(g.Q, mtol);
    if (k == 1) continue;
    if (k == 3) c3++;
    if (k == 4) c4++;
    if (k == 5) c5++;
    const Vec3 ax = axis_angle(g.Q).first;
    bool found = false;
    for (auto& [a, o] : rot_axes)
      if (parallel(a, ax, 1e-6)) { o = std::max(o, k); found = true; }
    if (!found) rot_axes.push_back({ax, k});
    if (k > maxn) { maxn = k; principal = ax; }
  }
  G.principal_order = maxn;
  for (const auto& [a, o] : rot_axes) G.symmetry_elements.push_back({SymmetryElement::Kind::rotation, o, a});
  // Mirrors (normal = axis of the half turn −Q) and rotoreflections (axis of −Q, by its order).
  std::vector<std::pair<Vec3, int>> s_axes;
  for (const auto& g : G.elements) {
    if (g.proper()) continue;
    Mat3 minus = g.Q;
    for (auto& row : minus)
      for (auto& v : row) v = -v;
    if (max_abs_diff(minus, identity()) < mtol) continue;  // the inversion
    const Vec3 ax = axis_angle(minus).first;
    if (std::fabs(g.trace() - 1) < 1e-6) {
      G.symmetry_elements.push_back({SymmetryElement::Kind::mirror, 2, ax});
      continue;
    }
    const int k = element_order(g.Q, mtol);
    bool found = false;
    for (auto& [a, o] : s_axes)
      if (parallel(a, ax, 1e-6)) { o = std::max(o, k); found = true; }
    if (!found) s_axes.push_back({ax, k});
  }
  for (const auto& [a, o] : s_axes) G.symmetry_elements.push_back({SymmetryElement::Kind::rotoreflection, o, a});
  if (G.inversion) G.symmetry_elements.push_back({SymmetryElement::Kind::inversion, 2, {0, 0, 0}});

  // Schoenflies's symbol.
  const auto N = std::to_string(maxn);
  std::size_t expected = 0;
  bool sigma_h = false;
  std::size_t c2_perp = 0;
  for (const auto& e : G.symmetry_elements) {
    if (e.kind == SymmetryElement::Kind::mirror && parallel(e.axis, principal, 1e-6)) sigma_h = true;
    if (e.kind == SymmetryElement::Kind::rotation && e.order % 2 == 0 && std::fabs(dot(e.axis, principal)) < 1e-6) c2_perp++;
  }
  // A C₂ axis perpendicular to the principal one counts once per axis line; for D2 any axis serves.
  if (c3 >= 4) {
    const std::string base = c5 ? "I" : c4 ? "O" : "T";
    const std::size_t proper_order = c5 ? 60 : c4 ? 24 : 12;
    if (G.inversion) G.schoenflies = base + "h", expected = 2 * proper_order;
    else if (G.mirrors) G.schoenflies = base + "d", expected = 2 * proper_order;
    else G.schoenflies = base, expected = proper_order;
  } else if (maxn == 1) {
    if (G.mirrors) G.schoenflies = "Cs", expected = 2;
    else if (G.inversion) G.schoenflies = "Ci", expected = 2;
    else G.schoenflies = "C1", expected = 1;
  } else if (c2_perp >= static_cast<std::size_t>(maxn)) {
    if (sigma_h) G.schoenflies = "D" + N + "h", expected = 4 * maxn;
    else if (G.mirrors >= static_cast<std::size_t>(maxn)) G.schoenflies = "D" + N + "d", expected = 4 * maxn;
    else G.schoenflies = "D" + N, expected = 2 * maxn;
  } else {
    if (sigma_h) G.schoenflies = "C" + N + "h", expected = 2 * maxn;
    else if (G.mirrors >= static_cast<std::size_t>(maxn)) G.schoenflies = "C" + N + "v", expected = 2 * maxn;
    else if (n == static_cast<std::size_t>(2 * maxn)) G.schoenflies = "S" + std::to_string(2 * maxn), expected = 2 * maxn;
    else G.schoenflies = "C" + N, expected = maxn;
  }
  if (expected != n)
    throw std::runtime_error("physis::symmetry::point_group: " + std::to_string(n) + " isometries do not make the group " + G.schoenflies);
  G.group = morphe::algebra::FiniteGroup::from_table(n, std::move(table), G.schoenflies);
  return G;
}

// The closure of a set of orthogonal matrices under multiplication: the finite group they
// generate (throws past `limit` elements).
inline std::vector<Mat3> generate(const std::vector<Mat3>& generators, std::size_t limit = 240, double tol = 1e-9) {
  using namespace pg_detail;
  std::vector<Mat3> G{identity()};
  for (std::size_t i = 0; i < G.size(); i++)
    for (const auto& g : generators) {
      const Mat3 P = mul(g, G[i]);
      bool found = false;
      for (const auto& h : G) found = found || max_abs_diff(h, P) < tol;
      if (!found) {
        if (G.size() >= limit) throw std::runtime_error("physis::symmetry::generate: the group exceeds the limit");
        G.push_back(P);
      }
    }
  return G;
}

// The rotation by angle θ about a unit axis (Rodrigues's formula), and the reflection in the
// plane with the given unit normal.
inline Mat3 rotation(const Vec3& axis, double theta) {
  const double c = std::cos(theta), s = std::sin(theta), t = 1 - c;
  const double x = axis[0], y = axis[1], z = axis[2];
  return {{{t * x * x + c, t * x * y - s * z, t * x * z + s * y},
           {t * x * y + s * z, t * y * y + c, t * y * z - s * x},
           {t * x * z - s * y, t * y * z + s * x, t * z * z + c}}};
}
inline Mat3 reflection(const Vec3& normal) {
  Mat3 R = pg_detail::identity();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) R[i][j] -= 2 * normal[i] * normal[j];
  return R;
}

}  // namespace physis::symmetry
