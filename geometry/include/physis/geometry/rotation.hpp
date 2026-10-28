// Rotations of Euclidean 3-space: the group SO(3), its double cover SU(2), and the
// parametrisations by which the sciences write them.
//
// A rotation is stored as a unit quaternion q = w + x i + y j + z k, the rotation quaternion
// of Hamilton, acting on a vector v (a pure quaternion) by v ↦ q v q*: the active,
// right-handed rotation by the angle θ about the unit axis n has q = cos(θ/2) + sin(θ/2) n
// (Altmann 1986, Rotations, Quaternions and Double Groups, ch. 12). q and −q are the same
// rotation; the two are the two elements of SU(2) over it. physis `symmetry` writes SU(2) as
// U = a0 I + i a·σ, and its element over the rotation (w, x, y, z) is (a0, a) = (w, −x, −y, −z)
// (symmetry/su2.hpp, `adjoint`); `su2()` and `from_su2` pass between the two.
//
// Conventions (physis glossary, §1):
// - **Active and right-handed.** `matrix()` is the matrix R with R v the rotated vector. A
//   rotation of the *axes* of a frame by Q changes the coordinates of a fixed vector by Qᵀ;
//   that passive matrix is what the astronomical literature writes R_x(φ), R_y, R_z (SOFA's
//   iauRx: [[1, 0, 0], [0, cos φ, sin φ], [0, −sin φ, cos φ]]), here `frame_rotation`.
// - **Euler angles.** A sequence of three axes a1 a2 a3 with angles (α, β, γ). *Intrinsic*:
//   rotate by α about a1, then by β about the new a2, then by γ about the newest a3, so that
//   R = R_a1(α) R_a2(β) R_a3(γ) with active matrices. *Extrinsic*: the same rotations about
//   the fixed axes in the order a1, a2, a3, so R = R_a3(γ) R_a2(β) R_a1(α). The aerospace
//   attitude (yaw ψ, pitch θ, roll φ) is the intrinsic sequence Z-Y-X: the body's axes are
//   R_z(ψ) R_y(θ) R_x(φ) applied to the axes of the local frame (north, east, down). The six
//   Tait–Bryan sequences (three distinct axes) have a singularity at β = ±π/2, the six proper
//   Euler sequences (a1 = a3) at β = 0 and π; there the first angle carries the whole
//   rotation about the locked axis and the third is set to zero (Diebel 2006, Representing
//   Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors, Stanford, §5).
// - **The rotation vector** θ n is the exponential coordinate of SO(3): `exp` and `log`.
// - **From a matrix** by Shepperd's method, which divides by the largest of the four
//   squared components and so stays accurate near θ = π (Shepperd 1978, J. Guidance and
//   Control 1(3), 223–224).
// - **Interpolation** along the geodesic of S³, Shoemake's slerp on the shorter arc
//   (Shoemake 1985, SIGGRAPH '85, 245–254).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>

#include <physis/symmetry/su2.hpp>

namespace physis::geometry {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;

// ---- vectors and matrices of order 3 --------------------------------------------------------
//
// Vec3 and Mat3 are std::arrays, which argument-dependent lookup does not connect with this
// namespace; their operators live in the inline namespace `operators`, which code in other
// namespaces brings in with `using namespace physis::geometry::operators;`.

inline namespace operators {
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 operator-(const Vec3& a) { return {-a[0], -a[1], -a[2]}; }
inline Vec3 operator*(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
}  // namespace operators

inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(const Vec3& a) { return std::hypot(a[0], a[1], a[2]); }
inline Vec3 unit(const Vec3& a) {
  const double n = norm(a);
  if (n == 0) throw std::domain_error("physis::geometry::unit: the zero vector has no direction");
  return (1 / n) * a;
}

inline Mat3 identity3() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }
inline Mat3 zero3() { return {}; }
inline namespace operators {
inline Mat3 operator*(const Mat3& a, const Mat3& b) {
  Mat3 c{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
  return c;
}
inline Vec3 operator*(const Mat3& a, const Vec3& v) {
  return {a[0][0] * v[0] + a[0][1] * v[1] + a[0][2] * v[2], a[1][0] * v[0] + a[1][1] * v[1] + a[1][2] * v[2],
          a[2][0] * v[0] + a[2][1] * v[1] + a[2][2] * v[2]};
}
inline Mat3 operator+(const Mat3& a, const Mat3& b) {
  Mat3 c{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) c[i][j] = a[i][j] + b[i][j];
  return c;
}
inline Mat3 operator-(const Mat3& a, const Mat3& b) {
  Mat3 c{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) c[i][j] = a[i][j] - b[i][j];
  return c;
}
inline Mat3 operator*(double s, const Mat3& a) {
  Mat3 c{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) c[i][j] = s * a[i][j];
  return c;
}
}  // namespace operators

inline Mat3 transpose(const Mat3& a) {
  Mat3 t{};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) t[i][j] = a[j][i];
  return t;
}
// [ω]×, the matrix of v ↦ ω × v.
inline Mat3 skew(const Vec3& w) { return {{{0, -w[2], w[1]}, {w[2], 0, -w[0]}, {-w[1], w[0], 0}}}; }
// The vector of the antisymmetric part of a matrix: vee([ω]×) = ω.
inline Vec3 vee(const Mat3& a) { return {0.5 * (a[2][1] - a[1][2]), 0.5 * (a[0][2] - a[2][0]), 0.5 * (a[1][0] - a[0][1])}; }
// The largest entry of |a − b|.
inline double max_difference(const Mat3& a, const Mat3& b) {
  double m = 0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m = std::max(m, std::fabs(a[i][j] - b[i][j]));
  return m;
}
// ‖R Rᵀ − I‖ and det R − 1: how far a matrix is from SO(3).
inline double orthogonality_defect(const Mat3& r) { return max_difference(r * transpose(r), identity3()); }
inline double determinant(const Mat3& a) {
  return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
         a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

// ---- elementary rotations ------------------------------------------------------------------

enum class Axis { x = 0, y = 1, z = 2 };

// The active rotation by `angle` about a coordinate axis.
inline Mat3 rotation_matrix(Axis axis, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  switch (axis) {
    case Axis::x: return {{{1, 0, 0}, {0, c, -s}, {0, s, c}}};
    case Axis::y: return {{{c, 0, s}, {0, 1, 0}, {-s, 0, c}}};
    default: return {{{c, -s, 0}, {s, c, 0}, {0, 0, 1}}};
  }
}
// The change of coordinates when a frame's axes turn by `angle` about one of its axes: the
// transpose of the active rotation, SOFA's iauRx, iauRy, iauRz.
inline Mat3 frame_rotation(Axis axis, double angle) { return rotation_matrix(axis, -angle); }

// ---- Euler sequences -------------------------------------------------------------------------

enum class EulerConvention { intrinsic, extrinsic };

// A sequence of three axes, a1 ≠ a2 ≠ a3: Tait–Bryan when a1 ≠ a3, proper Euler when a1 = a3.
struct EulerSequence {
  Axis a1, a2, a3;
  EulerConvention convention = EulerConvention::intrinsic;

  bool proper() const { return a1 == a3; }
  // "ZYX" or "zyx"; a trailing "-extrinsic" or a lower-case spelling chooses the convention
  // only when written explicitly: parse("ZYX") is intrinsic, parse("ZYX", extrinsic) extrinsic.
  static EulerSequence parse(std::string_view axes, EulerConvention convention = EulerConvention::intrinsic) {
    if (axes.size() != 3) throw std::invalid_argument("an Euler sequence has three axes: " + std::string(axes));
    auto axis = [&](char c) {
      switch (c) {
        case 'x': case 'X': return Axis::x;
        case 'y': case 'Y': return Axis::y;
        case 'z': case 'Z': return Axis::z;
        default: throw std::invalid_argument("not an axis: " + std::string(1, c));
      }
    };
    EulerSequence s{axis(axes[0]), axis(axes[1]), axis(axes[2]), convention};
    if (s.a1 == s.a2 || s.a2 == s.a3) throw std::invalid_argument("consecutive Euler axes must differ: " + std::string(axes));
    return s;
  }
  std::string name() const {
    const char* letters = "XYZ";
    std::string s{letters[int(a1)], letters[int(a2)], letters[int(a3)]};
    return s + (convention == EulerConvention::intrinsic ? " (intrinsic)" : " (extrinsic)");
  }
};

// ---- the rotation ----------------------------------------------------------------------------

class Rotation {
 public:
  Rotation() = default;
  static Rotation identity() { return {}; }

  // From a rotation quaternion (normalised here; its sign is kept).
  static Rotation from_quaternion(double w, double x, double y, double z) {
    const double n = std::sqrt(w * w + x * x + y * y + z * z);
    if (!(n > 0)) throw std::domain_error("physis::geometry::Rotation: the zero quaternion");
    Rotation r;
    r.w_ = w / n, r.x_ = x / n, r.y_ = y / n, r.z_ = z / n;
    return r;
  }
  // The active rotation by `angle` about `axis` (normalised).
  static Rotation about(const Vec3& axis, double angle) {
    const Vec3 n = unit(axis);
    const double s = std::sin(angle / 2);
    return from_quaternion(std::cos(angle / 2), s * n[0], s * n[1], s * n[2]);
  }
  static Rotation about(Axis axis, double angle) {
    Vec3 n{0, 0, 0};
    n[int(axis)] = 1;
    return about(n, angle);
  }
  // The exponential of the rotation vector v = θ n.
  static Rotation exp(const Vec3& v) {
    const double t = norm(v);
    const double k = t < 1e-5 ? 0.5 - t * t / 48 : std::sin(t / 2) / t;  // sin(t/2)/t
    return from_quaternion(std::cos(t / 2), k * v[0], k * v[1], k * v[2]);
  }
  // Shepperd's method. The matrix must be a rotation to rounding; anything else is refused.
  static Rotation from_matrix(const Mat3& m, double tolerance = 1e-9) {
    if (orthogonality_defect(m) > tolerance || std::fabs(determinant(m) - 1) > tolerance)
      throw std::domain_error("physis::geometry::Rotation::from_matrix: not a rotation matrix");
    const double tr = m[0][0] + m[1][1] + m[2][2];
    const std::array<double, 4> d{tr, m[0][0], m[1][1], m[2][2]};
    int k = 0;
    for (int i = 1; i < 4; i++)
      if (d[i] > d[k]) k = i;
    double w, x, y, z;
    if (k == 0) {
      const double s = 2 * std::sqrt(1 + tr);
      w = s / 4, x = (m[2][1] - m[1][2]) / s, y = (m[0][2] - m[2][0]) / s, z = (m[1][0] - m[0][1]) / s;
    } else if (k == 1) {
      const double s = 2 * std::sqrt(1 + 2 * m[0][0] - tr);
      w = (m[2][1] - m[1][2]) / s, x = s / 4, y = (m[0][1] + m[1][0]) / s, z = (m[0][2] + m[2][0]) / s;
    } else if (k == 2) {
      const double s = 2 * std::sqrt(1 + 2 * m[1][1] - tr);
      w = (m[0][2] - m[2][0]) / s, x = (m[0][1] + m[1][0]) / s, y = s / 4, z = (m[1][2] + m[2][1]) / s;
    } else {
      const double s = 2 * std::sqrt(1 + 2 * m[2][2] - tr);
      w = (m[1][0] - m[0][1]) / s, x = (m[0][2] + m[2][0]) / s, y = (m[1][2] + m[2][1]) / s, z = s / 4;
    }
    Rotation r = from_quaternion(w, x, y, z);
    return r.w_ < 0 ? r.negated() : r;
  }
  // The rotation covered by an element of SU(2) (symmetry/su2.hpp).
  static Rotation from_su2(const symmetry::SU2& u) { return from_quaternion(u.a0, -u.a1, -u.a2, -u.a3); }
  symmetry::SU2 su2() const { return {w_, -x_, -y_, -z_}; }

  static Rotation from_euler(const EulerSequence& s, double alpha, double beta, double gamma) {
    const Rotation r1 = about(s.a1, alpha), r2 = about(s.a2, beta), r3 = about(s.a3, gamma);
    return s.convention == EulerConvention::intrinsic ? r1 * r2 * r3 : r3 * r2 * r1;
  }
  // The Euler angles (α, β, γ) of this rotation in a sequence. α, γ ∈ (−π, π]; β ∈ [−π/2, π/2]
  // for Tait–Bryan sequences and [0, π] for proper ones. At the singular β the third angle is 0.
  std::array<double, 3> euler(const EulerSequence& s) const {
    if (s.convention == EulerConvention::extrinsic) {
      const EulerSequence reversed{s.a3, s.a2, s.a1, EulerConvention::intrinsic};
      const auto a = euler(reversed);
      return {a[2], a[1], a[0]};
    }
    const Mat3 r = matrix();
    const int i = int(s.a1), j = int(s.a2);
    const double lock = 1 - 1e-12;
    if (!s.proper()) {
      const int k = int(s.a3);
      const double sigma = ((j - i + 3) % 3 == 1) ? 1.0 : -1.0;  // +1 when (i, j, k) is cyclic
      const double sb = sigma * r[i][k];
      // β by atan2, well conditioned near ±π/2 where asin is not.
      const double beta = std::atan2(sb, std::hypot(r[i][i], r[i][j]));
      if (std::fabs(sb) < lock) return {std::atan2(-sigma * r[j][k], r[k][k]), beta, std::atan2(-sigma * r[i][j], r[i][i])};
      const double sign = sb > 0 ? 1.0 : -1.0;
      return {std::atan2(sign * r[j][i], r[j][j]), sign * std::numbers::pi / 2, 0.0};
    }
    const int k = 3 - i - j;
    const double sigma = ((j - i + 3) % 3 == 1) ? 1.0 : -1.0;
    const double cb = r[i][i];
    const double beta = std::atan2(std::hypot(r[i][j], r[i][k]), cb);  // well conditioned near 0 and π
    if (std::fabs(cb) < lock) return {std::atan2(r[j][i], -sigma * r[k][i]), beta, std::atan2(r[i][j], sigma * r[i][k])};
    const double sign = cb > 0 ? 1.0 : -1.0;
    return {std::atan2(-sigma * sign * r[j][k], r[j][j]), sign > 0 ? 0.0 : std::numbers::pi, 0.0};
  }

  // The active rotation matrix.
  Mat3 matrix() const {
    const double w = w_, x = x_, y = y_, z = z_;
    return {{{1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)},
             {2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)},
             {2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)}}};
  }
  // R v, by the quaternion: v + 2w (q × v) + 2 q × (q × v), q the vector part.
  Vec3 apply(const Vec3& v) const {
    const Vec3 q{x_, y_, z_};
    const Vec3 t = 2.0 * cross(q, v);
    return v + w_ * t + cross(q, t);
  }
  Rotation inverse() const {
    Rotation r;
    r.w_ = w_, r.x_ = -x_, r.y_ = -y_, r.z_ = -z_;
    return r;
  }
  // a * b: first b, then a (the composition of maps a ∘ b).
  friend Rotation operator*(const Rotation& a, const Rotation& b) {
    Rotation r;
    r.w_ = a.w_ * b.w_ - a.x_ * b.x_ - a.y_ * b.y_ - a.z_ * b.z_;
    r.x_ = a.w_ * b.x_ + a.x_ * b.w_ + a.y_ * b.z_ - a.z_ * b.y_;
    r.y_ = a.w_ * b.y_ - a.x_ * b.z_ + a.y_ * b.w_ + a.z_ * b.x_;
    r.z_ = a.w_ * b.z_ + a.x_ * b.y_ - a.y_ * b.x_ + a.z_ * b.w_;
    return r;
  }
  // The angle θ ∈ [0, π] of the rotation, accurate near 0 and π (atan2 of |q|, |w|).
  double angle() const { return 2 * std::atan2(std::sqrt(x_ * x_ + y_ * y_ + z_ * z_), std::fabs(w_)); }
  // The rotation vector θ n with θ ∈ [0, π].
  Vec3 log() const {
    const double v = std::sqrt(x_ * x_ + y_ * y_ + z_ * z_);
    if (v == 0) return {0, 0, 0};
    const double sign = w_ < 0 ? -1.0 : 1.0;
    const double theta = 2 * std::atan2(v, std::fabs(w_));
    return (sign * theta / v) * Vec3{x_, y_, z_};
  }
  // The angle of a⁻¹ b: the geodesic distance on SO(3).
  friend double distance(const Rotation& a, const Rotation& b) { return (a.inverse() * b).angle(); }
  // Shoemake's slerp on the shorter arc: a at t = 0, b at t = 1.
  friend Rotation slerp(const Rotation& a, const Rotation& b, double t) {
    Rotation d = a.inverse() * b;
    if (d.w_ < 0) d = d.negated();
    return a * exp(t * d.log());
  }

  double w() const { return w_; }
  double x() const { return x_; }
  double y() const { return y_; }
  double z() const { return z_; }
  Rotation negated() const { return from_raw(-w_, -x_, -y_, -z_); }

 private:
  static Rotation from_raw(double w, double x, double y, double z) {
    Rotation r;
    r.w_ = w, r.x_ = x, r.y_ = y, r.z_ = z;
    return r;
  }
  double w_ = 1, x_ = 0, y_ = 0, z_ = 0;
};

// Re-orthonormalise a matrix that has drifted from SO(3) (the nearest rotation in the
// Frobenius norm is the polar factor; for a small drift one Newton step of the polar
// iteration, R ← (R + R⁻ᵀ)/2, is used twice).
inline Mat3 orthonormalised(Mat3 r) {
  for (int k = 0; k < 2; k++) {
    // R⁻ᵀ = cofactor(R) / det R.
    const double det = determinant(r);
    Mat3 cof{};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        const int i1 = (i + 1) % 3, i2 = (i + 2) % 3, j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        cof[i][j] = r[i1][j1] * r[i2][j2] - r[i1][j2] * r[i2][j1];
      }
    r = 0.5 * (r + (1 / det) * cof);
  }
  return r;
}

}  // namespace physis::geometry
