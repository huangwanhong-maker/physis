// Changes of frame in Euclidean 3-space and their kinematics: the group of rigid motions
// SE(3) = SO(3) ⋉ ℝ³ as it acts on the coordinates of points, and its first derivatives in
// time, which act on velocities.
//
// A *transformation* from a source frame to a target frame maps the coordinates of one
// point in the source to its coordinates in the target:
//   x_t = R x_s + b,    R ∈ SO(3), b ∈ ℝ³ (the source's origin in target coordinates).
// When the frames move relative to each other, R and b depend on time, and the velocity of
// the point, measured in each frame, transforms with the motion of the frames
// (Goldstein, Poole and Safko 2002, Classical Mechanics, §4.9; the formalism of
// Cosmopoiesis, Proposition 2.2):
//   v_t = R v_s + Ṙ x_s + ḃ,
//   a_t = R a_s + 2 Ṙ v_s + R̈ x_s + b̈.
// With the angular velocity ω of the source relative to the target, in target coordinates,
// Ṙ = [ω]× R, and the velocity law reads v_t = R v_s + ω × (R x_s) + ḃ. A free vector (a force,
// a direction) transforms by R alone: it has no position for b or Ṙ to act on.
//
// `StateTransform` carries (R, Ṙ, b, ḃ), the "state transformation" of NAIF's SPICE toolkit
// (Frames Required Reading, §"State transformations"); `KinematicTransform` adds (R̈, b̈) for
// accelerations. Both form groups under composition, with the derivatives following the
// product rule, and both are checked against finite differences of moving frames.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include <physis/geometry/rotation.hpp>

namespace physis::geometry {

// x_t = R x_s + b, at one instant.
struct Transform {
  Mat3 r = identity3();
  Vec3 b{0, 0, 0};

  Vec3 position(const Vec3& x) const { return r * x + b; }
  Vec3 direction(const Vec3& v) const { return r * v; }
  Transform inverse() const {
    const Mat3 rt = transpose(r);
    return {rt, -(rt * b)};
  }
  // (a ∘ b)(x) = a(b(x)): first b, then a.
  friend Transform operator*(const Transform& a, const Transform& b) { return {a.r * b.r, a.r * b.b + a.b}; }
};

// The transformation with its first derivatives in time.
struct StateTransform {
  Mat3 r = identity3(), r_dot = zero3();
  Vec3 b{0, 0, 0}, b_dot{0, 0, 0};

  static StateTransform fixed(const Mat3& r, const Vec3& b = {0, 0, 0}) { return {r, zero3(), b, {0, 0, 0}}; }
  // A frame turning with angular velocity ω (target coordinates) about the target's origin.
  static StateTransform rotating(const Mat3& r, const Vec3& omega, const Vec3& b = {0, 0, 0}, const Vec3& b_dot = {0, 0, 0}) {
    return {r, skew(omega) * r, b, b_dot};
  }

  Transform at() const { return {r, b}; }
  Vec3 position(const Vec3& x) const { return r * x + b; }
  Vec3 direction(const Vec3& v) const { return r * v; }
  // (x, v) in the source to (x, v) in the target.
  std::array<Vec3, 2> state(const Vec3& x, const Vec3& v) const { return {r * x + b, r * v + r_dot * x + b_dot}; }
  // The angular velocity of the source relative to the target, in target coordinates:
  // vee(Ṙ Rᵀ).
  Vec3 angular_velocity() const { return vee(r_dot * transpose(r)); }

  StateTransform inverse() const {
    const Mat3 rt = transpose(r), rt_dot = transpose(r_dot);
    return {rt, rt_dot, -(rt * b), -(rt_dot * b + rt * b_dot)};
  }
  friend StateTransform operator*(const StateTransform& a, const StateTransform& c) {
    return {a.r * c.r, a.r_dot * c.r + a.r * c.r_dot, a.r * c.b + a.b, a.r_dot * c.b + a.r * c.b_dot + a.b_dot};
  }
};

// The transformation with its first and second derivatives in time.
struct KinematicTransform {
  Mat3 r = identity3(), r_dot = zero3(), r_ddot = zero3();
  Vec3 b{0, 0, 0}, b_dot{0, 0, 0}, b_ddot{0, 0, 0};

  StateTransform state_transform() const { return {r, r_dot, b, b_dot}; }
  // (x, v, a) in the source to (x, v, a) in the target.
  std::array<Vec3, 3> state(const Vec3& x, const Vec3& v, const Vec3& a) const {
    return {r * x + b, r * v + r_dot * x + b_dot, r * a + 2.0 * (r_dot * v) + r_ddot * x + b_ddot};
  }
  KinematicTransform inverse() const {
    // S = Rᵀ: Ṡ = Ṙᵀ, S̈ = R̈ᵀ; c = −Rᵀ b, ċ = −(Ṙᵀ b + Rᵀ ḃ), c̈ = −(R̈ᵀ b + 2 Ṙᵀ ḃ + Rᵀ b̈).
    const Mat3 s = transpose(r), s_dot = transpose(r_dot), s_ddot = transpose(r_ddot);
    return {s, s_dot, s_ddot, -(s * b), -(s_dot * b + s * b_dot), -(s_ddot * b + 2.0 * (s_dot * b_dot) + s * b_ddot)};
  }
  friend KinematicTransform operator*(const KinematicTransform& a, const KinematicTransform& c) {
    return {a.r * c.r,
            a.r_dot * c.r + a.r * c.r_dot,
            a.r_ddot * c.r + 2.0 * (a.r_dot * c.r_dot) + a.r * c.r_ddot,
            a.r * c.b + a.b,
            a.r_dot * c.b + a.r * c.b_dot + a.b_dot,
            a.r_ddot * c.b + 2.0 * (a.r_dot * c.b_dot) + a.r * c.b_ddot + a.b_ddot};
  }
};

// How far a transformation is from the identity. By Chasles's theorem every rigid motion is a
// screw: a rotation by θ about an axis and a slide d along it (Selig 2005, Geometric
// Fundamentals of Robotics, §3.3). θ and d are invariant under conjugation, so they measure the
// holonomy of a loop of frames (frame_graph.hpp) whatever frames the loop is read from; the
// translation b is the displacement of the origin of the frame the loop is read in, and the
// derivatives' norms are its rates.
struct Deviation {
  double angle = 0;        // θ, rad
  double slide = 0;        // |d|, the slide along the screw axis (|b| for a pure translation)
  double translation = 0;  // |b|
  double angular_rate = 0, velocity = 0;
};
inline Deviation deviation(const StateTransform& t) {
  // The angle of R from its trace, made accurate near 0 by the antisymmetric part.
  const Vec3 a = vee(t.r);
  const double c = std::clamp((t.r[0][0] + t.r[1][1] + t.r[2][2] - 1) / 2, -1.0, 1.0);
  const double theta = std::atan2(norm(a), c);
  double slide = norm(t.b);
  if (theta > 1e-12) {
    // The axis: the antisymmetric part away from π; near π, the symmetric part R + Rᵀ + (1 − tr R) I = 2(1 − cos θ) n nᵀ.
    Vec3 n;
    if (theta < 3) {
      n = unit(a);
    } else {
      const Mat3 m = t.r + transpose(t.r);
      int k = 0;
      for (int i = 1; i < 3; i++)
        if (m[i][i] > m[k][k]) k = i;
      n = unit(Vec3{m[0][k], m[1][k], m[2][k]});
    }
    slide = std::fabs(dot(t.b, n));
  }
  return {theta, slide, norm(t.b), norm(t.angular_velocity()), norm(t.b_dot)};
}

// The derivative of a time-dependent transformation by central differences of step h, for
// edges whose derivatives are not written out: Ṙ ≈ (R(t + h) − R(t − h)) / 2h, with an error
// of order h² |R⃛|.
template <class F>
StateTransform differentiated(F&& at, double h) {
  const Transform plus = at(h), minus = at(-h), here = at(0.0);
  return {here.r, (1 / (2 * h)) * (plus.r - minus.r), here.b, (1 / (2 * h)) * (plus.b - minus.b)};
}

}  // namespace physis::geometry
