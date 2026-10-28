// Checks of rotations, the kinematics of moving frames, and the frame graph
// (rotation.hpp, kinematics.hpp, frame_graph.hpp).
//
// - Rotations: quaternions against their matrices, composition, inverses, the double cover of
//   physis symmetry, Shepperd's method near θ = π, the exponential and logarithm, slerp.
// - Euler angles: all twelve sequences, intrinsic and extrinsic, round trips at random angles
//   and at the singular angles; the aerospace convention on a known attitude.
// - Kinematics: the velocity and acceleration laws of a moving frame against finite
//   differences of a moving point's coordinates; composition and inverses with derivatives.
// - The frame graph: paths through the common ancestor, transformations along them against
//   compositions by hand, and the holonomy of a loop, zero for a consistent second route and
//   equal to the disagreement planted in an inconsistent one.

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/geometry/frame_graph.hpp>
#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>
#include <physis/symmetry/su2.hpp>

using namespace physis::geometry;
using calx::verification::Suite;

namespace {

constexpr double pi = std::numbers::pi;

Rotation random_rotation(calx::random::Stream& z) {
  return Rotation::from_quaternion(z.normal(), z.normal(), z.normal(), z.normal());
}
Vec3 random_vector(calx::random::Stream& z, double scale = 1) { return {scale * z.normal(), scale * z.normal(), scale * z.normal()}; }
double max_abs(const Vec3& v) { return std::max({std::fabs(v[0]), std::fabs(v[1]), std::fabs(v[2])}); }

// A frame moving in a known way: R(t) = R_z(ω t) R_x(φ0 + κ t), b(t) = b0 + u t + ½ c t².
struct Moving {
  double omega = 0.3, phi0 = 0.4, kappa = -0.17;
  Vec3 b0{1, -2, 0.5}, u{0.2, 0.1, -0.3}, c{0.05, -0.02, 0.01};

  Transform at(double t) const {
    const Mat3 r = rotation_matrix(Axis::z, omega * t) * rotation_matrix(Axis::x, phi0 + kappa * t);
    return {r, b0 + t * u + (0.5 * t * t) * c};
  }
  KinematicTransform kinematic(double t) const {
    const double a = omega * t, p = phi0 + kappa * t;
    const Mat3 z = rotation_matrix(Axis::z, a), x = rotation_matrix(Axis::x, p);
    // d/dθ R_axis(θ) = [e_axis]× R_axis(θ); second derivative [e]×² R.
    const Mat3 ez = skew({0, 0, 1}), ex = skew({1, 0, 0});
    const Mat3 z1 = omega * (ez * z), x1 = kappa * (ex * x);
    const Mat3 z2 = (omega * omega) * (ez * (ez * z)), x2 = (kappa * kappa) * (ex * (ex * x));
    KinematicTransform k;
    k.r = z * x;
    k.r_dot = z1 * x + z * x1;
    k.r_ddot = z2 * x + 2.0 * (z1 * x1) + z * x2;
    k.b = b0 + t * u + (0.5 * t * t) * c;
    k.b_dot = u + t * c;
    k.b_ddot = c;
    return k;
  }
};

}  // namespace

int main() {
  Suite s;
  calx::random::Stream z(20261028, 1, 2, 3);

  // ---- rotations ------------------------------------------------------------------------------
  {
    double matrix = 0, apply = 0, compose = 0, inverse = 0, cover = 0, round = 0, logexp = 0, defect = 0;
    for (int i = 0; i < 500; i++) {
      const Rotation a = random_rotation(z), b = random_rotation(z);
      const Vec3 v = random_vector(z);
      defect = std::max({defect, orthogonality_defect(a.matrix()), std::fabs(determinant(a.matrix()) - 1)});
      apply = std::max(apply, max_abs(a.apply(v) - a.matrix() * v));
      compose = std::max(compose, max_difference((a * b).matrix(), a.matrix() * b.matrix()));
      inverse = std::max(inverse, max_difference((a * a.inverse()).matrix(), identity3()));
      cover = std::max(cover, max_difference(physis::symmetry::adjoint(a.su2()), a.matrix()));
      cover = std::max(cover, distance(Rotation::from_su2(a.su2()), a));
      round = std::max(round, distance(Rotation::from_matrix(a.matrix()), a));
      logexp = std::max(logexp, distance(Rotation::exp(a.log()), a));
      matrix = std::max(matrix, std::fabs(a.angle() - norm(a.log())));
    }
    s.check("rotation matrices are orthogonal with determinant 1 (500 random)", defect, 1e-14);
    s.check("the quaternion acts as its matrix", apply, 1e-14);
    s.check("composition is the product of matrices", compose, 1e-14);
    s.check("a rotation times its inverse is the identity", inverse, 1e-14);
    s.check("the double cover: symmetry's adjoint of su2() is the matrix", cover, 1e-14);
    s.check("Shepperd's method recovers the rotation from its matrix", round, 1e-7, "rad");
    s.check("exp(log R) = R", logexp, 1e-7, "rad");
    s.check("the angle is the length of the rotation vector", matrix, 1e-14, "rad");
  }
  {
    // Shepperd's method at and near half turns, where the trace is −1.
    double worst = 0;
    for (const double theta : {pi, pi - 1e-9, pi - 1e-5, 1e-9, 0.0}) {
      for (const Vec3& n : std::vector<Vec3>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 2, 3}, {-3, 1, 0.2}}) {
        const Rotation r = Rotation::about(n, theta);
        worst = std::max(worst, max_difference(Rotation::from_matrix(r.matrix()).matrix(), r.matrix()));
      }
    }
    s.check("Shepperd's method at and near θ = π and 0", worst, 1e-14);
  }
  {
    const Rotation r = Rotation::about(Axis::z, pi / 2);
    s.check("the active rotation by 90° about z takes x to y", max_abs(r.apply({1, 0, 0}) - Vec3{0, 1, 0}), 1e-15);
    s.check("frame_rotation is SOFA's R_z: (cos φ, −sin φ, 0) for the x axis",
            max_abs(frame_rotation(Axis::z, 0.3) * Vec3{1, 0, 0} - Vec3{std::cos(0.3), -std::sin(0.3), 0}), 1e-16);
    const Rotation a = Rotation::about({1, 1, 0}, 0.2), b = Rotation::about({0, -1, 2}, 1.4);
    const Rotation m = slerp(a, b, 0.5);
    s.check("slerp's midpoint is halfway along the geodesic", distance(a, m) - distance(m, b), 1e-14, "rad");
    s.check("slerp ends at its endpoints", distance(slerp(a, b, 1.0), b) + distance(slerp(a, b, 0.0), a), 1e-7, "rad");
    Mat3 drifted = Rotation::about({1, 2, 3}, 0.7).matrix();
    drifted[0][1] += 3e-7;
    drifted[2][0] -= 2e-7;
    s.check("re-orthonormalisation removes a drift of 3e-7", orthogonality_defect(orthonormalised(drifted)), 1e-15);
  }

  // ---- Euler angles ---------------------------------------------------------------------------
  {
    const char* sequences[] = {"XYZ", "YZX", "ZXY", "XZY", "ZYX", "YXZ", "XYX", "XZX", "YXY", "YZY", "ZXZ", "ZYZ"};
    double worst = 0, singular = 0, ranges = 0;
    for (const char* name : sequences) {
      for (const auto convention : {EulerConvention::intrinsic, EulerConvention::extrinsic}) {
        const EulerSequence seq = EulerSequence::parse(name, convention);
        for (int i = 0; i < 200; i++) {
          const Rotation r = random_rotation(z);
          const auto a = r.euler(seq);
          worst = std::max(worst, max_difference(Rotation::from_euler(seq, a[0], a[1], a[2]).matrix(), r.matrix()));
          const double lo = seq.proper() ? 0 : -pi / 2, hi = seq.proper() ? pi : pi / 2;
          if (a[1] < lo - 1e-15 || a[1] > hi + 1e-15) ranges = std::max(ranges, 1.0);
        }
        // The singular middle angles, where the first and third axes coincide.
        for (const double beta : seq.proper() ? std::vector<double>{0.0, pi} : std::vector<double>{pi / 2, -pi / 2}) {
          for (int i = 0; i < 20; i++) {
            const Rotation r = Rotation::from_euler(seq, z.uniform(-pi, pi), beta, z.uniform(-pi, pi));
            const auto a = r.euler(seq);
            singular = std::max(singular, max_difference(Rotation::from_euler(seq, a[0], a[1], a[2]).matrix(), r.matrix()));
          }
        }
      }
    }
    s.check("Euler angles of 12 sequences × 2 conventions × 200 rotations rebuild the rotation", worst, 1e-12);
    s.check("Euler angles at the singular middle angle rebuild the rotation", singular, 1e-12);
    s.require("middle angles lie in [−π/2, π/2] (Tait–Bryan) or [0, π] (proper)", ranges == 0);
    // Extrinsic a1 a2 a3 with (α, β, γ) is intrinsic a3 a2 a1 with (γ, β, α).
    const Rotation e = Rotation::from_euler(EulerSequence::parse("XYZ", EulerConvention::extrinsic), 0.1, 0.2, 0.3);
    const Rotation i = Rotation::from_euler(EulerSequence::parse("ZYX"), 0.3, 0.2, 0.1);
    s.check("extrinsic XYZ (α, β, γ) is intrinsic ZYX (γ, β, α)", distance(e, i), 1e-7, "rad");
    // Aerospace: in north-east-down, a yaw of 90° points the nose east; a pitch of 30° then
    // raises it 30° above the horizon (−30° in "down").
    const Rotation yaw = Rotation::from_euler(EulerSequence::parse("ZYX"), pi / 2, 0, 0);
    s.check("yaw 90°: the nose (body x) points east in NED", max_abs(yaw.apply({1, 0, 0}) - Vec3{0, 1, 0}), 1e-15);
    const Rotation climb = Rotation::from_euler(EulerSequence::parse("ZYX"), pi / 2, pi / 6, 0);
    s.check("yaw 90°, pitch 30°: the nose is 30° above the horizon",
            max_abs(climb.apply({1, 0, 0}) - Vec3{0, std::cos(pi / 6), -std::sin(pi / 6)}), 1e-15);
  }

  // ---- kinematics ------------------------------------------------------------------------------
  {
    const Moving m;
    const double t = 1.7, h = 1e-4;
    const Vec3 x0{0.3, -1.2, 2.0}, u{0.5, 0.25, -0.4}, acc{0.1, -0.3, 0.2};
    auto point = [&](double tt) { return x0 + tt * u + (0.5 * tt * tt) * acc; };  // in the source frame
    auto image = [&](double tt) { return m.at(tt).position(point(tt)); };          // in the target frame
    const KinematicTransform k = m.kinematic(t);
    const auto st = k.state(point(t), u + t * acc, acc);
    const Vec3 v_fd = (1 / (2 * h)) * (image(t + h) - image(t - h));
    const Vec3 a_fd = (1 / (h * h)) * (image(t + h) - 2.0 * image(t) + image(t - h));
    s.check("velocity in a moving frame: R v + Ṙ x + ḃ against finite differences", max_abs(st[1] - v_fd), 1e-7);
    s.check("acceleration: R a + 2 Ṙ v + R̈ x + b̈ against finite differences", max_abs(st[2] - a_fd), 1e-5);
    const StateTransform sd = differentiated([&](double dt) { return m.at(t + dt); }, 1e-4);
    s.check("Ṙ by central differences", max_difference(sd.r_dot, k.r_dot), 1e-8);
    // Ω = Ṙ Rᵀ is antisymmetric, and for R = R_z(ωt) R_x(...) its vector is ω e_z + κ R_z e_x.
    const Vec3 expected = Vec3{0, 0, m.omega} + m.kappa * (rotation_matrix(Axis::z, m.omega * t) * Vec3{1, 0, 0});
    s.check("the angular velocity of the moving frame", max_abs(k.state_transform().angular_velocity() - expected), 1e-14, "rad/s");
    // Groups: T ∘ T⁻¹ = id with vanishing derivatives; (a ∘ b) ∘ c = a ∘ (b ∘ c).
    const KinematicTransform a = m.kinematic(0.4), b = m.kinematic(-1.1), c = m.kinematic(2.3);
    const KinematicTransform e = a * a.inverse();
    s.check("T ∘ T⁻¹ = identity, with its derivatives",
            std::max({max_difference(e.r, identity3()), max_difference(e.r_dot, zero3()), max_difference(e.r_ddot, zero3()),
                      max_abs(e.b), max_abs(e.b_dot), max_abs(e.b_ddot)}),
            1e-14);
    const KinematicTransform l = (a * b) * c, r = a * (b * c);
    s.check("composition is associative, with its derivatives",
            std::max({max_difference(l.r, r.r), max_difference(l.r_dot, r.r_dot), max_difference(l.r_ddot, r.r_ddot),
                      max_abs(l.b - r.b), max_abs(l.b_dot - r.b_dot), max_abs(l.b_ddot - r.b_ddot)}),
            1e-13);
    // The state transforms of a composite act as the composite of their actions.
    const StateTransform sa = a.state_transform(), sb = b.state_transform();
    const auto twice = sa.state(sb.state(x0, u)[0], sb.state(x0, u)[1]);
    const auto once = (sa * sb).state(x0, u);
    s.check("the composite acts on states as the composition of actions", std::max(max_abs(twice[0] - once[0]), max_abs(twice[1] - once[1])), 1e-13);
  }

  // ---- the frame graph -------------------------------------------------------------------------
  {
    using Edge = FrameGraph<double>::Edge;
    const Moving m;
    FrameGraph<double> g;
    const int I = g.root("I", "the inertial root");
    const int A = g.add("A", I, Edge([&](double t) { return m.kinematic(t).state_transform(); }), "a moving frame");
    const int B = g.add("B", A, Edge([](double t) { return StateTransform::rotating(rotation_matrix(Axis::y, 0.5 * t), {0, 0.5, 0}, {0, 0, 3}); }));
    const int C = g.add("C", I, Edge([](double) { return StateTransform::fixed(Rotation::about({1, 1, 1}, 0.8).matrix(), {10, 0, 0}); }));
    const int D = g.add("D", "C", Edge([](double t) { return StateTransform::rotating(rotation_matrix(Axis::x, -0.2 * t), {-0.2, 0, 0}); }));
    const double t = 0.9;
    const std::vector<int> expected{B, A, I, C, D};
    s.require("the path from B to D runs through the root", g.path(B, D) == expected);
    s.require("a path to itself is one frame", g.path(A, A) == std::vector<int>{A});
    const StateTransform byhand = g.frames()[D].to_parent(t).inverse() * g.frames()[C].to_parent(t).inverse() *
                                  g.frames()[A].to_parent(t) * g.frames()[B].to_parent(t);
    const StateTransform bd = g.transform(B, D, t);
    s.check("the transformation from B to D along the tree, against composition by hand",
            std::max({max_difference(bd.r, byhand.r), max_difference(bd.r_dot, byhand.r_dot), max_abs(bd.b - byhand.b), max_abs(bd.b_dot - byhand.b_dot)}),
            1e-13);
    const Deviation round = deviation(g.transform(D, B, t) * bd);
    s.check("B to D and back is the identity", std::max({round.angle, round.translation, round.angular_rate, round.velocity}), 1e-13);
    // A consistent second route from B to D closes a loop with zero holonomy; one that is off
    // by a screw of 1 µrad about z and 1 mm along it is measured as such, whatever frame the
    // loop is read in (the screw's angle and slide are invariant under conjugation).
    const int consistent = g.relate("B→D by hand", B, D, Edge([&](double tt) {
      return g.frames()[D].to_parent(tt).inverse() * g.frames()[C].to_parent(tt).inverse() * g.frames()[A].to_parent(tt) * g.frames()[B].to_parent(tt);
    }));
    const int planted = g.relate("B→D off by a screw of 1 µrad and 1 mm", B, D, Edge([&](double tt) {
      return StateTransform::fixed(Rotation::about({0, 0, 1}, 1e-6).matrix(), {0, 0, 1e-3}) * g.transform(B, D, tt);
    }));
    const Deviation h0 = deviation(g.holonomy(consistent, t)), h1 = deviation(g.holonomy(planted, t));
    s.check("the holonomy of a consistent loop vanishes (angle)", h0.angle, 1e-14, "rad");
    s.check("the holonomy of a consistent loop vanishes (translation)", h0.translation, 1e-13);
    s.check("the holonomy of a planted disagreement is its angle", h1.angle - 1e-6, 1e-14, "rad");
    s.check("the holonomy of a planted disagreement is its slide", h1.slide - 1e-3, 1e-12);
    bool refused = false;
    try {
      g.require("nowhere");
    } catch (const std::invalid_argument&) {
      refused = true;
    }
    s.require("an unknown frame is refused by name", refused);
  }
  return s.done();
}
