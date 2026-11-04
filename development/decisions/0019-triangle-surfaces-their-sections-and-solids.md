# 0019. Triangle Surfaces in Space: Their Sections, the Solids They Bound, and Close Pairs of Points

- **Status:** Accepted, provisionally. Welding is not transitive. A section of a surface that is not closed is returned as open chains, which have no area.
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** The user asked for a Body Atlas, the anatomy of the whole body drawn from BodyParts3D's meshes.
  - **The rule.** General methods are written in the foundations.
  - **What Claude did.** It wrote these functions, their checks and this record.
  - **The composition.** The infrastructure's [ADR 0066](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0066-the-anatomy-domain-and-the-body-atlas.md). Its anatomy domain measures each part of the body and cuts the body in the planes of imaging.

## 1. Context

A part of the body, as BodyParts3D draws it, is a closed triangle surface in millimetres. To say what it is, the anatomy domain needs four things that name no object in nature:
- its area;
- the solid it bounds: volume, centre, inertia and principal axes;
- its sections by a plane;
- which parts touch.

physis had the mass properties of boxes and point masses, composed by the parallel axis theorem (decision 0013), and the faces of a plane graph (decision 0018). It had nothing for a surface in space. The topology of such a surface (closed or not, its genus, whether its orientation is consistent) is morphe's, through its decision 0010.

## 2. Decision

### 2.1 Triangle Surfaces (`computational_geometry/surfaces.hpp`)

- **A triangle surface.** Points and triangles, each oriented by the order of its corners.
- **Measures.** Its bounds and its area.
- **Welding.** Points within a tolerance become one. Points are hashed to cells of that side and compared with the 27 neighbouring cells. Triangles that lose a corner are dropped.

### 2.2 Sections by a Plane

The section of a surface by the plane N·(x − o) = 0 is returned as chains of crossing points.
- **The perturbation.** A corner with N·(x − o) = 0 is counted on the positive side, as if the plane moved infinitesimally toward −N. Then every triangle crosses the plane in no segment or in one segment between two of its edges.
- **Joining without tolerance.** A crossing point is named by its edge, so segments are joined into chains exactly.
- **Direction.** Each segment runs along N × n. The direction is read from the signs of the corners and their order alone. So it stays exact where the crossing points coincide, as they do when a corner lies on the plane.
- **Orientation.** For a closed surface oriented outward, the solid lies to the left of every segment seen from +N. So the loops run counter-clockwise round the pieces of the section and clockwise round its holes.
- **Area.** The section's area is the sum of the loops' shoelace areas in a basis with u × v = N (Green's theorem).
- **Cavalieri's principle.** Areas integrate over the planes' offsets to the solid's volume: ∫ A(h) dh = V. The checks hold the two together.

The slicing follows Minetto et al. (2017, *Computer-Aided Design* 92, 1).

### 2.3 The Solid a Closed Surface Bounds (`mechanics/rigid_body.hpp`, `solid_polyhedron`)

The solid is the signed sum of the tetrahedra joining a reference point to the triangles (the divergence theorem). For each tetrahedron, with its vertices relative to the reference point:
- volume: det[a b c]/6;
- first moment: V(a + b + c)/4;
- second moment: (V/20)(Σ_k p_k p_kᵀ + s sᵀ), where s = a + b + c (Tonon 2004).

The reference point is the mean of the corners, which keeps the sums small far from the origin. The result is a `MassProperties`, so solids compose with boxes and point masses by `combine`. A surface that encloses no positive volume is refused: an open surface, or one oriented inward. Its caller decides closure, with morphe.

### 2.4 Close Pairs (`computational_geometry/surfaces.hpp`, `close_pairs`)

Given points with labels, it finds the pairs of labels some of whose points lie within a distance d, with the least distance between them.
- **Method.** The points are sorted into cells of side d. Each cell is compared with itself and its 13 forward neighbours, so every pair within d is examined once.
- **What it measures.** Contiguity, not shared boundary.

## 3. Alternatives Considered

- **A borrowed geometry library (CGAL, libigl).** Rejected. Parts of CGAL are copyleft (CLAUDE.md §7). These functions are short, exact in what they claim, and are how the methods are learnt.
- **Sections joined by snapping nearby points.** Rejected. A tolerance joins what it should not on thin parts, such as vessels a millimetre wide. Naming crossing points by their edges is exact.
- **Orienting segments by their crossing points.** Rejected. It fails where a corner lies on the plane, and the first check found it so. The sign pattern decides instead.
- **A bounding-volume hierarchy for close pairs.** Rejected. A grid of cells of the distance's side is simpler and linear for the uniform density of a body's surfaces.

## 4. Consequences

- **Checks.** physis gains 29 checks:
  - `check_surfaces`, 20: welding a soup of corners; the cube's sections, among them the regular hexagon and the planes through its faces; Cavalieri's principle on a sphere of 5,120 triangles and on a torus, to 10⁻⁷; a torus's two loops of opposite signs; an open chain; close pairs against every pair compared.
  - `check_polyhedron`, 9: a box as twelve triangles is `solid_box`, also far from the origin and turned by a rotation; a regular tetrahedron has the isotropic inertia m a²/20; an inward surface is refused.
- **The anatomy domain** measures each of BodyParts3D's 2,234 parts with these functions. Its volumes agree with the database's at a median ratio of 1.0004.
- **Planned.** The volumes of open surfaces closed by caps; sections by curved surfaces.
