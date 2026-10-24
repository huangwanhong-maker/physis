# computational_geometry: Computational Geometry

*Status: draft (2026-10-24): points on the sphere, with their spherical Delaunay triangulation, great-circle angles and spherical areas ([decision 0008](../development/decisions/0008-cochains-hodge-and-the-sphere.md)).* Convex hulls and triangulations in the plane and in space, Voronoi diagrams, meshes, spatial indices and robust predicates are planned.

## 1. Purpose

Geometric constructions that the sciences share. The first is a triangulation of places on the sphere, which gives a network of places its faces: the economics domain's map, and later any field on scattered points of the Earth or the sky.

## 2. Concepts

- **The spherical Delaunay triangulation** of points on S² is the boundary of their convex hull in ℝ³ (Brown 1979). When the points surround the centre, its triangles tile the sphere (F = 2V − 4, E = 3V − 6) and each circumcircle is empty.
- **The incremental hull.** A first, well-conditioned tetrahedron, then each point removes the triangles it sees and joins itself to their horizon (Clarkson and Shor 1989; de Berg et al. 2008, ch. 11). A deterministic perturbation of 10⁻⁹ resolves degenerate configurations in one valid way.
- **Spherical measures.** The great-circle angle by atan2 of |a × b| and a·b, accurate near 0 and π. The signed area of a triangle by Van Oosterom and Strackee's formula.

## 3. Content

| Header | Content |
|---|---|
| `sphere.hpp` | `Vec3`, `unit_vector(latitude, longitude)`, `dot`, `cross`, `angle`, `spherical_area`; `spherical_delaunay(points)` → `SphericalTriangulation{points, point_of, triangles, edges, covers_sphere}` |

## 4. Checks

`tests/check_sphere.cpp` (11 checks):

- a quarter of the equator is π/2, and the octant has area π/2;
- the octahedron's 8 faces;
- the cube's 12 triangles, a closed sphere despite four points on each face's circle;
- on 500 random points: F = 2V − 4, E = 3V − 6, every edge on two triangles, every triangle anticlockwise and Delaunay, the areas adding up to 4π;
- a hemisphere does not cover the sphere;
- a duplicate point is merged.

## 5. Dependencies

None in its headers. In the checks, calx `random` and `verification`.

## 6. References

- Brown, K. Q. (1979), Voronoi diagrams from convex hulls, *Inf. Process. Lett.* 9, 223–228.
- Clarkson, K. L. and Shor, P. W. (1989), Applications of random sampling in computational geometry, II, *Discrete Comput. Geom.* 4, 387–421.
- de Berg, M., Cheong, O., van Kreveld, M. and Overmars, M. (2008), *Computational Geometry*, 3rd ed., Springer, ch. 11.
- Renka, R. J. (1997), Algorithm 772: STRIPACK, *ACM TOMS* 23, 416–434.
- Van Oosterom, A. and Strackee, J. (1983), The solid angle of a plane triangle, *IEEE Trans. Biomed. Eng.* 30, 125–126.
- Edelsbrunner, H. and Mücke, E. P. (1990), Simulation of simplicity, *ACM Trans. Graph.* 9, 66–104.

## 7. Open Questions

1. **Exact predicates** (Shewchuk's adaptive arithmetic) in place of the perturbation, when inputs come from exact data.
2. **Planar and three-dimensional Delaunay**, and Voronoi cells as dual complexes (morphe).
