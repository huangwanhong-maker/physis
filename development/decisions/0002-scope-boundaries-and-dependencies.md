# 0002. Scope, Boundaries and Dependencies

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user, who proposed a theoretical-physics library in its own repository. Claude drew the boundaries ([laboratory ADR 0008](https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory/blob/main/development/decisions/0008-independent-libraries.md)).

## 1. Context

The laboratory first planned `physics/` and `mathematics/geometry/` as layers inside its own core. These layers describe general laws, not the Earth or its ocean and atmosphere. Other work can reuse them, and they change on a different rhythm from the laboratory's domain code.

## 2. Decision

1. **Scope.** physis covers:
   - units and physical constants;
   - geometry as physics uses it: Euclidean, spherical and ellipsoidal geometry, rotations, differential geometry, map projections;
   - structure-preserving discretisation;
   - mechanics, gravitation, thermodynamics, continuum mechanics (including geophysical fluid dynamics), chemistry and radiation.
2. **No particular bodies.** Parameters of particular bodies are always inputs.
3. **Dependencies.** physis depends on calx and morphe, and never on the laboratory. It does not nest them as submodules. Inside the laboratory it uses the sibling copies; built on its own, it fetches them.
4. **Language.** physis is written in C++, following laboratory ADR 0003.

## 3. Alternatives Considered

- **Geometry in morphe.** Rejected. Differential geometry as physics uses it is numerical and works in coordinates. Algebraic geometry and topology stay in morphe.
- **Units in the laboratory's metrology.** Rejected. Physical laws need units at their own interfaces. Only the astronomical time scales and reference frames stay in the laboratory.
- **A separate geometry library.** Rejected as premature. Geometry can be split out later if it gains users outside physics.

## 4. Consequences

- A change that spans libraries is committed in dependency order: calx and morphe first, then physis, then the laboratory.
- The laboratory's domains (oceanography, atmosphere, tides, …) apply the laws in physis to particular bodies.
