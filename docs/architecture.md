# Architecture

*Status: draft, 2026-10-21.* This document describes how physis is structured. The decisions behind it are recorded in [development/decisions/](../development/decisions/).

## 1. Principles

### 1.1 General Laws, Not Particular Bodies

physis knows potentials, fluids, radiation and geometry. It does not know about WGS 84, the Moon's orbit or ERA5. Parameters of particular bodies arrive as inputs from the caller, which is usually a domain of The Science Infrastructures.

### 1.2 Where Structures Meet Numbers

physis is the meeting point of the project's two kinds of mathematics:

- the exact structures of morphe;
- the numerical methods of calx.

Methods that need both live here. Examples:

- **Discrete exterior calculus**, which uses complexes from morphe and sparse linear algebra from calx.
- **Lie-group integrators**, which use groups and ODE solvers.

### 1.3 Explicit Units

Physical dimensions form the free abelian group ℤ⁷ over the seven SI base units. A quantity carries its dimension, and an operation is defined only when the dimensions agree. The representation is an open decision; it may use morphe's abelian groups and calx's typed IR annotations.

### 1.4 Intrinsic Geometry, Explicit Charts

- Geometric objects are defined intrinsically wherever that is practical.
- Coordinates are charts on them.
- Every coordinate states its chart or frame, and transformations are explicit.

### 1.5 Theory with Code

- Every implemented law cites its derivation.
- Theory notes for what physics implements live in `docs/theory/`, which is created with the first note.

## 2. Modules and Their Dependencies

| Module | Depends on (within physis) |
|---|---|
| `units` | — |
| `geometry` | `units` |
| `computational_geometry` | `geometry` |
| `analysis` | `units` |
| `dynamics` | `units` |
| `control` | `dynamics` |
| `quantum` | `geometry`, `dynamics` |
| `statistical_mechanics` | `thermodynamics`, `dynamics`, `quantum` |
| `discretisation` | `geometry`, `analysis`, `dynamics` |
| `mechanics` | `geometry`, `discretisation`, `dynamics` |
| `gravitation` | `geometry`, `mechanics` |
| `thermodynamics` | `units` |
| `continuum` | `geometry`, `discretisation`, `thermodynamics` |
| `chemistry` | `thermodynamics`, `continuum` |
| `radiation` | `geometry`, `thermodynamics` |

Every module may also use calx and morphe. This table is provisional, and is revised as modules are written.

## 3. Dependencies on Sibling Libraries

```
physis ──► calx     numerical methods, computation stack
   └─────► morphe   structures: groups, complexes, exact arithmetic
```

- **Inside the infrastructure.** All three libraries are siblings under `foundations/`, and physis uses the existing calx and morphe build targets.
- **Built on its own.** physis fetches calx and morphe at pinned versions.
- **Never nest submodules.** Nesting would produce two copies of the same library.

## 4. Module Layout

```
<module>/
├── README.md                 purpose, status, dependencies, references, open questions
├── CMakeLists.txt            the target physis::<module>, and its checks
├── include/physis/<module>/  public headers
├── src/                      implementation
├── tests/                    unit tests
└── bindings/python/          Python bindings (optional)
```
