# CLAUDE.md

This file guides Claude Code in the physis repository.

physis is one of the foundation libraries of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures), which includes it as the git submodule `foundations/physis`. It is a project in its own right and must make sense without the infrastructure.

## 1. Purpose and Scope

physis implements general physical laws, and the geometry they are written in:

- units and constants;
- geometry, and computational geometry;
- analysis;
- dynamical systems and control;
- quantum theory and statistical mechanics;
- symmetry: Lie groups and their representations;
- spacetime, continuous and discrete;
- field theory, gauge theory and quantum gravity (spin foams, group field theory);
- electromagnetism;
- structure-preserving discretisation;
- mechanics;
- gravitation;
- thermodynamics;
- continuum mechanics, including geophysical fluid dynamics;
- chemistry;
- radiation.

It builds on two sibling libraries:

- **morphe**, for structures such as groups and cell complexes;
- **calx**, for numerical methods and the computation stack.

Out of scope:

- **Particular bodies and datasets.** These belong to the domains of The Science Infrastructures. Examples: the Earth's ellipsoid parameters, lunar ephemerides, seawater formulations such as TEOS-10, and astronomical time scales.
- **Exact algebraic and topological structures.** These belong to morphe.
- **Domain-free numerics and devices.** These belong to calx.

## 2. Invariants

1. **Dependencies.** physis depends only on calx and morphe among project libraries, and never on the infrastructure's domains or applications.
2. **No particular bodies.** Parameters of particular bodies are inputs, never built-in defaults.
3. **Explicit units.**
   - Every quantity at a public interface carries its units.
   - Internally, use SI units and radians.
4. **Explicit charts and frames.** Coordinates are always relative to an explicit chart or frame. Transformations between them are explicit functions.
5. **Cited sources.** Every implemented law, equation or scheme cites its derivation or source, with equation numbers.

## 3. Design Philosophy

physis shares its governing principle with its sibling repositories: **modular, extensible, revisable**.

- **Status vocabulary:**
  - `planned`: named, no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, with a pointer to its replacement.
- **Intrinsic formulations.** Prefer intrinsic, coordinate-free formulations where they are practical.
- **Structure-preserving schemes.** Prefer schemes that preserve the structure of the continuous problem, such as conservation laws, symplecticity and d∘d = 0.

## 4. Repository Layout

```
physis/
├── js/            a provisional JavaScript stand-in (ADR 0003)
├── docs/          architecture
└── development/   decisions · ideas · release notes
```

Each module directory is created when its first content is written. The modules are `units/`, `symmetry/`, `geometry/`, `computational_geometry/`, `analysis/`, `dynamics/`, `control/`, `quantum/`, `statistical_mechanics/`, `discretisation/`, `spacetime/`, `field_theory/`, `gauge/`, `quantum_gravity/`, `electromagnetism/`, `mechanics/`, `gravitation/`, `thermodynamics/`, `continuum/`, `chemistry/` and `radiation/` (ADR 0004 added six, decision 0005 six more). Until a module has content, it is listed as `planned` in the README. Written so far: `units/`, `dynamics/`, `quantum/`, `symmetry/`, `geometry/` (in part), `spacetime/`, `field_theory/`, `gauge/` and `quantum_gravity/`, each with `CMakeLists.txt`, `include/physis/<module>/` and `tests/`; the root `CMakeLists.txt` adds them and, inside the infrastructure, the siblings calx and morphe. Conventions are in `docs/glossary.md`.

## 5. Working Procedures

### 5.1 Module Contract

Each module has:

- a `README.md` giving its purpose, status, dependencies, references and open questions;
- public headers in `include/physis/<module>/`;
- the implementation in `src/`;
- unit tests in `tests/`;
- optionally, Python bindings in `bindings/python/`.

Module dependencies follow [architecture §2](docs/architecture.md).

### 5.2 Decisions and Ideas

- Decisions are numbered records in `development/decisions/`. Records are never rewritten; a later record supersedes an earlier one.
- Ideas are files named `development/ideas/YYYY-MM-DD-slug.md`, with status `seed`, `promoted` or `abandoned`.
- **Release notes.** At each commit, add to `development/release-notes/YYYY-MM-DD.md` (one file per day, one short section per commit) a concise summary in plain prose of what is new, what one can now do and what works better, not a list of changes (the git log keeps those); not between commits. This convention is shared with the sibling repositories.

### 5.3 Change Discipline

- Before restructuring or changing an invariant (§2), propose the change and wait for the user's agreement.
- Keep each change to one concept.

### 5.4 Shared Conventions

calx, morphe, physis and the infrastructure keep these conventions identical:

- writing style;
- status vocabulary;
- decision, idea and release-note conventions.

Change them in every repository of [The Science Infrastructures](https://github.com/huangwanhong-maker/The-Science-Infrastructures) together: the foundation libraries and the infrastructure's repository.

### 5.5 Working Inside the Infrastructure

When you edit physis inside the infrastructure's `foundations/physis`:

1. commit in physis first;
2. then commit the updated pointer in the infrastructure.

When a change needs calx or morphe to change as well, commit them in dependency order: calx and morphe first, then physis.

## 6. Physical and Mathematical Conventions

- **Units.** Use SI internally and radians for angles. Convert only at the input and output boundaries.
- **Sign and orientation conventions.** State each of these once in `docs/glossary.md`, which is created with the first module, and follow it everywhere:
  - handedness;
  - active versus passive rotations;
  - metric signature;
  - Fourier transform normalisation;
  - spherical-harmonic normalisation and the Condon–Shortley phase.
- **Rigour.** Write notes with the precision of a mathematics text, using LaTeX for mathematics.
- **Writing style.**
  - Headings are declarative, numbered and academic, never questions.
  - Concepts come before structure.
  - Label asides with *Rationale.*, *Remark.* or *Caveat.*
  - Use English throughout.

## 7. Technology

- **Language: C++.** This follows ADR 0003 of The Science Infrastructures and decision 0002 here.
- **Provisional:** C++20, CMake 3.20 or later (in use since ADR 0004; checks are programs on calx's `verification`, registered with CTest), and nanobind for Python bindings. On the development machine the compiler is GCC in WSL.
- **Dependency resolution.** In standalone builds, calx and morphe are fetched through the build system. Inside the infrastructure, physis uses the build targets of `foundations/calx` and `foundations/morphe` that already exist. Never nest submodules.
- **Licences.** physis is MIT-licensed. Check every dependency before adopting it.

## 8. Open Decisions

1. **How to represent units.** Run-time dimension vectors are in use since ADR 0004, checked where quantities are read; still open:
   - compile-time C++ types beside them;
   - The relation to morphe's free abelian groups.
   - The relation to calx's typed IR annotations.
2. **Sign and normalisation conventions** (§6). The signature, curvature, SU(2), recoupling and lattice conventions are fixed in `docs/glossary.md` (decision 0005); those of Fourier transforms and spherical harmonics wait for `analysis`.
3. **The first modules.** The infrastructure's first vertical slice needs:
   - `units`;
   - rotations in `geometry`;
   - spherical-harmonic potential theory in `gravitation`.
4. **The scope of general relativity.** Answered by decision 0005: a general metric toolkit (`geometry`, `spacetime`); post-Newtonian dynamics and numerical relativity stay planned in `gravitation`.

## 9. Current State

As of 2026-10-23, design documents exist; a provisional JavaScript stand-in for `discretisation` (the spectral transform on the sphere, and values on a longitude-latitude grid), with its tests ([ADR 0003](development/decisions/0003-provisional-javascript-stand-in.md)); the first C++, `units`, `dynamics` and `quantum` ([ADR 0004](development/decisions/0004-the-continuum-units-and-dynamics.md)); and the theoretical physics of [decision 0005](development/decisions/0005-symmetry-spacetime-fields-and-gauge.md) on morphe's complexes: `symmetry`, `geometry` (metrics, curvature, geodesics), `spacetime`, `field_theory`, `gauge` and `quantum_gravity`. With [decision 0006](development/decisions/0006-the-generation-of-spacetime-by-field-theories.md), the generation of spacetime by field theories: invariant theories of every rank and their large-N solutions (`field_theory`), group field theories, spin foams, random geometries and condensates (`quantum_gravity`), emergent dimensions (`spacetime`). In all, 19 check programs and 316 checks. The repository is published at `https://github.com/huangwanhong-maker/physis`.
