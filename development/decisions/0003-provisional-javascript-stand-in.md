# 0003. A Provisional JavaScript Stand-In

- **Status:** Accepted, provisionally
- **Date:** 2026-10-22
- **Decided by:** the user agreed that the domains of The Science Infrastructures become independent libraries, with the general numerics they need moved down into the foundation libraries ([ADR 0026 of the infrastructure](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0026-domains-as-libraries.md)). Claude wrote this record.

## 1. Context

The infrastructure's web applications compute in the browser, in provisional JavaScript, until the C++ libraries exist and can be compiled for it. Two of their modules are general physics and numerics, not knowledge of any particular body:

- the spectral transform method on the sphere (spherical harmonics, Gaussian quadrature, the fast Fourier transform), on which the infrastructure's meteorology builds its barotropic and shallow-water models;
- values on a regular longitude-latitude grid: where a place falls, interpolation in space and time, and the area-weighted mean.

Both declared themselves stand-ins for physis's `discretisation` module. Left inside a domain, they would make that domain the owner of general numerics; left inside an application, the domains could not use them.

## 2. Decision

1. **A `js/` directory** holds provisional ES modules that stand in for physis's planned modules. Its public interface is `js/index.js`; code outside physis imports nothing else.
2. **The first content** is `spectral.js` and `grid.js`, for `discretisation`.
3. **No particular bodies.** The radius of the sphere is an argument of every operation that depends on it, with no default; the Earth's radius and rotation rate, as Williamson et al. (1992) give them, belong to the infrastructure's meteorology. A call without a radius is refused.
4. **Charts in degrees.** The grid functions take longitude and latitude in degrees, the chart in which gridded data are written. This is an input and output boundary, where §6 of CLAUDE.md allows the conversion.
5. **Tests** are in `js/tests/run-tests.mjs`, run with Node, with no dependency outside physis.
6. **Replacement.** When the C++ module exists, it must pass the same checks. The JavaScript is then marked `deprecated`, or becomes the C++ compiled for the browser.

## 3. Alternatives Considered

- **Keep the modules in the infrastructure's meteorology.** Its models would own the spectral transform, which an ocean model or a model of another planet needs equally.
- **Wait for the C++.** The applications run today, and their numbers are checked today; the stand-in keeps that, with its tests, until the C++ takes over.

## 4. Consequences

- physis has its first code, in a language other than its decided one, marked provisional.
- The checks that were the Earth Atlas's for these modules are physis's now, with one more: that an operation without a radius is refused.
- *Caveat.* The stand-in follows the planned C++ interfaces only loosely. When the C++ interfaces are designed, the JavaScript follows them, not the other way round.
