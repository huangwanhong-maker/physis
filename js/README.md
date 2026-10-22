# js: A Provisional JavaScript Stand-In

*Status: draft (2026-10-22).* See [ADR 0003](../development/decisions/0003-provisional-javascript-stand-in.md).

## 1. Purpose

These ES modules stand in for physis's planned C++ modules wherever programs run in a browser today. They are tested, and their numbers are relied on; they are provisional because the C++ is to replace them.

## 2. Content

| File | Stands in for | Content |
|---|---|---|
| `spectral.js` | `discretisation` | Gauss–Legendre quadrature; the fast Fourier transform; spherical harmonics on a Gaussian grid for triangular truncation T_M, with analysis, synthesis, divergence, curl, the wind of a stream function and the inverse Laplacian. The radius of the sphere is an argument, never a default. |
| `grid.js` | `discretisation` | Values on a regular longitude-latitude grid: the position of a place, bilinear interpolation that renormalises around missing values, the bracketing times of an instant, the area-weighted mean |
| `index.js` | | The public interface: everything above, and nothing else |

Sources are cited in each file: Durran (2010), Hack and Jakob (1992), Press et al. (2007).

## 3. Tests

```
node js/tests/run-tests.mjs
```

The checks need nothing outside physis: Gauss–Legendre weights and moments, the normalisation of the associated Legendre functions, the round trip of synthesis and analysis, the wind and vorticity of a solid-body rotation, the curl and divergence of a random non-divergent flow, the fast Fourier transform against the direct sum, and bilinear interpolation of linear fields.
