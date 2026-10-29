# analysis

*Status: draft.* Harmonic analysis as computation needs it. Its first part is potential fields outside a sphere: internal expansions in spherical harmonics, evaluated at points and on grids, their spectra, and random fields with a given spectrum ([decision 0010](../development/decisions/0010-potential-fields.md)); and the covariance of such a random field between any two points, the law of a field known in part ([decision 0013](../development/decisions/0013-discrete-mechanics-and-the-law-of-a-field-known-in-part.md)).

## 1. Concepts

A field $B=-\nabla V$ with no sources outside a sphere of radius $a$ is harmonic there. It is the gradient of an internal potential
$$
V(r,\psi,\lambda)=a\sum_{n=1}^{N}\Big(\frac ar\Big)^{n+1}\sum_{m=0}^{n}\big(g_n^m\cos m\lambda+h_n^m\sin m\lambda\big)\,\breve P_n^m(\sin\psi),
$$
with $\psi$ the planetocentric latitude, $\lambda$ the east longitude, and $\breve P_n^m$ Schmidt's semi-normalised associated Legendre functions without the Condon–Shortley phase (calx `special_functions`). This is the convention of geomagnetism: Gauss, the IGRF, the WMM, and the crustal fields of the Moon and Mars. The field's components in the sphere's local frame are north $X$, east $Y$ and down $Z$. The power of degree $n$ is the Lowes–Mauersberger spectrum, $R_n(r)=(n+1)(a/r)^{2n+4}\sum_m[(g_n^m)^2+(h_n^m)^2]$, whose sum is the mean of $\lvert B\rvert^2$ over the sphere.

An isotropic Gaussian field has independent coefficients of variance $\sigma_n^2$, and its expected spectrum is $(n+1)(2n+1)\sigma_n^2$. A crustal field's spectrum is close to $A\,n^\beta\rho^{2n}$: the statistics of its sources, attenuated over their depth $(1-\rho)a$. An estimate from data is a conditional mean, so the part it leaves unresolved has the law's power less its own.

The law at two points is the covariance $\operatorname{Cov}(V(x),V(y)) = a^2\sum_n\sigma_n^2\rho^{n+1}P_n(\hat x\cdot\hat y)$, $\rho = a^2/(|x||y|)$, by Schmidt's addition theorem. That of $B$ is its mixed Hessian, by the chain rule over five sums of $P_n$, $P_n'$ and $P_n''$. A mission that conditions the field on a magnetometer's readings (tyche `spatial`) needs it at every pair of points.

## 2. Contents

| Header | Content |
|---|---|
| `potential_field.hpp` | `HarmonicPotential` (degree, radius, $g$, $h$; truncation, sums); `read_shtools` (SHTOOLS's text format with its header); `field_at` (the potential and $X$, $Y$, $Z$ at a point); `field_grid` (on a grid of latitudes and longitudes, by rows: $O(LN^2)+O(LMN)$); `lowes_spectrum` |
| `covariance.hpp` | `PotentialCovariance` (the covariance of $V$ and of $B$ at two points, from a spectrum; the five sums, their assembly by the chain rule), `TabulatedCovariance` (at one radius, the sums tabulated in the angle and read by Lagrange's cubic) |
| `random_field.hpp` | `SpectralLaw` and `fit_spectral_law` (log-linear least squares on a band of degrees; its depth); `random_potential` (an isotropic Gaussian field of a given spectrum from a counter-based stream); `unresolved_spectrum` and `unresolved_part` (the power an estimate leaves to the law, and a realisation of it) |

## 3. Accuracy

The Legendre functions are calx's, in extended range, so the expansions hold at every latitude, the poles included, to the degrees of the crustal models and beyond: the sum rule holds to $n^2\varepsilon$ up to degree 2700.

The checks (`tests/check_analysis.cpp`) are:

- a dipole against its closed form;
- an expansion of degree 30 and one of 449 against pyshtools at 48 points on and above the sphere;
- near the poles, against the same expansion evaluated in 40-digit arithmetic. There pyshtools loses some $10^{-9}$ of its horizontal components, because its derivatives divide by $\sin\theta$; physis is right to $10^{-13}$;
- the gradient by finite differences;
- the grid against the points;
- SHTOOLS's format read back;
- random fields' spectra within their sampling error;
- a law fitted exactly to its own spectrum;
- the unresolved part.

The fixtures are made by `tests/make_fixtures.py` from a closed formula for the coefficients, so that no file of coefficients is kept.

`tests/check_covariance.cpp` (6 checks):
- the covariances of $V$ and of $B$ at five pairs of points (near, one point, far, a pole, antipodes) against the sum over the coefficients, $\sum_k\sigma_k^2 b_k(x)b_k(y)^{\mathsf T}$, each $b_k$ by `field_at`, to rounding;
- at degree 720, the variance at the sphere is $\sum_n E[R_n]$, and the field's covariance is the mixed second difference of the potential's;
- the table against the exact sums from a metre to 16 km, to $4\times10^{-12}$ of the variance, and its fallback to the exact sums.

## 4. Dependencies

calx `special_functions` (Schmidt's Legendre functions), `linalg` (the normal equations of the fit), `random` (counter-based streams).

## 5. References

- Langel, R. A. (1987), The main field, in Jacobs, J. A. (ed.), *Geomagnetism*, vol. 1, Academic Press, 249–512.
- Lowes, F. J. (1966), Mean-square values on sphere of spherical harmonic vector fields, *Journal of Geophysical Research* 71, 2179.
- Chulliat, A. et al. (2025), *The US/UK World Magnetic Model for 2025–2030: Technical Report*, NOAA NCEI (eqs. 10–12).
- Jackson, A. (1994), Statistical treatment of crustal magnetization, *Geophysical Journal International* 119, 991–998.
- Voorhies, C. V., Sabaka, T. J. and Purucker, M. (2002), On magnetic spectra of Earth and Mars, *Journal of Geophysical Research* 107, 5034.
- Wieczorek, M. A. and Meschede, M. (2018), SHTools: tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems* 19, 2574–2592.

## 6. Open Questions

- **External potentials.** Only internal expansions are written. A field with sources outside (the magnetosphere's, a tidal potential) needs $(r/a)^n$ terms.
- **Grids by the FFT.** The sum over the orders on a row is a discrete Fourier transform, which calx's FFT would do in $O(M\log M)$, against $O(MN)$ now.
- **Gravity.** `gravitation` will give fully normalised coefficients ($\bar C$, $\bar S$, scaled by $GM/a$) on these functions.
- **Fourier normalisation**, the other half of the open decision on conventions, waits for the transforms.
