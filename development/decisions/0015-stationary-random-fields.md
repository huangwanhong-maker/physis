# 0015. Stationary Random Fields Read at Any Point

- **Status:** Accepted, provisionally (one covariance, the squared exponential)
- **Date:** 2026-10-30
- **Decided by:** the user, who asked that a simulated world presume the physical properties of a place where no data give them; the infrastructure's [ADR 0059](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0059-robots-in-their-environments.md). Claude wrote this module and record.

## 1. Context

A property presumed at a place varies over it: the friction of regolith from one step to the next. A simulated world needs a realisation of that variation, the truth a robot's feet press on, at any point of a body's surface, at a metre's scale, the same on every machine. `analysis` held random potential fields on a sphere by their spectrum (decision 0010), which reach a metre on the Moon only at degree 10⁷.

## 2. Decision

`analysis/stationary_field.hpp`: a stationary random field on ℝᴰ with mean 0 and the covariance $C(r) = \sigma^2 e^{-|r|^2/2\ell^2}$, by the spectral representation (Shinozuka and Deodatis 1991):

$$g(x) = \sigma\sqrt{2/M}\sum_{m=1}^M \cos\big(k_m\cdot(x - x_0) + \varphi_m\big),$$

with the wave vectors drawn from the spectral density (normal, covariance $I/\ell^2$) and the phases uniform, from one of calx's counter-based streams. For every $M$ its covariance over the realisations is $C$ exactly; its law tends to the Gaussian as $M$ grows (256 by default). It is read at any point with no grid, and a seed fixes it.

## 3. Alternatives Considered

- **Spherical harmonics** (decision 0010). Rejected at a metre's scale: the degree needed is far beyond any expansion.
- **A grid and its Cholesky factor, or circulant embedding.** Rejected: a grid fixes the region and the resolution in advance; the spectral representation needs neither.
- **Other covariances** (Matérn, exponential). Not yet: their spectral densities are Student's t and Cauchy laws, drawn the same way, when a property needs them.

## 4. Consequences

- `check_physis_stationary` (9): over 4000 realisations, the mean, the variance and the covariance at ℓ and 3ℓ along two axes within four standard errors; the kurtosis near 3; the seed and the stream; the law kept 1.7 × 10⁶ m from the coordinates' origin.
- The planetary-science domain's regolith and the robotics domain's checks use it.
