# 0010. Potential Fields outside a Sphere, Their Spectra, and Random Fields with a Given Spectrum

- **Status:** Accepted, provisionally (internal expansions only; grids by direct sums over the orders, not yet by the FFT)
- **Date:** 2026-10-28
- **Decided by:** the user asked Cosmopoiesis for "a generative model … from historical data of moon magnet field etc. And input to virtual sensor instances", and set the rule that a domain or a work writes no general method of its own. Claude wrote the general methods here and this record. Their use is Cosmopoiesis's [decision 0005](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0005-fields-known-in-part-and-declared-instruments.md) and the infrastructure's [ADR 0052](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0052-fields-known-in-part.md).

## 1. Context

The Moon's crustal field is published as a spherical-harmonic expansion of degree 449 (T2015_449). The geophysics domain's WMM2025 evaluated its own expansion of degree 12, with Legendre functions that hold only to modest degrees. Nothing in the foundations evaluated an expansion of high degree at the poles, gave its spectrum, or drew a random field of a given spectrum. These name no object in nature, and gravity, the main field and crustal fields all need them.

## 2. Decision

### 2.1 The Module

`analysis` is written, its first part potential fields outside a sphere ([analysis/README.md](../../analysis/README.md)). It sits in `analysis` rather than `gravitation` because the methods are harmonic analysis on the sphere, shared by magnetic and gravitational potentials. `gravitation` will give a gravity field's normalisation and scaling on top of them.

### 2.2 The Convention

This fixes the spherical-harmonic half of open decision 2:

- **The functions.** Expansions use Schmidt's semi-normalised associated Legendre functions without the Condon–Shortley phase, $\breve P_n^m$, the convention of geomagnetism (Langel 1987), in calx `special_functions`.
- **The internal potential.** It is $V=a\sum_n(a/r)^{n+1}\sum_m(g_n^m\cos m\lambda+h_n^m\sin m\lambda)\breve P_n^m(\sin\psi)$, with the coefficients in the field's unit.
- **The components.** They are north, east and down in the sphere's local frame (the WMM2025 report's eqs. 10–12).
- **The spectrum.** It is Lowes's and Mauersberger's.

The fully normalised functions of gravity are $\bar P_n^m=\sqrt{2n+1}\,\breve P_n^m$. They are recorded in the glossary, and `gravitation` will use them.

### 2.3 Evaluation

`field_at` evaluates at a point. `field_grid` evaluates on a grid by rows: the functions once per latitude, the sums over the degrees once per order, so that a grid of $L$ latitudes and $M$ longitudes costs $O(LN^2)+O(LMN)$. The orders of one longitude are stored side by side, so that the inner sum runs along memory; laid out by order instead, it ran at a stride of the grid's width, which the caches serve poorly.

The functions are calx's in extended range (calx decision 0010), so the expansions hold at the poles and to the degrees of every crustal model.

### 2.4 Random Fields

A spectral law $A\,n^\beta\rho^{2n}$ is fitted by log-linear least squares on a band of degrees. Its $\rho$ gives an effective depth $(1-\rho)a$ of the sources and of an estimate's regularisation together.

An isotropic Gaussian field of a given spectrum is drawn from a counter-based stream, so that a seed gives the same field on every machine.

The part an estimate leaves unresolved has, by the law of total variance, the law's power less the estimate's, and the law's whole power above the estimate's degree. `unresolved_part` draws it.

## 3. Alternatives Considered

- **pyshtools behind the interface.** Rejected for the engines' C++: it is Python, and near the poles its horizontal components lose some $10^{-9}$ of their value. It is the reference of the checks instead, and an exact evaluation in 40-digit arithmetic arbitrates at the poles.
- **The Condon–Shortley phase**, as quantum mechanics writes it. Rejected for potentials: every published geomagnetic and crustal model omits it. Where physics needs it, `symmetry`'s recoupling keeps it.
- **A learned generator of fields.** Rejected: it would not be a potential field, and would invent where the data do not reach. The Gaussian law is the least committal choice for a given spectrum.

## 4. Consequences

- physis has 25 check programs. `check_physis_analysis` holds the module:
  - a dipole against its closed form;
  - pyshtools at degree 30 and 449;
  - exact references at the poles;
  - the gradient by finite differences;
  - the grid against the points;
  - SHTOOLS's format;
  - the statistics of random fields;
  - the law fitted exactly, and from a random field within four standard errors;
  - the unresolved part.
- An out-of-bounds read in the first `field_grid` (the powers of $a/r$ to degree $N+2$) was found and corrected while it was written.
- Open: external potentials, grids by calx's FFT, and gravity's normalisation in `gravitation`.
