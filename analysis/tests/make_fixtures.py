"""Reference values for physis's analysis checks, from pyshtools (Wieczorek and Meschede 2018).

The coefficients are a closed formula, which tests/check_analysis.cpp computes in the same way, so
that no file of coefficients needs to be kept: to degree N = 449, Schmidt semi-normalised without
the Condon–Shortley phase, on a sphere of radius a = 1 737 400 m,

    g(n, m) = 1e-9 sin(1.3 n + 0.7 m + 0.1) / (n + 1)^1.2,
    h(n, m) = 1e-9 cos(0.9 n − 1.7 m) / (n + 1)^1.2   (h(n, 0) = 0),     g(0, 0) = 0,

in tesla. pyshtools evaluates the field (B_r, B_θ, B_φ) at points on and above the sphere, and the
Lowes–Mauersberger spectrum.

Within a hundredth of a degree of the poles pyshtools's horizontal components lose some 10⁻⁹ of
their value (its derivatives divide by sin θ, about 10⁻⁴ there), so the reference at those points
is exact: the degree-30 expansion evaluated in 40-digit arithmetic (mpmath), the Legendre
polynomials as exact rationals differentiated m times, the gradient by mpmath's differentiation.
The check compares pyshtools there by B_r alone. (Established on 2026-10-28: at the six polar
points pyshtools differs from the exact value by what physis's implementation differs from
pyshtools, to 10⁻¹⁴.)

    python -I foundations/physis/analysis/tests/make_fixtures.py        (needs pyshtools 4.13 or later and mpmath)
"""

import datetime
import math
from pathlib import Path

from fractions import Fraction

import mpmath as mp
import numpy as np
import pyshtools as pysh

N, A = 449, 1737400.0
OUT = Path(__file__).resolve().parent / "fixtures" / "potential-reference.txt"


def coefficients(n_max: int) -> np.ndarray:
    c = np.zeros((2, n_max + 1, n_max + 1))
    for n in range(1, n_max + 1):
        for m in range(n + 1):
            c[0, n, m] = 1e-9 * math.sin(1.3 * n + 0.7 * m + 0.1) / (n + 1) ** 1.2
            if m:
                c[1, n, m] = 1e-9 * math.cos(0.9 * n - 1.7 * m) / (n + 1) ** 1.2
    return c


def exact_field(n_max: int, lat: float, lon: float, r: float) -> tuple:
    """(B_r, B_θ, B_φ) in 40-digit arithmetic, from the potential."""
    mp.mp.dps = 40
    a = mp.mpf(A)
    coef = {}
    for n in range(1, n_max + 1):
        for m in range(n + 1):
            g = mp.mpf("1e-9") * mp.sin(mp.mpf("1.3") * n + mp.mpf("0.7") * m + mp.mpf("0.1")) / mp.power(n + 1, mp.mpf("1.2"))
            h = mp.mpf("1e-9") * mp.cos(mp.mpf("0.9") * n - mp.mpf("1.7") * m) / mp.power(n + 1, mp.mpf("1.2")) if m else mp.mpf(0)
            coef[(n, m)] = (g, h)
    # Legendre polynomials as exact rationals (Bonnet's recurrence), their m-th derivatives exact.
    poly = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for n in range(1, n_max):
        c = [Fraction(0)] * (n + 2)
        for k, v in enumerate(poly[n]):
            c[k + 1] += Fraction(2 * n + 1, n + 1) * v
        for k, v in enumerate(poly[n - 1]):
            c[k] -= Fraction(n, n + 1) * v
        poly.append(c)
    der = {}
    for n in range(n_max + 1):
        d = poly[n]
        for m in range(n + 1):
            der[(n, m)] = d
            d = [k * d[k] for k in range(1, len(d))]

    def schmidt(n, m, x):
        v = mp.mpf(0)
        for q in reversed(der[(n, m)]):
            v = v * x + mp.mpf(q.numerator) / q.denominator
        return (1 - x * x) ** (mp.mpf(m) / 2) * v * (mp.sqrt(2 * mp.factorial(n - m) / mp.factorial(n + m)) if m else 1)

    def potential(rr, psi, lam):
        x, total = mp.sin(psi), mp.mpf(0)
        for n in range(1, n_max + 1):
            t = sum((coef[(n, m)][0] * mp.cos(m * lam) + coef[(n, m)][1] * mp.sin(m * lam)) * schmidt(n, m, x) for m in range(n + 1))
            total += a * (a / rr) ** (n + 1) * t
        return total

    psi, lam, rr = mp.radians(mp.mpf(lat)), mp.radians(mp.mpf(lon)), mp.mpf(r)
    br = -mp.diff(lambda x: potential(x, psi, lam), rr)
    bt = mp.diff(lambda x: potential(rr, x, lam), psi) / rr
    bp = -mp.diff(lambda x: potential(rr, psi, x), lam) / (rr * mp.cos(psi))
    return float(br), float(bt), float(bp)


def main() -> None:
    lines = [
        f"# The field of the closed-formula expansion (make_fixtures.py), from pyshtools {pysh.__version__}, {datetime.date.today().isoformat()}.",
        "# field: degree latitude_deg longitude_deg radius_m B_r B_theta B_phi (tesla)",
        "# spectrum: degree R_n (tesla², Lowes–Mauersberger at the reference sphere)",
        "# exact: degree latitude_deg longitude_deg radius_m B_r B_theta B_phi, in 40-digit arithmetic (mpmath), near the poles",
    ]
    for n_max in (30, N):
        c = pysh.SHMagCoeffs.from_array(coefficients(n_max), r0=A, normalization="schmidt", csphase=1, units="T")
        points = [(89.99, 10.0), (-89.995, 200.0), (80.0, -45.0), (45.0, 123.0), (0.0, 0.0), (-30.0, 300.0), (-75.5, 60.25), (7.5, -59.0)]
        for lat, lon in points:
            for r in (A, A + 30e3, 3 * A):
                b = c.expand(lat=lat, lon=lon, r=r)
                lines.append(f"field {n_max} {lat:.6f} {lon:.6f} {r:.3f} {b[0]:.17e} {b[1]:.17e} {b[2]:.17e}")
        if n_max == N:
            s = c.spectrum(function="total", unit="per_l")
            for n in (1, 2, 10, 100, 300, 449):
                lines.append(f"spectrum {n} {s[n]:.17e}")
    for lat, lon in [(89.99, 10.0), (-89.995, 200.0)]:
        for r in (A, A + 30e3, 3 * A):
            b = exact_field(30, lat, lon, r)
            lines.append(f"exact 30 {lat:.6f} {lon:.6f} {r:.3f} {b[0]:.17e} {b[1]:.17e} {b[2]:.17e}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
