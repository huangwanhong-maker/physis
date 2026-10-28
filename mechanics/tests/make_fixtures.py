"""Reference values for physis's two-body problem: Vallado's own RV2COE.

Vallado's rv2coe (the code accompanying Fundamentals of Astrodynamics and Applications, as python-sgp4
carries it in sgp4.ext, MIT licence) gives the classical elements of a state. This script evaluates it
at the state of Vallado's Example 2-5 and at seeded random states of eccentric, inclined orbits (the
kind whose six elements are all defined), and writes tests/fixtures/rv2coe-reference.txt, which the
C++ checks read. It needs sgp4:

    python foundations/physis/mechanics/tests/make_fixtures.py
"""

from __future__ import annotations

import datetime
import math
from pathlib import Path

import numpy as np
import sgp4
from sgp4.ext import rv2coe

OUT = Path(__file__).resolve().parent / "fixtures" / "rv2coe-reference.txt"
MU = 398600.4418  # km³/s², as Vallado uses


def main() -> None:
    rng = np.random.default_rng(20261028)
    states = [((6524.834, 6862.875, 6448.296), (4.901327, 5.533756, -1.976341))]  # Vallado's Example 2-5
    while len(states) < 120:
        r = rng.normal(0, 1, 3)
        r = r / np.linalg.norm(r) * rng.uniform(6600, 45000)
        speed = math.sqrt(MU / np.linalg.norm(r)) * rng.uniform(0.75, 1.3)
        v = rng.normal(0, 1, 3)
        v = v / np.linalg.norm(v) * speed
        states.append((tuple(r), tuple(v)))
    rows = [f"# sgp4 {sgp4.__version__}: sgp4.ext.rv2coe (Vallado's RV2COE), mu = {MU} km^3/s^2; computed {datetime.date.today().isoformat()}",
            "# by foundations/physis/mechanics/tests/make_fixtures.py. rx ry rz (km) vx vy vz (km/s) p a (km) e i raan argp nu (rad)"]
    for r, v in states:
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(list(r), list(v), MU)
        if ecc < 1e-6 or incl < 1e-6 or abs(incl - math.pi) < 1e-6 or ecc >= 1:
            continue
        rows.append(" ".join(repr(float(x)) for x in (*r, *v, p, a, ecc, incl, omega, argp, nu)))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {OUT}: {len(rows) - 2} states")


if __name__ == "__main__":
    main()
