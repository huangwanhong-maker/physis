"""Reference values for the H3 grid, from the reference library itself (the `h3` package, on the reference C library).

Random points over the sphere, the poles and the antimeridian, at every resolution; for each its cell, the cell's
centre, boundary and parent; the twelve pentagons of three resolutions; and the areas of some cells on the reference's
sphere. Run where `h3` is installed (it is not a dependency of the infrastructure):

    python foundations/physis/computational_geometry/tests/make_h3_fixtures.py
"""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

import h3

OUT = Path(__file__).resolve().parent / "fixtures" / "h3.txt"


def entry(lat: float, lon: float, res: int) -> dict:
    c = h3.latlng_to_cell(lat, lon, res)
    return {
        "lat": lat, "lon": lon, "res": res, "cell": c,
        "centre": list(h3.cell_to_latlng(c)),
        "boundary": [list(p) for p in h3.cell_to_boundary(c)],
        "parent": h3.cell_to_parent(c, max(0, res - 1)) if res else c,
        "pentagon": bool(h3.is_pentagon(c)), "base": h3.get_base_cell_number(c),
    }


def main() -> None:
    rng = random.Random(20261107)
    points = []
    for _ in range(800):
        # Uniform on the sphere.
        points.append(entry(math.degrees(math.asin(rng.uniform(-1, 1))), rng.uniform(-180, 180), rng.randint(0, 15)))
    for lat, lon in [(89.9999, 10.0), (-89.9999, -120.0), (0.0, 179.9999), (0.0, -179.9999), (-37.8136, 144.9631), (40.7128, -74.006), (8.4657, -13.2317)]:
        for res in range(16):
            points.append(entry(lat, lon, res))
    lines = [f"# H3 reference values: h3 {h3.__version__} (Python, on the reference C library), made by make_h3_fixtures.py",
             "# P lat lon res cell parent pentagon base centre_lat centre_lon n (lat lon)*n ; A cell area_m2 ; all in degrees"]
    number = lambda v: f"{v:.15g}"
    for e in points:
        b = e["boundary"]
        lines.append(" ".join(["P", number(e["lat"]), number(e["lon"]), str(e["res"]), e["cell"], e["parent"], str(int(e["pentagon"])), str(e["base"]), number(e["centre"][0]), number(e["centre"][1]), str(len(b)), *[number(x) for p in b for x in p]]))
    for res in (0, 3, 6, 9, 10):
        for c in h3.get_pentagons(res):
            la, lo = h3.cell_to_latlng(c)
            b = h3.cell_to_boundary(c)
            lines.append(" ".join(["P", number(la), number(lo), str(res), c, h3.cell_to_parent(c, max(0, res - 1)) if res else c, "1", str(h3.get_base_cell_number(c)), number(la), number(lo), str(len(b)), *[number(x) for p in b for x in p]]))
    for _ in range(40):
        c = h3.latlng_to_cell(rng.uniform(-80, 80), rng.uniform(-180, 180), rng.randint(0, 9))
        lines.append(f"A {c} {number(h3.cell_area(c, unit='m^2'))}")
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(chr(10).join(lines) + chr(10), encoding="utf-8", newline=chr(10))
    print(f"{OUT}: {OUT.stat().st_size} bytes, {len(lines)} lines", file=sys.stderr)


if __name__ == "__main__":
    main()
