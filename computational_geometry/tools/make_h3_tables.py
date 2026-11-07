"""The tables of the H3 grid, copied from the reference implementation at a pinned tag into h3_tables.hpp.

H3's grid is defined by tables: where the icosahedron's faces stand, how they adjoin, which base cell lies at
each position of a face and how it is turned. They are data to be copied, not results to be derived, and this script
copies them, so that a reader can see where each number came from. The reference is Apache-2.0; its notice is kept in
the header written.

    python foundations/physis/computational_geometry/tools/make_h3_tables.py
"""

from __future__ import annotations

import re
import urllib.request
from pathlib import Path

TAG = "v4.3.0"
BASE = f"https://raw.githubusercontent.com/uber/h3/{TAG}/src/h3lib/lib/"
OUT = Path(__file__).resolve().parents[1] / "include" / "physis" / "computational_geometry" / "h3_tables.hpp"
AGENT = {"User-Agent": "science-infrastructures/1.0"}


def source(name: str) -> str:
    with urllib.request.urlopen(urllib.request.Request(BASE + name, headers=AGENT), timeout=60) as r:
        text = r.read().decode("utf-8")
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)


def initialiser(text: str, name: str) -> str:
    """The braces after `name[...] = `, whole, without comments."""
    m = re.search(re.escape(name) + r"(\s*\[[^\]]*\])+\s*=\s*\{", text)
    if not m:
        raise SystemExit(f"{name} not found in the reference")
    start = text.index("{", m.start())
    depth, i = 0, start
    while True:
        c = text[i]
        depth += c == "{"
        depth -= c == "}"
        i += 1
        if depth == 0:
            break
    body = re.sub(r"//[^\n]*", "", text[start:i])
    return re.sub(r"\s+", " ", body).strip()


def main() -> None:
    face, base = source("faceijk.c"), source("baseCells.c")
    tables = [
        ("face_center_geo", "constexpr GeoRad", "[20]", initialiser(face, "faceCenterGeo")),
        ("face_center_point", "constexpr double", "[20][3]", initialiser(face, "faceCenterPoint")),
        ("face_axes_az", "constexpr double", "[20][3]", initialiser(face, "faceAxesAzRadsCII")),
        ("face_neighbors", "constexpr FaceOrient", "[20][4]", initialiser(face, "faceNeighbors")),
        ("adjacent_face_dir", "constexpr int", "[20][20]", initialiser(face, "adjacentFaceDir")),
        ("max_dim", "constexpr int", "[17]", initialiser(face, "maxDimByCIIres")),
        ("unit_scale", "constexpr int", "[17]", initialiser(face, "unitScaleByCIIres")),
        ("face_ijk_base_cells", "constexpr BaseCellRotation", "[20][3][3][3]", initialiser(base, "faceIjkBaseCells")),
        ("base_cell_data", "constexpr BaseCell", "[122]", initialiser(base, "baseCellData")),
    ]
    out = [
        f"// The tables of the H3 grid, copied by tools/make_h3_tables.py from the reference implementation (uber/h3, tag {TAG}). Do not edit.",
        "//",
        "// Copyright 2016-2021 Uber Technologies, Inc. Licensed under the Apache License, Version 2.0",
        "// (http://www.apache.org/licenses/LICENSE-2.0); distributed on an \"AS IS\" BASIS, without warranties or conditions of any kind.",
        "",
        "#pragma once",
        "",
        "namespace physis::computational_geometry::h3::tables {",
        "",
        "struct GeoRad { double lat, lng; };",
        "struct Ijk { int i, j, k; };",
        "struct FaceOrient { int face; Ijk translate; int ccw_rot60; };",
        "struct BaseCellRotation { int base_cell, ccw_rot60; };",
        "struct FaceIjkBase { int face; Ijk coord; };",
        "struct BaseCell { FaceIjkBase home; int is_pentagon; int cw_offset_pent[2]; };",
        "",
        "// Directions from a face to its neighbours.",
        "inline constexpr int IJ = 1, KI = 2, JK = 3, INVALID_FACE = -1;",
        "",
    ]
    for name, kind, dims, body in tables:
        # The reference's names for these constants, as written in its tables.
        out.append(f"inline {kind} {name}{dims} = {body};")
        out.append("")
    out.append("}  // namespace physis::computational_geometry::h3::tables")
    OUT.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    print(f"{OUT}: {OUT.stat().st_size} bytes")


if __name__ == "__main__":
    main()
