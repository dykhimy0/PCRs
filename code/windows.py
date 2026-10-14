"""
windows.py -- single source of truth for the analysis windows.

Every window is one contiguous 3072 x 3072 block of the ESA WorldCover v200
(2021) 10 m map, read from the public COG on S3.  A window is named by the
tile it sits in and the (row, col) offset of its north-west corner.

The three ORIGINAL windows all lie in the temperate-forest belt of central
Europe and in the same product, so they vary the site but not the biome.
The three ADDED windows give two further biomes and, for the Mediterranean
biome, two sites, so that a biome effect can be told apart from a site effect
within one biome.

Row index increases southwards, column index increases eastwards; the tile is
3 x 3 degrees = 36000 x 36000 cells at 0.0000833 deg.
"""

from __future__ import annotations

# tile,       row,   col,   label,                    biome,          group
WINDOWS = [
    ("N48E006", 26000, 26000, "N48E006 (26k,26k)",     "temperate", "original"),
    ("N48E006", 28000, 20000, "N48E006 (28k,20k)",     "temperate", "original"),
    ("N48E009", 26000, 26000, "N48E009 (26k,26k)",     "temperate", "original"),
    ("N39W009", 18000, 18000, "N39W009 (18k,18k)",     "mediterranean", "added"),
    ("N36W006", 15000, 15000, "N36W006 (15k,15k)",     "mediterranean", "added"),
    ("N66E024", 26000, 20000, "N66E024 (26k,20k)",     "boreal", "added"),
]

ORIGINAL_TILES = [w[0] for w in WINDOWS if w[5] == "original"]

# place on Earth, only for describing the sites in prose
PLACE = {
    "N48E006 (26k,26k)": "Black Forest, southwest Germany",
    "N48E006 (28k,20k)": "Kraichgau, southwest Germany",
    "N48E009 (26k,26k)": "Swabian Jura, southern Germany",
    "N39W009 (18k,18k)": "Alentejo montado, Portugal",
    "N36W006 (15k,15k)": "Andalusian dehesa, Spain",
    "N66E024 (26k,20k)": "Finnish Lapland",
}

# latitude of the window centre, needed because the map is in EPSG:4326 and a
# cell is therefore not square on the ground: at 48 N a degree of longitude is
# 0.67 of a degree of latitude, at 66 N it is 0.41.  A "30 m" block built from
# 3 x 3 cells is 27.8 m (N-S) x 18.6 m (E-W) in the temperate windows and
# 27.8 m x 11.3 m in the boreal one.  This is reported, not corrected, so that
# the reader can see it.
LAT_CENTRE = {
    "N48E006 (26k,26k)": 48.9, "N48E006 (28k,20k)": 48.9,
    "N48E009 (26k,26k)": 48.9, "N39W009 (18k,18k)": 39.9,
    "N36W006 (15k,15k)": 36.9, "N66E024 (26k,20k)": 66.9,
}


def tile_of(label: str) -> str:
    return label.split(" ")[0]


def biome_of(label: str) -> str:
    for t, _r, _c, lb, biome, _g in WINDOWS:
        if lb == label:
            return biome
    return "unknown"
