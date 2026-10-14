"""
windows_v2.py -- single source of truth for the analysis windows (revision 2).

Revision 2 widens the biome coverage from three biomes / six windows to seven
biomes / eleven windows.  Five windows are new:

  N63E027 (14k, 6k)  boreal          Finnish Kainuu           (2nd boreal site)
  N21E099 (24k,12k)  tropical        Xishuangbanna, Yunnan    (new biome)
  N24E117 (14k, 6k)  subtropical     Wuyi Shan, Fujian        (new biome)
  N27E093 (10k,16k)  alpine          Bomi, Tibet              (new biome)
  N45E126 ( 2k,25k)  cool-temperate  Xiaoxing'anling, Heilongjiang (new biome)

The point of the additions is that four of them lie in mainland China, which
supplies the tropical, subtropical and alpine belts that the original design
lacked outright plus a cool-temperate forest belt, and that the second boreal
window turns the boreal biome from a single site into a biome with an
estimable within-biome spread.

Latitude is computed from the tile transform rather than assumed: an ESA
WorldCover tile is named after its south-west corner and runs 3 degrees north
and east of that, with row 0 at the NORTH edge, so a window centred at
(row + 1536, col + 1536) sits at

    lat = tile_N + 3 - (row + 1536) / 12000
    lon = tile_E     + (col + 1536) / 12000

(12000 cells per degree, 36000 cells per tile).  Revision 1 of this file used
the north EDGE of the window and therefore overstated the centre latitude by
about 0.2 degrees; the ground dimensions of a 3 x 3 block change by < 1 %, but
the stated latitudes did not match the product grid, so they are corrected here.
"""

from __future__ import annotations

CELLS_PER_DEG = 12000.0
TILE_SPAN = 36000
WINDOW = 3072

# tile,        row,   col,   biome,            group,      place
WINDOWS = [
    ("N48E006", 26000, 26000, "temperate",       "original", "Black Forest, southwest Germany"),
    ("N48E006", 28000, 20000, "temperate",       "original", "Kraichgau, southwest Germany"),
    ("N48E009", 26000, 26000, "temperate",       "original", "Swabian Jura, southern Germany"),
    ("N39W009", 18000, 18000, "mediterranean",   "original", "Alentejo montado, Portugal"),
    ("N36W006", 15000, 15000, "mediterranean",   "original", "Andalusian dehesa, Spain"),
    ("N66E024", 26000, 20000, "boreal",          "original", "Finnish Lapland"),
    ("N63E027", 14000,  6000, "boreal",          "added",    "Kainuu, Finland"),
    ("N21E099", 24000, 12000, "tropical",        "added",    "Xishuangbanna, Yunnan, China"),
    ("N24E117", 14000,  6000, "subtropical",     "added",    "Wuyi Shan, Fujian, China"),
    ("N27E093", 10000, 16000, "alpine",          "added",    "Bomi, Tibet, China"),
    ("N45E126",  2000, 25000, "cool-temperate",  "added",    "Xiaoxing'anling, Heilongjiang, China"),
]


def label_of(tile: str, row: int, col: int) -> str:
    """Stable, human-readable window name; also the join key in the CSVs."""
    return f"{tile} ({row // 1000}k,{col // 1000}k)"


def _tile_origin(tile: str) -> tuple[float, float]:
    """South-west corner of a WorldCover tile, e.g. N48E006 -> (48.0, 6.0)."""
    ns, lat = tile[0], float(tile[1:3])
    ew, lon = tile[3], float(tile[4:])
    return (lat if ns == "N" else -lat), (lon if ew == "E" else -lon)


def centre(tile: str, row: int, col: int) -> tuple[float, float]:
    """(lat, lon) of the window centre, in degrees."""
    lat0, lon0 = _tile_origin(tile)
    off = (row + WINDOW / 2) / CELLS_PER_DEG
    offx = (col + WINDOW / 2) / CELLS_PER_DEG
    return lat0 + 3.0 - off, lon0 + offx


# label -> (lat, lon)
CENTRE = {label_of(t, r, c): centre(t, r, c) for t, r, c, *_ in WINDOWS}
LABEL = {label_of(t, r, c): f"{t} ({r // 1000}k,{c // 1000}k)" for t, r, c, *_ in WINDOWS}
BIOME = {label_of(t, r, c): b for t, r, c, b, *_ in WINDOWS}
PLACE = {label_of(t, r, c): p for t, r, c, _b, _g, p in WINDOWS}
LABELS = [label_of(t, r, c) for t, r, c, *_ in WINDOWS]

# biome display order and the accent colour used by the figures
BIOME_ORDER = ["temperate", "cool-temperate", "mediterranean", "boreal",
               "alpine", "subtropical", "tropical"]
BIOME_SHORT = {"temperate": "temp", "cool-temperate": "cool-temp",
               "mediterranean": "medit", "boreal": "boreal",
               "alpine": "alpine", "subtropical": "subtrop",
               "tropical": "trop"}


def biome_of(label: str) -> str:
    return BIOME.get(label, "unknown")


def lat_of(label: str) -> float:
    return CENTRE[label][0]


def lon_of(label: str) -> float:
    return CENTRE[label][1]


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print(f"{'label':18s} {'biome':15s} {'lat':>7s} {'lon':>8s}  place")
    for t, r, c, b, g, p in WINDOWS:
        lb = label_of(t, r, c)
        la, lo = centre(t, r, c)
        print(f"{lb:18s} {b:15s} {la:7.3f} {lo:8.3f}  {p}")
    print(f"\n{len(WINDOWS)} windows, {len(set(w[3] for w in WINDOWS))} biomes")
    from collections import Counter
    print(dict(Counter(w[3] for w in WINDOWS)))
