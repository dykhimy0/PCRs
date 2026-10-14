# -*- coding: utf-8 -*-
"""
compare_metric_geographic.py -- what the move to metric support actually changes.

The reprojection was not performed to make the numbers prettier; it was
performed because a fitted range in "pixels" is not a distance when the pixel
is a latitude-dependent rectangle.  This script lines the two sets of numbers
up and answers the only question that matters for the revision: which of the
manuscript's claims survive the rebuild, which have to be restated, and are
there any that the rebuild overturns?

The central column is the RANGE AMBIGUITY.  A geographic range of `a` pixels
corresponds to a ground distance that depends on which way you read it:

    read north-south   a * 27.643 m          (the same at every latitude)
    read east-west     a * 27.643 * cos(lat) m

because a 3 x 3 block of 1/12000-degree cells is 27.64 m tall at every
latitude but only 27.64 * cos(lat) m wide.  The ratio between the two readings
is 1/cos(lat) and reaches 2.51 at 66.7 N.  The metric range is a single number
and must lie in that bracket if the two analyses are describing the same
landscape -- so the bracket is also a consistency check, not just a criticism.

Outputs
    _metric_vs_geographic_impact.csv   one row per window, every comparable
    a printed report grouped by what the change means for the text
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "code"))
sys.stdout.reconfigure(encoding="utf-8")

RES_G = ROOT / "results_v2"
RES_M = ROOT / "results_metric"
MET = ROOT / "data" / "metric" / "_metric_fields.csv"
OUT = ROOT / "_metric_vs_geographic_impact.csv"

PIXEL_M = 30.0


def _pct(new, old):
    return 100.0 * (new - old) / old if old else np.nan


def load() -> pd.DataFrame:
    meta = pd.read_csv(MET)
    gb = pd.read_csv(RES_G / "real_blocks.csv")
    mb = pd.read_csv(RES_M / "metric_blocks.csv")
    gf = pd.read_csv(RES_G / "baseline_fits.csv")
    mf = pd.read_csv(RES_M / "metric_fits.csv")
    gv = pd.read_csv(RES_G / "baseline_variogram.csv")
    mv = pd.read_csv(RES_M / "metric_variogram.csv")
    gr = pd.read_csv(RES_G / "rho_scalar.csv")
    mr = pd.read_csv(RES_M / "metric_rho_scalar.csv")
    ma = pd.read_csv(RES_M / "metric_anisotropy.csv")

    def pick(df, model, col):
        s = df[df["model"] == model].set_index("window")[col]
        return s

    rows = []
    for _, r in meta.iterrows():
        lab = r["window"]
        g = gb[gb["window"] == lab].set_index("m")
        m = mb[mb["window"] == lab].set_index("m")
        gvw = gv[gv["window"] == lab].set_index("m")
        mvw = mv[mv["window"] == lab].set_index("m")
        grw = gr[gr["window"] == lab]
        mrw = mr[mr["window"] == lab]
        arow = ma[ma["window"] == lab].iloc[0]

        # ---- the range ambiguity -------------------------------------
        a_geo_px = float(pick(gf, "wide", "a1")[lab])
        a_met_m = float(pick(mf, "wide", "a1_m")[lab])
        ns_m = float(r["geo_pixel_ns_m"])
        ew_m = float(r["geo_pixel_ew_m"])
        rec = dict(
            window=lab, biome=r["biome"], place=r["place"],
            lat=r["lat"],
            geo_pixel_ns_m=ns_m, geo_pixel_ew_m=ew_m,
            geo_pixel_aspect=r["geo_pixel_aspect"],
            met_pixel_m=PIXEL_M,

            sigma2_geo=float(g["pixel_marginal_var"].iloc[0]),
            sigma2_met=float(m["sigma2"].iloc[0]),

            tree_geo=r["tree_geo"], tree_met=r["tree_met"],
            tree_shift_pp=r["tree_shift_pp"],
            met_pixels=r["met_pixels"],

            range_geo_px=a_geo_px,
            range_geo_if_ns_m=a_geo_px * ns_m,
            range_geo_if_ew_m=a_geo_px * ew_m,
            range_geo_ambiguity=ns_m / ew_m,
            range_met_m=a_met_m,
            range_met_px=a_met_m / PIXEL_M,
            range_inside_bracket=bool(
                min(a_geo_px * ew_m, a_geo_px * ns_m) <= a_met_m
                <= max(a_geo_px * ew_m, a_geo_px * ns_m)),

            range_trunc_geo_px=float(pick(gf, "truncated", "a1")[lab]),
            range_trunc_met_m=float(pick(mf, "truncated", "a1_m")[lab]),

            aniso_ratio_ew_over_ns=float(arow["range_ratio_ew_over_ns"]),
            range_ew_m=float(arow["range_ew_m"]),
            range_ns_m=float(arow["range_ns_m"]),
        )
        rec["range_change_pct"] = _pct(a_met_m, a_geo_px * ns_m)
        rec["range_trunc_change_pct"] = _pct(
            rec["range_trunc_met_m"], rec["range_trunc_geo_px"] * ns_m)

        # ---- change of support ---------------------------------------
        for mm in (4, 9, 16, 25):
            if mm not in g.index or mm not in m.index:
                continue
            rec[f"m{mm}_vtrue_geo"] = float(g.loc[mm, "truth_block_var"])
            rec[f"m{mm}_vtrue_met"] = float(m.loc[mm, "truth_block_var"])
            rec[f"m{mm}_vtrue_pct"] = _pct(float(m.loc[mm, "truth_block_var"]),
                                           float(g.loc[mm, "truth_block_var"]))
            rec[f"m{mm}_under_geo"] = 1.0 / float(
                g.loc[mm, "independent_understates_by"])
            rec[f"m{mm}_under_met"] = 1.0 / float(
                m.loc[mm, "independent_understates_by"])
            rec[f"m{mm}_aci_geo"] = float(g.loc[mm, "aci"])
            rec[f"m{mm}_aci_met"] = float(m.loc[mm, "aci"])
            rec[f"m{mm}_inside_geo"] = bool(g.loc[mm, "truth_inside"])
            rec[f"m{mm}_inside_met"] = bool(m.loc[mm, "truth_inside"])

        # ---- the geostatistical route's error against truth ----------
        for mm in (4, 25):
            if mm in gvw.index and mm in mvw.index:
                rec[f"m{mm}_errwide_geo"] = float(gvw.loc[mm, "err_wide_pct"])
                rec[f"m{mm}_errwide_met"] = float(mvw.loc[mm, "err_wide_pct"])

        # ---- one stored scalar ---------------------------------------
        rec["rho_spread_geo"] = float(grw["rho_exact_spread"].iloc[0])
        rec["rho_spread_met"] = float(mrw["rho_exact_spread"].iloc[0])
        rec["err_lag1_geo_p25"] = float(
            grw[grw["m"] == 25]["err_lag1_pct"].iloc[0]) if (grw["m"] == 25).any() else np.nan
        rec["err_lag1_met_p25"] = float(
            mrw[mrw["m"] == 25]["err_lag1_pct"].iloc[0]) if (mrw["m"] == 25).any() else np.nan
        rows.append(rec)

    return pd.DataFrame(rows)


def report(df: pd.DataFrame) -> None:
    pd.set_option("display.width", 200)
    print("=" * 100)
    print("1.  THE RANGE IS NOT A DISTANCE ON THE GEOGRAPHIC GRID")
    print("=" * 100)
    print(f"{'window':18s}{'lat':>7}{'geo px':>8}{'if N-S':>9}{'if E-W':>9}"
          f"{'metric':>9}{'in bracket':>12}{'ambiguity':>11}")
    print("-" * 100)
    for _, r in df.sort_values("lat", ascending=False).iterrows():
        print(f"{r['window']:18s}{r['lat']:>7.1f}{r['range_geo_px']:>8.1f}"
              f"{r['range_geo_if_ns_m']:>9.0f}{r['range_geo_if_ew_m']:>9.0f}"
              f"{r['range_met_m']:>9.0f}"
              f"{('yes' if r['range_inside_bracket'] else 'NO'):>12}"
              f"{r['range_geo_ambiguity']:>11.2f}")
    print("-" * 100)
    print(f"the same fitted range can mean anything from "
          f"{df['range_geo_if_ew_m'].min():.0f} m to "
          f"{df['range_geo_if_ns_m'].max():.0f} m depending on the direction it "
          f"is read in;")
    print(f"the north-south / east-west ambiguity of one pixel runs from "
          f"{df['range_geo_ambiguity'].min():.2f} to "
          f"{df['range_geo_ambiguity'].max():.2f}")
    print(f"the metric range falls inside the geographic bracket in "
          f"{int(df['range_inside_bracket'].sum())} of {len(df)} windows")

    print()
    print("=" * 100)
    print("2.  WHAT THE REBUILD DOES TO THE CHANGE-OF-SUPPORT NUMBERS")
    print("=" * 100)
    print(f"{'window':18s}{'sigma2 %':>10}{'m4 V %':>9}{'m25 V %':>9}"
          f"{'m25 under g':>12}{'m25 under m':>13}"
          f"{'aci25 g':>9}{'aci25 m':>9}")
    print("-" * 100)
    for _, r in df.iterrows():
        print(f"{r['window']:18s}{_pct(r['sigma2_met'], r['sigma2_geo']):>+10.2f}"
              f"{r['m4_vtrue_pct']:>+9.2f}{r['m25_vtrue_pct']:>+9.2f}"
              f"{r['m25_under_geo']:>12.2f}{r['m25_under_met']:>13.2f}"
              f"{r['m25_aci_geo']:>9.3f}{r['m25_aci_met']:>9.3f}")
    print("-" * 100)
    for m in (4, 25):
        c = df[f"m{m}_vtrue_pct"]
        print(f"m={m:2d}  V_true change {c.min():+.1f}% to {c.max():+.1f}%  "
              f"mean {c.mean():+.1f}%")
    print(f"truth inside the identified band : geographic "
          f"{bool(df[[c for c in df.columns if c.endswith('_inside_geo')]].all().all())}"
          f" | metric "
          f"{bool(df[[c for c in df.columns if c.endswith('_inside_met')]].all().all())}")
    ug = df["m25_under_geo"]
    um = df["m25_under_met"]
    print(f"m=25 independence understates the truth by "
          f"geographic {ug.min():.2f}-{ug.max():.2f}x | "
          f"metric {um.min():.2f}-{um.max():.2f}x")
    print(f"m=25 ACI : geographic {df['m25_aci_geo'].min():.3f}-"
          f"{df['m25_aci_geo'].max():.3f} | metric {df['m25_aci_met'].min():.3f}-"
          f"{df['m25_aci_met'].max():.3f}")

    print()
    print("=" * 100)
    print("3.  ANISOTROPY -- only measurable once the grid is metric")
    print("=" * 100)
    print(f"{'window':18s}{'biome':15s}{'range E-W':>11}{'range N-S':>11}"
          f"{'E-W / N-S':>11}")
    print("-" * 100)
    for _, r in df.sort_values("aniso_ratio_ew_over_ns").iterrows():
        print(f"{r['window']:18s}{r['biome']:15s}{r['range_ew_m']:>11.0f}"
              f"{r['range_ns_m']:>11.0f}{r['aniso_ratio_ew_over_ns']:>11.2f}")
    print("-" * 100)
    rt = df["aniso_ratio_ew_over_ns"]
    print(f"ratio {rt.min():.2f} to {rt.max():.2f}; "
          f"{int(((rt < 0.8) | (rt > 1.25)).sum())} of {len(df)} windows are "
          f"outside +/-25 % of isotropic")
    print("this is a property of the landscape, not of the grid: on the metric")
    print("grid one cell is 10 m in both directions, so a difference between the")
    print("two directions cannot be a grid artefact.")

    print()
    print("=" * 100)
    print("4.  ONE STORED SCALAR -- does the support dependence survive?")
    print("=" * 100)
    print(f"{'window':18s}{'rho spread g':>14}{'rho spread m':>14}"
          f"{'lag1 err g %':>14}{'lag1 err m %':>14}")
    print("-" * 100)
    for _, r in df.iterrows():
        print(f"{r['window']:18s}{r['rho_spread_geo']:>14.4f}"
              f"{r['rho_spread_met']:>14.4f}"
              f"{r['err_lag1_geo_p25']:>14.2f}{r['err_lag1_met_p25']:>14.2f}")
    print("-" * 100)
    print(f"rho_exact spread across m : geographic "
          f"{df['rho_spread_geo'].min():.4f}-{df['rho_spread_geo'].max():.4f} | "
          f"metric {df['rho_spread_met'].min():.4f}-{df['rho_spread_met'].max():.4f}")

    print()
    print("=" * 100)
    print("5.  WHAT THIS MEANS FOR THE MANUSCRIPT")
    print("=" * 100)
    claims = [
        ("The truth lies inside the identified band in every case",
         bool(df[[c for c in df.columns if c.endswith('_inside_met')]].all().all()),
         "unchanged -- the central claim is support-independent"),
        ("Independence understates the block variance by 3-4x",
         bool(um.min() >= 3.0),
         f"metric says {um.min():.1f}-{um.max():.1f}x at m=25; restate, do not "
         f"re-derive"),
        ("ACI stays close to its worst case 1 - 1/m",
         bool(df["m25_aci_met"].min() > 0.9),
         f"metric {df['m25_aci_met'].min():.3f}-{df['m25_aci_met'].max():.3f}"),
        ("The fitted range is 10-88 pixels",
         False,
         "RETRACT as a distance: it is "
         f"{df['range_geo_if_ew_m'].min():.0f}-"
         f"{df['range_geo_if_ns_m'].max():.0f} m depending on direction; "
         f"metric {df['range_met_m'].min():.0f}-{df['range_met_m'].max():.0f} m"),
        ("A truncated lag window biases the range low",
         bool((df["range_trunc_met_m"] < df["range_met_m"]).all()),
         "survives on the metric grid, and is now a statement about metres"),
        ("The fields are treated as isotropic",
         not bool(((rt < 0.8) | (rt > 1.25)).any()),
         f"{int(((rt < 0.8) | (rt > 1.25)).sum())} windows are outside "
         f"+/-25 % -- this is now a measured limitation, not an assumption"),
    ]
    for text, ok, note in claims:
        print(f"  [{'holds  ' if ok else 'CHANGED'}]  {text}")
        print(f"            {note}")


def main() -> None:
    df = load()
    df.to_csv(OUT, index=False)
    report(df)
    print(f"\nCSV: {OUT}")


if __name__ == "__main__":
    main()
