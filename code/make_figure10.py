"""
make_figure10.py -- the mechanism behind the biome grading of Section 6.13.

Section 6.13 reports that the truncated-window bias is largest in the
Mediterranean windows and smallest in the boreal one, and attributes that order
to how much of the sill a 60-pixel window resolves.  An order is a weak claim;
this figure turns it into a relation that can be extrapolated.

Panel (a) the mean absolute error of the 60-pixel fit within each window against
          the fraction of the sill the empirical variogram has reached at 60
          pixels.  If the attribution is right the six points fall on a monotone
          line; the line does not use the biome label, so the biome ordering is
          a consequence and not an input.
Panel (b) the same resolved fraction against the fitted range.  A window with a
          short range has already reached its sill at 60 pixels, so a truncated
          fit has little left to miss.  The guide line is an empirical fit in
          1/a, not the exponential model evaluated at 60, because the empirical
          variogram rises faster near the origin than a single exponential does
          -- which is also why the wide exponential fit remains biased low.

Reads results/baseline_variogram.csv and results/baseline_fits.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 120,
    "savefig.dpi": 320,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "legend.frameon": True,
    "legend.fontsize": 7.2,
})

BIOME_COLOR = {"temperate": "#2F5D8A", "mediterranean": "#C4762A",
               "boreal": "#5B4B8A"}
MARKER = {"temperate": "o", "mediterranean": "s", "boreal": "^"}


def r2_of(y: np.ndarray, pred: np.ndarray) -> float:
    return float(1.0 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum())


def fig10() -> None:
    bv = pd.read_csv(RES / "baseline_variogram.csv")
    ft = pd.read_csv(RES / "baseline_fits.csv")

    per = (bv.assign(a=bv["err_truncated_pct"].abs())
             .groupby(["window", "biome"], as_index=False)
             .agg(err=("a", "mean")))
    tr = ft[ft["model"] == "truncated"][["window", "gamma_at_trunc_over_sill"]]
    wi = (ft[ft["model"] == "wide"][["window", "a1"]]
          .rename(columns={"a1": "range_wide"}))
    d = per.merge(tr, on="window").merge(wi, on="window")

    fig, axs = plt.subplots(1, 2, figsize=(9.0, 3.7))

    # ---------------- (a) error against the resolved fraction ------------
    ax = axs[0]
    x = d["gamma_at_trunc_over_sill"].to_numpy()
    y = d["err"].to_numpy()
    slope, icpt = np.polyfit(x, y, 1)
    r2a = r2_of(y, slope * x + icpt)
    xs = np.linspace(x.min() - 0.008, x.max() + 0.008, 50)
    ax.plot(xs, slope * xs + icpt, "-", color="#8A8A83", lw=1.5, zorder=0,
            label=f"OLS on $x$, $R^2={r2a:.2f}$")
    for b, sub in d.groupby("biome"):
        ax.plot(sub["gamma_at_trunc_over_sill"], sub["err"], MARKER[b], ms=8.5,
                color=BIOME_COLOR[b], mec="white", mew=0.9, ls="none", label=b)
    for r in d.itertuples():
        ax.annotate(r.window.split(" ")[0], (r.gamma_at_trunc_over_sill, r.err),
                    textcoords="offset points", xytext=(7, -3), fontsize=6.2,
                    color="#555555")
    ax.set_xlabel(r"$\gamma(60\,\mathrm{px})\;/\;\mathrm{sill}$  (empirical)")
    ax.set_ylabel("mean |bias| of the 60-px fit (%)")
    ax.set_xlim(0.745, 0.975)
    ax.set_ylim(8, 50)
    ax.legend(loc="upper right", framealpha=0.95, edgecolor="none")
    ax.set_title("(a) The truncated-window bias is set by\n"
                 "how much of the sill 60 px resolves")

    # ---------------- (b) the cause: range versus resolved fraction ------
    ax = axs[1]
    a = d["range_wide"].to_numpy()
    inv = 1.0 / a
    k, c0 = np.polyfit(inv, x, 1)
    r2b = r2_of(x, k * inv + c0)
    aa = np.linspace(a.min() - 1.5, a.max() + 3.0, 200)
    ax.plot(aa, c0 + k / aa, "-", color="#8A8A83", lw=1.5, zorder=0,
            label=f"empirical fit in $1/a$, $R^2={r2b:.2f}$")
    for b, sub in d.groupby("biome"):
        ax.plot(sub["range_wide"], sub["gamma_at_trunc_over_sill"], MARKER[b],
                ms=8.5, color=BIOME_COLOR[b], mec="white", mew=0.9, ls="none",
                label=b)
    for r in d.itertuples():
        ax.annotate(r.window.split(" ")[0],
                    (r.range_wide, r.gamma_at_trunc_over_sill),
                    textcoords="offset points", xytext=(7, -3), fontsize=6.2,
                    color="#555555")
    ax.set_xlabel("fitted range over the wide window, $a$ (pixels)")
    ax.set_ylabel(r"$\gamma(60\,\mathrm{px})\;/\;\mathrm{sill}$")
    ax.set_xlim(10, 76)
    ax.set_ylim(0.74, 0.99)
    ax.legend(loc="upper right", framealpha=0.95, edgecolor="none")
    ax.set_title("(b) A short range saturates inside 60 px,\n"
                 "a long one does not")

    fig.tight_layout()
    fig.savefig(FIG / "fig10_truncation_mechanism.png")
    plt.close(fig)
    print(f"wrote fig10_truncation_mechanism.png  (n={len(d)})")
    print(f"  (a) |bias| ~ gamma60/sill : slope {slope:.1f}, R2 {r2a:.4f}")
    print(f"  (b) gamma60/sill ~ 1/a    : R2 {r2b:.4f}")
    print(d.sort_values("gamma_at_trunc_over_sill")
          .to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    fig10()
