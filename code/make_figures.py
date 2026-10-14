"""
make_figures.py -- publication figures for the scale-consistency study.

All figures are generated from the CSV outputs of exp_theory.py,
exp_spatial.py and exp_real.py; nothing is hand-entered.
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
DATA = ROOT / "data"
FIG.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 120,
    "savefig.dpi": 320,
    "savefig.bbox": "tight",
    "legend.frameon": False,
    "legend.fontsize": 7.6,
})

C_TRUTH = "#1A1A1A"
C_IND = "#D1823C"
C_BAND = "#9DBBD3"
C_UPPER = "#2F5D8A"
C_LOWER = "#7FA8C9"
C_COMO = "#3E7C59"
C_CTR = "#A8474B"
C_LEVB = "#5B4B8A"
C_ACC = "#B03A5B"
L_TICKS = [0, 1, 2, 4, 8, 16, 32]


def _save(fig, name):
    fig.savefig(FIG / name)
    plt.close(fig)
    print("wrote", name)


def _lx(ax):
    ax.set_xlim(-1.2, 34)
    ax.set_xticks(L_TICKS)
    ax.set_xticklabels([str(t) for t in L_TICKS])
    ax.set_xlabel("correlation length (fine cells)")


# --------------------------------------------------------------------------
# Figure 1 -- scale amplification
# --------------------------------------------------------------------------
def fig1():
    df = pd.read_csv(RES / "scale_amplification.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0))

    ax = axs[0]
    for p, col, lab in ((0.5, C_UPPER, "$p=0.5$"), (0.3, "#7E9DBA", "$p=0.3$")):
        d = df[df["p"] == p].sort_values("m")
        ax.fill_between(d["m"], d["v_lower_exact"], d["v_upper_exact"],
                        color=col, alpha=0.22 if p == 0.5 else 0.13)
        ax.plot(d["m"], d["v_upper_exact"], "-", color=col, lw=1.5)
        ax.plot(d["m"], d["v_lower_exact"], "-", color=col, lw=1.2, ls="-.")
    d = df[df["p"] == 0.5].sort_values("m")
    ax.plot(d["m"], d["v_independent"], "--o", color=C_IND, ms=3.4, lw=1.4,
            label="independence coupling")
    ax.fill_between([], [], [], color=C_UPPER, alpha=0.22,
                    label="identified set of Var($\\nu_C$)")
    ax.annotate("", xy=(16, 0.246), xytext=(16, 0.0156),
                arrowprops=dict(arrowstyle="<->", color="#555555", lw=0.9))
    ax.text(17.0, 0.10, "factor $m$", fontsize=7.6, color="#555555", rotation=90,
            va="center")
    ax.text(2.1, 0.235, "$V^+=V_1$", fontsize=7.6, color=C_UPPER)
    ax.text(2.05, 0.135, "$p=0.3$", fontsize=7.2, color="#5C7A96")
    ax.set_xscale("log", base=2)
    ax.set_xticks([2, 4, 8, 16]); ax.set_xticklabels(["2", "4", "8", "16"])
    ax.set_ylim(-0.012, 0.275)
    ax.set_xlabel("number of pixels aggregated, $m$")
    ax.set_ylabel("Var($\\nu_C$)")
    ax.set_title("(a) Aggregation widens the identified set")
    ax.legend(loc="upper right")

    ax = axs[1]
    d = df[df["p"] == 0.5].sort_values("m")
    ax.plot(d["m"], d["aci_exact"], "o-", color=C_ACC, lw=1.7, ms=4.2,
            label="$\\mathrm{ACI}(m)$, exact LP, $p=0.5$")
    ms = np.array(sorted(d["m"]))
    ax.plot(ms, 1 - 1 / ms, "--", color="#7A7A7A", lw=1.3,
            label="worst-case bound $1-1/m$")
    d3 = df[df["p"] == 0.3].sort_values("m")
    ax.plot(d3["m"], d3["aci_exact"], "s:", color="#B9846A", lw=1.3, ms=3.2,
            label="$\\mathrm{ACI}(m)$, $p=0.3$")
    ax.set_xscale("log", base=2)
    ax.set_xticks([2, 4, 8, 16]); ax.set_xticklabels(["2", "4", "8", "16"])
    ax.set_ylim(0.4, 1.06)
    ax.set_xlabel("number of pixels aggregated, $m$")
    ax.set_ylabel("aggregation-ambiguity index")
    ax.set_title("(b) Representational ambiguity grows with $m$")
    ax.legend(loc="lower left")
    fig.tight_layout()
    _save(fig, "fig1_scale_amplification.png")


# --------------------------------------------------------------------------
# Figure 2 -- exceedance envelope and the information ladder
# --------------------------------------------------------------------------
def fig2():
    env = pd.read_csv(RES / "exceedance_envelope.csv")
    lad = pd.read_csv(RES / "identification_ladder.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.7))

    ax = axs[0]
    ax.fill_between(env["tau"], env["lower_exact"], env["upper_exact"],
                    color=C_BAND, alpha=0.55, label="sharp envelope (marginals only)")
    ax.plot(env["tau"], env["upper_exact"], color=C_UPPER, lw=1.2)
    ax.plot(env["tau"], env["lower_exact"], color=C_LOWER, lw=1.2)
    ax.plot(env["tau"], env["comonotonic"], color=C_COMO, lw=1.5,
            label="comonotonic coupling")
    ax.plot(env["tau"], env["countermonotonic"], color=C_CTR, lw=1.5,
            label="countermonotonic coupling")
    ax.plot(env["tau"], env["independent"], color=C_IND, lw=1.5, ls="--",
            label="independent coupling")
    ax.axvline(0.5, color="#999999", lw=0.7, ls=":")
    ax.annotate("the two couplings\nswap roles here", xy=(0.53, 0.30),
                xytext=(0.60, 0.10), fontsize=7.2, color="#555555",
                arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8))
    ax.set_xlabel("threshold $\\tau$")
    ax.set_ylabel("P($\\nu_C \\geq \\tau$),  $m=4$")
    ax.set_title("(a) No coupling is extremal at all thresholds")
    ax.set_ylim(-0.04, 1.06)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
              columnspacing=1.0, handlelength=1.5)

    ax = axs[1]
    x = np.arange(len(lad))
    w = 0.28
    ax.bar(x - w, lad["level_A_width"], w, color=C_BAND, edgecolor=C_UPPER,
           lw=0.6, label="A: marginals only")
    ax.bar(x, lad["level_B_true_width"], w, color=C_LEVB, edgecolor="#3A2F5C",
           lw=0.6, alpha=0.9, label="B: + pairwise covariance matrix")
    ax.bar(x + w, np.full(len(lad), 0.012), w, color="#BBBBBB",
           edgecolor="#666666", lw=0.6, label="C: full joint law (identified)")
    ax.plot(x, lad["truth_block_tail"], "D", color=C_TRUTH, ms=5,
            label="true value on the landscape", zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"$\\tau={t:g}$" for t in lad["tau"]])
    ax.set_ylabel("width of the identified interval")
    ax.set_ylim(0, 1.12)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
    ax.set_title("(b) The identification ladder for exceedance")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
              columnspacing=1.0, handlelength=1.5)
    for xi, (a, b) in enumerate(zip(lad["level_A_width"], lad["level_B_true_width"])):
        ax.text(xi - w, a + 0.018, f"{a:.2f}", ha="center", fontsize=6.8, color=C_UPPER)
        ax.text(xi, b + 0.018, f"{b:.2f}", ha="center", fontsize=6.8, color=C_LEVB)
    fig.tight_layout()
    _save(fig, "fig2_envelope_and_ladder.png")


# --------------------------------------------------------------------------
# Figure 3 -- landscape sweep with known ground truth
# --------------------------------------------------------------------------
def fig3():
    df = pd.read_csv(RES / "landscape_sweep.csv")
    d = df[df["proportion"] == 0.5]
    g = d.groupby(["correlation_length", "m"]).mean(numeric_only=True).reset_index()
    fig, axs = plt.subplots(1, 3, figsize=(10.8, 3.5))
    ms = [4, 9, 16]
    cols = {4: "#1F4E79", 9: "#2E8B7A", 16: "#B03A5B"}

    ax = axs[0]
    s4 = g[g["m"] == 4].sort_values("correlation_length")
    ax.fill_between(s4["correlation_length"], s4["v_lower_exact"],
                    s4["v_upper_exact"], color=C_BAND, alpha=0.45,
                    label="identified set, $m=4$")
    for m in ms:
        s = g[g["m"] == m].sort_values("correlation_length")
        ax.plot(s["correlation_length"], s["v_independent"], ":", color="#B08A5A",
                lw=1.2, label="independence" if m == 4 else None)
        ax.plot(s["correlation_length"], s["truth_block_var"], "-o", ms=3,
                color=cols[m], lw=1.6, label=f"truth, $m={m}$")
    _lx(ax)
    ax.set_xlim(-1.2, 34)
    ax.set_ylim(-0.012, 0.28)
    ax.set_ylabel("Var($\\nu_C$)")
    ax.set_title("(a) Truth inside the identified bounds")
    ax.legend(loc="upper left")

    ax = axs[1]
    for m in ms:
        s = g[g["m"] == m].sort_values("correlation_length")
        ax.plot(s["correlation_length"], s["implied_reduction_truth"], "-o", ms=3,
                color=cols[m], lw=1.6, label=f"achieved, $m={m}$")
        ax.plot([-1.2, 34], [m, m], ls=":", color=cols[m], lw=1.0)
        ax.text(33, m * 1.06, f"$m={m}$", fontsize=6.8, color=cols[m], ha="right")
    _lx(ax)
    ax.set_ylim(0.9, 22)
    ax.set_yscale("log")
    ax.set_xlim(-1.2, 34)
    ax.set_yticks([1, 2, 4, 8, 16])
    ax.set_yticklabels(["1", "2", "4", "8", "16"])
    ax.set_ylabel("variance-reduction factor $R$")
    ax.set_title("(b) Independence promises $R=m$; reality delivers $\\approx 1$")
    ax.legend(loc="lower left", ncol=3, columnspacing=0.7, handlelength=1.2)

    ax = axs[2]
    s = g[g["m"] == 4].sort_values("correlation_length")
    ax.fill_between(s["correlation_length"], s["lower_tail_0.5"],
                    s["upper_tail_0.5"], color=C_BAND, alpha=0.55,
                    label="sharp envelope")
    ax.plot(s["correlation_length"], s["truth_tail_0.5"], "-o", ms=3.4,
            color="#1F4E79", lw=1.8, label="truth")
    ax.plot(s["correlation_length"], s["independent_tail_0.5"], "--",
            color=C_IND, lw=1.5, label="independence")
    ax.plot(s["correlation_length"], s["comonotonic_tail_0.5"], "-.",
            color=C_COMO, lw=1.3, label="comonotonic")
    _lx(ax)
    ax.set_xlim(-1.2, 34)
    ax.set_ylim(0.25, 1.02)
    ax.set_ylabel("P($\\nu_C \\geq 0.5$)")
    ax.set_title("(c) Exceedance, $m=4$")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2,
              columnspacing=0.8, handlelength=1.4)
    fig.tight_layout()
    _save(fig, "fig3_landscape_sweep.png")


# --------------------------------------------------------------------------
# Figure 4 -- calibration of the independence assumption
# --------------------------------------------------------------------------
def fig4():
    df = pd.read_csv(RES / "coverage_calibration.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0))
    cols = {4: "#1F4E79", 9: "#2E8B7A", 16: "#B03A5B"}

    ax = axs[0]
    for p, ls in ((0.3, "-"), (0.5, "--")):
        d = df[df["proportion"] == p]
        for m in (4, 9, 16):
            s = d[d["m"] == m].groupby("correlation_length").empirical_coverage.mean()
            ax.plot(s.index, s.values, ls, marker="o", ms=2.8, color=cols[m], lw=1.4,
                    label=f"$m={m}$, $p={p}$")
    ax.axhline(0.90, color="#444444", lw=1.0, ls=":")
    ax.text(0.2, 0.915, "nominal 90%", fontsize=7, color="#444444")
    _lx(ax)
    ax.set_ylim(-0.03, 1.08)
    ax.set_xticks(L_TICKS)
    ax.set_ylabel("empirical coverage")
    ax.set_title("(a) Coverage of a nominal 90% interval")
    ax.legend(loc="lower left", ncol=2, columnspacing=0.8, handlelength=1.3)

    ax = axs[1]
    d = df[df["proportion"] == 0.3]
    for m in (4, 9, 16):
        s = d[d["m"] == m].groupby("correlation_length")
        ax.plot(s.nominal_width.mean().index,
                (s.nominal_width.mean() / s.truth_width.mean()).values,
                "-o", ms=3, color=cols[m], lw=1.4, label=f"$m={m}$")
    ax.axhline(1.0, color="#444444", lw=1.0, ls=":")
    ax.text(0.2, 1.04, "widths agree", fontsize=7, color="#444444")
    _lx(ax)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 1.3)
    ax.set_yticks([0.3, 0.5, 0.7, 1.0])
    ax.set_yticklabels(["0.3", "0.5", "0.7", "1.0"])
    ax.set_ylabel("independence width / truth width")
    ax.set_title("(b) The interval is also the wrong shape")
    ax.legend(loc="lower left", title="$p=0.3$", ncol=3, columnspacing=0.7,
              handlelength=1.2)
    fig.tight_layout()
    _save(fig, "fig4_calibration.png")


# --------------------------------------------------------------------------
# Figure 5 -- real land cover
# --------------------------------------------------------------------------
def fig5():
    blk = pd.read_csv(RES / "real_blocks.csv")
    marg = pd.read_csv(RES / "real_marginal.csv")

    win = blk["window"].unique()[0]
    b = blk[blk["window"] == win].sort_values("m")
    mrow = marg[marg["window"] == win].iloc[0]

    fig = plt.figure(figsize=(10.6, 3.9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1, 1.15],
                          hspace=0.95, wspace=0.34)

    ax = fig.add_subplot(gs[:, 0])
    arr = None
    for f in sorted(DATA.glob("worldcover_*.npy")):
        a = np.load(f)
        if win.startswith("N48E006") and f.name.startswith("worldcover_N48E006_26000_26000"):
            arr = a
            break
    if arr is None:
        arr = np.load(sorted(DATA.glob("worldcover_*.npy"))[0])
    sub = arr[1000:1000 + 384, 1000:1000 + 384]
    sub = np.where(sub == 10, 2, np.where(sub == 50, 1, 0))
    cmap = mpl.colors.ListedColormap(["#EDE7D9", "#B0623C", "#2F6B45"])
    ax.imshow(sub, cmap=cmap, interpolation="nearest")
    ax.set_xticks(np.arange(-0.5, 384, 30), minor=True)
    ax.set_yticks(np.arange(-0.5, 384, 30), minor=True)
    ax.grid(which="minor", color="#4A4A4A", lw=0.25, alpha=0.75)
    ax.tick_params(which="both", length=0, labelbottom=False, labelleft=False)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color("#888888")
    ax.set_title("(a) 10 m reference, 30 m pixel grid", fontsize=9)
    ax.text(0.5, -0.055, "tree  /  built-up  /  other", transform=ax.transAxes,
            ha="center", fontsize=7, color="#555555")

    ax = fig.add_subplot(gs[0, 1])
    pk = [k for k in mrow.index if k.startswith("marg_p_k")]
    vals = np.array([int(k.replace("marg_p_k", "")) for k in pk]) / 9.0
    pr = np.array([float(mrow[k]) for k in pk])
    ax.bar(vals, pr, width=0.082, color="#7FA8C9", edgecolor="#2F5D8A", lw=0.6)
    ax.set_xlim(-0.06, 1.06)
    ax.set_ylim(0, 0.92)
    ax.set_xlabel("within-pixel tree fraction", labelpad=1)
    ax.set_ylabel("stored $\\Lambda$", labelpad=1)
    ax.set_title("(b) Stored marginal", fontsize=9)

    ax = fig.add_subplot(gs[1, 1])
    ax.fill_between(b["m"], b["v_lower"].clip(lower=0), b["v_upper"],
                    color=C_BAND, alpha=0.35)
    ax.plot(b["m"], b["truth_block_var"], "-o", ms=3.4, color=C_TRUTH, lw=1.7,
            label="truth")
    ax.plot(b["m"], b["v_independent"], "--s", ms=3.2, color=C_IND, lw=1.5,
            label="independence")
    ax.plot(b["m"], b["v_upper"], "--", color="#6E8296", lw=1.0)
    ax.plot(b["m"], b["v_lower"].clip(lower=0), "--", color="#6E8296", lw=1.0,
            label="bounds")
    ax.set_xscale("log", base=2)
    ax.set_xticks([4, 9, 16, 25]); ax.set_xticklabels(["4", "9", "16", "25"])
    ax.set_ylim(-0.008, 0.158)
    ax.set_xlabel("pixels per block, $m$", labelpad=1)
    ax.set_ylabel("Var($\\nu_C$)", labelpad=1)
    ax.set_title("(c) Real block variance", fontsize=9)
    ax.annotate("identified bounds", xy=(17, 0.108), xytext=(7.2, 0.142),
                fontsize=6.8, color="#5C7284", ha="center",
                arrowprops=dict(arrowstyle="-", color="#8FA3B4", lw=0.6))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=3,
              columnspacing=0.7, handlelength=1.3, fontsize=7.0)

    ax = fig.add_subplot(gs[:, 2])
    ax.plot(b["m"], b["implied_reduction_truth"], "-o", ms=4.2, color=C_ACC, lw=1.8,
            label="reduction actually achieved")
    ax.plot(b["m"], b["m"], "--", color=C_IND, lw=1.5,
            label="reduction implied by independence")
    for _, r in b.iterrows():
        ax.annotate(f"$\\times${1.0 / r['independent_understates_by']:.0f}",
                    xy=(r["m"], r["implied_reduction_truth"]),
                    xytext=(0, -15), textcoords="offset points",
                    ha="center", fontsize=7.2, color="#555555")
    ax.set_xscale("log", base=2)
    ax.set_xticks([4, 9, 16, 25]); ax.set_xticklabels(["4", "9", "16", "25"])
    ax.set_yscale("log")
    ax.set_ylim(0.85, 60)
    ax.set_yticks([1, 2, 4, 8, 16, 32])
    ax.set_yticklabels(["1", "2", "4", "8", "16", "32"])
    ax.set_xlabel("pixels per block, $m$")
    ax.set_ylabel("variance-reduction factor")
    ax.set_title("(d) The gap on real land cover", fontsize=9)
    ax.text(9.5, 1.10, "actually achieved", fontsize=7.0, color=C_ACC)
    ax.text(6.2, 13.5, "implied by independence", fontsize=7.0, color=C_IND,
            rotation=31)
    ax.text(0.03, 0.88, "labels: factor by which\nindependence understates\nthe block variance",
            transform=ax.transAxes, ha="left", fontsize=6.6, color="#777777")
    fig.suptitle(f"{win}: ESA WorldCover 10 m, tree cover as the class of interest",
                 fontsize=9.5, y=1.00)
    _save(fig, "fig5_real_landcover.png")


# --------------------------------------------------------------------------
# Figure 6 -- the value of one extra dependence parameter
# --------------------------------------------------------------------------
def fig6():
    cf = pd.read_csv(RES / "copula_family.csv")
    lad = pd.read_csv(RES / "identification_ladder.csv")
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.0))

    ax = axs[0]
    cols = {0.25: "#2F5D8A", 0.5: "#B03A5B", 0.75: "#3E7C59"}
    for tau in (0.25, 0.5, 0.75):
        s = cf[cf["tau"] == tau].sort_values("rho")
        ax.plot(s["rho"], s["exceedance"], "-", lw=1.7, color=cols[tau],
                label=f"$\\tau={tau:g}$")
    ax.axvline(-1 / 3, color="#AAAAAA", lw=0.8, ls=":")
    ax.text(-0.325, 0.06, "PSD limit\n$\\rho=-1/(m-1)$", fontsize=6.6, color="#888888")
    ax.set_xlabel("dependence parameter $\\rho$")
    ax.set_ylabel("P($\\nu_C \\geq \\tau$), $m=4$")
    ax.set_title("(a) One stored scalar spans a band")
    ax.legend(loc="center right")

    ax = axs[1]
    x = np.arange(len(lad))
    w = 0.36
    rho_cols = [c for c in lad.columns if c.startswith("level_B_rho") and
                c.endswith("_width")]
    ref = rho_cols[2] if len(rho_cols) > 2 else rho_cols[0]
    ax.bar(x - w / 2, lad["level_A_width"], w, color=C_BAND, edgecolor=C_UPPER,
           lw=0.6, label="marginals only")
    ax.bar(x + w / 2, lad[ref], w, color=C_LEVB, edgecolor="#3A2F5C", lw=0.6,
           alpha=0.9, label="+ the true $\\rho$")
    ax.set_xticks(x)
    ax.set_xticklabels([f"$\\tau={t:g}$" for t in lad["tau"]])
    ax.set_ylim(0, 0.92)
    ax.set_ylabel("width of the identified interval")
    ax.set_title("(b) One scalar shrinks, never closes, the interval")
    ax.legend(loc="upper right")
    for xi, (a, c) in enumerate(zip(lad["level_A_width"], lad[ref])):
        ax.text(xi - w / 2, a + 0.02, f"{a:.2f}", ha="center", fontsize=7)
        ax.text(xi + w / 2, c + 0.02, f"{c:.2f}", ha="center", fontsize=7)
        ax.text(xi, 0.045, f"-{100*(1-c/a):.0f}%", ha="center", fontsize=7.2,
                color="#FFFFFF")
    fig.tight_layout()
    _save(fig, "fig6_dependence_parameter.png")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
