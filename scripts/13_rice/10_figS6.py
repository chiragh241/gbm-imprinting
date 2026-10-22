#!/usr/bin/env python
"""Supplemental Figure S6: rice maternal against paternal gene body CG methylation, drawn as the
arabidopsis supplement is.

  rice_matpat_metagene  a-d: maternal (red) and paternal (blue) CG methylation on Figure 2a-d's
                        axis, one panel per set, shading a 95% bootstrap interval over genes
  rice_matpat_violins   e: per gene, paternal minus maternal CG methylation (a maternal deficit
                        is positive, as in Fig. 2e), one violin per set with a box inside
  rice_matpat_legend_numbers.txt   n, medians, IQRs and p values per set

Input   $GI_RESULTS/tables/rice_matpat_per_gene.csv, rice_matpat_summary.csv,
        rice_matpat_tests.csv (07), rice_metagene_source.csv (08)
Output  $GI_RESULTS/figures/rice_matpat_violins.{png,pdf}, rice_matpat_metagene.{png,pdf},
        rice_matpat_legend_numbers.txt
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import FIG, TAB as RES
BASE = 8
if "Arial" not in {f.name for f in fm.fontManager.ttflist}:
    raise SystemExit("Arial is not installed; matplotlib would fall back to DejaVu Sans")
matplotlib.rcParams.update({
    "font.family": "Arial", "font.size": BASE,
    "axes.labelsize": BASE, "axes.titlesize": BASE,
    "xtick.labelsize": BASE, "ytick.labelsize": BASE, "legend.fontsize": BASE,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
MAT, PAT = "#c0392b", "#2c7fb8"
SETS = [("gbM-PEG", "gbM-PEG", "#c0392b"), ("gbM-biparental", "gbM-bip", "#f4a582"),
        ("UM-PEG", "UM-PEG", "#2c7fb8"), ("UM-biparental", "UM-bip", "#a6cee3")]
NUP, NIN, NDOWN = 20, 20, 20
NBIN = NUP + 2 * NIN + NDOWN
MID = NUP + NIN


def save(fig, name):
    """Write both formats after checking nothing was drawn off the sheet."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    arts = list(fig.texts)
    for ax in fig.axes:
        arts += (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                 + list(ax.texts) + [ax.xaxis.label, ax.yaxis.label, ax.title])
    for a in arts:
        if not a.get_text():
            continue
        bb = a.get_window_extent(r)
        if bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6 or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6:
            raise SystemExit(f"!!! {name}: {a.get_text().splitlines()[0][:40]!r} runs off the sheet")
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w:.2f} x {h:.2f} in")
    plt.close(fig)


def violins():
    t = pd.read_csv(RES / "rice_matpat_per_gene.csv")
    t = t[t.n_comparisons > 0]
    groups = [(i, t[t["class"] == c].pat_minus_mat.to_numpy(), c, s, col)
              for i, (c, s, col) in enumerate(SETS)]
    lo = min(v.min() for _, v, *_ in groups)
    hi = max(v.max() for _, v, *_ in groups)
    ylim = (np.floor(lo * 4) / 4, np.ceil(hi * 4) / 4)
    W, LEFT, RIGHT, PAD, XLAB, AX, TITLE = 1.90, 0.48, 0.04, 0.06, 0.62, 1.59, 0.35
    H = PAD + XLAB + AX + TITLE
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([LEFT / W, (PAD + XLAB) / H, (W - LEFT - RIGHT) / W, AX / H])
    parts = ax.violinplot([v for _, v, *_ in groups], positions=[x for x, *_ in groups],
                          widths=0.8, showextrema=False)
    for body, (_, _, _, _, col) in zip(parts["bodies"], groups):
        body.set_facecolor(col); body.set_alpha(0.8); body.set_edgecolor("none")
    ax.boxplot([v for _, v, *_ in groups], positions=[x for x, *_ in groups], widths=0.2,
               patch_artist=True, showfliers=False, manage_ticks=False,
               boxprops=dict(facecolor="white", edgecolor="black", linewidth=0.6),
               medianprops=dict(color="black", linewidth=1.1),
               whiskerprops=dict(color="black", linewidth=0.6),
               capprops=dict(color="black", linewidth=0.6))
    ax.axhline(0, color="0.5", ls=":", lw=0.6, zorder=0)
    for x, _, _, s, _ in groups:
        ax.text(x, -0.05, s, ha="center", va="top", rotation=90, transform=ax.get_xaxis_transform())
    ax.set_xticks([]); ax.set_xlim(-0.6, len(groups) - 0.4); ax.set_ylim(*ylim)
    ticks = np.arange(np.ceil(ylim[0] * 4) / 4, ylim[1] + 1e-9, 0.25)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{v:g}".replace("-", "−") for v in ticks])
    ax.grid(axis="y", alpha=0.25, lw=0.5); ax.set_axisbelow(True)
    ax.set_ylabel("pat − mat", labelpad=1)
    crosses = [name for alt, name in (("KAS", "Nip x Kas"), ("9311", "Nip x 93-11"))
               if f"NIP_in_{alt}_mat_m" in t.columns]
    ax.set_title("rice endosperm, 5 DAP\n" + ", ".join(crosses), pad=4)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    save(fig, "rice_matpat_violins")


def metagene():
    src = pd.read_csv(RES / "rice_metagene_source.csv")
    W, L_LAB, R_PAD, GAP = 6.30, 0.40, 0.04, 0.12
    PAD, TITLE, AX, KEY = 0.04, 0.16, 0.95, 0.20
    H = PAD + TITLE + AX + KEY
    pw = (W - L_LAB - R_PAD - 3 * GAP) / 4
    fig = plt.figure(figsize=(W, H))
    for k, (name, short, _) in enumerate(SETS):
        d = src[src.set == name].sort_values("bin")
        ax = fig.add_axes([(L_LAB + k * (pw + GAP)) / W, (H - PAD - TITLE - AX) / H, pw / W, AX / H])
        x = np.arange(NBIN)
        for side, colour in (("mat", MAT), ("pat", PAT)):
            m, lo, hi = (d[c].to_numpy() for c in (side, f"{side}_lo", f"{side}_hi"))
            for a, b in ((0, MID), (MID, NBIN)):
                ax.fill_between(x[a:b], lo[a:b], hi[a:b], color=colour, alpha=0.18, lw=0)
                ax.plot(x[a:b], m[a:b], color=colour, lw=0.9)
        ax.axvline(NUP - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(NBIN - NDOWN - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(MID - 0.5, color="0.75", ls=":", lw=0.6)
        ax.set_xticks([]); ax.set_xlim(0, NBIN - 1)
        ax.set_ylim(0, 1); ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"] if k == 0 else [])
        if k == 0:
            ax.set_ylabel("mCG", labelpad=1)
        ax.set_title(f"{short}, n = {int(d.genes.iloc[0]):,}", pad=2, loc="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.legend(handles=[Line2D([], [], color=MAT, lw=1.2, label="mat"),
                        Line2D([], [], color=PAT, lw=1.2, label="pat")],
               ncol=2, frameon=False, loc="lower left", bbox_to_anchor=(L_LAB / W, 0.0),
               handletextpad=0.35, columnspacing=0.8, handlelength=1.2, borderaxespad=0,
               borderpad=0.1)
    save(fig, "rice_matpat_metagene")


def numbers():
    s = pd.read_csv(RES / "rice_matpat_summary.csv")
    t = pd.read_csv(RES / "rice_matpat_tests.csv")
    lines = ["Rice 5 DAP endosperm (Tonosaki et al. 2024 WGBS), per gene pat - mat CG methylation,",
             "mean of the same-sequence comparisons; points = percentage points.", ""]
    for ms, d in s.groupby("min_sites", sort=False):
        lines.append(f">= {ms} covered CG cytosines on both sides of a comparison")
        for r in d.itertuples():
            lines.append(f"  {r.cls:15s} n {r.n:5d} of {r.genes_in_set:6d}  median {r.median_pat_minus_mat_points:+.1f} "
                         f"(IQR {r.iqr_low:+.1f} to {r.iqr_high:+.1f})  pat > mat {r.pct_pat_above_mat:.0f}%  "
                         f"Wilcoxon p {r.p_wilcoxon_mat_vs_pat:.2g}")
        for r in t[t.min_sites == ms].itertuples():
            lines.append(f"  {r.a} vs {r.b}: medians {r.median_a:+.1f} vs {r.median_b:+.1f}, "
                         f"Mann-Whitney p {r.p_mannwhitney:.2g} (n {r.n_a}, {r.n_b})")
        lines.append("")
    (FIG / "rice_matpat_legend_numbers.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    violins()
    metagene()
    numbers()
