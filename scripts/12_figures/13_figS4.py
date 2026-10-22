#!/usr/bin/env python
"""Supplemental Figure 5: maternal expression ratio of biparentally expressed genes by
methylation class, at 11, 14, 17 and 21 DAP, parallel to Figure 1a.

Are gbM genes as a group biased toward maternal repression in endosperm, though not strongly
enough to be called PEGs? Each gene's ratio is the mean of its B73 x W22 and W22 x B73 maternal
ratios, for one-to-one B73/W22 syntelogs with at least 5 RPM in each direction
(lib/timecourse.py), as in Figure 1a, for the biparental genes of the master gene sets in their
gbM, UM and teM classes. Mann-Whitney U, gbM vs UM, at each stage.

Drawn in the house style of the Figure 2 parts: Arial 8 pt at full page width, 6.3 x 2.0 in;
gbM-bip and UM-bip in Figure 2's colours, teM grey; the 2/3 line solid grey and named at the
right edge, as in Figure 1a; no panel letters and no overall title, which the legend carries.

Input   $GI_DATA/maize_tables/B73_W22_maternalpreference_with_avg_rpm.csv   (00_references/02)
        $GI_WORK/02_gene_sets/b73v5_to_w22.tsv                               (02_gene_sets/02)
        $GI_WORK/02_gene_sets/master_table.csv                               (02_gene_sets/03)
Output  $GI_RESULTS/figures/figS_matpref_by_methylation.{png,pdf}
        $GI_RESULTS/tables/figS_matpref_source.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.transforms as mtransforms

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB
from lib.timecourse import averaged_ratios

BASE = 8
if "Arial" not in {f.name for f in fm.fontManager.ttflist}:
    raise SystemExit("Arial is not installed; matplotlib would silently fall back to "
                     "DejaVu Sans and the figure would not match the published one")
matplotlib.rcParams.update({
    "font.family": "Arial", "font.size": BASE,
    "axes.labelsize": BASE, "axes.titlesize": BASE,
    "xtick.labelsize": BASE, "ytick.labelsize": BASE,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

FIG.mkdir(parents=True, exist_ok=True)
RES = TAB; RES.mkdir(parents=True, exist_ok=True)

STAGES = ["11DAP", "14DAP", "17DAP", "21DAP"]     # the source table's keys; drawn as "11 DAP"
CATS = ["gbM", "Unmethylated", "TE-like"]
# gbM-bip and UM-bip as in Figure 2e-f and 3a; no main figure draws teM-bip, so a neutral grey
COLORS = {"gbM": "#f4a582", "Unmethylated": "#a6cee3", "TE-like": "#bdbdbd"}
LABEL = {"gbM": "gbM", "Unmethylated": "UM", "TE-like": "teM"}
NULL_C = "#7f8c8d"                               # Figure 1a's 2/3 line

# Full page width, so pasted at 100% the text is 8 pt. Inches: title band, axes, tick-label
# band; the y label and ticks on the left, the "2/3" name on the right.
W, PAD, TITLE, AX_H, TICKS = 6.3, 0.04, 0.30, 1.30, 0.36
LEFT, RIGHT, GAP = 0.40, 0.22, 0.10
H = PAD + TITLE + AX_H + TICKS
AX_W = (W - LEFT - RIGHT - (len(STAGES) - 1) * GAP) / len(STAGES)

RPM_MIN = 5
both, _ = averaged_ratios(RPM_MIN)
master = pd.read_csv(WORK / "02_gene_sets" / "master_table.csv")
bip = master[master.imp_class == "biparental"]
meth = dict(zip(bip.b73_gene, bip.meth_class.map(
    {"gbM": "gbM", "UM": "Unmethylated", "teM": "TE-like"})))
df = pd.DataFrame({"dataset": "B73xW22 and W22xB73, averaged",
                   "group": both.stage.astype(str) + "DAP",
                   "methylation": both.b73_gene.map(meth),
                   "mat_pref": both.mat_ratio}).dropna(subset=["methylation"])
print("stages:", sorted(df.group.unique()))
print(df.groupby(["group", "methylation"]).size().unstack().to_string())

rows = []
fig = plt.figure(figsize=(W, H))
axes = [fig.add_axes([(LEFT + k * (AX_W + GAP)) / W, TICKS / H, AX_W / W, AX_H / H])
        for k in range(len(STAGES))]
for ax, stage in zip(axes, STAGES):
    sub = df[df.group == stage]
    series = []
    for i, cat in enumerate(CATS):
        v = sub[sub.methylation == cat].mat_pref.dropna().to_numpy()
        series.append(v)
        parts = ax.violinplot([v], positions=[i], widths=0.75,
                              showmedians=True, showextrema=False)
        parts["bodies"][0].set_facecolor(COLORS[cat])
        parts["bodies"][0].set_alpha(0.75); parts["bodies"][0].set_edgecolor("none")
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.1)
        rows.append(dict(stage=stage, category=LABEL[cat], n=len(v),
                         median=float(np.median(v)) if len(v) else np.nan,
                         q1=float(np.percentile(v, 25)) if len(v) else np.nan,
                         q3=float(np.percentile(v, 75)) if len(v) else np.nan))
    # the claim is gbM ~ UM, so report that comparison explicitly. On the B73-mother rows alone
    # the sign of (gbM - UM) flipped between 14 and 17 DAP; averaged over both directions gbM is
    # below UM at all four stages.
    u = mannwhitneyu(series[0], series[1]).pvalue if len(series[0]) and len(series[1]) else np.nan
    for r in rows[-len(CATS):]:
        r["gbM_vs_UM_p"] = u
        r["gbM_minus_UM_median"] = float(np.median(series[0]) - np.median(series[1]))
    ax.axhline(2 / 3, color=NULL_C, lw=0.7, zorder=0)
    ax.set_xticks(range(len(CATS)))
    ax.set_xticklabels([f"{LABEL[c]}\nn={len(s):,}" for c, s in zip(CATS, series)],
                       linespacing=1.1)
    ax.set_xlim(-0.45, len(CATS) - 0.55)
    ax.set_title(f"{stage[:-3]} DAP\ngbM vs UM p = {u:.1e}", pad=2, linespacing=1.1)
    # The whole argument is read against the 2/3 dosage null, so zooming to the data makes
    # the direction of each shift visible; at 0-1 every violin is a thin band near 0.65.
    ax.set_ylim(0.35, 0.95)
    ax.set_yticks([0.4, 0.6, 0.8])
    if ax is not axes[0]:
        ax.set_yticklabels([])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

axes[0].set_ylabel("maternal expression ratio", labelpad=2)
axes[-1].text(1.03, 2 / 3, "2/3", va="center", ha="left",
              transform=mtransforms.blended_transform_factory(axes[-1].transAxes,
                                                              axes[-1].transData))

# At these widths neighbouring tick labels collide first; and nothing may leave the sheet.
fig.canvas.draw()
_r = fig.canvas.get_renderer()
for ax, stage in zip(axes, STAGES):
    bbs = sorted((t.get_window_extent(_r) for t in ax.get_xticklabels() if t.get_text()),
                 key=lambda b: b.x0)
    for a, b in zip(bbs, bbs[1:]):
        if (b.x0 - a.x1) / fig.dpi < 0.03:
            raise SystemExit(f"{stage}: x tick labels {(b.x0 - a.x1) / fig.dpi:.3f} in apart")
for ax in axes:
    for a in (list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + list(ax.texts)
              + [ax.yaxis.label, ax.title]):
        if not a.get_text() or not a.get_visible():
            continue
        bb = a.get_window_extent(_r)
        if bb.x0 < -0.5 or bb.x1 / fig.dpi > W + 1e-6 or bb.y0 < -0.5 or bb.y1 / fig.dpi > H + 1e-6:
            raise SystemExit(f"{a.get_text()[:40]!r} runs off the {W} x {H:.2f} in sheet")
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"figS_matpref_by_methylation.{ext}", dpi=600)
print(f"figure {W} x {H:.2f} in")

src = pd.DataFrame(rows)
src.to_csv(RES / "figS_matpref_source.csv", index=False)
print(src.to_string(index=False))
print(f"\nwrote {FIG}/figS_matpref_by_methylation.png")
