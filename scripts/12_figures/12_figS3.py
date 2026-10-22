#!/usr/bin/env python
"""Supplemental Figure 4: mdr1 (a-c) and dng102 (d-f) in mature endosperm, parallel to the
15 DAP mutant comparisons of the main text.

One mutant per row, each panel comparing WT mature endosperm against that mutant: gbM-PEGs,
gbM-biparental genes and the shuffled control regions (c, f). The y label is the quantity only,
mCG; the panel titles say genes or shuffled regions. The gene panels are paired per gene from
Supplemental Dataset 1, with module 06's filter (at least 5 covered CG sites and 10 reads in
both genotypes), as Figure 3a; the control panels are paired per region from the files behind
control_deltas.csv. Every n, median and IQR is asserted against mutant_deltas_per_gene.csv and
control_deltas.csv, the tables the text quotes. Paired two-sided Wilcoxon signed-rank tests.
The shuffled control also moves, about a third as far as the gene sets and significantly at n
over 3,000, so the control panels say "smaller shift", never "flat", and the numbers are printed.

mathtext draws digits from mathtext.rm whatever they are wrapped in, so mathtext.rm is Arial
italic, which makes the 1 of mdr1 italic too; mathtext here holds gene names only.

Input   $GI_RESULTS/tables/Supplemental_Dataset_1_per_gene.csv   (06_per_gene_methylation/01)
        $GI_WORK/06_per_gene_methylation/ctrl_<genotype>.tsv     (06_per_gene_methylation/03)
        $GI_RESULTS/tables/mutant_deltas_per_gene.csv, control_deltas.csv   (06/02, 06/04)
Output  $GI_RESULTS/figures/figS3_mature.{png,pdf}
        $GI_RESULTS/tables/figS3_source_data.csv
"""
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch

# ---------------------------------------------------------------- house style
FULL_WIDTH = 6.3
BASE = 8
LETTER = 10
if "Arial" not in {f.name for f in fm.fontManager.ttflist}:
    raise SystemExit("Arial is not installed; matplotlib would silently fall back to "
                     "DejaVu Sans and the figure would not match the published one")
matplotlib.rcParams.update({
    "font.family": "Arial", "font.size": BASE,
    "axes.labelsize": BASE, "axes.titlesize": BASE,
    "xtick.labelsize": BASE, "ytick.labelsize": BASE, "legend.fontsize": BASE,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial:italic",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
IT = "$\\mathit{%s}$"

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, FIG, TAB

WT = "wt_mature"
ROWS = [("mdr1_mature", "mdr1", "#31a354"), ("dng102_mature", "dng102", "#8856a7")]
COLS = [("PEG", "gbM-PEG", "per gene"),
        ("biparental", "gbM-biparental", "per gene"),
        ("control", "shuffled control", "per region")]
WT_C = "#c0392b"

S = pd.read_csv(TAB / "Supplemental_Dataset_1_per_gene.csv")


def usable(tag):
    # 06_per_gene_methylation/02_mutant_deltas.py's filter, as in 05_fig3.py
    return (S[f"nCG_{tag}"] >= 5) & (S[f"cov_{tag}"] >= 10) & S[f"mCG_{tag}"].notna()


def ctrl(genotype):
    return pd.read_csv(WORK / "06_per_gene_methylation" / f"ctrl_{genotype}.tsv", sep="\t",
                       header=None, names=["region", "meth"])


def paired(mut, key):
    """WT and mutant mCG, paired per gene (per region for the control)."""
    if key == "control":
        c = ctrl(WT).merge(ctrl(mut), on="region", suffixes=("_wt", "_mut"))
        return c.meth_wt.to_numpy(), c.meth_mut.to_numpy()
    s = S[usable(WT) & usable(mut) & (S.imp_class == key) & (S.meth_class == "gbM")]
    return s[f"mCG_{WT}"].to_numpy(), s[f"mCG_{mut}"].to_numpy()


TEXT = pd.concat([pd.read_csv(TAB / "mutant_deltas_per_gene.csv"),
                  pd.read_csv(TAB / "control_deltas.csv")])


def letter_at(fig, ax, ch):
    """Bold panel letter at the top-left of this panel, clear of its title.

    Measured, not hand-placed: the x offset is in inches off the axes box so the gap is the
    same whatever the panel width, and the y comes from the rendered title extent. All three
    title artists are checked, since a loc="left" or loc="center" title lives in a different
    one. Call it after the layout is final; get_position() returns the default box until then.
    """
    bb = ax.get_position()
    w, h = fig.get_size_inches()
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    top = bb.y1
    for t in (ax.title, ax._left_title, ax._right_title):
        if t.get_text():
            top = max(top, inv.transform((0, t.get_window_extent(r).y1))[1])
    fig.text(bb.x0 - 0.28 / w, top + 0.03 / h, ch,
             fontsize=LETTER, fontweight="bold", va="bottom", ha="left")


fig, axes = plt.subplots(2, 3, figsize=(FULL_WIDTH, 5.20), sharey=True)
fig.subplots_adjust(left=0.092, right=0.990, top=0.848, bottom=0.235,
                    hspace=0.52, wspace=0.12)
rows = []
for ri, (mut, mlab, mcol) in enumerate(ROWS):
    for ci, (key, slab, unit) in enumerate(COLS):
        ax = axes[ri, ci]
        x, y = paired(mut, key)
        delta = (y - x) * 100.0
        p = float(wilcoxon(x, y).pvalue)
        # the text quotes these numbers citing this figure, so it has to give them back
        t = TEXT[TEXT.comparison == f"{mlab} vs WT, mature endosperm"].set_index("gene_set").loc[slab]
        q1, med, q3 = np.percentile(delta, [25, 50, 75])
        assert len(x) == t.n and round(med, 1) == t.delta_points, (mlab, slab, len(x), med, t)
        assert round(q1, 1) == t.iqr_lo and round(q3, 1) == t.iqr_hi, (mlab, slab, q1, q3, t)

        parts = ax.violinplot([x, y], positions=[0, 1], widths=0.8,
                              showmedians=True, showextrema=False)
        for body, colour in zip(parts["bodies"], (WT_C, mcol)):
            body.set_facecolor(colour); body.set_alpha(0.68); body.set_edgecolor("none")
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.2)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["WT", IT % mlab])
        ax.set_ylim(0, 1)
        ax.grid(axis="y", alpha=0.25, lw=0.5)
        # "smaller shift", never "flat": the control rises too, and at n over 3,000 that rise
        # is significant. Calling it unchanged would contradict the numbers below it.
        note = "smaller shift" if key == "control" else "elevated"
        ax.set_title(f"{slab}, {unit}\nn = {len(x):,},  {np.median(delta):+.1f} pts, "
                     f"{note}", pad=3, loc="center", linespacing=1.3)
        if ci == 0:
            ax.set_ylabel("mCG")
        rows.append(dict(panel="abcdef"[ri * 3 + ci], mutant=mlab, gene_set=slab,
                         n=len(x), median_wt=float(np.median(x)),
                         median_mutant=float(np.median(y)),
                         median_delta_points=float(np.median(delta)),
                         pct_regions_up=float(100 * (delta > 0).mean()), wilcoxon_p=p))

fig.legend(handles=[Patch(facecolor=WT_C, alpha=0.68, label="WT mature endosperm"),
                    Patch(facecolor=ROWS[0][2], alpha=0.68, label=(IT % "mdr1") + " mature"),
                    Patch(facecolor=ROWS[1][2], alpha=0.68, label=(IT % "dng102") + " mature")],
           ncol=3, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.130),
           handletextpad=0.5, columnspacing=2.0)

_notes = ("Paired per gene (per region for the shuffled control), two-sided Wilcoxon signed-rank "
          "on the weighted sum(mc)/sum(total) of each gene or region, needing at least 5 covered "
          "CG sites and 10 reads in both genotypes. pts is the median of the paired "
          "differences in percentage points. The shuffled control rises too, about a third as far as the gene sets, and at this n "
          "that rise is significant: it is a difference of degree, not a flat control.")
fig.text(0.012, 0.010, textwrap.fill(_notes, 118), va="bottom", linespacing=1.4)

# fontsize is set explicitly: matplotlib's default suptitle is 1.2x the base, which would put a
# third size on a figure that otherwise uses only 8 pt and the 10 pt panel letters.
fig.suptitle("mature endosperm: elevated gene-body methylation in the " + (IT % "mdr1") +
             " (a-c) and " + (IT % "dng102") + " (d-f) single mutants", y=0.978, fontsize=BASE)
for ri in range(len(ROWS)):
    for ci in range(len(COLS)):
        letter_at(fig, axes[ri, ci], "abcdef"[ri * 3 + ci])

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"figS3_mature.{ext}", dpi=600)
out = pd.DataFrame(rows)
out.to_csv(TAB / "figS3_source_data.csv", index=False)
print(out.to_string(index=False))
print(f"\nwrote {FIG / 'figS3_mature'}.png / .pdf at {FULL_WIDTH} in wide, Arial {BASE} pt")
