#!/usr/bin/env python
"""Supplemental Figure 1: H3K27me3 outside endosperm, in B73 developing ear and leaf (Ricci et
al. 2019).

Six gene sets as ROWS (gbM-PEG, gbM-biparental, UM-PEG, UM-biparental, UM-MEG, teM-MEG) by two
tissues as COLUMNS. Each cell is a profile of mean ChIP/input over the metagene (5 kb flanks,
scaled gene body), on a linear axis from 0 with the reference line at 1 (input level), above a
heatmap of the per-gene log2 enrichment, genes sorted by gene-body enrichment. The y limits are
shared across the six sets within a tissue, so the sets can be compared within a library; each
tissue column carries its own y numbers (ear 0 to 55, leaf 0 to 17). The heatmaps keep a log2
colour scale: a linear one would leave every gbM row one flat colour.

Genes need at least MIN_INPUT input reads. Each profile title carries the genes enriched at least
twofold over input (the number the text quotes); the source table also gives the count at a
looser 1.25-fold.

Input   $GI_WORK/09_chip_leaf_ear/per_gene_bins.csv    per-gene per-bin enrichment (09_chip_leaf_ear/08)
        $GI_WORK/09_chip_leaf_ear/figS1_per_gene.csv   per-gene pooled ChIP and input reads (09_chip_leaf_ear/05)
Output  $GI_RESULTS/figures/figS1_h3k27me3.{png,pdf}
        $GI_RESULTS/tables/figS1_source_data.csv
"""
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ---------------------------------------------------------------- house style
FULL_WIDTH = 6.3
BASE = 8
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
})

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB


MIN_INPUT = 10
NBIN, NUP, NBODY = 100, 20, 60
MAIN_LOG2 = 1.0        # twofold over input, the number the text quotes
LOOSE_LOG2 = np.log2(1.25)
# PEG next to its own methylation-matched control, then the MEGs
SETS = [("gbM-PEG", "#c0392b"), ("gbM-biparental", "#e8a5a0"),
        ("UM-PEG", "#e67e22"), ("UM-biparental", "#f2c9a0"),
        ("UM-MEG", "#2c7fb8"), ("teM-MEG", "#7b3294")]
TISSUES = [("ear", "B73 developing ear"), ("leaf", "B73 leaf")]
LETTERS = "abcdefghijkl"       # the source table's panel key; not drawn


pg = pd.read_csv(WORK / "09_chip_leaf_ear" / "per_gene_bins.csv")
pg = pg[pg.scheme == "wide"]
depth = pd.read_csv(WORK / "09_chip_leaf_ear" / "figS1_per_gene.csv")
keep = {(r.gene, r.tissue) for r in depth.itertuples() if r.input >= MIN_INPUT}

# Build every panel's matrix first: the shared y limit per tissue cannot be known until all six
# profiles in that tissue exist.
panels = {}
for t, _ in TISSUES:
    for s, _ in SETS:
        d = pg[(pg.tissue == t) & (pg["set"] == s)]
        genes = sorted({g for g in d.gene.unique() if (g, t) in keep})
        if not genes:
            panels[(t, s)] = None
            continue
        M = (d[d.gene.isin(genes)]
             .pivot_table(index="gene", columns="bin", values="enrichment")
             .reindex(columns=range(NBIN)))
        body = M.iloc[:, NUP:NUP + NBODY].mean(axis=1)
        panels[(t, s)] = M.loc[body.sort_values(ascending=False).index]

ylim = {}
for t, _ in TISSUES:
    hi = []
    for s, _ in SETS:
        M = panels[(t, s)]
        if M is None:
            continue
        prof = M.mean(axis=0).to_numpy()
        prof = prof[np.isfinite(prof)]
        if len(prof):
            hi.append(prof.max())
    # linear, so the axis starts at 0 and ChIP/input = 1 is input level
    ylim[t] = (0, 1.08 * max(hi))

# Geometry, every number in inches and the fractions derived at the end. The profiles are
# 0.70 in tall, to show the difference between UM and gbM genes, above 0.50 in heatmaps; the
# figure runs past the A4 text height, which a supplementary figure can.
# Above the top row the page has to fit, in this order upward: a two-line title (0.26 in), the
# panel letter (0.14 in) and the column header (0.14 in). At a 0.28 in top margin the letters
# for a and b were pushed clean off the canvas.
PROF_IN, HEAT_IN, INNER_GAP = 0.70, 0.50, 0.06
ROW_GAP = 0.44          # holds the next row's two-line title and its letter
TOP_IN = 0.62           # column header, then the top row's letter and title
BOTTOM_IN = 1.45        # bottom row's x labels, then the colourbar band, then the footnote
CELL_IN = PROF_IN + HEAT_IN + INNER_GAP
H = TOP_IN + len(SETS) * CELL_IN + (len(SETS) - 1) * ROW_GAP + BOTTOM_IN
TOP = 1 - TOP_IN / H
fig = plt.figure(figsize=(FULL_WIDTH, H))
# right 0.972, not 0.995: the last x tick sits on the final bin, so "+5 Kb" is centred on the
# axes edge and half of it fell off the canvas.
# hspace is a fraction of the average axes height, so both gaps convert the same way.
gs = fig.add_gridspec(len(SETS), len(TISSUES), left=0.105, right=0.972, top=TOP,
                      bottom=BOTTOM_IN / H, hspace=ROW_GAP / CELL_IN, wspace=0.13)

rows, im = [], None
for si, (s, colour) in enumerate(SETS):
    for ti, (t, tlab) in enumerate(TISSUES):
        M = panels[(t, s)]
        letter = LETTERS[si * len(TISSUES) + ti]
        inner = gs[si, ti].subgridspec(2, 1, height_ratios=[PROF_IN, HEAT_IN],
                                      hspace=INNER_GAP / ((PROF_IN + HEAT_IN) / 2))
        axp = fig.add_subplot(inner[0])
        axh = fig.add_subplot(inner[1])
        if M is None:
            for ax in (axp, axh):
                ax.set_xticks([]); ax.set_yticks([])
                ax.text(0.5, 0.5, "no genes pass\nthe input filter", ha="center",
                        va="center", transform=ax.transAxes)
            continue

        # The enriched count uses the gene-level log2(ChIP RPKM / input RPKM) from
        # figS1_per_gene.csv, not the mean of the per-bin ratios drawn in the heatmap. A mean
        # of ratios is unweighted by coverage and gave 2/60 where the text says 0/60.
        dd = depth[(depth.tissue == t) & (depth["set"] == s) & depth.gene.isin(M.index)]
        n_main = int((dd.log2_enrich > MAIN_LOG2).sum())
        n_loose = int((dd.log2_enrich > LOOSE_LOG2).sum())

        prof = M.mean(axis=0).to_numpy()
        axp.plot(range(NBIN), prof, color=colour, lw=1.4)
        axp.axhline(1, color="black", lw=0.6, alpha=0.5)
        for b in (NUP - 0.5, NUP + NBODY - 0.5):
            axp.axvline(b, color="0.55", ls="--", lw=0.8)
        axp.set_xlim(0, NBIN - 1)
        axp.set_ylim(*ylim[t])           # shared within the tissue
        axp.set_xticks([])
        # The counts ride on the title, on two lines: a rotated label on the heatmap's left edge
        # needs more height than a 0.6 in heatmap has, and on one line the title is wider than
        # a 2.5 in panel.
        axp.set_title(f"{s}\nn = {len(M):,},  {n_main} enriched",
                      loc="center", pad=3, linespacing=1.3)
        for sp in ("top", "right"):
            axp.spines[sp].set_visible(False)

        with np.errstate(divide="ignore", invalid="ignore"):
            im = axh.imshow(np.log2(np.clip(M.to_numpy(), 0.05, None)), aspect="auto",
                            cmap="RdYlBu_r", vmin=-2, vmax=4, interpolation="nearest")
        axh.set_yticks([])
        for b in (NUP - 0.5, NUP + NBODY - 0.5):
            axh.axvline(b, color="0.55", ls="--", lw=0.8)
        axh.set_xticks([0, NUP - 0.5, NUP + NBODY - 0.5, NBIN - 1])
        if si == len(SETS) - 1:
            axh.set_xticklabels(["-5 Kb", "TSS", "polyA", "+5 Kb"])
        else:
            axh.set_xticklabels([])

        rows.append(dict(panel=letter, tissue=t, gene_set=s, n_plotted=len(M),
                         n_enriched_2fold=n_main,
                         pct_enriched_2fold=round(100 * n_main / len(M), 2),
                         n_enriched_1p25fold=n_loose,
                         pct_enriched_1p25fold=round(100 * n_loose / len(M), 2),
                         median_body_enrichment=float(
                             M.iloc[:, NUP:NUP + NBODY].mean(axis=1).median())))

# One y label for the whole figure. Per-axes, the rotated "log2(ChIP/input)" is 0.95 in against
# a 0.42 in profile, so it overflowed its own panel top and bottom and ran into the panel letter.
fig.text(0.014, (TOP + BOTTOM_IN / H) / 2, "ChIP/input", rotation=90,
         va="center", ha="left")

# Column headers, once each, rather than the tissue repeated in twelve panel titles.
for ti, (t, tlab) in enumerate(TISSUES):
    box = gs[0, ti].get_position(fig)
    fig.text((box.x0 + box.x1) / 2, TOP + 0.46 / H, tlab, ha="center", va="bottom",
             fontweight="bold")

# One horizontal colourbar for the whole figure. Twelve panels do not each need their own, and
# one per tissue row put its label through the next panel's title.
# Its own band between the bottom row's x labels and the footnote; at 0.082 it sat on the text.
cax = fig.add_axes([0.36, 1.02 / H, 0.28, 0.075 / H])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("log2(ChIP/input), per gene", labelpad=2)
cb.outline.set_linewidth(0.6)

_notes = (
    f"Ricci 2019, GSE120304. Heatmap rows are genes, sorted by gene-body mean. A gene enters "
    f"only if it has at least {MIN_INPUT} pooled input reads, since the value plotted is a ratio "
    "to input. Profile y limits are shared across the six sets within a tissue, so panels can be "
    "compared down a column. Enriched means a gene-level ChIP/input above twofold.")
fig.text(0.012, 0.10 / H, textwrap.fill(_notes, 118), va="bottom", linespacing=1.4)

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"figS1_h3k27me3.{ext}", dpi=600)
src = pd.DataFrame(rows)
src.to_csv(TAB / "figS1_source_data.csv", index=False)
print(src.to_string(index=False))
print(f"\nwrote {FIG / 'figS1_h3k27me3'}.png / .pdf at {FULL_WIDTH} x {H:.2f} in; "
      f"profile {PROF_IN} in, heatmap {HEAT_IN} in")
