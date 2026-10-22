#!/usr/bin/env python
"""Supplemental Figure 3: CHG (a-d) and CHH (e-h) methylation metagenes, parallel to Figure 2a-d.

Do DNGs also demethylate CHG and CHH where present in gene bodies? The estimator, the axis, the
gene sets and the tissues are Figure 2a-d's, and only the cytosine context changes: the pooled
weighted ratio per bin (per-gene CHH is mostly zero, so a median trace would sit on the axis).

The panels are meant to be read against Figure 2a-d bin for bin, so they are drawn as Figure
2a-d's 2x2 twice, CHG above CHH, and every panel is the same size as its CG counterpart. There
is no x tick text; the dashed verticals are named in the legend, as in Figure 2a-d.

The y scale is not shared between the CHG and CHH blocks, CHH being roughly an order of
magnitude lower, and UM-PEG in c and g is noisy because n = 34; both belong in the legend.

Input   $GI_WORK/04_emseq_methylation/metagene/<tissue>_metagene_bycontext_pergene.tsv
        (04_emseq_methylation/09_metagene_by_context.py)
Output  $GI_RESULTS/figures/figS2_chg_chh.{png,pdf}
        $GI_RESULTS/tables/figS2_source_data.csv
"""
from pathlib import Path
import textwrap

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D

# ---------------------------------------------------------------- house style
FULL_WIDTH = 6.3
BASE = 11      # Jonathan 2026-10-20: "increase font size"; pasted at 74%, so 11 pt reads as about 8 pt
LETTER = 13
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
from lib.paths import WORK as GI_WORK, FIG, TAB
WORK = GI_WORK / "04_emseq_methylation" / "metagene"

# axis identical to Figure 2a-d
NUP, NIN, NDOWN = 20, 20, 20
NBIN = NUP + 2 * NIN + NDOWN
HALF = NUP + NIN
BREAK = HALF - 0.5

PANELS = [("PEG__gbM", "gbM-PEG"), ("biparental__gbM", "gbM-biparental"),
          ("PEG__UM", "UM-PEG"), ("biparental__UM", "UM-biparental")]
TISSUES = [("wt_15dap", "15 DAP endosperm", "#c0392b"),
           ("embryo", "15 DAP embryo", "#2c7fb8"),
           ("leaf", "2nd leaf", "#31a354"),
           ("tassel", "premeiotic tassel", "#7b3294")]
CONTEXTS = ["CHG", "CHH"]


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
    fig.text(bb.x0 - 0.34 / w, top + 0.04 / h, ch,
             fontsize=LETTER, fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------- data
frames = []
for t, _, _ in TISSUES:
    f = WORK / f"{t}_metagene_bycontext_pergene.tsv"
    if not f.exists():
        print(f"  !! {t}: no by-context file yet, skipping")
        continue
    d = pd.read_csv(f, sep="\t", header=None,
                    names=["context", "set", "gene", "bin", "mc", "total"])
    d["tissue"] = t
    frames.append(d)
if not frames:
    raise SystemExit("no by-context files: run 04_emseq_methylation/09_metagene_by_context.py first")
D = pd.concat(frames, ignore_index=True)

agg = (D.groupby(["context", "tissue", "set", "bin"])
         .agg(mc=("mc", "sum"), total=("total", "sum"), n_genes=("gene", "nunique"))
         .reset_index())
agg["ratio"] = agg.mc / agg.total.replace(0, np.nan)

# ---------------------------------------------------------------- figure
# Four rows of two, which is Figure 2a-d's 2x2 stacked twice. Each panel is 2.7 in, the same
# width as its CG counterpart, which is the whole point of "parallel to Figure 2a-d".
H = 8.55
fig = plt.figure(figsize=(FULL_WIDTH, H))
gs = fig.add_gridspec(4, 2, left=0.118, right=0.985, top=0.893, bottom=0.222,
                      hspace=0.46, wspace=0.10)

rows, axmap = [], {}
for ci, ctx in enumerate(CONTEXTS):
    # One y limit per context block, so the four sets are comparable within a context. The two
    # blocks deliberately do NOT share it: CHH is about an order of magnitude below CHG and a
    # shared axis would flatten the CHH block onto the baseline.
    ymax = float(np.nanpercentile(agg.loc[agg.context == ctx, "ratio"], 99.5)) * 1.15
    if ctx == "CHH":
        # Jonathan 2026-10-20: "Panel g endosperm runs off the top of the axis. Increase value
        # to 1 for panels e-h?" Done as asked. g's endosperm peak (bin 16, n = 34) is 0.112,
        # so one bin of g still runs off the top at 0.10; 0.12 would hold every CHH value.
        ymax = 0.10
    for k, (key, label) in enumerate(PANELS):
        r, c = ci * 2 + k // 2, k % 2
        ax = fig.add_subplot(gs[r, c])
        axmap[(ci, k)] = ax
        nmax = 0
        for t, tlab, colour in TISSUES:
            d = (agg[(agg.context == ctx) & (agg.tissue == t) & (agg["set"] == key)]
                 .set_index("bin").reindex(range(NBIN)))
            if d.ratio.isna().all():
                continue
            nmax = max(nmax, int(np.nanmax(d.n_genes.to_numpy())))
            x = np.arange(NBIN)
            for a, b in ((0, HALF), (HALF, NBIN)):
                ax.plot(x[a:b], d.ratio.to_numpy()[a:b], color=colour, lw=1.3)
            rows.append(dict(context=ctx, tissue=t, gene_set=label,
                             n_genes=int(np.nanmax(d.n_genes.to_numpy())),
                             body_mean=float(np.nanmean(d.ratio.to_numpy()[NUP:NBIN - NDOWN]))))
        ax.axvline(NUP - 0.5, color="0.55", ls="--", lw=0.8)
        ax.axvline(NBIN - NDOWN - 0.5, color="0.55", ls="--", lw=0.8)
        ax.axvline(BREAK, color="0.75", ls=":", lw=0.8)
        ax.set_xticks([0, NUP - 0.5, BREAK, NBIN - NDOWN - 0.5, NBIN - 1])
        ax.set_xticklabels([])          # the legend names the landmarks, as in Figure 2a-d
        ax.set_xlim(0, NBIN - 1)
        ax.set_ylim(0, ymax)
        # ticks pinned: at 11 pt the auto locator switched CHG to 0.25 steps
        ax.set_yticks(np.arange(0, ymax + 1e-9, 0.02 if ctx == "CHH" else 0.2))
        ax.set_title(f"{label}, n = {nmax:,}", pad=4, loc="center")
        if c == 0:
            ax.set_ylabel(f"m{ctx}")
        else:
            ax.set_yticklabels([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

handles = [Line2D([], [], color=c, lw=1.3, label=lab) for _, lab, c in TISSUES]
handles += [Line2D([], [], color="0.55", ls="--", lw=0.8, label="TSS and polyA"),
            Line2D([], [], color="0.75", ls=":", lw=0.8, label="mid-gene break")]
# Three columns, not six: six entries in one row runs past both edges of a 6.3 in figure.
fig.legend(handles=handles, ncol=3, frameon=False, loc="lower center",
           bbox_to_anchor=(0.5, 0.118), handletextpad=0.4, columnspacing=1.6,
           labelspacing=0.45)

_notes = ("Pooled weighted ratio, sum(mc) / sum(total) over the genes in the set, the same "
          "estimator Figure 2a-d draws. The y scale is shared within a context block but not "
          "between them: CHH runs about an order of magnitude below CHG. UM-PEG, c and g, is "
          "noisy because n = 34, not because the signal is unstable.")
fig.text(0.012, 0.010, textwrap.fill(_notes, 84), va="bottom", linespacing=1.4)

# fontsize is set explicitly: matplotlib's default suptitle is 1.2x the base, which would put a
# third size on a figure that otherwise uses only 8 pt and the 10 pt panel letters.
fig.suptitle("CHG (a-d) and CHH (e-h) methylation across tissues (W22 inbred),\n"
             "same axis and gene sets as Figure 2a-d", y=0.992, fontsize=BASE)
# No CHG/CHH block headers: they landed on the a and e panel letters, and the context is
# already on every left-hand y label and in the suptitle.

for ci in range(len(CONTEXTS)):
    for k in range(len(PANELS)):
        letter_at(fig, axmap[(ci, k)], "abcdefgh"[ci * 4 + k])

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"figS2_chg_chh.{ext}", dpi=600)
src = pd.DataFrame(rows)
src.to_csv(TAB / "figS2_source_data.csv", index=False)
print(src.pivot_table(index=["context", "gene_set"], columns="tissue",
                      values="body_mean").round(4).to_string())
print(f"\nwrote {FIG / 'figS2_chg_chh'}.png / .pdf at {FULL_WIDTH} in wide, Arial {BASE} pt")
