#!/usr/bin/env python
"""Figure 1: maternal expression ratio over endosperm development (a), tissue specificity by
methylation class (b), and expression of the imprinted genes across ten tissues (c).

a  Per gene and stage, the maternal expression ratio averaged over the two reciprocal B73 x W22
   crosses, for one-to-one B73/W22 syntelogs with at least RPM_MIN in each direction
   (lib/timecourse.py). PEGs and MEGs are dotted; lines mark 2/3 (the 2:1 endosperm null) and
   the 2/5 and 6/7 cutoffs.
b  Tau over the ten tissues for the genes tested for imprinting, by methylation class, with the
   PEGs and MEGs of each class dotted. Kruskal-Wallis across the three classes, and Mann-Whitney
   U for gbM-PEG vs UM-PEG.
c  log2(TPM + 1) in the ten tissues for the imprinted genes, in blocks by methylation and
   imprinting class, ordered by position within a block.

The row is budgeted in inches across the 6.3 in page width, and the script stops if any label
runs off the page or into a neighbouring panel. The numbers for the legend are printed and
written to fig1_legend_numbers.txt.

Input   $GI_DATA/maize_tables/B73_W22_maternalpreference_with_avg_rpm.csv   (00_references/02)
        $GI_WORK/02_gene_sets/b73v5_to_w22.tsv                               (02_gene_sets/02)
        $GI_WORK/01_imprinting_calls/b73_imprinting_calls_all.csv            (01_imprinting_calls)
        $GI_DATA/maize_tables/b73.all.6_38.txt       expression in ten tissues (00_references/02)
        $GI_WORK/02_gene_sets/master_table.csv       imprinting and methylation classes
Output  $GI_RESULTS/figures/fig1_compact.{png,pdf}, fig1_legend_numbers.txt
        $GI_RESULTS/tables/fig1_a_source_data.csv, fig1_tau_source_data.csv,
                           fig1_heatmap_source_data.csv
"""
import io
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from scipy.stats import kruskal, mannwhitneyu

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB
from lib.timecourse import averaged_ratios

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

FULL_WIDTH = 6.3     # A4 text width: 21.0 cm less 2.5 cm margins each side
BASE = 8             # smallest font size, pt
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
})


# PEG red and MEG blue: a red/green pair is unreadable for many colour-blind readers.
PEG_C, MEG_C = "#c0392b", "#2166ac"
STAGES = [11, 14, 17, 21]
RPM_MIN = 5
NULL_RATIO, MEG_CUT, PEG_CUT = 2 / 3, 6 / 7, 2 / 5
TISSUES = ["embryo", "endosperm", "root", "shoot", "anther",
           "base", "middle", "tip", "ear", "tassel"]
# Display labels for c. TISSUES itself is the set of column names in b73.all.6_38.txt and
# cannot be abbreviated. "endosperm" is also the longest rotated label, so shortening it is
# what keeps the bottom band at 0.50 in.
TISSUE_LABELS = ["endo" if t == "endosperm" else t for t in TISSUES]
CATS = ["gbM", "UM", "teM"]
# One neutral violin fill in a and b, so the two violin panels read as a pair and the only
# colours in the row are PEG red and MEG blue.
VIOLIN_FACE, VIOLIN_EDGE = "#c9d6e3", "#5a6b7d"
GROUPS = [("gbM", "PEG", "gbM-PEG"), ("UM", "PEG", "UM-PEG"),
          ("UM", "MEG", "UM-MEG"), ("teM", "MEG", "teM-MEG")]


def tau(values):
    v = np.clip(np.asarray(values, float), 0, None)
    if v.max() <= 0:
        return np.nan
    return float(np.sum(1.0 - v / v.max()) / (len(v) - 1))


# ---------------------------------------------------------------- data for a
# Both directions, averaged per gene; the helper also checks each direction's recomputed ratio
# against the file's own column.
filt, n_pairs = averaged_ratios(RPM_MIN)
print(f"a: {n_pairs:,} one-to-one B73/W22 pairs present in both directions")
call = pd.read_csv(WORK / "01_imprinting_calls" / "b73_imprinting_calls_all.csv").set_index("b73_gene").final_call
filt["imprint_class"] = filt.b73_gene.map(call).fillna("not_called")
filt.loc[~filt.imprint_class.isin(["MEG", "PEG"]), "imprint_class"] = "other"

# ---------------------------------------------------------------- data for b and c
txt = (DATA / "maize_tables" / "b73.all.6_38.txt").read_text()
for a_, b_ in [("Near-Core Gene", "Near-Core_Gene"), ("Core Gene", "Core_Gene"),
               ("Dispensable Gene", "Dispensable_Gene"), ("Private Gene", "Private_Gene")]:
    txt = txt.replace(a_, b_)
ref = pd.read_csv(io.StringIO(txt), sep=r"\s+", low_memory=False)
# Zm00001eb007510 and Zm00001eb390120 each carry two pan-gene IDs with identical expression.
ndup = int(ref.gene.duplicated().sum())
ref = ref.drop_duplicates(subset="gene", keep="first").reset_index(drop=True)
print(f"dropped {ndup} duplicate gene rows (two pan-gene IDs, same expression)")
for c in TISSUES:
    ref[c] = pd.to_numeric(ref[c], errors="coerce")
ref["tissue_tau"] = ref[TISSUES].apply(tau, axis=1)
master = pd.read_csv(WORK / "02_gene_sets" / "master_table.csv")[
    ["b73_gene", "imp_class", "meth_class"]].rename(columns={"b73_gene": "gene"})
ref = ref.merge(master, on="gene", how="left")
scoped = ref[ref.imp_class.notna() & ref.meth_class.isin(CATS)]

# ---------------------------------------------------------------- layout
# Everything below is budgeted in inches and converted at the end, because the three panels
# have very different label overheads and equal width_ratios would starve c.
# No titles or panel letters, so the top margin is a hairline, and the bottom band holds tick
# labels only: c's rotated tissue names are the tallest and set it.
AX_H = 1.85                          # the axes band shared by a, b and c
TOP_PAD, XTICK_H = 0.10, 0.50
H = TOP_PAD + AX_H + XTICK_H
L_MARGIN = 0.62                      # a's y label and tick numbers; 0.42 clipped the label
# 0.18 in of white space on each panel boundary, between a's 2/5 2/3 6/7 labels and b's y
# ticks, and between b's axes and c's slanted labels.
A_AX, A_RIGHT = 1.10, 0.30           # a's axes, then its 2/5 2/3 6/7 twin labels
GAP_AB = GAP_BC = 0.18
B_LAB = 0.38
# 0.93, not 0.94: with gbM at mCG >= 0.2 b's count blocks are wider ("102 PEG") and sat
# 0.079 in apart; the 0.01 in goes to b.
C_AX = 0.93
C_CBAR = 0.46                        # 0.44 pushed 'log2(TPM + 1)' 0.01 in off the page

# c's block labels, slanted at 45 degrees. Their width is MEASURED rather than guessed: a label
# lying across the neighbouring violins is not something a page-bounds check can see.
SLANT = 45
_cnt = {f"{lab}": int(((ref.meth_class == m) & (ref.imp_class == i)).sum())
        for m, i, lab in GROUPS}
C_LABELS = {lab: f"{lab} {n}" for lab, n in _cnt.items()}
_probe = plt.figure(figsize=(1, 1))
_r = _probe.canvas.get_renderer()
_w = max(_probe.text(0, 0, t).get_window_extent(_r).width / _probe.dpi
         for t in C_LABELS.values())
plt.close(_probe)
# plus 0.13 in between the labels and the heatmap
C_LAB = _w * math.cos(math.radians(SLANT)) + 0.13
print(f"c block labels: widest {_w:.2f} in flat, {C_LAB:.2f} in at {SLANT} degrees")

# b takes whatever is left, so the row always fills the page exactly.
B_AX = FULL_WIDTH - (L_MARGIN + A_AX + A_RIGHT + GAP_AB + B_LAB
                     + GAP_BC + C_LAB + C_AX + C_CBAR)
assert B_AX > 1.0, f"b's axes squeezed to {B_AX:.2f} in"
AX_BOTTOM, AX_TOP = XTICK_H, XTICK_H + AX_H

xa0 = L_MARGIN
xb0 = xa0 + A_AX + A_RIGHT + GAP_AB + B_LAB
xc0 = xb0 + B_AX + GAP_BC + C_LAB
assert xc0 + C_AX + C_CBAR <= FULL_WIDTH + 1e-6, \
    f"row is {xc0 + C_AX + C_CBAR:.3f} in, wider than the page"

fig = plt.figure(figsize=(FULL_WIDTH, H))


def box(x0, w):
    return [x0 / FULL_WIDTH, AX_BOTTOM / H, w / FULL_WIDTH, (AX_TOP - AX_BOTTOM) / H]


axA = fig.add_axes(box(xa0, A_AX))
axB = fig.add_axes(box(xb0, B_AX))
axC = fig.add_axes(box(xc0, C_AX))

# ---------------------------------------------------------------- a
parts = axA.violinplot([filt.loc[filt.stage == st, "mat_ratio"].values for st in STAGES],
                       positions=range(len(STAGES)), widths=0.84,
                       showextrema=False, showmedians=False)
for body in parts["bodies"]:
    body.set_facecolor(VIOLIN_FACE); body.set_edgecolor(VIOLIN_EDGE)
    body.set_alpha(0.9); body.set_linewidth(0.6)
for x, st in zip(range(len(STAGES)), STAGES):
    s = filt[filt.stage == st]
    axA.hlines(s.mat_ratio.median(), x - 0.21, x + 0.21, color="black", lw=1.2, zorder=3)
rng = np.random.default_rng(42)
for cls, colour, dx in (("PEG", PEG_C, -0.21), ("MEG", MEG_C, 0.21)):
    for x, st in zip(range(len(STAGES)), STAGES):
        y = filt.loc[(filt.stage == st) & (filt.imprint_class == cls), "mat_ratio"].values
        if len(y):
            axA.scatter(x + dx + rng.uniform(-0.055, 0.055, len(y)), y, s=4.0, color=colour,
                        edgecolor="white", linewidth=0.15, alpha=0.85, zorder=4)
for yv, style in ((NULL_RATIO, "-"), (MEG_CUT, ":"), (PEG_CUT, ":")):
    axA.axhline(yv, color="#7f8c8d", ls=style, lw=0.7, zorder=1)

# The named cutoffs sit on a twin axis so they cannot collide with the last violin.
refax = axA.twinx()
refax.set_position(axA.get_position())
refax.set_ylim(-0.03, 1.03)
refax.set_yticks([PEG_CUT, NULL_RATIO, MEG_CUT])
refax.set_yticklabels(["2/5", "2/3", "6/7"], color="black")
refax.tick_params(length=0, pad=2)
for side in ("top", "right", "left"):
    refax.spines[side].set_visible(False)

axA.set_xticks(range(len(STAGES)))
# The stage numbers alone, with DAP as the axis label: at 0.40 in per tick "11 DAP" and
# "14 DAP" touch each other.
axA.set_xticklabels([str(st) for st in STAGES])
axA.set_xlabel("DAP", labelpad=2)
axA.tick_params(axis="x", pad=2)
axA.set_ylabel("maternal expression ratio")
axA.set_ylim(-0.03, 1.03)
axA.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
axA.set_xlim(-0.62, len(STAGES) - 0.38)
for side in ("top", "right"):
    axA.spines[side].set_visible(False)

# ---------------------------------------------------------------- b
pools, rowsB = {}, []
dropped = int(scoped.tissue_tau.isna().sum())
for i, cat in enumerate(CATS):
    v = scoped.loc[scoped.meth_class == cat, "tissue_tau"].dropna().to_numpy()
    pools[cat] = v
    parts = axB.violinplot([v], positions=[i], widths=0.80,
                           showmedians=False, showextrema=False)
    parts["bodies"][0].set_facecolor(VIOLIN_FACE)
    parts["bodies"][0].set_alpha(0.9); parts["bodies"][0].set_edgecolor(VIOLIN_EDGE)
    parts["bodies"][0].set_linewidth(0.6)
    axB.hlines(np.median(v), i - 0.27, i + 0.27, color="black", lw=1.2, zorder=4)
    # No printed median above the violin, where it would compete with the PEG and MEG dots;
    # the numbers are in fig1_tau_source_data.csv.
    rowsB.append(dict(axis="tissue_tau", category=cat, group="all", n=len(v),
                      median_tau=float(np.median(v))))
    for impc, colour in (("PEG", PEG_C), ("MEG", MEG_C)):
        d = ref.loc[(ref.meth_class == cat) & (ref.imp_class == impc),
                    "tissue_tau"].dropna().to_numpy()
        if not len(d):
            continue
        x = i + (0.21 if impc == "PEG" else -0.21)
        axB.scatter(np.random.default_rng(0).normal(x, 0.045, len(d)), d, s=4.0,
                    color=colour, edgecolor="white", linewidth=0.15, zorder=3, alpha=0.9)
        rowsB.append(dict(axis="tissue_tau", category=cat, group=impc, n=len(d),
                          median_tau=float(np.median(d))))
kw_p = kruskal(*[pools[c] for c in CATS])[1]
gp = ref.loc[(ref.meth_class == "gbM") & (ref.imp_class == "PEG"), "tissue_tau"].dropna()
up = ref.loc[(ref.meth_class == "UM") & (ref.imp_class == "PEG"), "tissue_tau"].dropna()
pp = mannwhitneyu(gp, up)[1]
rowsB.append(dict(axis="tissue_tau", category="gbM-PEG vs UM-PEG", group="contrast",
                  n=len(gp) + len(up), median_tau=float(np.median(gp) - np.median(up)),
                  mwu_p=float(pp)))

axB.set_xticks(range(len(CATS)))
# The class name alone; the three count lines are drawn inside the panel below.
axB.set_xticklabels(CATS)
axB.set_xlim(-0.5, len(CATS) - 0.5)
axB.tick_params(axis="x", pad=2)
axB.set_ylabel("Tau specificity index")
# The floor is -0.08 rather than 0 because the count block inside the panel is a fixed 0.33 in
# of text, and at AX_H 1.85 that is 0.19 of a 0 to 1.05 range, which would reach above the UM
# violin's tail at 0.172. The assertion below checks it.
B_YLIM = (-0.08, 1.05)
axB.set_ylim(*B_YLIM)
axB.set_yticks(np.arange(0, 1.01, 0.25))
axB.grid(axis="y", alpha=0.25, lw=0.4)
axB.set_axisbelow(True)

# The counts sit inside the panel, anchored to the floor. Three 8 pt lines at linespacing 1.0
# are 0.33 in, and the lowest violin tail has to clear the top of the block. That is asserted
# below rather than eyeballed, because a violin tail behind a number is unreadable and silent.
COUNT_TOP = B_YLIM[0] + 3 * (BASE / 72.0) / AX_H * (B_YLIM[1] - B_YLIM[0])
b_counts = []
for i, cat in enumerate(CATS):
    npeg = int((ref.meth_class.eq(cat) & ref.imp_class.eq("PEG")).sum())
    nmeg = int((ref.meth_class.eq(cat) & ref.imp_class.eq("MEG")).sum())
    b_counts.append(axB.text(i, B_YLIM[0] + 0.005, f"{len(pools[cat]):,}\n{npeg} PEG\n{nmeg} MEG",
                             ha="center", va="bottom", linespacing=1.0, zorder=5))
# Centring each block on its own violin does not space them evenly, because the blocks are
# different widths: 0.39, 0.39 and 0.35 in, so the two wide ones end up closer to each other
# than the narrow one is to either. They are laid out on equal EDGE gaps instead, which moves
# each block at most 0.04 in off its violin.
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_w = [t.get_window_extent(_r).width / fig.dpi for t in b_counts]
_space = (B_AX - sum(_w)) / (len(_w) + 1)
_x, _run = [], _space
for w in _w:
    _x.append(_run + w / 2)
    _run += w + _space
_span = len(CATS)                      # xlim is -0.5 to len(CATS) - 0.5
for t, xin in zip(b_counts, _x):
    t.set_x(-0.5 + xin / B_AX * _span)
print(f"   b  even spacing: blocks {', '.join(f'{w:.2f}' for w in _w)} in, "
      f"gaps {_space:.2f} in each")

for cat in CATS:
    lo = float(np.min(pools[cat]))
    assert lo > COUNT_TOP, (f"b: the {cat} violin reaches {lo:.3f}, below the count block top "
                            f"at {COUNT_TOP:.3f}; the numbers would sit on the violin")
print(f"   b  count block top {COUNT_TOP:.3f}; lowest violin tail "
      + ", ".join(f"{c} {float(np.min(pools[c])):.3f}" for c in CATS))
for side in ("top", "right"):
    axB.spines[side].set_visible(False)

kw_txt = "< 1e-300" if kw_p == 0 else f"= {kw_p:.2g}"

# ---------------------------------------------------------------- c
imp = ref[ref.imp_class.notna()].copy()
blocks, labels, bounds = [], [], []
for meth, impc, label in GROUPS:
    g = imp[(imp.meth_class == meth) & (imp.imp_class == impc)].copy()
    g["_chr"] = pd.to_numeric(g.chr, errors="coerce")
    g = g.sort_values(["_chr", "start"])
    assert len(g) == _cnt[label], f"{label}: {len(g)} here, {_cnt[label]} when C_LAB was sized"
    blocks.append(g); labels.append(C_LABELS[label]); bounds.append(len(g))
genes = pd.concat(blocks)
M = np.log2(genes[TISSUES].to_numpy(float) + 1.0)
im = axC.imshow(M, aspect="auto", cmap="viridis", interpolation="nearest",
                vmin=0, vmax=np.nanmax(M))
axC.set_xticks(range(len(TISSUES)))
axC.set_xticklabels(TISSUE_LABELS, rotation=90)
axC.tick_params(axis="x", pad=2)
y, centres = 0, []
for n in bounds:
    centres.append(y + n / 2 - 0.5)
    y += n
    if y < len(M):
        axC.axhline(y - 0.5, color="white", lw=1.0)
# Label positions are nudged apart from the block centres where they have to be. UM-MEG and
# teM-MEG are 9 and 7 rows of 113, so their centres are 8 rows apart; slanted at 45 degrees
# that leaves 0.105 in between the two lines of text, which is less than the 0.111 in a line
# of 8 pt is tall, and they touch. The separator lines and the blocks themselves do not move.
_rows = len(M)
_in_per_row = AX_H / _rows
_need = (BASE / 72.0 + 0.02) / math.cos(math.radians(SLANT)) / _in_per_row   # in rows
lab_y = list(centres)
for _ in range(50):
    moved = False
    for i in range(len(lab_y) - 1):
        gap = lab_y[i + 1] - lab_y[i]
        if gap < _need:
            push = (_need - gap) / 2
            lab_y[i] -= push
            lab_y[i + 1] += push
            moved = True
    if not moved:
        break
for lab, c0, c1 in zip(labels, centres, lab_y):
    if abs(c1 - c0) > 0.5:
        print(f"   c  label {lab!r} nudged {c1 - c0:+.1f} rows off its block centre "
              f"({abs(c1 - c0) * _in_per_row:.3f} in) to clear its neighbour")
axC.set_yticks(lab_y)
# Slanted at 45 degrees, the longest label needs 0.44 in of width instead of 0.72, and the
# difference goes to the heatmap.
axC.set_yticklabels(labels, rotation=SLANT, ha="right", rotation_mode="anchor")
axC.tick_params(axis="y", length=0, pad=4)

# A vertical colourbar to the right of the heatmap, so no colourbar band is needed under it.
# Inset 0.12 in top and bottom: a vertical bar flush with the axes puts its end tick labels
# half above the canvas, and matplotlib clips them without saying so.
cax = fig.add_axes([(xc0 + C_AX + 0.07) / FULL_WIDTH, (AX_BOTTOM + 0.12) / H,
                    0.085 / FULL_WIDTH, (AX_H - 0.24) / H])
cb = fig.colorbar(im, cax=cax, orientation="vertical", ticks=[0, 5, 10])
# TPM, not RPM: the ten tissue columns are Zeng et al. 2023's TPM, as the Methods say.
cb.set_label("log2(TPM + 1)", labelpad=2)
cb.ax.tick_params(length=2, pad=1)
cb.outline.set_linewidth(0.6)

# ---------------------------------------------------------------- the PEG / MEG key
# In space that already exists rather than new height. The bottom band is as tall as c's
# rotated tissue names need; a's stage numbers and DAP label and b's class names use only its
# upper half, so the lower half is empty under both. The key sits low in it, centred on the gap
# between a's axes and b's, clear of a's DAP label and of b's class names.
fig.canvas.draw()
_r0 = fig.canvas.get_renderer()
_tick_bottom = min(t.get_window_extent(_r0).y0 for t in axB.get_xticklabels()
                   if t.get_text()) / fig.dpi
key = fig.legend(
    handles=[Line2D([], [], marker="o", ls="", color=PEG_C, markersize=3, label="PEG"),
             Line2D([], [], marker="o", ls="", color=MEG_C, markersize=3, label="MEG")],
    frameon=False, ncol=2, loc="upper center",
    bbox_to_anchor=(((xa0 + A_AX) + xb0) / 2 / FULL_WIDTH, 0.17 / H),
    handletextpad=0.25, columnspacing=1.2, borderpad=0, borderaxespad=0)

# ---------------------------------------------------------------- numbers for the legend
_n = "  ".join(f"{len(filt[filt.stage == st]):,}" for st in STAGES)
_m = "  ".join(str(int((filt[filt.stage == st].imprint_class == "MEG").sum())) for st in STAGES)
_p = "  ".join(str(int((filt[filt.stage == st].imprint_class == "PEG").sum())) for st in STAGES)
_dots_B = int((imp.imp_class.isin(["PEG", "MEG"]) & imp.meth_class.isin(CATS)).sum())

# The numbers the figure legend quotes, written to a file rather than drawn on the canvas.
NOTES = [
    (f"a, one-to-one B73/W22 syntelogs at {RPM_MIN} RPM or more in each cross direction, "
     f"ratio averaged over the two directions:  {_n}   (11, 14, 17, 21 DAP)", False),
    (f"a, of them MEG  {_m}  (of 27);   PEG  {_p}  (of 153)", False),
    (f"b, gbM-PEG {np.median(gp):.2f} vs UM-PEG {np.median(up):.2f}, MWU p = {pp:.1g};   "
     f"Kruskal-Wallis across the three classes p {kw_txt}", True),
    (f"b, the {len(scoped):,} genes evaluable for imprinting that carry a methylation class; "
     f"{dropped:,} excluded for Tau undefined", False),
    (f"c, {_dots_B - sum(bounds)} of the {_dots_B} genes dotted in b have no block here "
     "(gbM-MEG and teM-PEG are one gene each)", False),
]
(FIG / "fig1_legend_numbers.txt").write_text(
    "\n".join(t for t, _ in NOTES) + "\n", encoding="utf-8")
print("\nnumbers for the legend text:")
for t, _ in NOTES:
    print("   " + t)

# Nothing may be drawn outside the canvas, which has no slack anywhere. Every tick label and
# axis label is measured against the page.
fig.canvas.draw()
r = fig.canvas.get_renderer()
over = []
for ax in (axA, axB, axC, cb.ax):
    arts = list(ax.get_xticklabels()) + list(ax.get_yticklabels())
    arts += [ax.xaxis.label, ax.yaxis.label]
    arts += [t for t in b_counts] if ax is axB else []
    for a in arts:
        if not a.get_text():
            continue
        bb = a.get_window_extent(r)
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > FULL_WIDTH + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > H + 1e-6):
            over.append((a.get_text()[:40],
                         round(bb.x0 / fig.dpi, 2), round(bb.x1 / fig.dpi, 2),
                         round(bb.y0 / fig.dpi, 2), round(bb.y1 / fig.dpi, 2)))
# A page-bounds check cannot see a label lying across the neighbouring panel, so c's slanted
# block labels are checked against b's right edge as well.
# Two more collisions that a page-bounds check cannot see: c's slanted labels against each
# other (UM-MEG and teM-MEG are 9 and 7 rows apart out of 113, so their centres are close),
# and b's three count blocks against each other.
def clashes(arts, name):
    """Axis-aligned box overlap. Correct for upright text only."""
    bad = []
    for i in range(len(arts)):
        for j in range(i + 1, len(arts)):
            a, b_ = arts[i].get_window_extent(r), arts[j].get_window_extent(r)
            if a.overlaps(b_):
                bad.append(f"{name}: {arts[i].get_text()!r} overlaps {arts[j].get_text()!r}")
    return bad


# get_window_extent on rotated text returns the axis-aligned box around it, which for a 0.65 in
# label at 45 degrees is 0.53 in square. Two labels a comfortable 0.28 in apart still report as
# overlapping, so a box test here raises on figures that are fine. What matters for parallel
# slanted lines is the perpendicular distance between them.
clash = list(clashes(b_counts, "b count block"))
_bw = [a.get_window_extent(r).width / fig.dpi for a in b_counts]
_bx = [a.get_window_extent(r).x0 / fig.dpi for a in b_counts]
_gaps = [(_bx[i + 1]) - (_bx[i] + _bw[i]) for i in range(len(_bw) - 1)]
print("   b  count blocks " + ", ".join(f"{w:.3f}" for w in _bw)
      + " in wide; gaps " + ", ".join(f"{g:.3f}" for g in _gaps) + " in")
if min(_gaps) < 0.08:
    clash.append(f"b count blocks are only {min(_gaps):.2f} in apart")
_perp_min = BASE / 72.0
for i in range(len(lab_y) - 1):
    perp = (lab_y[i + 1] - lab_y[i]) * _in_per_row * math.cos(math.radians(SLANT))
    if perp < _perp_min:
        clash.append(f"c block label: {labels[i]!r} and {labels[i + 1]!r} are {perp:.3f} in "
                     f"apart, less than the {_perp_min:.3f} in a line of {BASE} pt is tall")
    else:
        print(f"   c  {labels[i]:>12} to {labels[i + 1]:<12} {perp:.3f} in apart")
for c_ in clash:
    print("!!! " + c_)

kb = key.get_window_extent(r)
_k = (kb.x0 / fig.dpi, kb.x1 / fig.dpi, kb.y0 / fig.dpi, kb.y1 / fig.dpi)
print(f"   key  x {_k[0]:.2f}-{_k[1]:.2f} in, y {_k[2]:.2f}-{_k[3]:.2f} in "
      f"(b's ticks end at {_tick_bottom:.2f}, canvas floor 0.00)")
if _k[2] < 0:
    clash.append(f"PEG/MEG key runs {-_k[2]:.3f} in below the canvas")
if _k[3] > _tick_bottom + 1e-6:
    clash.append("PEG/MEG key overlaps b's class labels")
for a in axB.get_xticklabels():
    if a.get_text() and kb.overlaps(a.get_window_extent(r)):
        clash.append(f"PEG/MEG key overlaps b tick {a.get_text()!r}")
for a in list(axA.get_xticklabels()) + [axA.xaxis.label]:
    if a.get_text() and kb.overlaps(a.get_window_extent(r)):
        clash.append(f"PEG/MEG key overlaps a's {a.get_text()!r}")
for a in axC.get_xticklabels():
    if a.get_text() and kb.overlaps(a.get_window_extent(r)):
        clash.append(f"PEG/MEG key overlaps c tick {a.get_text()!r}")

b_right = xb0 + B_AX + GAP_BC
for a in axC.get_yticklabels():
    x0 = a.get_window_extent(r).x0 / fig.dpi
    if x0 < b_right - 1e-6:
        over.append((a.get_text()[:40], round(x0, 2), round(b_right, 2), 0.0, 0.0))
        print(f"!!! c label {a.get_text()!r} starts at {x0:.2f} in, inside b which ends "
              f"at {b_right:.2f} in")
if over or clash:
    for t, x0, x1, y0, y1 in over:
        print(f"!!! runs off the {FULL_WIDTH} x {H:.2f} in page: {t!r} x {x0}-{x1} y {y0}-{y1}")
    raise SystemExit("compact layout does not hold; see the list above")

for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig1_compact.{ext}", dpi=600)
assert abs(H - (TOP_PAD + AX_H + XTICK_H)) < 1e-9, "the key changed the figure height"
print(f"wrote {FIG / 'fig1_compact.png'} at {FULL_WIDTH} x {H:.2f} in, Arial {BASE} pt")

filt[["b73_gene", "w22_gene", "stage", "rpm_b73_mother", "ratio_b73_mother", "rpm_w22_mother",
      "ratio_w22_mother", "mat_ratio", "imprint_class"]].to_csv(
    TAB / "fig1_a_source_data.csv", index=False)
pd.DataFrame(rowsB).to_csv(TAB / "fig1_tau_source_data.csv", index=False)
genes[["gene", "chr", "start", "imp_class", "meth_class"] + TISSUES].to_csv(
    TAB / "fig1_heatmap_source_data.csv", index=False)

for st in STAGES:
    s = filt[filt.stage == st]
    print(f"   a  {st:>2} DAP  n={len(s):>6}  MEG={int((s.imprint_class=='MEG').sum()):>3}/27  "
          f"PEG={int((s.imprint_class=='PEG').sum()):>3}/153  median={s.mat_ratio.median():.3f}")
# Imprinted genes outside their own cutoff at each stage.
for st in STAGES:
    s = filt[filt.stage == st]
    meg, peg = s[s.imprint_class == "MEG"], s[s.imprint_class == "PEG"]
    print(f"   a  {st:>2} DAP  MEGs below 6/7: {int((meg.mat_ratio < MEG_CUT).sum())} of {len(meg)};"
          f"  PEGs above 2/5: {int((peg.mat_ratio > PEG_CUT).sum())} of {len(peg)}")
print("   b  " + "   ".join(f"{c} {np.median(pools[c]):.2f} (n={len(pools[c]):,})" for c in CATS))
print(f"      gbM-PEG {np.median(gp):.2f} vs UM-PEG {np.median(up):.2f}, MWU p = {pp:.2g}; "
      f"Kruskal-Wallis p {kw_txt};  {dropped:,} undefined Tau")
print("   c  block sizes", {l: b for l, b in zip(labels, bounds)},
      f"\n      colour scale 0 to {np.nanmax(M):.2f}, nothing clipped")
