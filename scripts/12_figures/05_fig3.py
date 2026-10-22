#!/usr/bin/env python
"""Figure 3 at half-page width, as three parts: a, b and c.

a  Per-gene gene-body mCG in W22 15 DAP endosperm, mdr1 minus WT, for gbM-PEGs, gbM-biparental
   genes and the shuffled control regions (genes with at least 5 covered CG sites and 10 reads
   in both genotypes). Each B73 gene is counted once. The medians and IQRs are asserted equal to
   the per-set tables of module 06.
b  The reciprocal cross design of the mdr1 allele-specific RNA-seq (Xu et al. 2022): the mutant
   as mother (upper pair) and as father (lower pair), each against its wild-type cross.
c  Maternal expression ratio per gene class, wild type vs mutant, with the mutant as mother
   (left) and as father (right, the control). Black bars mark each violin's interquartile range;
   the axis runs to 1.05 so medians at 1 are not cut by the top of the axes. Paired two-sided
   Wilcoxon signed-rank tests.

Sized so that, placed at 100%, the parts assemble to 3.15 in wide with 8 pt text:
  fig3p_a  W_A x H_TOP    top row, left
  fig3p_b  W_B x H_TOP    top row, right, 0.05 in gutter
  fig3p_c  3.15 x H_C     bottom row, 0.06 in below

Input   $GI_RESULTS/tables/Supplemental_Dataset_1_per_gene.csv     a, genes    (06_per_gene_methylation/01)
        $GI_WORK/06_per_gene_methylation/ctrl_{wt,mdr1}_15dap.tsv  a, control  (06_per_gene_methylation/03)
        $GI_RESULTS/tables/mutant_deltas_per_gene.csv, control_deltas.csv
                                                        a, the per-set numbers (06_per_gene_methylation/02, 04)
        $GI_WORK/08_rnaseq_allele/per_gene_maternal_fraction.tsv   c           (08_rnaseq_allele/05)
        $GI_WORK/02_gene_sets/master_table.csv                     c, gene classes
Output  $GI_RESULTS/figures/fig3p_{a,b,c}.{png,pdf}, fig3_legend_numbers.txt
        $GI_RESULTS/tables/fig3a_source_data.csv, fig3c_source_data.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Polygon

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
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
IT = "$\\mathit{%s}$"      # digits stay upright in mathtext, a known and accepted limitation

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB


WIDTH, GUTTER, ROW_GAP = 3.15, 0.05, 0.06
PAD = 0.03
CUR_W = WIDTH

# a: Figure 2e's set colours; the control grey
A_SETS = [("PEG", "gbM-PEG", "#c0392b"), ("biparental", "gbM-bip", "#f4a582")]
CTRL_C = "#bdbdbd"
A_YLIM = (-0.1, 0.2)
# b: genome colours and the mutant ring
B73_C, MO17_C, MUT_C, LINE_C = "#2f6fb0", "#c8791a", "#c0392b", "#5b6570"
MUT_RING_LW = 2.2
# c: gene sets, abbreviated as in Figure 2
WT_C = "#7f8c8d"
C_SETS = [("PEG", "gbM", "gbM-PEG"), ("PEG", "UM", "UM-PEG"),
          ("MEG", "UM", "UM-MEG"), ("MEG", "teM", "teM-MEG"),
          ("biparental", "gbM", "gbM-bip"), ("biparental", "UM", "UM-bip")]
C_SPECS = [("BM_frac", "bM_frac", "mutant as mother", "♀"),
           ("MB_frac", "Mb_frac", "mutant as father", "♂")]
C_YTOP = 1.05


def new_fig(w, h):
    global CUR_W
    CUR_W = w
    return plt.figure(figsize=(w, h))


def box(fig, x0, w, top_in, h_in):
    h = fig.get_size_inches()[1]
    return fig.add_axes([x0 / CUR_W, (h - top_in - h_in) / h, w / CUR_W, h_in / h])


def text_w(s, **kw):
    """Rendered width of a string in inches, so the cartoon is laid out on measured text."""
    fig = plt.figure(figsize=(2, 1))
    t = fig.text(0, 0, s, **kw)
    fig.canvas.draw()
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
    plt.close(fig)
    return w


def minus(v):
    return f"{v:g}".replace("-", "−")


def save(fig, name):
    """Write both formats after checking no text runs off the sheet or under a key."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    arts = list(fig.texts)
    for ax in fig.axes:
        arts += list(ax.texts)
        if ax.axison:
            arts += (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                     + [ax.xaxis.label, ax.yaxis.label, ax.title, ax._left_title])
    arts = [a for a in arts if a.get_text() and a.get_visible()]
    for a in arts:
        bb = a.get_window_extent(r)
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
            raise SystemExit(
                f"!!! {name}: {a.get_text().splitlines()[0][:44]!r} runs off the "
                f"{w:.2f} x {h:.2f} in sheet, x {bb.x0 / fig.dpi:.2f}-{bb.x1 / fig.dpi:.2f} "
                f"y {bb.y0 / fig.dpi:.2f}-{bb.y1 / fig.dpi:.2f}")
    for lg in [a.get_legend() for a in fig.axes if a.get_legend()]:
        lb = lg.get_window_extent(r)
        for a in arts:
            if a not in lg.get_texts() and lb.overlaps(a.get_window_extent(r)):
                raise SystemExit(f"{name}: the key runs into {a.get_text()[:30]!r}")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w:.2f} x {h:.2f} in")
    plt.close(fig)


def pfmt(p):
    if not np.isfinite(p):
        return "n/a"
    return "< 1e-300" if p <= 0 else (f"{p:.2g}" if p >= 0.001 else f"{p:.1e}")


notes = []

# ---------------------------------------------------------------- a, data
S = pd.read_csv(TAB / "Supplemental_Dataset_1_per_gene.csv")


def usable(tag):
    # the filter of 06_per_gene_methylation/02_mutant_deltas.py
    return (S[f"nCG_{tag}"] >= 5) & (S[f"cov_{tag}"] >= 10) & S[f"mCG_{tag}"].notna()


ok = S[usable("wt_15dap") & usable("mdr1_15dap")]
a_data, a_rows = [], []
for imp, name, colour in A_SETS:
    s = ok[(ok.imp_class == imp) & (ok.meth_class == "gbM")]
    a_data.append((name, colour, (s.mCG_mdr1_15dap - s.mCG_wt_15dap).to_numpy(),
                   s.mCG_wt_15dap.to_numpy(), s.mCG_mdr1_15dap.to_numpy()))
    a_rows += [dict(set=name, id=g, wt=w_, mdr1=m_, mdr1_minus_wt=m_ - w_) for g, w_, m_ in
               zip(s.b73_gene, s.mCG_wt_15dap, s.mCG_mdr1_15dap)]
cw = pd.read_csv(WORK / "06_per_gene_methylation" / "ctrl_wt_15dap.tsv", sep="\t", header=None,
                 names=["region", "wt"])
cm = pd.read_csv(WORK / "06_per_gene_methylation" / "ctrl_mdr1_15dap.tsv", sep="\t", header=None,
                 names=["region", "mdr1"])
c = cw.merge(cm, on="region")
a_data.append(("control", CTRL_C, (c.mdr1 - c.wt).to_numpy(), c.wt.to_numpy(),
               c.mdr1.to_numpy()))
a_rows += [dict(set="control", id=g, wt=w_, mdr1=m_, mdr1_minus_wt=m_ - w_) for g, w_, m_ in
           zip(c.region, c.wt, c.mdr1)]

# The per-set tables report mdr1 - WT (06_per_gene_methylation/02 and 04). The same genes must
# give them back.
text = pd.concat([pd.read_csv(TAB / "mutant_deltas_per_gene.csv"),
                  pd.read_csv(TAB / "control_deltas.csv")])
text = text[text.comparison == "mdr1 vs WT, 15 DAP endosperm"].set_index("gene_set")
for (name, _, v, _, _), key in zip(a_data, ("gbM-PEG", "gbM-biparental", "shuffled control")):
    q1, med, q3 = np.percentile(v * 100, [25, 50, 75])
    t = text.loc[key]
    assert len(v) == t.n and round(med, 1) == t.delta_points, (name, len(v), med, t)
    assert round(q1, 1) == t.iqr_lo and round(q3, 1) == t.iqr_hi, (name, q1, q3, t)

# ---------------------------------------------------------------- b, layout on measured text
B_HDR, B_BOX_H, B_PAIR_GAP, B_GROUP_GAP = 0.15, 0.28, 0.06, 0.12
LBL_W = max(text_w(s, fontweight="bold") for s in ("♀WT", "♀mut", "♂WT", "♂mut")) + 0.05
BOX_W = max(text_w(s, fontweight="bold") for s in ("B73", "Mo17")) + 0.08
X_W, ARR_W, END_W, END_H = 0.12, 0.16, 0.23, B_BOX_H
CMP_TXT = text_w("compare", fontweight="bold")
X_M = 0.02 + LBL_W
X_F = X_M + BOX_W + X_W
X_E = X_F + BOX_W + ARR_W
X_C = X_E + END_W + 0.10                     # centre of the compare arrow
W_B = X_C + 0.07 + 0.03 + 0.13 + 0.01       # arrow half-head, gap, turned text, margin
H_TOP = PAD + B_HDR + 2 * (2 * B_BOX_H + B_PAIR_GAP) + B_GROUP_GAP + PAD
W_A = WIDTH - GUTTER - W_B
print(f"b: label {LBL_W:.3f} box {BOX_W:.3f} compare text {CMP_TXT:.3f}; W_B {W_B:.2f} W_A {W_A:.2f}")

# ---------------------------------------------------------------- a, figure
A_LAB, A_TICK = 0.45, 0.56
fig = new_fig(W_A, H_TOP)
ax = box(fig, A_LAB, W_A - A_LAB - 0.03, PAD + 0.05, H_TOP - PAD - 0.05 - A_TICK)
parts = ax.violinplot([v for _, _, v, _, _ in a_data], positions=range(len(a_data)),
                      widths=0.75, showmedians=True, showextrema=False)
for body, (_, colour, _, _, _) in zip(parts["bodies"], a_data):
    body.set_facecolor(colour); body.set_alpha(0.75); body.set_edgecolor("none")
parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.1)
ax.axhline(0, color="0.5", ls=":", lw=0.6, zorder=0)
ax.set_xticks(range(len(a_data)))
ax.set_xticklabels([n for n, _, _, _, _ in a_data], rotation=90)
ax.tick_params(axis="x", pad=1.5)
ax.set_xlim(-0.55, len(a_data) - 0.45)
ax.set_ylim(*A_YLIM)
ticks = [-0.1, 0, 0.1, 0.2]
ax.set_yticks(ticks)
ax.set_yticklabels([minus(t) for t in ticks])
ax.set_ylabel("mCG, " + IT % "mdr1" + " − WT", labelpad=1)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
save(fig, "fig3p_a")
for name, _, v, w_, m_ in a_data:
    q1, med, q3 = np.percentile(v * 100, [25, 50, 75])
    p = float(wilcoxon(w_, m_).pvalue)
    lo, hi = int((v < A_YLIM[0]).sum()), int((v > A_YLIM[1]).sum())
    notes.append(f"a {name}: n = {len(v):,}, median mdr1 - WT {med:+.1f} points (IQR {q1:+.1f} "
                 f"to {q3:+.1f}), paired two-sided Wilcoxon p = {pfmt(p)}; "
                 f"{lo + hi} values beyond the {A_YLIM[0]} to {A_YLIM[1]} axis "
                 f"({lo} below, {hi} above).")
notes.append("a: per-gene gene-body mCG in W22 15 DAP endosperm, mdr1 minus WT (the running "
             "text's sign), genes with at least 5 covered CG sites and 10 reads in both. Control: "
             "shuffled regions with the gbM-bip length distribution, same filter.")

# ---------------------------------------------------------------- b, the cross design
ROWS = [[("♀WT", "B73", False, "Mo17", False), ("♀mut", "B73", True, "Mo17", False)],
        [("♂WT", "Mo17", False, "B73", False), ("♂mut", "Mo17", False, "B73", True)]]
fig = new_fig(W_B, H_TOP)
axb = fig.add_axes([0, 0, 1, 1])
axb.set_xlim(0, W_B)
axb.set_ylim(0, H_TOP)
axb.axis("off")


def parent(x, yc, name, is_mut):
    col = B73_C if name == "B73" else MO17_C
    axb.add_patch(FancyBboxPatch((x, yc - B_BOX_H / 2), BOX_W, B_BOX_H,
                                 boxstyle="round,pad=0,rounding_size=0.035",
                                 facecolor=col, edgecolor=MUT_C if is_mut else col,
                                 linewidth=MUT_RING_LW if is_mut else 0.6, zorder=2))
    axb.text(x + BOX_W / 2, yc + 0.058, name, ha="center", va="center", color="white",
             fontweight="bold", zorder=3)
    axb.text(x + BOX_W / 2, yc - 0.062, "mut" if is_mut else "WT", ha="center", va="center",
             color="white", zorder=3)


def endosperm(x, yc, mother, father):
    mcol = B73_C if mother == "B73" else MO17_C
    pcol = B73_C if father == "B73" else MO17_C
    axb.add_patch(FancyBboxPatch((x, yc - END_H / 2), END_W, END_H,
                                 boxstyle="round,pad=0,rounding_size=0.035",
                                 facecolor="#f6f7f9", edgecolor=LINE_C, linewidth=0.6,
                                 zorder=2))
    for i, col in enumerate((mcol, mcol, pcol)):
        axb.add_patch(FancyBboxPatch((x + 0.035, yc + 0.07 - i * 0.07 - 0.024),
                                     END_W - 0.07, 0.048,
                                     boxstyle="round,pad=0,rounding_size=0.012",
                                     facecolor=col, edgecolor="none", zorder=3))


def compare(x, y_hi, y_lo, shaft=0.05, head_len=0.075, head_w=0.13):
    """A fat double-headed block arrow, upright."""
    t, hw = shaft / 2, head_w / 2
    axb.add_patch(Polygon(
        [(x, y_hi), (x + hw, y_hi - head_len), (x + t, y_hi - head_len),
         (x + t, y_lo + head_len), (x + hw, y_lo + head_len), (x, y_lo),
         (x - hw, y_lo + head_len), (x - t, y_lo + head_len),
         (x - t, y_hi - head_len), (x - hw, y_hi - head_len)],
        closed=True, facecolor=MUT_C, edgecolor="none", zorder=4))
    axb.text(x + hw + 0.03, (y_hi + y_lo) / 2, "compare", rotation=90, ha="left",
             va="center", color=MUT_C, fontweight="bold")


y = H_TOP - PAD
hdr = y - B_HDR / 2
axb.text(X_M + BOX_W / 2, hdr, "♀", ha="center", va="center")
axb.text(X_F + BOX_W / 2, hdr, "♂", ha="center", va="center")
axb.text(X_E + END_W / 2, hdr, "endosperm", ha="center", va="center")
y -= B_HDR
for g, crosses in enumerate(ROWS):
    centres = []
    for k, (label, mo, mo_mut, fa, fa_mut) in enumerate(crosses):
        yc = y - B_BOX_H / 2
        centres.append(yc)
        axb.text(0.02, yc, label, ha="left", va="center", fontweight="bold")
        parent(X_M, yc, mo, mo_mut)
        axb.text(X_M + BOX_W + X_W / 2, yc, "×", ha="center", va="center", fontsize=9)
        parent(X_F, yc, fa, fa_mut)
        axb.add_patch(FancyArrowPatch((X_F + BOX_W + 0.03, yc), (X_E - 0.03, yc),
                                      arrowstyle="-|>", mutation_scale=6, color=LINE_C,
                                      linewidth=0.7, zorder=3))
        endosperm(X_E, yc, mo, fa)
        y -= B_BOX_H + (B_PAIR_GAP if k == 0 else 0)
    compare(X_C, centres[0] + 0.02, centres[1] - 0.02)
    y -= B_GROUP_GAP
save(fig, "fig3p_b")

# ---------------------------------------------------------------- c, data and figure
gm = pd.read_csv(WORK / "08_rnaseq_allele" / "per_gene_maternal_fraction.tsv", sep="\t")
master = pd.read_csv(WORK / "02_gene_sets" / "master_table.csv")[
    ["b73_gene", "imp_class", "meth_class"]]
gm = gm.merge(master, left_on="gene", right_on="b73_gene", how="left")

C_LAB, C_GAP, C_TITLE, C_AX, C_TICK = 0.42, 0.08, 0.17, 1.30, 0.57
C_AXW = (WIDTH - C_LAB - C_GAP - PAD) / 2
H_C = PAD + C_TITLE + C_AX + C_TICK
STEP, OFF, VW = 2.5, 1.0, 0.9
fig = new_fig(WIDTH, H_C)
c_rows = []
for j, (a_col, b_col, title, sym) in enumerate(C_SPECS):
    ax = box(fig, C_LAB + j * (C_AXW + C_GAP), C_AXW, PAD + C_TITLE, C_AX)
    ticks = []
    for i, (imp, meth, label) in enumerate(C_SETS):
        d = gm[(gm.imp_class == imp) & (gm.meth_class == meth)]
        d = d[d[a_col].notna() & d[b_col].notna()]
        x, yv = d[a_col].to_numpy(), d[b_col].to_numpy()
        pos = [i * STEP, i * STEP + OFF]
        parts = ax.violinplot([x, yv], positions=pos, widths=VW, showmedians=True,
                              showextrema=False)
        for body, colour in zip(parts["bodies"], (WT_C, MUT_C)):
            body.set_facecolor(colour); body.set_alpha(0.65); body.set_edgecolor("none")
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.1)
        iqr = []
        for p_, v in zip(pos, (x, yv)):
            q1, q3 = np.percentile(v, [25, 75])
            ax.vlines(p_, q1, q3, color="black", lw=1.6, zorder=3)
            iqr.append((q1, q3))
        try:
            pv = float(wilcoxon(x, yv).pvalue)
        except ValueError:                   # every paired difference is zero
            pv = np.nan
        ticks.append((i * STEP + OFF / 2, label))
        c_rows.append(dict(panel=title, gene_set=label, n=len(x),
                           median_wt=float(np.median(x)), iqr_wt_lo=iqr[0][0],
                           iqr_wt_hi=iqr[0][1], median_mut=float(np.median(yv)),
                           iqr_mut_lo=iqr[1][0], iqr_mut_hi=iqr[1][1],
                           median_paired_delta_pts=float(np.median((yv - x) * 100)),
                           wilcoxon_p=pv))
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks], rotation=90)
    ax.tick_params(axis="x", pad=1.5)
    ax.set_xlim(-0.75, (len(C_SETS) - 1) * STEP + OFF + 0.75)
    ax.set_ylim(0, C_YTOP)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(axis="y", alpha=0.25, lw=0.4)
    ax.set_axisbelow(True)
    ax.axhline(2 / 3, color="black", ls=":", lw=0.8)
    if j == 0:
        ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.set_ylabel("maternal expression ratio", labelpad=1)
    else:
        ax.set_yticklabels([])
    ax.set_title(title, pad=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    leg = ax.legend(handles=[Patch(facecolor=WT_C, alpha=0.65, label=f"{sym}WT"),
                       Patch(facecolor=MUT_C, alpha=0.65, label=f"{sym}mut")],
              # a white ground so the 0.75 gridline does not run through the key
              frameon=True, facecolor="white", edgecolor="none", framealpha=1,
              loc="upper left", bbox_to_anchor=(0.03, 0.87),
              handlelength=0.9, handleheight=0.7, handletextpad=0.35, labelspacing=0.2,
              borderaxespad=0, borderpad=0)
    leg.set_zorder(2)          # over the gridlines, under the spine and its ticks
save(fig, "fig3p_c")

src_c = pd.DataFrame(c_rows)
for _, r in src_c.iterrows():
    notes.append(f"c {r['panel']:<16s} {r['gene_set']:<8s} n = {int(r['n']):>5,}   "
                 f"WT {r['median_wt']:.3f} (IQR {r['iqr_wt_lo']:.3f}-{r['iqr_wt_hi']:.3f})   "
                 f"mut {r['median_mut']:.3f} (IQR {r['iqr_mut_lo']:.3f}-{r['iqr_mut_hi']:.3f})   "
                 f"paired median {r['median_paired_delta_pts']:+.2f} points, "
                 f"p = {pfmt(r['wilcoxon_p'])}")
notes.append("c: black bars are the interquartile range of each violin, black segments the "
             "median; the dotted line is 2/3. The mutant-as-father panel is the control.")

TAB.mkdir(parents=True, exist_ok=True)
pd.DataFrame(a_rows).to_csv(TAB / "fig3a_source_data.csv", index=False)
src_c.to_csv(TAB / "fig3c_source_data.csv", index=False)
(FIG / "fig3_legend_numbers.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
print("\n".join(notes))
print(f"\nassembled: {W_A:.2f} + {GUTTER} + {W_B:.2f} = {W_A + GUTTER + W_B:.2f} in wide, "
      f"top row {H_TOP:.2f} in, then c at {WIDTH} x {H_C:.2f} in; "
      f"{H_TOP + ROW_GAP + H_C:.2f} in tall with {ROW_GAP} in between. Paste at 100%.")
