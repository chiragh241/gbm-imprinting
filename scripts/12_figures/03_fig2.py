#!/usr/bin/env python
"""Figure 2, drawn as four parts at half-page width (a-d, e, f, g), and the supplementary figure of
maternal and paternal violins.

a-d  CG methylation metagenes of gbM-PEG, gbM-biparental, UM-PEG and UM-biparental genes in four
     W22 tissues, pooled over the genes of each set (02_fig2_metagene_table.py).
e    Per gene, paternal minus maternal gene-body mCG in B73/W22 hybrid mature endosperm
     (allele-resolved EM-seq, PRJNA759188), counts pooled over the libraries of both cross
     directions, for genes with at least 3 CG sites and 10 reads on each allele. Paired two-sided
     Wilcoxon signed-rank test, maternal vs paternal.
f    Cumulative distribution of the distance from each gene to the nearest endosperm-vs-leaf DMR,
     for five gene groups.
g    Where that DMR lies: more than 2 kb or within 2 kb upstream, overlapping the gene, within 2 kb
     or more than 2 kb downstream.
S    Maternal and paternal mCG violins of the two gene sets in e (figS_matpat_violins).

The parts are sized so that, placed at 100%, they assemble to a 3.15 x 3.96 in half-page figure
with 8 pt text:

  fig2p_ad  1.62 x 3.96    left column
  fig2p_e   1.48 x 0.95   \
  fig2p_f   1.48 x 1.00    >  right column, 0.05 in gutter, parts 0.06 in apart
  fig2p_g   1.48 x 1.44   /

f's key is a separate image, fig2p_f_key, one column on a transparent background, laid over the
bottom right of f.

Before drawing e, the per-gene pooling is checked against the per-set n and medians in
12_class_summary.csv (05_allele_methylation/11_gene_summary.py); the script stops if they differ.

Input   $GI_RESULTS/tables/fig2ad_source_data.csv          a-d  (02_fig2_metagene_table.py)
        $GI_WORK/05_allele_methylation/gene_summary/12_per_gene_allele_methylation.csv
                                                            e, S (05_allele_methylation/11)
        $GI_WORK/05_allele_methylation/gene_summary/12_class_summary.csv
                                                            check on e (05_allele_methylation/11)
        $GI_WORK/02_gene_sets/w22/*.bed                    gene sets (02_gene_sets/03)
        $GI_RESULTS/tables/fig_dmr_distance_source_data.csv f, g (07_dmr_distances/03)
Output  $GI_RESULTS/figures/fig2p_{ad,e,f,f_key,g}.{png,pdf}, figS_matpat_violins.{png,pdf}
        $GI_RESULTS/figures/fig2_legend_numbers.txt         numbers the legend quotes
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import wilcoxon

# The two columns and the height of the assembled figure. new_fig() sets CUR_W; box() and the
# figure-level legends divide by it.
AD_W, COL_W, GUTTER, TOTAL_H = 1.62, 1.48, 0.05, 3.96
CUR_W = AD_W
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
})

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB
FIG.mkdir(parents=True, exist_ok=True)

NUP, NIN, NDOWN = 20, 20, 20
NBIN = NUP + 2 * NIN + NDOWN
MID = NUP + NIN
BREAK = MID - 0.5
PROMOTER = 2000

PANELS_AD = [("PEG__gbM", "gbM-PEG"), ("biparental__gbM", "gbM-bip"),
             ("PEG__UM", "UM-PEG"), ("biparental__UM", "UM-bip")]
# Stages go in the legend: endo and embryo are 15 DAP, leaf is the second leaf, tassel premeiotic.
TISSUES = [("wt_15dap", "endo", "#c0392b"), ("embryo", "embryo", "#2c7fb8"),
           ("leaf", "leaf", "#31a354"), ("tassel", "tassel", "#7b3294")]
MAT, PAT = "#c0392b", "#2c7fb8"
# Order for f and g: the two biparental sets, then gbM-PEG, then the two groups it is
# contrasted with. MEG in Figure 1's MEG blue, UM-PEG purple.
DMR_CLASSES = [("UM-bip", "#a6cee3"), ("gbM-bip", "#f4a582"), ("gbM-PEG", "#c0392b"),
               ("UM-PEG", "#7b3294"), ("MEG", "#2166ac")]
DMR_SRC = {"gbM-PEG": "gbM-PEG", "gbM-bip": "gbM-biparental",
           "UM-PEG": "UM-PEG", "UM-bip": "UM-biparental", "MEG": "MEG"}
# Stacked bottom to top, so downstream is on top. Dark for the overlap and the two 2 Kb
# flanks, light beyond them.
ZONES = [("upstream >2 Kb", "5' >2 Kb", "#9ecae1"),
         ("promoter (<=2 Kb up)", "5' ≤2 Kb", "#2c7fb8"),
         ("gene body (overlap)", "overlap", "#000000"),
         ("3' (<=2 Kb down)", "3' ≤2 Kb", "#238b45"),
         ("downstream >2 Kb", "3' >2 Kb", "#a1d99b")]
EF_SETS = [("PEG__gbM", "gbM-PEG", "#c0392b"), ("biparental__gbM", "gbM-bip", "#f4a582")]
PAD, PART_GAP = 0.03, 0.06


def pfmt(p):
    return "< 1e-300" if p <= 0 else (f"{p:.2f}" if p >= 0.01
                                      else f"{p:.0e}".replace("e-0", "e-"))


def new_fig(w, h):
    global CUR_W
    CUR_W = w
    return plt.figure(figsize=(w, h))


def box(fig, x0, w, top_in, h_in):
    h = fig.get_size_inches()[1]
    return fig.add_axes([x0 / CUR_W, (h - top_in - h_in) / h, w / CUR_W, h_in / h])


def ticks_apart(ax, name, axis="x", gap=0.03):
    """At these widths neighbouring tick labels are the first thing to collide."""
    r = ax.figure.canvas.get_renderer()
    labs = [t for t in (ax.get_xticklabels() if axis == "x" else ax.get_yticklabels())
            if t.get_text() and t.get_visible()]
    bbs = sorted((t.get_window_extent(r) for t in labs), key=lambda b: b.x0)
    dpi = ax.figure.dpi
    for a, b in zip(bbs, bbs[1:]):
        if axis == "x" and (b.x0 - a.x1) / dpi < gap:
            raise SystemExit(f"{name}: x tick labels {(b.x0 - a.x1) / dpi:.3f} in apart")


def save(fig, name):
    """Write both formats, after checking nothing was drawn off the sheet."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    over = []
    arts = list(fig.texts)
    for ax in fig.axes:
        arts += list(ax.get_xticklabels()) + list(ax.get_yticklabels()) + list(ax.texts)
        arts += [ax.xaxis.label, ax.yaxis.label, ax.title]
    for a in arts:
        if not a.get_text() or not a.get_visible():
            continue
        bb = a.get_window_extent(r)
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
            over.append((name, a.get_text()[:44], round(bb.x0 / fig.dpi, 2),
                         round(bb.x1 / fig.dpi, 2), round(bb.y0 / fig.dpi, 2),
                         round(bb.y1 / fig.dpi, 2)))
    for lg in list(fig.legends) + [a.get_legend() for a in fig.axes if a.get_legend()]:
        bb = lg.get_window_extent(r)
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
            over.append((name, f"legend, {len(lg.get_texts())} entries",
                         round(bb.x0 / fig.dpi, 2), round(bb.x1 / fig.dpi, 2),
                         round(bb.y0 / fig.dpi, 2), round(bb.y1 / fig.dpi, 2)))
        # at half-page width the keys sit right against the tick labels, so check contact too
        for a in arts:
            if a.get_text() and a.get_visible() and bb.overlaps(a.get_window_extent(r)):
                raise SystemExit(f"{name}: the key runs into {a.get_text()[:30]!r}")
    if over:
        for n, t, x0, x1, y0, y1 in over:
            print(f"!!! {n}: {t!r} runs off the {w} x {h:.2f} in sheet, "
                  f"x {x0}-{x1} y {y0}-{y1}")
        raise SystemExit(f"{name} does not fit")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w} x {h:.2f} in")
    plt.close(fig)


# ---------------------------------------------------------------- data
st = pd.read_csv(TAB / "fig2ad_source_data.csv")

# ---- per-gene maternal and paternal CG methylation, pooled over the six libraries
MIN_SITES, MIN_READS = 3, 10

# ---------------------------------------------------------------- gene panels
per_gene = pd.read_csv(WORK / "05_allele_methylation" / "gene_summary"
                       / "12_per_gene_allele_methylation.csv")
cg = per_gene[per_gene.context == "CG"]

pooled = (cg.groupby("gene")[["mat_met", "mat_un", "pat_met", "pat_un", "n_sites"]]
            .sum().reset_index())
pooled["mat_n"] = pooled.mat_met + pooled.mat_un
pooled["pat_n"] = pooled.pat_met + pooled.pat_un
pooled = pooled[(pooled.n_sites >= MIN_SITES) &
                (pooled.mat_n >= MIN_READS) & (pooled.pat_n >= MIN_READS)]
pooled["mat"] = pooled.mat_met / pooled.mat_n
pooled["pat"] = pooled.pat_met / pooled.pat_n

gene2set, dup = {}, 0
# Sort on the NAME STRING, not the Path. Path ordering is case-insensitive on Windows and
# case-sensitive on Linux, so "biparental__gbM.bed" sorts before "biparental__UM.bed" on one
# and after it on the other. With first-wins over the 18 genes that appear in two class files,
# the other order moves seven genes into gbM-biparental (n 4,854 instead of 4,847).
for bed in sorted((WORK / "02_gene_sets" / "w22").glob("*.bed"), key=lambda b: b.name):
    for line in bed.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 4:
            # col 4 is a composite "W22|B73" id (Zm00004b000082|Zm00001eb000870); the
            # allele tables are keyed on the W22 half alone.
            g = f[3].split("|")[0]
            # A handful of W22 ids appear in more than one class file. The allele gene
            # summary (05_allele_methylation/11_gene_summary.py) resolves this with
            # drop_duplicates(keep="first") over the same sorted file order, so first wins
            # here too. Last-wins moves 12 genes (gbM-biparental n 4,851).
            if g in gene2set:
                dup += 1
            else:
                gene2set[g] = bed.stem
pooled["set"] = pooled["gene"].map(gene2set)
print(f"gene ids in more than one class file: {dup} (first wins, as in 05_allele_methylation/11_gene_summary.py)")
print(f"genes mapped to a set: {pooled['set'].notna().sum():,} of {len(pooled):,}")

# ---- reproduce 05_allele_methylation/11_gene_summary.py's per-set n and medians before drawing
# anything. Read from its table rather than typed in, so the check follows any rerun.
_s12 = pd.read_csv(WORK / "05_allele_methylation" / "gene_summary" / "12_class_summary.csv")
_s12 = _s12[(_s12.context == "CG") & (_s12.scope == "pooled")].set_index(
    ["imprint_class", "meth_class"])
EXPECT = {f"{i}__{m}": (int(_s12.loc[(i, m), "n_genes"]),
                        *_s12.loc[(i, m), ["med_mat", "med_pat", "med_delta"]].astype(float))
          for i, m in (("PEG", "gbM"), ("biparental", "gbM"))}
print("=== check against 12_class_summary.csv ===")
ok = True
for s, (n_e, mat_e, pat_e, dlt_e) in EXPECT.items():
    d = pooled[pooled["set"] == s]
    n, mat, pat = len(d), d.mat.median(), d.pat.median()
    dlt = float(np.median(d.mat.to_numpy() - d.pat.to_numpy()))
    good = (n == n_e and abs(mat - mat_e) < 5e-4 and abs(pat - pat_e) < 5e-4
            and abs(dlt - dlt_e) < 5e-4)
    ok &= good
    print(f"  {s:<18} n {n:5d} (exp {n_e:5d})  mat {mat:.4f} (exp {mat_e:.4f})  "
          f"pat {pat:.4f} (exp {pat_e:.4f})  med_delta {dlt:+.4f} (exp {dlt_e:+.4f})  "
          f"{'ok' if good else 'MISMATCH'}")
if not ok:
    sys.exit("!!! does not reproduce 12_class_summary.csv; refusing to draw a figure that "
             "disagrees with it")

pooled28 = pooled
EF = [(pooled28[pooled28["set"] == key], name, colour) for key, name, colour in EF_SETS]

dmr = pd.read_csv(TAB / "fig_dmr_distance_source_data.csv")
notes = []

# ---------------------------------------------------------------- a to d
L_LAB, MG_TITLE, MG_GAP, MG_LEG = 0.36, 0.14, 0.05, 0.34
MG_AX = (TOTAL_H - PAD - 4 * MG_TITLE - 3 * MG_GAP - MG_LEG) / 4
fig = new_fig(AD_W, TOTAL_H)
y = PAD
for k, (key, label) in enumerate(PANELS_AD):
    y += MG_TITLE
    ax = box(fig, L_LAB, AD_W - L_LAB - 0.03, y, MG_AX)
    nmax = 0
    for t, _lab, colour in TISSUES:
        d = st[(st.tissue == t) & (st.set == key)].set_index("bin").reindex(range(NBIN))
        if d.pooled.isna().all():
            continue
        nmax = max(nmax, int(np.nanmax(d.n_genes.to_numpy())))
        x = np.arange(NBIN)
        for a, b in ((0, MID), (MID, NBIN)):
            ax.plot(x[a:b], d.pooled.to_numpy()[a:b], color=colour, lw=0.8)
    ax.axvline(NUP - 0.5, color="0.55", ls="--", lw=0.6)
    ax.axvline(NBIN - NDOWN - 0.5, color="0.55", ls="--", lw=0.6)
    ax.axvline(BREAK, color="0.75", ls=":", lw=0.6)
    ax.set_xticks([])
    ax.set_xlim(0, NBIN - 1)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_ylabel("mCG", labelpad=1)
    ax.set_title(f"{label}, n = {nmax:,}", pad=1.5, loc="left")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    y += MG_AX + (MG_GAP if k < 3 else 0)
fig.legend(handles=[Line2D([], [], color=c, lw=1.2, label=lab) for _, lab, c in TISSUES],
           ncol=2, frameon=False, loc="upper left",
           # below the half of the "0" tick label that hangs under the last axis
           bbox_to_anchor=(0.02, (TOTAL_H - y - 0.07) / TOTAL_H),
           handletextpad=0.35, labelspacing=0.2, columnspacing=0.8, handlelength=1.2,
           borderaxespad=0, borderpad=0)
save(fig, "fig2p_ad")
notes.append("a-d key: endo = 15 DAP endosperm, embryo = 15 DAP embryo, leaf = second leaf, "
             "tassel = premeiotic tassel (W22).")

# ---------------------------------------------------------------- e
# One panel, both sets, per-gene paternal minus maternal. The p-values are on the panel as a
# second tick line; n and the medians go to the legend numbers.
E_LAB, E_AX, E_TICK = 0.40, 0.55, 0.34
h_e = PAD + 0.03 + E_AX + E_TICK
fig = new_fig(COL_W, h_e)
ax = box(fig, E_LAB, COL_W - E_LAB - 0.03, PAD + 0.03, E_AX)
vals = [(d.pat - d.mat).to_numpy() for d, _, _ in EF]
parts = ax.violinplot(vals, positions=[0, 1], widths=0.75, showmedians=True,
                      showextrema=False)
for body, (_, _, colour) in zip(parts["bodies"], EF):
    body.set_facecolor(colour); body.set_alpha(0.75); body.set_edgecolor("none")
parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.1)
ax.axhline(0, color="0.5", ls=":", lw=0.6, zorder=0)
labels = []
for (d, name, _), v in zip(EF, vals):
    pv = float(wilcoxon(d.mat.to_numpy(), d.pat.to_numpy()).pvalue)
    q1, med, q3 = np.percentile(v, [25, 50, 75]) * 100
    labels.append(f"{name}\np {pfmt(pv)}")
    notes.append(f"e {name}: n = {len(v):,}, median pat - mat {med:+.1f} points "
                 f"(IQR {q1:+.1f} to {q3:+.1f}), paired two-sided Wilcoxon signed-rank "
                 # an underflowed p (gbM-bip, n = 4,847) would print as "p = 0.0e+00"
                 + (f"p = {pv:.1e}." if pv > 0 else "p < 1e-300."))
ax.set_xticks([0, 1])
ax.set_xticklabels(labels, linespacing=1.1)
# -0.45 to 1.45, not -0.55 to 1.55: the gbM-bip p underflows (n = 4,847) and its label reads
# "p < 1e-300", which meets the gbM-PEG label at the wider range.
ax.set_xlim(-0.45, 1.45)
ax.set_ylim(-1, 1)
ax.set_yticks([-1, 0, 1])
ax.set_yticklabels(["−1", "0", "1"])
ax.set_ylabel("pat − mat", labelpad=1)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
fig.canvas.draw()
ticks_apart(ax, "fig2p_e")
save(fig, "fig2p_e")
notes.append("e: per-gene paternal minus maternal gene-body mCG in B73/W22 hybrid MATURE "
             "endosperm (PRJNA759188), counts pooled over both cross directions.")

# ---------------------------------------------------------------- f
F_LAB, F_TOP, F_AX, F_TICK = 0.40, 0.06, 0.63, 0.31
h_f = F_TOP + F_AX + F_TICK
fig = new_fig(COL_W, h_f)
axf = box(fig, F_LAB, COL_W - F_LAB - 0.05, F_TOP, F_AX)
handles = []
for c, colour in DMR_CLASSES:
    v = np.sort(dmr.loc[dmr.gene_class == DMR_SRC[c], "abs_dist"].to_numpy())
    yy = np.arange(1, len(v) + 1) / len(v)
    lw = 1.2 if c.endswith("PEG") else 0.9
    axf.step(np.concatenate([[0], v]), np.concatenate([[0], yy]), where="post",
             color=colour, lw=lw)
    handles.append(Line2D([], [], color=colour, lw=lw, label=c))
    notes.append(f"f, g {c}: n = {len(v):,}")
axf.set_xscale("symlog", linthresh=1000, linscale=0.45)
axf.set_xlim(0, 1e6)
axf.set_ylim(0, 1)
# 1,000 Kb at the right end cannot be labelled at this width without touching 100; the axis
# still runs to 1 Mb, where every curve has reached 1.
axf.set_xticks([0, 1e3, 1e4, 1e5])
axf.set_xticklabels(["0", "1", "10", "100"])
axf.set_yticks([0, 0.5, 1])
axf.set_yticklabels(["0", "0.5", "1"])
axf.set_xlabel("distance to DMR (Kb)", labelpad=1)
axf.set_ylabel("cum. fraction", labelpad=1)
axf.axvline(PROMOTER, color="0.5", ls=":", lw=0.7)
axf.grid(alpha=0.25, lw=0.4)
for side in ("top", "right"):
    axf.spines[side].set_visible(False)
fig.canvas.draw()
ticks_apart(axf, "fig2p_f")
save(fig, "fig2p_f")
notes.append("f: the dotted line is 2 Kb.")

# ---------------------------------------------------------------- f key
# Laid over the bottom right of f. One column, since the corner under the curves is narrower
# than a two-column key is wide; transparent, so only the entries show.
fig = new_fig(COL_W, 1.0)
fig.legend(handles=handles, ncol=1, frameon=False, loc="center", handlelength=1.1,
           handletextpad=0.4, labelspacing=0.2, borderpad=0)
for ext in ("png", "pdf"):
    fig.savefig(FIG / f"fig2p_f_key.{ext}", dpi=600, transparent=True, bbox_inches="tight",
                pad_inches=0.02)
print("wrote fig2p_f_key.png / .pdf")
plt.close(fig)

# ---------------------------------------------------------------- g
# G_TICK is the turned class names: "gbM-PEG" is 0.49 in at 8 pt, plus the tick and its pad.
G_LAB, G_TOP, G_AX, G_TICK, G_KEY = 0.40, 0.06, 0.46, 0.55, 0.41
h_g = G_TOP + G_AX + G_TICK + G_KEY
fig = new_fig(COL_W, h_g)
axg = box(fig, G_LAB, COL_W - G_LAB - 0.03, G_TOP, G_AX)
keys = [c for c, _ in DMR_CLASSES]
bottom = np.zeros(len(keys))
for zone, _short, col in ZONES:
    v = np.array([100 * (dmr.loc[dmr.gene_class == DMR_SRC[c], "zone"] == zone).mean()
                  for c in keys])
    axg.bar(range(len(keys)), v, bottom=bottom, color=col, width=0.62,
            edgecolor="white", lw=0.4)
    bottom += v
assert np.allclose(bottom, 100), "a DMR zone is missing from ZONES"
axg.set_xticks(range(len(keys)))
axg.set_xticklabels(keys, rotation=90)
axg.tick_params(axis="x", pad=1.5)
axg.set_xlim(-0.5, len(keys) - 0.5)
axg.set_ylim(0, 100)
axg.set_yticks([0, 50, 100])
axg.set_ylabel("% genes", labelpad=1)
for side in ("top", "right"):
    axg.spines[side].set_visible(False)
fig.legend(handles=[Patch(facecolor=c, label=s) for _, s, c in ZONES], ncol=2,
           frameon=False, loc="lower left", bbox_to_anchor=(0.02, 0.0),
           handlelength=1.0, handletextpad=0.4, labelspacing=0.2, columnspacing=0.8,
           borderaxespad=0, borderpad=0.1)
save(fig, "fig2p_g")
notes.append("g zones: 5' >2 Kb = upstream beyond 2 Kb, 5' ≤2 Kb = within 2 Kb upstream, "
             "overlap = DMR overlaps the gene, 3' ≤2 Kb = within 2 Kb downstream, "
             "3' >2 Kb = downstream beyond 2 Kb. DMRs are endosperm vs leaf "
             "(B73_endo-leaf_hyperDMRs.bed).")

# ---------------------------------------------------------------- supplement
# Maternal and paternal violins of the two gene sets in e, at half-page width like Figure 2.
S_W, S_LAB, S_TITLE, S_AX, S_TICK = 3.15, 0.40, 0.30, 1.05, 0.18
h_s = PAD + S_TITLE + S_AX + S_TICK
fig = new_fig(S_W, h_s)
s_w = (S_W - S_LAB - 0.05 - 0.15) / 2
for i, (d, name, _) in enumerate(EF):
    ax = box(fig, S_LAB + i * (s_w + 0.15), s_w, PAD + S_TITLE, S_AX)
    parts = ax.violinplot([d.mat.to_numpy(), d.pat.to_numpy()], positions=[0, 1],
                          widths=0.8, showmedians=True, showextrema=False)
    for body, colour in zip(parts["bodies"], (MAT, PAT)):
        body.set_facecolor(colour); body.set_alpha(0.65); body.set_edgecolor("none")
    parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.2)
    delta = float(np.median(d.pat.to_numpy() - d.mat.to_numpy()) * 100)
    pv = float(wilcoxon(d.mat.to_numpy(), d.pat.to_numpy()).pvalue)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["mat", "pat"])
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.grid(axis="y", alpha=0.25, lw=0.4)
    ax.set_axisbelow(True)
    # points are pat - mat, the same quantity as Figure 2e
    ax.set_title(f"{name}, n = {len(d):,}\n{delta:+.1f} points, p {pfmt(pv)}",
                 pad=2, loc="left", linespacing=1.25)
    if i == 0:
        ax.set_ylabel("mCG", labelpad=1)
    else:
        ax.set_yticklabels([])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
save(fig, "figS_matpat_violins")

(FIG / "fig2_legend_numbers.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
print("\n".join(notes))
right = h_e + h_f + h_g + 2 * PART_GAP
print(f"\nassembled: {AD_W} + {GUTTER} + {COL_W} = {AD_W + GUTTER + COL_W:.2f} in wide; "
      f"left column {TOTAL_H:.2f} in, right column {right:.2f} in with {PART_GAP} in between "
      "parts. paste at 100% and 8 pt stays 8 pt.")
if right > TOTAL_H + 1e-6:
    raise SystemExit(f"right column {right:.2f} in is taller than the left, {TOTAL_H} in")
