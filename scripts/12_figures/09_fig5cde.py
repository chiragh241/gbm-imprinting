#!/usr/bin/env python
"""Former Figure 5 c, d and e as difference violins; c and d are Supplemental Figure S7e and f.

The Tau panel 5a (its script and the arabidopsis Tau were removed), the metagene 5b
(08_fig5b.py) and the fie panel 5e are no longer in the paper; the notes below keep the former
Figure 5 lettering.

The bidirectional, dme and fie panels are 5c, 5d and 5e; the metagene is 5b (08_fig5b.py,
its own full-width row).

WHAT IS PLOTTED. Per gene, paternal minus maternal CG methylation, one violin per group, so a
maternal deficit is positive, as in Figure 2e. Each violin carries a box (quartiles, median
line, whiskers to 1.5 IQR, no outliers), as Figure 4's per-gene panels do.

GENE SETS: the consensus PEG-36, as in the former 5a (see gene_classes()).

THE ROW, adding to 6.30 in, all four 2.62 in tall:

  a  2.06    Tau violins (removed)
  c  1.50    carries the y label; bidirectional Col-0 x C24, GSE280598
  d  1.37    Col-0 (dme het) x Ler, GSE38935
  e  1.37    Ler (fie het) x Col-0, GSE38935

c has its own y range: GSE280598 is deep and its per-gene differences are narrow, while the
single Ibarra libraries spread across -1 to 1, which would flatten c to a line. n, the median,
its IQR and the p-values are in fig5_legend_numbers.txt.

Input   $GI_WORK/11_arabidopsis/10_paternal_allele_endosperm_mCG.csv  c  (11_arabidopsis/08)
        $GI_WORK/11_arabidopsis/11_ibarra_per_gene.csv   d, e  (11_arabidopsis/09)
        $GI_WORK/11_arabidopsis/06_master_table.csv                       (11_arabidopsis)
        $GI_RESULTS/tables/consensus_pegs.csv, listed_imprinted_genes.csv (11_arabidopsis/03)
Output  $GI_RESULTS/figures/fig5p_{c,d,e}.{png,pdf}
        $GI_RESULTS/figures/fig5_legend_text.txt, fig5_legend_numbers.txt
        $GI_RESULTS/tables/fig5cde_source.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import matplotlib.font_manager as fm
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB

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
    # Gene names are italicised, and mathtext defaults to a different family,
    # which would put DejaVu glyphs on an Arial figure.
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
IT = "$\\mathit{%s}$"

RES = WORK / "11_arabidopsis"
MIN_SITES_A = 10        # GSE280598 is deep: covered CG sites per allele
MIN_SITES_BC = 5        # Ibarra is shallow; as 11_arabidopsis/09_dme_fie_methylation.py
MAT, PAT = "#c0392b", "#2c7fb8"
BODY_ALPHA = 0.65       # violins and the legend swatches have to read as the same colour
CLASSES = [("gbM-PEG", 1.00), ("gbM-biparental", 0.42)]   # alpha separates the two sets


def gene_classes():
    """imprinting x methylation label per gene.

    PEGs are the consensus PEG-36, the paper's canonical set and the one 5a draws, not the
    union of the three screens in 09_mechanism_per_gene.csv (135 PEGs, 30 gbM). Methylation
    class is Bewick's, from 06_master_table.csv, as in 5a. The biparental control excludes
    every gene on any of the three screens' PEG or MEG lists, so a gene one screen called
    imprinted is never a control.
    """
    master = pd.read_csv(RES / "06_master_table.csv")
    meth = dict(zip(master.AGI, master.meth_class))
    cons = set(pd.read_csv(TAB / "consensus_pegs.csv").gene)
    cls = {g: f"{meth[g]}-PEG" for g in cons if meth.get(g) in ("gbM", "UM")}
    imprinted = set(pd.read_csv(TAB / "listed_imprinted_genes.csv").gene)
    for c in ("gbM", "UM", "teM"):
        for g in set(master[master.meth_class == c].AGI) - imprinted:
            cls[g] = f"{c}-biparental"
    return cls



W_A, W_C, W_DE = 2.06, 1.50, 1.37
LEFT_C, LEFT_DE, RIGHT = 0.48, 0.30, 0.04
# Budgeted from the bottom: pad, the genotype heading, the turned class names, axes, title.
# GENO_H is reserved in c as well, which has no heading, so the four images sit level.
PAD, GENO_H, XLAB_H, AX_H, TITLE_H = 0.06, 0.17, 0.62, 1.42, 0.35
H = PAD + GENO_H + XLAB_H + AX_H + TITLE_H
GENO_Y = -0.05 - (XLAB_H / AX_H) - 0.02
# Figure 2's class colours, so a gbM-PEG is the same red in both figures
COL = {"gbM-PEG": "#c0392b", "gbM-biparental": "#f4a582", "UM-biparental": "#a6cee3"}
SHORT = {"gbM-PEG": "gbM-PEG", "gbM-biparental": "gbM-bip", "UM-biparental": "UM-bip"}
YLIM_C, YTICKS_C = (-0.4, 0.7), [-0.25, 0, 0.25, 0.5]
YLIM_DE, YTICKS_DE = (-1, 1), [-1, -0.5, 0, 0.5, 1]


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
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
            raise SystemExit(
                f"!!! {name}: {a.get_text().splitlines()[0][:44]!r} runs off the "
                f"{w:.2f} x {h:.2f} in sheet, x {bb.x0 / fig.dpi:.2f}-{bb.x1 / fig.dpi:.2f} "
                f"y {bb.y0 / fig.dpi:.2f}-{bb.y1 / fig.dpi:.2f}")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w:.2f} x {h:.2f} in")
    plt.close(fig)


def sheet(w, left):
    fig = plt.figure(figsize=(w, H))
    ax = fig.add_axes([left / w, (PAD + GENO_H + XLAB_H) / H, (w - left - RIGHT) / w, AX_H / H])
    return fig, ax


def violins(ax, groups):
    """groups: (x position, per-gene pat - mat, class). One violin each, a box inside it."""
    vals, pos = [v for _, v, _ in groups], [x for x, _, _ in groups]
    parts = ax.violinplot(vals, positions=pos, widths=0.8, showextrema=False)
    for body, (_, _, c) in zip(parts["bodies"], groups):
        body.set_facecolor(COL[c]); body.set_alpha(0.8); body.set_edgecolor("none")
    ax.boxplot(vals, positions=pos, widths=0.2, patch_artist=True, showfliers=False,
               manage_ticks=False,
               boxprops=dict(facecolor="white", edgecolor="black", linewidth=0.6),
               medianprops=dict(color="black", linewidth=1.1),
               whiskerprops=dict(color="black", linewidth=0.6),
               capprops=dict(color="black", linewidth=0.6))
    ax.axhline(0, color="0.5", ls=":", lw=0.6, zorder=0)


def dress(ax, title, ylim, yticks, ylabel):
    ax.set_xticks([])
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{t:g}".replace("-", "−") for t in yticks])
    ax.grid(axis="y", alpha=0.25, lw=0.5)
    ax.set_axisbelow(True)
    if ylabel:
        ax.set_ylabel("pat − mat", labelpad=1)
    ax.set_title(title, loc="center", pad=4)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def stats(v, m, p):
    return dict(n=len(v), median_pat_minus_mat_points=float(np.median(v) * 100),
                iqr_low=float(np.percentile(v, 25) * 100),
                iqr_high=float(np.percentile(v, 75) * 100),
                pct_positive=float(100 * (v > 0).mean()),
                p_paired_mat_vs_pat=float(wilcoxon(m, p).pvalue))


def main():
    cls = gene_classes()
    rows, between = [], []

    # ---------------------------------------------------------------- c, bidirectional
    pa = pd.read_csv(RES / "10_paternal_allele_endosperm_mCG.csv")
    pa["cls"] = pa.AGI.map(cls)
    pa = pa[(pa.pat_nCG >= MIN_SITES_A) & (pa.mat_nCG >= MIN_SITES_A)
            & pa.pat_mCG.notna() & pa.mat_mCG.notna()]
    fig, ax = sheet(W_C, LEFT_C)
    groups = []
    for i, c in enumerate(("gbM-PEG", "gbM-biparental", "UM-biparental")):
        d = pa[pa.cls == c]
        v = (d.pat_mCG - d.mat_mCG).to_numpy()
        groups.append((i, v, c))
        rows.append(dict(panel="c", source="GSE280598 Col x C24 reciprocal", genotype="WT",
                         cls=c, **stats(v, d.mat_mCG.to_numpy(), d.pat_mCG.to_numpy())))
        ax.text(i, -0.05, SHORT[c], ha="center", va="top", rotation=90,
                transform=ax.get_xaxis_transform())
    lo_c = min(v.min() for _, v, _ in groups)
    hi_c = max(v.max() for _, v, _ in groups)
    if lo_c < YLIM_C[0] or hi_c > YLIM_C[1]:
        raise SystemExit(f"c data run {lo_c:.2f} to {hi_c:.2f}, outside the y range {YLIM_C}")
    violins(ax, groups)
    ax.set_xlim(-0.6, 2.6)
    dress(ax, "bidirectional\nCol-0 x C24", YLIM_C, YTICKS_C, ylabel=True)
    save(fig, "fig5p_c")

    # ---------------------------------------------------------------- d and e, Ibarra
    ib = pd.read_csv(RES / "11_ibarra_per_gene.csv")
    ib["cls"] = ib.AGI.map(cls)
    for stem, letter, pairs, title, source in [
        ("fig5p_d", "d",
         [("matWT_col", "patWT_ler", "WT", "maternal WT"),
          ("matDME_col", "patDME_ler", IT % "dme", "maternal dme")],
         "Col-0 (" + IT % "dme" + " het)\nx Ler (WT)", "GSE38935 Ibarra Col x Ler"),
        ("fig5p_e", "e",
         [("matWT_ler", "patWT_col", "WT", "maternal WT"),
          ("matFIE_ler", "patFIE_col", IT % "fie", "maternal fie")],
         "Ler (" + IT % "fie" + " het)\nx Col-0", "GSE38935 Ibarra Ler x Col"),
    ]:
        fig, ax = sheet(W_DE, LEFT_DE)
        groups, per = [], {}
        for i, (mat, pat, lab, plain) in enumerate(pairs):
            for j, c in enumerate(("gbM-PEG", "gbM-biparental")):
                d = ib[ib.cls == c]
                d = d[(d[mat + "_nsite"] >= MIN_SITES_BC)
                      & (d[pat + "_nsite"] >= MIN_SITES_BC)
                      & d[mat].notna() & d[pat].notna()]
                v = (d[pat] - d[mat]).to_numpy()
                x = i * 2.4 + j
                groups.append((x, v, c))
                per[(plain, c)] = d.assign(diff=d[pat] - d[mat]).set_index("AGI")["diff"]
                rows.append(dict(panel=letter, source=source, genotype=plain, cls=c,
                                 **stats(v, d[mat].to_numpy(), d[pat].to_numpy())))
                ax.text(x, -0.05, SHORT[c], ha="center", va="top", rotation=90,
                        transform=ax.get_xaxis_transform())
            ax.text(i * 2.4 + 0.5, GENO_Y, lab, ha="center", va="top",
                    fontweight="bold" if lab == "WT" else "normal",
                    transform=ax.get_xaxis_transform())
        violins(ax, groups)
        ax.set_xlim(-0.6, 2.4 + 1.6)
        dress(ax, title, YLIM_DE, YTICKS_DE, ylabel=False)
        save(fig, stem)
        # WT against the mutant, two ways: unpaired over every gene in each group, and paired
        # over the genes covered in both.
        wt_lab, mut_lab = pairs[0][3], pairs[1][3]
        for c in ("gbM-PEG", "gbM-biparental"):
            a, b = per[(wt_lab, c)], per[(mut_lab, c)]
            both = a.index.intersection(b.index)
            between.append(dict(panel=letter, cls=c, n_wt=len(a), n_mut=len(b),
                                p_mannwhitney=float(mannwhitneyu(a, b).pvalue),
                                n_paired=len(both),
                                p_paired_wilcoxon=float(wilcoxon(a[both], b[both]).pvalue)
                                if len(both) >= 5 else np.nan))

    t = pd.DataFrame(rows)
    TAB.mkdir(parents=True, exist_ok=True)
    t.to_csv(TAB / "fig5cde_source.csv", index=False)
    lines = [f"{r.panel}  {r.genotype:12s} {SHORT[r.cls]:8s} n = {r.n:>6,}   "
             f"pat - mat {r.median_pat_minus_mat_points:+.1f} points "
             f"(IQR {r.iqr_low:+.1f} to {r.iqr_high:+.1f})   paired p = {r.p_paired_mat_vs_pat:.2g}"
             for r in t.itertuples()]
    lines.append("")
    lines += [f"{b['panel']}  WT vs mutant, {SHORT[b['cls']]}: Mann-Whitney p = "
              f"{b['p_mannwhitney']:.2g} (n = {b['n_wt']} vs {b['n_mut']}); paired Wilcoxon "
              f"over the {b['n_paired']} genes covered in both, p = {b['p_paired_wilcoxon']:.2g}"
              for b in between]
    (FIG / "fig5_legend_numbers.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    notes = (
        "c is GSE280598, 3 replicates per direction; d and e are GSE38935 (Ibarra 2012), one "
        "library per cross. Each value is one gene's paternal minus maternal CG methylation, "
        "measured inside one library, so depth and conversion cancel; in c both alleles are Col, "
        "compared across the two cross directions, so parent of origin is not confounded with "
        "ecotype. The box in each violin spans the interquartile range, the line in it is the "
        "median, and the whiskers reach the furthest gene within 1.5 IQR of the box. c has its "
        "own y range because its deep libraries give narrow per-gene "
        "differences. PEGs are the consensus set (called in at least two of three screens); "
        "the biparental controls exclude every gene on any imprinted list. Class names are "
        "abbreviated on the panels: gbM-bip is gbM-biparental and UM-bip the UM-biparental "
        "control.")
    (FIG / "fig5_legend_text.txt").write_text(notes + "\n", encoding="utf-8")
    print(f"\nsingle row: a {W_A} + c {W_C} + d {W_DE} + e {W_DE} = "
          f"{W_A + W_C + 2 * W_DE:.2f} in wide, all four {H:.2f} in tall; "
          "b is its own 6.30 in row. paste at 100% and 8 pt stays 8 pt.")


if __name__ == "__main__":
    main()
