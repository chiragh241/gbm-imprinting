#!/usr/bin/env python
"""Supplemental Figure S2: per-gene CG methylation, 15-DAP endosperm minus embryo, one violin for
each gene set of Figure 2a-d and one for the shuffled control regions.

These are the values the text quotes with Figure 2a-d (06_per_gene_methylation/05_tissue_deltas.py):
per B73 gene from Supplemental Dataset 1, genes with at least MIN_SITES covered CG sites and
MIN_COV reads in both tissues, and the shuffled control regions (the gbM-biparental length
distribution) paired the same way with the same filter in each tissue. So the medians and IQRs
here are the text's. Figure 2a-d draw the metagenes on W22 coordinates and title each with the
most W22 genes in any bin, so their n differ slightly.

Drawing: the set colours of Figure 2f and g, the grey of Figure 3a's control, black median
segments and interquartile bars, a dotted line at no difference; Arial 8 pt, half page width.

Input   $GI_RESULTS/tables/Supplemental_Dataset_1_per_gene.csv   (06_per_gene_methylation/01)
        $GI_WORK/04_emseq_methylation/violin/{embryo,wt_15dap}_hits.tsv
                                                            (04_emseq_methylation/06_region_hits.sh)
Output  $GI_RESULTS/figures/figS_endo_embryo_mcg.{png,pdf}, figS_endo_embryo_mcg_legend_numbers.txt
        $GI_RESULTS/tables/figS_endo_embryo_mcg_source_data.csv   one row per gene or region
"""
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

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

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, FIG, TAB

V = WORK / "04_emseq_methylation" / "violin"
NAME = "figS_endo_embryo_mcg"

# Figure 2a-d's order, Figure 2f and g's colours
SETS = [("PEG", "gbM", "gbM-PEG", "#c0392b"), ("biparental", "gbM", "gbM-bip", "#f4a582"),
        ("PEG", "UM", "UM-PEG", "#2c7fb8"), ("biparental", "UM", "UM-bip", "#a6cee3")]
CTRL_C = "#bdbdbd"                  # Figure 3a's control grey
MIN_COV, MIN_SITES = 10, 5          # 05_tissue_deltas.py's filter, which is Figure 3a's
YLIM, YTICKS = (-0.4, 0.2), [-0.4, -0.2, 0, 0.2]
WIDTH, HEIGHT = 3.15, 1.95
LEFT, RIGHT, TOP, BOTTOM = 0.57, 0.03, 0.05, 0.34
MIN_LABEL_GAP = 0.03                # inches between neighbouring x labels, about one space

# 05_tissue_deltas.py's control aggregation: CG cytosines (second base G) summed per region
CTRL_AWK = ('$1=="shuffled_control" && substr($3,2,1)=="G" '
            '{mc[$2]+=$4; tot[$2]+=$5; n[$2]++} '
            'END{for(k in tot) print k"\\t"mc[k]"\\t"tot[k]"\\t"n[k]}')


def minus(v):
    return f"{v:g}".replace("-", "−")


def pfmt(p):
    return "p < 1e-300" if p <= 0 else (f"p = {p:.2g}" if p >= 0.001 else f"p = {p:.1e}")


def ctrl_regions(tag):
    """Per-region CG methylation of the shuffled control regions in one tissue, filtered."""
    path = V / f"{tag}_hits.tsv"
    out = subprocess.run(["bash", "-c", f"awk -F'\\t' '{CTRL_AWK}' '{path}'"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"awk failed on {path}: {out.stderr[:400]}")
    t = pd.read_csv(io.StringIO(out.stdout), sep="\t", header=None,
                    names=["region", "mc", "tot", "nsite"])
    t = t[(t.tot >= MIN_COV) & (t.nsite >= MIN_SITES)]
    return (t.mc / t.tot).set_axis(t.region).rename(tag)


def save(fig, name):
    """Write both formats after checking no text runs off the sheet."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    for ax in fig.axes:
        for a in (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                  + [ax.xaxis.label, ax.yaxis.label]):
            if not a.get_text():
                continue
            bb = a.get_window_extent(r)
            if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                    or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
                raise SystemExit(f"!!! {name}: {a.get_text().splitlines()[0]!r} runs off the "
                                 f"{w:.2f} x {h:.2f} in sheet")
        # five labels on a half-page axis: neighbours must not touch
        boxes = sorted((a.get_window_extent(r) for a in ax.get_xticklabels() if a.get_text()),
                       key=lambda b: b.x0)
        gaps = [(b1.x0 - b0.x1) / fig.dpi for b0, b1 in zip(boxes, boxes[1:])]
        if gaps and min(gaps) < MIN_LABEL_GAP:
            raise SystemExit(f"!!! {name}: x labels {min(gaps):.3f} in apart")
        if gaps:
            print(f"   x labels at least {min(gaps):.3f} in apart")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w:.2f} x {h:.2f} in")
    plt.close(fig)


def main():
    S = pd.read_csv(TAB / "Supplemental_Dataset_1_per_gene.csv")
    ok = np.ones(len(S), bool)
    for tag in ("wt_15dap", "embryo"):
        ok &= (S[f"nCG_{tag}"] >= MIN_SITES) & (S[f"cov_{tag}"] >= MIN_COV) & S[f"mCG_{tag}"].notna()

    data, rows, notes = [], [], []
    for imp, meth, name, colour in SETS:
        s = S[ok & (S.imp_class == imp) & (S.meth_class == meth)]
        d = (s.mCG_wt_15dap - s.mCG_embryo).to_numpy()
        data.append((name, colour, d))
        rows += [dict(set=name, id=g, mCG_embryo=e, mCG_endosperm=n_, endo_minus_embryo=n_ - e)
                 for g, e, n_ in zip(s.b73_gene, s.mCG_embryo, s.mCG_wt_15dap)]
        q1, med, q3 = np.percentile(d * 100, [25, 50, 75])
        lo, hi = int((d < YLIM[0]).sum()), int((d > YLIM[1]).sum())
        total = int(((S.imp_class == imp) & (S.meth_class == meth)).sum())
        notes.append(f"{name}: n = {len(d):,} of {total:,} in the set; median mCG embryo "
                     f"{s.mCG_embryo.median() * 100:.1f}%, endosperm "
                     f"{s.mCG_wt_15dap.median() * 100:.1f}%; median endosperm - embryo "
                     f"{med:+.1f} points (IQR {q1:+.1f} to {q3:+.1f}), paired two-sided Wilcoxon "
                     f"{pfmt(float(wilcoxon(s.mCG_wt_15dap, s.mCG_embryo).pvalue))}; "
                     f"{lo + hi} values beyond the {YLIM[0]} to {YLIM[1]} axis ({lo} below, "
                     f"{hi} above).")

    # the control regions, filtered the same way in each tissue and paired on the regions in both
    c = pd.concat([ctrl_regions("embryo"), ctrl_regions("wt_15dap")], axis=1, join="inner")
    c = c.loc[sorted(c.index, key=lambda r: int(str(r).rsplit("_", 1)[-1]))]
    d = (c.wt_15dap - c.embryo).to_numpy()
    data.append(("control", CTRL_C, d))
    rows += [dict(set="control", id=g, mCG_embryo=e, mCG_endosperm=n_, endo_minus_embryo=n_ - e)
             for g, e, n_ in zip(c.index, c.embryo, c.wt_15dap)]
    q1, med, q3 = np.percentile(d * 100, [25, 50, 75])
    lo, hi = int((d < YLIM[0]).sum()), int((d > YLIM[1]).sum())
    notes.append(f"control: n = {len(d):,} shuffled regions; median mCG embryo "
                 f"{c.embryo.median() * 100:.1f}%, endosperm {c.wt_15dap.median() * 100:.1f}%; "
                 f"median endosperm - embryo {med:+.1f} points (IQR {q1:+.1f} to {q3:+.1f}), "
                 f"paired two-sided Wilcoxon {pfmt(float(wilcoxon(c.wt_15dap, c.embryo).pvalue))}; "
                 f"{lo + hi} values beyond the {YLIM[0]} to {YLIM[1]} axis ({lo} below, "
                 f"{hi} above).")

    fig = plt.figure(figsize=(WIDTH, HEIGHT))
    ax = fig.add_axes([LEFT / WIDTH, BOTTOM / HEIGHT, (WIDTH - LEFT - RIGHT) / WIDTH,
                       (HEIGHT - TOP - BOTTOM) / HEIGHT])
    parts = ax.violinplot([d for _, _, d in data], positions=range(len(data)), widths=0.75,
                          showmedians=True, showextrema=False)
    for body, (_, colour, _) in zip(parts["bodies"], data):
        body.set_facecolor(colour); body.set_alpha(0.75); body.set_edgecolor("none")
    parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.1)
    for i, (_, _, d) in enumerate(data):
        q1, q3 = np.percentile(d, [25, 75])
        ax.vlines(i, q1, q3, color="black", lw=1.6, zorder=3)
    ax.axhline(0, color="0.5", ls=":", lw=0.6, zorder=0)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels([f"{name}\nn = {len(d):,}" for name, _, d in data])
    ax.tick_params(axis="x", pad=1.5)
    ax.set_xlim(-0.47, len(data) - 0.58)   # five labels at half page width need the room
    ax.set_ylim(*YLIM)
    ax.set_yticks(YTICKS)
    ax.set_yticklabels([minus(t) for t in YTICKS])
    ax.set_ylabel("gene body mCG\nendo − embryo", labelpad=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    FIG.mkdir(parents=True, exist_ok=True)
    save(fig, NAME)

    notes.append("Per-gene CG methylation (methylated / total CG reads over the gene, TSS to "
                 "polyA, summed over the gene's W22 syntelogs), W22 wild-type 15 DAP endosperm "
                 "minus 15 DAP embryo, genes with at least 5 covered CG sites and 10 reads in both "
                 "tissues (Supplemental Dataset 1). Control: regions shuffled from the "
                 "gbM-biparental coordinates with the same length distribution, same filter in "
                 "both tissues. Black segments are medians, black bars interquartile ranges, the "
                 "dotted line is no difference.")
    TAB.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TAB / f"{NAME}_source_data.csv", index=False)
    (FIG / f"{NAME}_legend_numbers.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    print("\n".join(notes))


if __name__ == "__main__":
    main()
