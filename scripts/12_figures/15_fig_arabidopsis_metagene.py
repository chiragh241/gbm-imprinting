#!/usr/bin/env python
"""Supplemental Figure S7a-d: CG methylation metagenes in arabidopsis endosperm and leaf, the
arabidopsis counterpart of Figure 2a-d. S7e and f are 09_fig5cde.py's fig5p_c and fig5p_d.

DATA. GSE280598 (Kirkbride/Gehring): wild-type Col-0 endosperm at 7 DAP and young leaf, three
replicates each, EM-seq, per-cytosine CG beds on TAIR10 (11_arabidopsis/01_download.sh). The same
series supplies S7e. No embryo methylome exists in this series; the arabidopsis embryo methylomes
(Hsieh 2009, Ibarra 2012) are older, shallower and on TAIR8, so the tissue contrast here is
endosperm against leaf.

GENE SETS, AXIS AND ESTIMATOR as 08_fig5b.py (imported from there): consensus PEG-36 split by the
Bewick 2016 class; biparental = every gbM or UM gene on no imprinted list; Figure 2a-d's axis of
100 bp bins, 2 kb upstream, 2 kb inward from the TSS and from the TES, 2 kb downstream, a gene
reaching only the inward bins wholly inside it and on its own side of the midpoint; each bin
sum(mC) / sum(C) over genes and replicates, drawn with at least MIN_SITES covered CGs.
Arabidopsis genes are short and packed, so fewer genes reach the deep inward bins and the 2 kb
flanks run into neighbouring genes; the axis is kept identical to Figure 2a-d so the two can be
read side by side.

DRAWING as Figure 2a-d (03_fig2.py): one line per tissue, split at the mid-gene break, dashed TSS
and polyA, y 0 to 1, each panel titled with its set and the most genes behind any bin, the key
under the last panel. Same 1.62 x 3.96 in as fig2p_ad, so it pastes at the same scale.

Input   $GI_DATA/arabidopsis/gse280598_tissues/GSM86018{82,83,84,94,95,96}_*_CG.bed.gz
                                                            (11_arabidopsis/01_download.sh)
        $GI_WORK/11_arabidopsis/gene_targets.csv, 06_master_table.csv   (11_arabidopsis)
        $GI_RESULTS/tables/consensus_pegs.csv, listed_imprinted_genes.csv   (11_arabidopsis/03)
Output  $GI_RESULTS/figures/figS_ara_ad.{png,pdf}, figS_ara_ad_legend_numbers.txt
        $GI_RESULTS/tables/figS_ara_ad_source_data.csv
"""
import gzip
import importlib.util
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("_fig5b", HERE / "08_fig5b.py")
m08 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m08)           # house style, gene sets, bins; its main() is not run

FIG, TAB = m08.FIG, m08.TAB
BEDS = m08.DATA / "arabidopsis" / "gse280598_tissues"
TISSUES = [("endosperm", "endo", "#c0392b",
            ["GSM8601882_wt_1_Col_spiked_CG.bed.gz", "GSM8601883_wt_2_Col_spiked_CG.bed.gz",
             "GSM8601884_wt_3_Col_spiked_CG.bed.gz"]),
           # blue, not Figure 2a-d's leaf green, as the second line of 2a-d, 4 and 5b
           ("leaf", "leaf", "#2c7fb8",
            ["GSM8601894_wt_1_leaf_Col_spiked_CG.bed.gz",
             "GSM8601895_wt_2_leaf_Col_spiked_CG.bed.gz",
             "GSM8601896_wt_3_leaf_Col_spiked_CG.bed.gz"])]
NUP, NIN, NDOWN, NBIN = m08.NUP, m08.NIN, m08.NDOWN, m08.NBIN
MID = NUP + NIN
BREAK = MID - 0.5
MIN_SITES = m08.MIN_SITES
SETS = m08.SETS                         # gbM-PEG, gbM-bip, UM-PEG, UM-bip: Figure 2a-d's order

# Figure 2a-d's part, 03_fig2.py
AD_W, TOTAL_H, PAD = 1.62, 3.96, 0.03
L_LAB, MG_TITLE, MG_GAP, MG_LEG = 0.36, 0.14, 0.05, 0.34
MG_AX = (TOTAL_H - PAD - 4 * MG_TITLE - 3 * MG_GAP - MG_LEG) / 4


def load_bed(fn):
    with gzip.open(BEDS / fn, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 3, 4],
                        names=["chrom", "s0", "um", "m"],
                        dtype={"chrom": str, "s0": np.int64, "um": np.int32, "m": np.int32})
    return d[(d.um + d.m) > 0]


def save(fig, name):
    """Write both formats after checking no text runs off the sheet or into the key."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    arts = list(fig.texts)
    for ax in fig.axes:
        arts += (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                 + [ax.xaxis.label, ax.yaxis.label, ax.title, ax._left_title])
    arts = [a for a in arts if a.get_text()]
    for a in arts:
        bb = a.get_window_extent(r)
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
            raise SystemExit(f"!!! {name}: {a.get_text()[:40]!r} runs off the {w:.2f} x "
                             f"{h:.2f} in sheet")
    for lg in fig.legends:
        lb = lg.get_window_extent(r)
        for a in arts:
            if a not in lg.get_texts() and lb.overlaps(a.get_window_extent(r)):
                raise SystemExit(f"!!! {name}: the key runs into {a.get_text()[:30]!r}")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w:.2f} x {h:.2f} in")
    plt.close(fig)


def main():
    sets = m08.gene_sets()
    coords = pd.read_csv(m08.GENES, dtype={"chr": str})
    wanted = set().union(*sets.values())
    bins = m08.make_bins(coords[coords.AGI.isin(wanted)]).reset_index(drop=True)
    per_tissue = {}
    for tissue, _lab, _c, files in TISSUES:
        tot = np.zeros((len(bins), 3))
        for fn in files:
            bed = load_bed(fn)
            tot += m08.bin_sums(bed, bins)
            print(f"  {fn}: {len(bed):,} covered CGs, mCG {bed.m.sum() / (bed.m + bed.um).sum():.3f}")
        per_tissue[tissue] = tot

    rows, prof, notes = [], {}, []
    for name in SETS:
        b = bins[bins.gene.isin(sets[name])]
        genes = sorted(b.gene.unique())
        gi = {g: i for i, g in enumerate(genes)}
        prof[name] = {}
        for tissue, _lab, _c, _f in TISSUES:
            a = np.zeros((3, len(genes), NBIN))
            v = per_tissue[tissue][b.index.to_numpy()]
            for k in range(3):
                np.add.at(a[k], (b.gene.map(gi).to_numpy(), b.bin.to_numpy()), v[:, k])
            with np.errstate(invalid="ignore", divide="ignore"):
                pooled = a[0].sum(0) / a[1].sum(0)
            ok = a[2].sum(0) >= MIN_SITES
            pooled = np.where(ok, pooled, np.nan)
            n_genes = (a[2] > 0).sum(0)
            prof[name][tissue] = (pooled, n_genes)
            for k in range(NBIN):
                rows.append(dict(set=name, tissue=tissue, bin=k, pooled=pooled[k],
                                 n_genes=int(n_genes[k]), cg_sites=int(a[2][:, k].sum())))
        body = slice(NUP, NUP + 2 * NIN)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            e_, l_ = (np.nanmean(prof[name][t][0][body]) for t in ("endosperm", "leaf"))
        nmax = max(int(prof[name][t][1].max()) for t in ("endosperm", "leaf"))
        notes.append(f"{name}: {len(sets[name]):,} genes in the set, up to {nmax:,} behind a bin; "
                     f"mean of the gene-body bins, endosperm {e_:.3f}, leaf {l_:.3f} "
                     f"(endosperm - leaf {e_ - l_:+.3f})")

    fig = plt.figure(figsize=(AD_W, TOTAL_H))

    def box(x0, w, top_in, h_in):
        return fig.add_axes([x0 / AD_W, (TOTAL_H - top_in - h_in) / TOTAL_H, w / AD_W,
                             h_in / TOTAL_H])

    y = PAD
    for k, name in enumerate(SETS):
        y += MG_TITLE
        ax = box(L_LAB, AD_W - L_LAB - 0.03, y, MG_AX)
        nmax = 0
        for tissue, _lab, colour, _f in TISSUES:
            pooled, n_genes = prof[name][tissue]
            nmax = max(nmax, int(n_genes.max()))
            x = np.arange(NBIN)
            for a_, b_ in ((0, MID), (MID, NBIN)):
                ax.plot(x[a_:b_], pooled[a_:b_], color=colour, lw=0.8)
        ax.axvline(NUP - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(NBIN - NDOWN - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(BREAK, color="0.75", ls=":", lw=0.6)
        ax.set_xticks([])
        ax.set_xlim(0, NBIN - 1)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"])
        ax.set_ylabel("mCG", labelpad=1)
        ax.set_title(f"{name}, n = {nmax:,}", pad=1.5, loc="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        y += MG_AX + (MG_GAP if k < 3 else 0)
    fig.legend(handles=[Line2D([], [], color=c, lw=1.2, label=lab) for _, lab, c, _f in TISSUES],
               ncol=2, frameon=False, loc="upper left",
               bbox_to_anchor=(0.02, (TOTAL_H - y - 0.07) / TOTAL_H),
               handletextpad=0.35, labelspacing=0.2, columnspacing=0.8, handlelength=1.2,
               borderaxespad=0, borderpad=0)
    FIG.mkdir(parents=True, exist_ok=True)
    save(fig, "figS_ara_ad")

    notes.append("Key: endo = Col-0 endosperm, 7 DAP; leaf = Col-0 young leaf (GSE280598, "
                 "EM-seq, three replicates each, summed). Bins as Figure 2a-d: 100 bp, 2 kb "
                 "upstream, 2 kb in from the TSS, 2 kb in from the polyA site, 2 kb downstream; "
                 f"a bin is drawn with at least {MIN_SITES} covered CGs. n is the most genes "
                 "behind any bin. PEGs: consensus of three screens; biparental: every gbM or UM "
                 "gene on no imprinted list.")
    TAB.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TAB / "figS_ara_ad_source_data.csv", index=False)
    (FIG / "figS_ara_ad_legend_numbers.txt").write_text("\n".join(notes) + "\n",
                                                         encoding="utf-8")
    print("\n".join(notes))


if __name__ == "__main__":
    main()
