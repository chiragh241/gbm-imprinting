#!/usr/bin/env python
"""Former Figure 5b, not in the paper: maternal and paternal CG methylation metagenes in
arabidopsis endosperm. Supplemental Figure S7a-d (15_fig_arabidopsis_metagene.py) uses its
gene sets, bins and estimator.

Parallel to Figure 2a-d.

DATA. GSE280598 (Kirkbride/Gehring), reciprocal Col-0 x C24 F1 WT endosperm, 7 DAP, three
replicates per direction, the same files panel c uses. Both tracks are the COL allele: maternal
from Col x C24 (Col is the mother), paternal from C24 x Col (Col is the father). The sequence
is the same on both sides, so parent of origin is not confounded with ecotype.
Per-cytosine CG beds, TAIR10: chrom start0 end1 unmethylated methylated ratio.

GENE SETS, the paper's canonical ones, as in 5a and 5c-e:
  PEGs     consensus PEG-36 (called in at least two of Pignatta, Gehring, Wolff)
  class    Bewick 2016 via 06_master_table.csv, gbM or UM
  bip      every gbM or UM gene on none of the three screens' PEG or MEG lists
           (listed_imprinted_genes.csv), the same control rule as 09_fig5cde.py

AXIS, Figure 2a-d's: 100 bp bins, 2 kb upstream, 2 kb inward from the TSS, 2 kb inward from
the TES, 2 kb downstream, a gene only reaching the inward bins that sit wholly inside it and
on its own side of the midpoint (04_emseq_methylation/07_metagene_bins.sh). Arabidopsis
genes are short, so
fewer genes reach the deep inward bins and the flanks often run into neighbouring genes.

ESTIMATOR, as 04_fig2_supp_matpat_metagene.py: counts summed over replicates and genes, each bin
sum(mC) / sum(C) per allele, a bin drawn only with at least MIN_SITES covered CGs per allele.
Shading is a 95% bootstrap interval over genes (1,000 resamples, multinomial weights).

Input   $GI_DATA/arabidopsis/gse280598/{mat,pat}_r{1,2,3}.CG.bed.gz    (11_arabidopsis/01)
        $GI_WORK/11_arabidopsis/06_master_table.csv                    (11_arabidopsis)
        $GI_RESULTS/tables/consensus_pegs.csv, listed_imprinted_genes.csv  (11_arabidopsis/03)
Output  $GI_RESULTS/figures/fig5p_b.{png,pdf}
        $GI_RESULTS/tables/fig5b_source_data.csv
"""
from pathlib import Path
import gzip
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D

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
RES = WORK / "11_arabidopsis"
BEDS = DATA / "arabidopsis" / "gse280598"
GENES = RES / "gene_targets.csv"
FILES = {"mat": ["mat_r1.CG.bed.gz", "mat_r2.CG.bed.gz", "mat_r3.CG.bed.gz"],
         "pat": ["pat_r1.CG.bed.gz", "pat_r2.CG.bed.gz", "pat_r3.CG.bed.gz"]}

NUP, NIN, NDOWN, BINW = 20, 20, 20, 100
NBIN = NUP + 2 * NIN + NDOWN
MID = NUP + NIN
MIN_SITES, BOOT = 20, 1000
MAT, PAT = "#c0392b", "#2c7fb8"
RNG = np.random.default_rng(0)
SETS = ["gbM-PEG", "gbM-bip", "UM-PEG", "UM-bip"]


def gene_sets():
    master = pd.read_csv(RES / "06_master_table.csv")
    mc = dict(zip(master.AGI, master.meth_class))
    cons = set(pd.read_csv(TAB / "consensus_pegs.csv").gene)
    imprinted = set(pd.read_csv(TAB / "listed_imprinted_genes.csv").gene)
    return {"gbM-PEG": {g for g in cons if mc.get(g) == "gbM"},
            "UM-PEG": {g for g in cons if mc.get(g) == "UM"},
            "gbM-bip": set(master[master.meth_class == "gbM"].AGI) - imprinted,
            "UM-bip": set(master[master.meth_class == "UM"].AGI) - imprinted}


def make_bins(genes):
    """One row per (gene, bin): chrom, start0, end, bin index 5' to 3'."""
    rows = []
    for g in genes.itertuples():
        s, e = int(g.start) - 1, int(g.end)       # TAIR10 1-based inclusive -> half-open
        L = e - s
        if L < 2 * BINW:
            continue
        mid = s + L // 2
        out = []
        for i in range(NUP):
            a = s - (NUP - i) * BINW
            if a >= 0:
                out.append((a, a + BINW, i))
        for i in range(NIN):
            a = s + i * BINW
            if a + BINW > mid:
                break
            out.append((a, a + BINW, NUP + i))
        for i in range(NIN):
            a = e - (NIN - i) * BINW
            if a < mid:
                continue
            out.append((a, a + BINW, NUP + NIN + i))
        for i in range(NDOWN):
            a = e + i * BINW
            out.append((a, a + BINW, NUP + 2 * NIN + i))
        for a, b, k in out:
            rows.append((f"Chr{g.chr}", a, b, k if g.strand == "+" else NBIN - 1 - k, g.AGI))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "bin", "gene"])


def load_bed(fn):
    with gzip.open(BEDS / fn, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 3, 4],
                        names=["chrom", "s0", "um", "m"],
                        dtype={"chrom": str, "s0": np.int64, "um": np.int32, "m": np.int32})
    d = d[(d.um + d.m) > 0]
    return d


def bin_sums(bed, bins):
    """(methylated, total, covered sites) inside each bin, via cumulative sums per chromosome."""
    out = np.zeros((len(bins), 3))
    for chrom, b in bins.groupby("chrom"):
        c = bed[bed.chrom == chrom].sort_values("s0")
        if c.empty:
            continue
        p = c.s0.to_numpy()
        vals = np.column_stack([c.m.to_numpy(float), (c.m + c.um).to_numpy(float),
                                np.ones(len(c))])
        cs = np.vstack([np.zeros(3), np.cumsum(vals, axis=0)])
        lo = np.searchsorted(p, b.start.to_numpy(), "left")
        hi = np.searchsorted(p, b.end.to_numpy(), "left")
        out[b.index.to_numpy()] = cs[hi] - cs[lo]
    return out


def main():
    sets = gene_sets()
    coords = pd.read_csv(GENES, dtype={"chr": str})
    wanted = set().union(*sets.values())
    bins = make_bins(coords[coords.AGI.isin(wanted)]).reset_index(drop=True)
    per_allele = {}
    for allele, files in FILES.items():
        tot = np.zeros((len(bins), 3))
        for fn in files:
            bed = load_bed(fn)
            tot += bin_sums(bed, bins)
            print(f"  {fn}: {len(bed):,} covered CGs")
        per_allele[allele] = tot

    rows, prof = [], {}
    for name in SETS:
        b = bins[bins.gene.isin(sets[name])]
        genes = sorted(b.gene.unique())
        gi = {g: i for i, g in enumerate(genes)}
        W = RNG.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=BOOT)
        est, n_linked = {}, 0
        for allele in FILES:
            a = np.zeros((3, len(genes), NBIN))
            v = per_allele[allele][b.index.to_numpy()]
            for k in range(3):
                np.add.at(a[k], (b.gene.map(gi).to_numpy(), b.bin.to_numpy()), v[:, k])
            ok = a[2].sum(0) >= MIN_SITES
            with np.errstate(invalid="ignore", divide="ignore"):
                pt = a[0].sum(0) / a[1].sum(0)
                bt = (W @ a[0]) / (W @ a[1])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                lo, hi = np.nanpercentile(bt, [2.5, 97.5], axis=0)
            est[allele] = tuple(np.where(ok, x, np.nan) for x in (pt, lo, hi))
            n_linked = max(n_linked, int((a[2].sum(1) > 0).sum()))
        prof[name] = (est, n_linked)
        for k in range(NBIN):
            rows.append({"set": name, "bin": k, "mat": est["mat"][0][k],
                         "mat_lo": est["mat"][1][k], "mat_hi": est["mat"][2][k],
                         "pat": est["pat"][0][k], "pat_lo": est["pat"][1][k],
                         "pat_hi": est["pat"][2][k]})
        body = slice(NUP, NUP + 2 * NIN)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            print(f"{name:8s} {len(sets[name]):6,} in set, {n_linked:6,} with covered CGs  "
                  f"body mat {np.nanmean(est['mat'][0][body]):.3f} "
                  f"pat {np.nanmean(est['pat'][0][body]):.3f}  bins drawn "
                  f"{int(np.isfinite(est['mat'][0]).sum())}/{NBIN}")
    TAB.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TAB / "fig5b_source_data.csv", index=False)
    draw(prof)


def draw(prof):
    # A second full-width row under a, c, d and e: four panels across 6.30 in, parallel to the
    # four stacked in Figure 2a-d.
    W, L_LAB, R_PAD, GAP = 6.30, 0.40, 0.04, 0.12
    PAD, TITLE, AX, KEY = 0.04, 0.16, 0.95, 0.20
    H = PAD + TITLE + AX + KEY
    pw = (W - L_LAB - R_PAD - 3 * GAP) / 4
    fig = plt.figure(figsize=(W, H))
    for k, (name, (est, n)) in enumerate(prof.items()):
        x0 = L_LAB + k * (pw + GAP)
        ax = fig.add_axes([x0 / W, (H - PAD - TITLE - AX) / H, pw / W, AX / H])
        x = np.arange(NBIN)
        for allele, colour in (("mat", MAT), ("pat", PAT)):
            m, lo, hi = est[allele]
            for a, b in ((0, MID), (MID, NBIN)):
                ax.fill_between(x[a:b], lo[a:b], hi[a:b], color=colour, alpha=0.18, lw=0)
                ax.plot(x[a:b], m[a:b], color=colour, lw=0.9)
        ax.axvline(NUP - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(NBIN - NDOWN - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(MID - 0.5, color="0.75", ls=":", lw=0.6)
        ax.set_xticks([])
        ax.set_xlim(0, NBIN - 1)
        # 0 to 1 as in Figure 2a-d, so the two figures read on one scale
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"])
        if k == 0:
            ax.set_ylabel("mCG", labelpad=1)
        else:
            ax.set_yticklabels([])
        ax.set_title(f"{name}, n = {n:,}", pad=2, loc="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.legend(handles=[Line2D([], [], color=MAT, lw=1.2, label="mat"),
                        Line2D([], [], color=PAT, lw=1.2, label="pat")],
               ncol=2, frameon=False, loc="lower left", bbox_to_anchor=(L_LAB / W, 0.0),
               handletextpad=0.35, columnspacing=0.8, handlelength=1.2, borderaxespad=0,
               borderpad=0.1)
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"fig5p_b.{ext}", dpi=600)
    print(f"wrote fig5p_b.png / .pdf at {W} x {H:.2f} in")
    plt.close(fig)


if __name__ == "__main__":
    main()
