#!/usr/bin/env python
"""Maternal and paternal CG methylation metagenes in rice 5 DAP endosperm (Supplemental Fig.
S6a-d), the profile form of Tonosaki et al.'s allelic metagenes (their Fig. 3d,e).

AXIS, Figure 2a-d's: 100 bp bins, 2 kb upstream, 2 kb inward from the TSS, 2 kb inward from the
TES, 2 kb downstream; a gene only reaches the inward bins that sit wholly inside it and on its
own side of its midpoint.

ESTIMATOR. The four same-sequence comparisons of 07_gene_matpat.py (Nipponbare allele in each
cross, Kasalath, 93-11). Per comparison and bin, counts summed over replicates and genes,
sum(mC) / sum(C) on each side; a comparison enters a bin when both of its sides have at least
MIN_SITES covered CG cytosines there. The drawn maternal and paternal values are means over the
comparisons that enter, so each parental genome weighs the same. Shading: 95% bootstrap over
genes (1,000 multinomial resamples, the same resample for every comparison).

Input   $GI_WORK/13_rice/meth/<run>.genome{1,2}.bismark.cov.gz, the gene sets and coordinates
        of 06 and 07
Output  $GI_RESULTS/tables/rice_metagene_source.csv  set, bin, mat, mat_lo, mat_hi, pat, pat_lo,
                                                     pat_hi, comparisons, genes
"""
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("m07", HERE / "07_gene_matpat.py")
m07 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m07)

NUP, NIN, NDOWN, BINW = 20, 20, 20, 100
NBIN = NUP + 2 * NIN + NDOWN
MIN_SITES, BOOT = 20, 1000
SETS = ["gbM-PEG", "gbM-biparental", "UM-PEG", "UM-biparental"]
RNG = np.random.default_rng(0)


def make_bins(genes):
    """One row per (gene, bin): chrom, start0, end, bin index 5' to 3'."""
    rows = []
    for g in genes.itertuples():
        s, e = int(g.start) - 1, int(g.end)       # MSU7 1-based inclusive -> half-open
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
            rows.append((g.chrom, a, b, k if g.strand == "+" else NBIN - 1 - k, g.gene))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "bin", "gene"])


def bin_sums(cov_file, bins):
    """(methylated, total, covered CG cytosines) inside each bin, by cumulative sums."""
    d = pd.read_csv(cov_file, sep="\t", header=None, usecols=[0, 1, 4, 5],
                    names=["chrom", "pos", "m", "u"],
                    dtype={"chrom": str, "pos": np.int64, "m": np.int64, "u": np.int64})
    out = np.zeros((len(bins), 3))
    for chrom, b in bins.groupby("chrom"):
        c = d[d.chrom == chrom].sort_values("pos")
        if c.empty:
            continue
        p = c.pos.to_numpy() - 1                   # Bismark coverage is 1-based
        vals = np.column_stack([c.m.to_numpy(float), (c.m + c.u).to_numpy(float), np.ones(len(c))])
        cs = np.vstack([np.zeros(3), np.cumsum(vals, axis=0)])
        lo = np.searchsorted(p, b.start.to_numpy(), "left")
        hi = np.searchsorted(p, b.end.to_numpy(), "left")
        out[b.index.to_numpy()] = cs[hi] - cs[lo]
    return out


def main():
    cls = m07.gene_sets()
    coords = m07.gene_coords()
    coords = coords[coords.gene.map(cls).isin(SETS)]
    bins = make_bins(coords).reset_index(drop=True)
    bins["set"] = bins.gene.map(cls)

    sums = {}           # (comparison, side) -> (nbins, 3)
    cache = {}
    for comp, sides in m07.SIDES.items():
        for side, files in sides.items():
            acc = np.zeros((len(bins), 3))
            for run, genome in files:
                f = m07.METH / f"{run}.{genome}.bismark.cov.gz"
                if f not in cache:
                    cache[f] = bin_sums(f, bins)
                    print(f"  {f.name}")
                acc += cache[f]
            sums[(comp, side)] = acc

    rows = []
    for name in SETS:
        b = bins[bins.set == name]
        genes = sorted(b.gene.unique())
        gi = b.gene.map({g: i for i, g in enumerate(genes)}).to_numpy()
        W = RNG.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=BOOT)
        est = {s: [] for s in ("mat", "pat")}
        boot = {s: [] for s in ("mat", "pat")}
        oks, with_data = [], set()
        for comp in m07.COMPARISONS:
            arr = {}
            for side in ("mat", "pat"):
                a = np.zeros((3, len(genes), NBIN))
                v = sums[(comp, side)][b.index.to_numpy()]
                for k in range(3):
                    np.add.at(a[k], (gi, b.bin.to_numpy()), v[:, k])
                arr[side] = a
            ok = (arr["mat"][2].sum(0) >= MIN_SITES) & (arr["pat"][2].sum(0) >= MIN_SITES)
            oks.append(ok)
            for side in ("mat", "pat"):
                a = arr[side]
                with np.errstate(invalid="ignore", divide="ignore"):
                    est[side].append(np.where(ok, a[0].sum(0) / a[1].sum(0), np.nan))
                    boot[side].append(np.where(ok, (W @ a[0]) / (W @ a[1]), np.nan))
                with_data |= {genes[i] for i in np.flatnonzero(a[2].sum(1) > 0)}
        ncomp = np.sum(oks, axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean = {s: np.nanmean(est[s], axis=0) for s in est}
            ci = {s: np.nanpercentile(np.nanmean(boot[s], axis=0), [2.5, 97.5], axis=0) for s in boot}
        for k in range(NBIN):
            rows.append(dict(set=name, bin=k, mat=mean["mat"][k], mat_lo=ci["mat"][0][k],
                             mat_hi=ci["mat"][1][k], pat=mean["pat"][k], pat_lo=ci["pat"][0][k],
                             pat_hi=ci["pat"][1][k], comparisons=int(ncomp[k]),
                             genes=len(with_data)))
        body = slice(NUP, NUP + 2 * NIN)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            print(f"{name:15s} {len(genes):6,} with bins, {len(with_data):6,} with covered CGs; "
                  f"body mat {np.nanmean(mean['mat'][body]):.3f} pat {np.nanmean(mean['pat'][body]):.3f}")
    m07.TAB.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(m07.TAB / "rice_metagene_source.csv", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
