#!/usr/bin/env python
"""Per-gene maternal and paternal CG methylation in rice 5 DAP endosperm (Tonosaki et al. 2024
WGBS, DRA015356), and the gene-set summaries of the rice Results paragraph (Supplemental Fig. S6e).

DESIGN. Two crosses, Nipponbare x Kasalath and Nipponbare x 93-11, each reciprocal, two
replicates per direction. SNPsplit gives each allele-informative read pair to Nipponbare
(genome1) or to the other parent (genome2). Parent of origin is compared on the SAME sequence,
as the arabidopsis analysis does:
    Nipponbare allele   maternal from NIP x alt, paternal from alt x NIP
    alt allele          maternal from alt x NIP, paternal from NIP x alt
which gives four comparisons (Nipponbare in each cross, Kasalath, 93-11). An allele's own
methylation and its mapping behaviour are then the same on both sides of a comparison and cancel.
Counts are summed over the two replicates of a direction.

PER GENE. Weighted CG methylation over the gene body (MSU7 locus, 5' to 3' end), sum(mC) / sum(C)
on each side of a comparison. A comparison counts for a gene when both sides have at least
MIN_SITES covered CG cytosines (10; 5 reported as a sensitivity). The gene's maternal and
paternal values, and pat - mat, are means over its counting comparisons, so each parental genome
weighs the same whatever its read depth.

Input   $GI_WORK/13_rice/meth/<run>.genome{1,2}.bismark.cov.gz   (05)
        $GI_DATA/rice/ref/osa1_r7.all_models.gff3.gz              (01)
        $GI_RESULTS/tables/rice_gene_sets.csv                     (06)
        data/rice/runs.tsv
Output  $GI_RESULTS/tables/rice_matpat_per_gene.csv       counts per comparison, per-gene means
                                                         (MIN_SITES 10)
        $GI_RESULTS/tables/rice_matpat_summary.csv        per class and threshold: n, median,
                                                         IQR, % pat > mat, Wilcoxon mat vs pat
        $GI_RESULTS/tables/rice_matpat_tests.csv          between classes, Mann-Whitney
        $GI_RESULTS/tables/rice_matpat_by_comparison.csv  per class, each comparison on its own
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, SHIPPED, TAB

METH = WORK / "13_rice" / "meth"
GFF = DATA / "rice" / "ref" / "osa1_r7.all_models.gff3.gz"
MIN_SITES = 10
SENSITIVITY = 5
CLASSES = ["gbM-PEG", "gbM-biparental", "UM-PEG", "UM-biparental", "teM-PEG",
           "gbM-MEG", "UM-MEG", "teM-MEG", "teM-biparental"]
TESTS = [("gbM-PEG", "gbM-biparental"), ("UM-PEG", "UM-biparental"),
         ("gbM-PEG", "UM-PEG"), ("gbM-biparental", "UM-biparental")]


def gene_sets():
    gs = pd.read_csv(TAB / "rice_gene_sets.csv")
    return dict(zip(gs.gene, gs["class"]))


def gene_coords():
    rows = []
    with gzip.open(GFF, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) == 9 and f[2] == "gene":
                gid = f[8].split("ID=")[1].split(";")[0]
                rows.append((gid, f[0], int(f[3]), int(f[4]), f[6]))
    return pd.DataFrame(rows, columns=["gene", "chrom", "start", "end", "strand"])


def gene_sums(cov_file, genes):
    """(methylated, total, covered CG cytosines) per gene body, by cumulative sums per chromosome."""
    d = pd.read_csv(cov_file, sep="\t", header=None, usecols=[0, 1, 4, 5],
                    names=["chrom", "pos", "m", "u"],
                    dtype={"chrom": str, "pos": np.int64, "m": np.int64, "u": np.int64})
    out = np.zeros((len(genes), 3))
    for chrom, g in genes.groupby("chrom"):
        c = d[d.chrom == chrom].sort_values("pos")
        if c.empty:
            continue
        p = c.pos.to_numpy()
        vals = np.column_stack([c.m.to_numpy(float), (c.m + c.u).to_numpy(float), np.ones(len(c))])
        cs = np.vstack([np.zeros(3), np.cumsum(vals, axis=0)])
        lo = np.searchsorted(p, g.start.to_numpy(), "left")
        hi = np.searchsorted(p, g.end.to_numpy(), "right")
        out[g.index.to_numpy()] = cs[hi] - cs[lo]
    return out


def summarise(tab, min_sites):
    """Per-gene means over the comparisons where both sides reach min_sites."""
    mats, pats = [], []
    for c in COMPARISONS:
        ok = (tab[f"{c}_mat_n"] >= min_sites) & (tab[f"{c}_pat_n"] >= min_sites)
        mats.append((tab[f"{c}_mat_m"] / tab[f"{c}_mat_t"]).where(ok))
        pats.append((tab[f"{c}_pat_m"] / tab[f"{c}_pat_t"]).where(ok))
    mat = pd.concat(mats, axis=1)
    pat = pd.concat(pats, axis=1)
    return pd.DataFrame({"n_comparisons": mat.notna().sum(axis=1),
                         "mat_mCG": mat.mean(axis=1), "pat_mCG": pat.mean(axis=1),
                         "pat_minus_mat": (pat - mat).mean(axis=1)}, index=tab.index)


def stats(d):
    v = d.pat_minus_mat.to_numpy()
    return dict(n=len(v), median_pat_minus_mat_points=np.median(v) * 100,
                iqr_low=np.percentile(v, 25) * 100, iqr_high=np.percentile(v, 75) * 100,
                pct_pat_above_mat=100 * (v > 0).mean(),
                p_wilcoxon_mat_vs_pat=wilcoxon(d.mat_mCG, d.pat_mCG).pvalue if len(v) > 1 else np.nan)


runs = pd.read_csv(SHIPPED / "rice" / "runs.tsv", sep="\t")
runs = runs[runs.type == "wgbs"]
COMPARISONS = []
SIDES = {}          # comparison -> {"mat": [(run, genome)], "pat": [...]}
for alt in ("KAS", "9311"):
    nip_mother = runs[(runs.alt == alt) & (runs.mother == "NIP")].run.tolist()
    alt_mother = runs[(runs.alt == alt) & (runs.mother == alt)].run.tolist()
    SIDES[f"NIP_in_{alt}"] = {"mat": [(r, "genome1") for r in nip_mother],
                              "pat": [(r, "genome1") for r in alt_mother]}
    SIDES[alt] = {"mat": [(r, "genome2") for r in alt_mother],
                  "pat": [(r, "genome2") for r in nip_mother]}
    COMPARISONS += [f"NIP_in_{alt}", alt]


def main():
    cls = gene_sets()
    genes = gene_coords()
    genes = genes[genes.gene.isin(cls)].reset_index(drop=True)
    genes["class"] = genes.gene.map(cls)
    missing = set(cls) - set(genes.gene)
    print(f"genes in sets {len(cls):,}; with MSU7 coordinates {len(genes):,}; missing {len(missing)}")

    cache = {}
    for comp, sides in SIDES.items():
        for side, files in sides.items():
            acc = np.zeros((len(genes), 3))
            for run, genome in files:
                f = METH / f"{run}.{genome}.bismark.cov.gz"
                if not f.exists():
                    sys.exit(f"missing {f}")
                if f not in cache:
                    cache[f] = gene_sums(f, genes)
                    print(f"  {f.name}: {int(cache[f][:, 2].sum()):,} covered CGs in gene bodies")
                acc += cache[f]
            for k, s in enumerate(("m", "t", "n")):
                genes[f"{comp}_{side}_{s}"] = acc[:, k]

    main_t = summarise(genes, MIN_SITES)
    out = pd.concat([genes, main_t], axis=1)
    TAB.mkdir(parents=True, exist_ok=True)
    out.to_csv(TAB / "rice_matpat_per_gene.csv", index=False)

    rows, tests, bycomp = [], [], []
    for ms in (MIN_SITES, SENSITIVITY):
        t = pd.concat([genes[["gene", "class"]], summarise(genes, ms)], axis=1)
        t = t[t.n_comparisons > 0]
        for c in CLASSES:
            d = t[t["class"] == c]
            if len(d):
                rows.append(dict(min_sites=ms, cls=c, genes_in_set=int((genes["class"] == c).sum()),
                                 **stats(d)))
        for a, b in TESTS:
            x = t[t["class"] == a].pat_minus_mat
            y = t[t["class"] == b].pat_minus_mat
            if len(x) and len(y):
                tests.append(dict(min_sites=ms, a=a, b=b, n_a=len(x), n_b=len(y),
                                  median_a=np.median(x) * 100, median_b=np.median(y) * 100,
                                  p_mannwhitney=mannwhitneyu(x, y).pvalue))
    for comp in COMPARISONS:
        ok = (genes[f"{comp}_mat_n"] >= MIN_SITES) & (genes[f"{comp}_pat_n"] >= MIN_SITES)
        diff = (genes[f"{comp}_pat_m"] / genes[f"{comp}_pat_t"]
                - genes[f"{comp}_mat_m"] / genes[f"{comp}_mat_t"])[ok]
        for c in CLASSES:
            v = diff[genes.loc[ok, "class"] == c]
            if len(v):
                bycomp.append(dict(comparison=comp, cls=c, n=len(v),
                                   median_pat_minus_mat_points=np.median(v) * 100,
                                   pct_pat_above_mat=100 * (v > 0).mean()))
    pd.DataFrame(rows).to_csv(TAB / "rice_matpat_summary.csv", index=False, float_format="%.4g")
    pd.DataFrame(tests).to_csv(TAB / "rice_matpat_tests.csv", index=False, float_format="%.4g")
    pd.DataFrame(bycomp).to_csv(TAB / "rice_matpat_by_comparison.csv", index=False,
                                float_format="%.4g")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(pd.DataFrame(rows).round(3).to_string(index=False))
        print(pd.DataFrame(tests).to_string(index=False))
        print(pd.DataFrame(bycomp).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
