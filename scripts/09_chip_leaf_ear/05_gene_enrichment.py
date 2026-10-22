#!/usr/bin/env python
"""H3K27me3 enrichment over gene bodies in B73 leaf and developing ear (Supplemental Fig. S1).

Data: Ricci et al. 2019 Nature Plants, GSE120304, AGPv4 mapped reads.
  leaf  H3K27me3 rep1+rep2 vs input rep1+rep2   (seedling second leaf)
  ear   H3K27me3 rep1+rep2 vs input rep1        (developing ear)

Enrichment is log2(ChIP RPKM / input RPKM) per gene, replicates pooled within tissue, with a
pseudocount of one read. Input normalisation matters here: gbM and UM genes differ in sequence
composition and mappability, so a raw ChIP density difference between them could be an artefact
of the very variable being compared. Genes need at least 10 input reads so the ratio has a
denominator (27 of 34 UM-PEGs pass in leaf; the gbM-PEG and UM-PEG results hold at 5, 20 and 40).
A gene is called enriched above twofold when log2 enrichment exceeds 1. UM-PEGs and gbM-PEGs are
compared by two-sided Mann-Whitney U test.

Input   $GI_WORK/09_chip_leaf_ear/counts/<lib>.counts.tsv, <lib>.total.txt   (04_count.sh)
        $GI_WORK/09_chip_leaf_ear/fig1e_gene_intervals.bed                  (03_gene_intervals.py)
Output  $GI_WORK/09_chip_leaf_ear/figS1_per_gene.csv    per gene and tissue: reads, enrichment
        $GI_RESULTS/tables/figS1_gene_enrichment.csv     per set and tissue: n, median, IQR,
                                                        percent enriched, and the test
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

D = WORK / "09_chip_leaf_ear"
CNT = D / "counts"
TAB.mkdir(parents=True, exist_ok=True)
MIN_INPUT = 10
PSEUDO = 1.0

TISSUES = {
    "leaf": (["leaf_H3K27me3_rep1", "leaf_H3K27me3_rep2"],
             ["leaf_input_rep1", "leaf_input_rep2"]),
    "ear":  (["ear_H3K27me3_rep1", "ear_H3K27me3_rep2"],
             ["ear_input_rep1"]),
}
SETS = ["gbM-PEG", "UM-PEG", "gbM-biparental", "UM-biparental"]

intervals = pd.read_csv(D / "fig1e_gene_intervals.bed", sep="\t", header=None,
                        names=["chrom", "start", "end", "gene", "set"])
intervals["kb"] = (intervals.end - intervals.start) / 1000.0


def load(lib):
    c = pd.read_csv(CNT / f"{lib}.counts.tsv", sep="\t", header=None,
                    names=["gene", "set", "n"])
    total = int(open(CNT / f"{lib}.total.txt").read().strip())
    return c.set_index("gene").n, total


def pooled(libs):
    """summed reads per gene and summed library size"""
    s, tot = None, 0
    for lib in libs:
        n, t = load(lib)
        s = n if s is None else s.add(n, fill_value=0)
        tot += t
    return s, tot


rows, per_gene = [], []
for key, (chip_libs, input_libs) in TISSUES.items():
    chip, chip_tot = pooled(chip_libs)
    inp, inp_tot = pooled(input_libs)

    d = intervals.copy()
    d["chip"] = d.gene.map(chip).fillna(0.0)
    d["input"] = d.gene.map(inp).fillna(0.0)
    d = d[d.input >= MIN_INPUT]
    d["chip_rpkm"] = (d.chip + PSEUDO) / d.kb / (chip_tot / 1e6)
    d["input_rpkm"] = (d.input + PSEUDO) / d.kb / (inp_tot / 1e6)
    d["log2_enrich"] = np.log2(d.chip_rpkm / d.input_rpkm)
    d["tissue"] = key
    per_gene.append(d[["gene", "set", "tissue", "chip", "input", "log2_enrich"]])

    series = []
    for s in SETS:
        v = d[d.set == s].log2_enrich.to_numpy()
        series.append(v)
        rows.append(dict(tissue=key, set=s, n=len(v),
                         median_log2_enrichment=float(np.median(v)),
                         q1=float(np.percentile(v, 25)), q3=float(np.percentile(v, 75)),
                         pct_enriched_above_1=float(100 * np.mean(v > 1))))
    # the claim is UM-PEG > gbM-PEG
    rows[-len(SETS)]["UMPEG_vs_gbMPEG_p"] = mannwhitneyu(series[1], series[0]).pvalue

src = pd.DataFrame(rows)
src.to_csv(TAB / "figS1_gene_enrichment.csv", index=False)
pd.concat(per_gene).to_csv(D / "figS1_per_gene.csv", index=False)
print(src.to_string(index=False))
