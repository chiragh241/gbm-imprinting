#!/usr/bin/env python
"""Gene and DMR BED files for the distance from each gene to its nearest demethylated region.

Genes on chromosomes 1 to 10 in B73 v5, each in one category: MEG, gbM-PEG, UM-PEG,
undefined-PEG, or biparental (the single teM-PEG and genes without a call are left out). DMRs are
the 27,402 200 bp B73 v5 regions in which endosperm CG methylation is lower than leaf CG
methylation (data/maize/B73_endo-leaf_hyperDMRs.bed; see data/README.md for its provenance).

Chromosome names: the gene table uses chr1..chr10, the DMR file 1..10, so "chr" is removed.

Input   $GI_WORK/02_gene_sets/master_table.csv   (02_gene_sets/03_build_sets.py)
        data/maize/B73_endo-leaf_hyperDMRs.bed
Output  $GI_WORK/07_dmr_distances/genes.sorted.bed    6 columns, name = gene id
        $GI_WORK/07_dmr_distances/gene_category.tsv   gene id -> category
        $GI_WORK/07_dmr_distances/dmrs.sorted.bed     3 columns, sorted the same way
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import SHIPPED, WORK as GI_WORK

WORK = GI_WORK / "07_dmr_distances"
WORK.mkdir(parents=True, exist_ok=True)

mt = pd.read_csv(GI_WORK / "02_gene_sets" / "master_table.csv")
mt["chrom"] = mt["b73_chr"].str.replace("chr", "", regex=False)
mt = mt[mt["chrom"].isin([str(i) for i in range(1, 11)])].copy()  # nuclear only


def category(r):
    ic, mc = r["imp_class"], r["meth_class"]
    if ic == "MEG":        return "MEG"
    if ic == "biparental": return "biparental"
    if ic == "PEG":
        if mc == "gbM":       return "gbM-PEG"
        if mc == "UM":        return "UM-PEG"
        if mc == "undefined": return "undefined-PEG"
        return None            # teM-PEG (n=1) excluded
    return None


mt["category"] = mt.apply(category, axis=1)
mt = mt.dropna(subset=["category"])
print("category counts (nuclear chr):")
print(mt["category"].value_counts().to_string())

# gene BED (6-col, 0-based half-open; b73_start0 is already 0-based)
genes = mt[["chrom", "b73_start0", "b73_end", "b73_gene", "b73_strand"]].copy()
genes.insert(4, "score", 0)
genes.columns = ["chrom", "start", "end", "name", "score", "strand"]
genes = genes.sort_values(["chrom", "start"])
genes.to_csv(WORK / "genes.sorted.bed", sep="\t", header=False, index=False)

mt[["b73_gene", "category"]].to_csv(WORK / "gene_category.tsv", sep="\t", index=False)

# DMR BED (3-col), sorted the same lexicographic way
dmr = pd.read_csv(SHIPPED / "maize" / "B73_endo-leaf_hyperDMRs.bed", sep="\t", header=None,
                  usecols=[0, 1, 2], names=["chrom", "start", "end"])
dmr["chrom"] = dmr["chrom"].astype(str)
dmr = dmr.sort_values(["chrom", "start"])
dmr.to_csv(WORK / "dmrs.sorted.bed", sep="\t", header=False, index=False)
print(f"\nwrote genes.sorted.bed ({len(genes)}), dmrs.sorted.bed ({len(dmr)})")
