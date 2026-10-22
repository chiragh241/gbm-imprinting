#!/usr/bin/env python3
"""Methylation category of every B73 gene: the rules of Zeng et al. 2023, with 10 covered sites.

The rules are those of Zeng's GeneMethylationType.sh (github.com/dawelab/Natural-methylation-
epialleles-correlate-with-gene-expression-in-maize, "patterns of gene methylation"), in its
order and with its inclusive bounds, applied to the per-gene CDS values of Zeng's own table:

    covered   cCG >= MIN_SITES and cCHG >= MIN_SITES
    UM        mCG <= 0.05 and mCHG <= 0.05
    gbM       mCG >= 0.2  and mCHG <= 0.05
    teM       mCG >= 0.4  and mCHG >= 0.4

Zeng required 40 covered CG and 40 covered CHG sites; this uses 10, which keeps short genes such
as four TE-like MEGs that 40 would drop. At 40 the rules reproduce the table's published epiallele
column for every gene (asserted), so the site minimum is the only departure from the paper.
A gene that is not covered or meets none of the three is left out here and becomes "undefined"
when the gene sets are built (03_build_sets.py).

Input   $GI_DATA/maize_tables/B73.all.csv   (00_references/02_download_maize_tables.sh)
Output  $GI_WORK/02_gene_sets/genes_all_categorized.csv   categorised genes only: every input
                                                          column plus methylation_category
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

OUT = WORK / "02_gene_sets"
OUT.mkdir(parents=True, exist_ok=True)
MIN_SITES = 10


def classify(d, sites):
    """GeneMethylationType.sh: the first rule a gene meets, or "ambiguous"."""
    cov = (d.cCG >= sites) & (d.cCHG >= sites)
    um = cov & (d.mCG <= 0.05) & (d.mCHG <= 0.05)
    gbm = cov & ~um & (d.mCG >= 0.2) & (d.mCHG <= 0.05)
    tem = cov & ~um & ~gbm & (d.mCG >= 0.4) & (d.mCHG >= 0.4)
    return pd.Series(np.select([um, gbm, tem], ["UM", "gbM", "teM"], "ambiguous"), index=d.index)


df = pd.read_csv(DATA / "maize_tables" / "B73.all.csv", low_memory=False)
# Two genes carry two pan-gene rows each, with identical methylation; one row per gene.
meth_cols = ["mCG", "cCG", "mCHG", "cCHG", "epiallele"]
dup = df[df.gene.duplicated(keep=False)]
assert (dup.groupby("gene")[meth_cols].nunique(dropna=False) <= 1).all().all()
df = df.drop_duplicates("gene").reset_index(drop=True)
print(f"genes in table: {len(df):,} ({len(dup) // 2} duplicated rows dropped)")

bad = int((classify(df, 40) != df.epiallele).sum())
assert bad == 0, f"the rules at 40 sites miss the published epiallele column for {bad} genes"
print("rules at 40 sites reproduce the published epiallele column for every gene")

df["methylation_category"] = classify(df, MIN_SITES)
print(f"at {MIN_SITES} sites: " + ", ".join(
    f"{k} {v:,}" for k, v in df.methylation_category.value_counts().items()))
out = df[df.methylation_category != "ambiguous"]
out.to_csv(OUT / "genes_all_categorized.csv", index=False)
print(f"wrote {OUT / 'genes_all_categorized.csv'} ({len(out):,} categorised genes)")
