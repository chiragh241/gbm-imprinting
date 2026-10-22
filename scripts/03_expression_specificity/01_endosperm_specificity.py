#!/usr/bin/env python
"""Endosperm specificity of imprinted genes: endosperm TPM over the highest TPM in the other nine
tissues, compared between gbM-PEGs and UM-PEGs.

A gene expressed in endosperm and in no other tissue has a zero denominator and an infinite
ratio. It is kept, ranked above every finite value: the median and the Mann-Whitney U test are
rank statistics, so +infinity enters at the top rank and needs no imputation. A gene with zero
endosperm expression as well is undefined (0/0) and is left out.

Tissues are the ten of Zeng et al. 2023 (TPM from the NAM RNA-seq of Hufford et al. 2021):
embryo and endosperm at 16 DAP, root and shoot, base, middle and tip of the tenth leaf, tassel,
ear and anther. Two genes appear twice in the table under different pan-gene identifiers with
identical values and are kept once.

Input   $GI_DATA/maize_tables/b73.all.6_38.txt      (00_references/02_download_maize_tables.sh)
        $GI_WORK/02_gene_sets/master_table.csv      (02_gene_sets/03_build_sets.py)
Output  $GI_RESULTS/tables/endosperm_specificity_per_gene.csv   imprinted genes with a category
        $GI_RESULTS/tables/endosperm_specificity_summary.csv    n and median per set, and the test
"""
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, TAB

TAB.mkdir(parents=True, exist_ok=True)
TISSUES = ["embryo", "endosperm", "root", "shoot", "anther",
           "base", "middle", "tip", "ear", "tassel"]
OTHERS = [t for t in TISSUES if t != "endosperm"]
SETS = [("PEG", "gbM", "gbM-PEG"), ("PEG", "UM", "UM-PEG"),
        ("MEG", "UM", "UM-MEG"), ("MEG", "teM", "teM-MEG")]

# The pan-gene class strings contain spaces in a whitespace-separated table.
txt = (DATA / "maize_tables" / "b73.all.6_38.txt").read_text()
for a, b in [("Near-Core Gene", "Near-Core_Gene"), ("Core Gene", "Core_Gene"),
             ("Dispensable Gene", "Dispensable_Gene"), ("Private Gene", "Private_Gene")]:
    txt = txt.replace(a, b)
e = pd.read_csv(io.StringIO(txt), sep=r"\s+", low_memory=False).drop_duplicates("gene")
for c in TISSUES:
    e[c] = pd.to_numeric(e[c], errors="coerce")
m = pd.read_csv(WORK / "02_gene_sets" / "master_table.csv")[
    ["b73_gene", "imp_class", "meth_class"]].rename(columns={"b73_gene": "gene"})
e = e.merge(m, on="gene", how="left")

e["max_nonendo"] = e[OTHERS].max(axis=1)
e["endo_specificity"] = np.where((e.max_nonendo == 0) & (e.endosperm > 0), np.inf,
                                 e.endosperm / e.max_nonendo.replace(0, np.nan))

rows, vals = [], {}
for imp, meth, name in SETS:
    s = e[(e.imp_class == imp) & (e.meth_class == meth)].endo_specificity.dropna().to_numpy()
    vals[name] = s
    rows.append(dict(set=name, n=len(s), median=float(np.median(s)),
                     n_infinite=int(np.isinf(s).sum())))
g, u = vals["gbM-PEG"], vals["UM-PEG"]
stat, p = mannwhitneyu(g, u, alternative="two-sided")
cliff = 2 * stat / (len(g) * len(u)) - 1
summary = pd.DataFrame(rows)
summary["test"] = ""
summary.loc[summary.set == "UM-PEG", "test"] = (
    f"gbM-PEG vs UM-PEG two-sided Mann-Whitney U p = {p:.3g}, Cliff's delta {cliff:+.3f}")
summary.to_csv(TAB / "endosperm_specificity_summary.csv", index=False)

cols = ["gene", "imp_class", "meth_class", "endosperm", "max_nonendo", "endo_specificity"]
e[e.imp_class.isin(["PEG", "MEG"]) & e.meth_class.isin(["gbM", "UM", "teM"])][cols].to_csv(
    TAB / "endosperm_specificity_per_gene.csv", index=False)
print(summary.to_string(index=False))
