#!/usr/bin/env python
"""Per-gene mutant-minus-WT methylation deltas, with interquartile ranges and paired tests.

The text quotes each comparison as a median per-gene change in points with its interquartile
range, for gbM-PEGs and gbM-biparental genes, in mdr1 at 15 DAP and in mdr1 and dng102 mature
endosperm. Figure 3a plots the 15-DAP per-gene differences (as wild type minus mdr1) and
Supplemental Fig. S4 the mature-endosperm methylation levels.

Genes need at least MIN_SITES covered CG sites and MIN_READS total reads in BOTH genotypes of
a pair, so the delta is paired per gene.

The shuffled control regions are not genes and are handled by 03_control_region_deltas.sh and
04_control_deltas_summary.py.

Input   $GI_RESULTS/tables/Supplemental_Dataset_1_per_gene.csv   (01_supplemental_dataset1.py)
Output  $GI_RESULTS/tables/mutant_deltas_per_gene.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import TAB

MIN_SITES, MIN_READS = 5, 10

PAIRS = [("wt_15dap", "mdr1_15dap", "mdr1 vs WT, 15 DAP endosperm"),
         ("wt_mature", "mdr1_mature", "mdr1 vs WT, mature endosperm"),
         ("wt_mature", "dng102_mature", "dng102 vs WT, mature endosperm")]
SETS = [("PEG", "gbM", "gbM-PEG"), ("biparental", "gbM", "gbM-biparental"),
        ("PEG", "UM", "UM-PEG"), ("biparental", "UM", "UM-biparental")]


def usable(d, tag):
    return (d[f"nCG_{tag}"] >= MIN_SITES) & (d[f"cov_{tag}"] >= MIN_READS) & d[f"mCG_{tag}"].notna()


def main():
    d = pd.read_csv(TAB / "Supplemental_Dataset_1_per_gene.csv")
    rows = []
    for wt, mut, label in PAIRS:
        ok = d[usable(d, wt) & usable(d, mut)].copy()
        ok["delta"] = 100 * (ok[f"mCG_{mut}"] - ok[f"mCG_{wt}"])
        for imp, meth, name in SETS:
            s = ok[(ok.imp_class == imp) & (ok.meth_class == meth)]
            if len(s) < 5:
                continue
            stat, p = wilcoxon(s[f"mCG_{mut}"], s[f"mCG_{wt}"])
            rows.append(dict(
                comparison=label, gene_set=name, n=len(s),
                wt_median_pct=round(100 * float(s[f"mCG_{wt}"].median()), 1),
                mutant_median_pct=round(100 * float(s[f"mCG_{mut}"].median()), 1),
                delta_points=round(float(s.delta.median()), 1),
                iqr_lo=round(float(s.delta.quantile(0.25)), 1),
                iqr_hi=round(float(s.delta.quantile(0.75)), 1),
                pct_genes_up=round(100 * float((s.delta > 0).mean()), 1),
                wilcoxon_p=float(p)))

    out = pd.DataFrame(rows)
    out.to_csv(TAB / "mutant_deltas_per_gene.csv", index=False)
    pd.set_option("display.width", 220)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
