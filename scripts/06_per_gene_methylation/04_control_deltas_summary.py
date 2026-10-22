#!/usr/bin/env python
"""Per-region shuffled-control deltas for every mutant comparison (Fig. 3a, Supplemental Fig. S4).

Regions are paired by region id between the two genotypes, so the delta is paired the same way
the per-gene deltas of 02_mutant_deltas.py are.

Input   $GI_WORK/06_per_gene_methylation/ctrl_<genotype>.tsv   (03_control_region_deltas.sh)
Output  $GI_RESULTS/tables/control_deltas.csv
"""
from pathlib import Path

import pandas as pd
from scipy.stats import wilcoxon

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

V = WORK / "06_per_gene_methylation"
PAIRS = [("wt_15dap", "mdr1_15dap", "mdr1 vs WT, 15 DAP endosperm"),
         ("wt_mature", "mdr1_mature", "mdr1 vs WT, mature endosperm"),
         ("wt_mature", "dng102_mature", "dng102 vs WT, mature endosperm")]

rows = []
for wt, mut, label in PAIRS:
    fw, fm = V / f"ctrl_{wt}.tsv", V / f"ctrl_{mut}.tsv"
    if not (fw.exists() and fm.exists()):
        print(f"### {label}: missing aggregation, skip")
        continue
    a = pd.read_csv(fw, sep="\t", header=None, names=["region", "wt"])
    b = pd.read_csv(fm, sep="\t", header=None, names=["region", "mut"])
    m = a.merge(b, on="region")
    d = 100 * (m["mut"] - m["wt"])
    stat, p = wilcoxon(m["mut"], m["wt"])
    rows.append(dict(comparison=label, gene_set="shuffled control", n=len(m),
                     wt_median_pct=round(100 * float(m["wt"].median()), 1),
                     mutant_median_pct=round(100 * float(m["mut"].median()), 1),
                     delta_points=round(float(d.median()), 1),
                     iqr_lo=round(float(d.quantile(0.25)), 1),
                     iqr_hi=round(float(d.quantile(0.75)), 1),
                     pct_regions_up=round(100 * float((d > 0).mean()), 1),
                     wilcoxon_p=float(p)))

out = pd.DataFrame(rows)
out.to_csv(TAB / "control_deltas.csv", index=False)
pd.set_option("display.width", 200)
print(out.to_string(index=False))
