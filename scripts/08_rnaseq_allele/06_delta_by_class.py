#!/usr/bin/env python
"""Per-gene change in maternal expression ratio in the mdr1 dng103 crosses, by gene class.

Delta is per gene: the maternal ratio with the double mutant as the mother minus that with the
wild-type mother (bM - BM), in percentage points; the paternal-side control is the same with
the mutant as the father (Mb - MB). Reported per class as the median delta with its
interquartile range and a two-sided Wilcoxon signed-rank test on the paired per-gene ratios
(scipy default zero_method), not corrected across classes. n is the genes scored in both
crosses of a pair (at least 20 allele-informative fragments in each; 05_maternal_fraction.py).

Input   $GI_WORK/08_rnaseq_allele/per_gene_maternal_fraction.tsv   (05_maternal_fraction.py)
        $GI_WORK/02_gene_sets/master_table.csv                       imprinting and methylation classes
Output  $GI_RESULTS/tables/xu_delta_by_methylation_class.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

SETS = [("PEG", "gbM", "gbM-PEG"), ("PEG", "UM", "UM-PEG"),
        ("MEG", "UM", "UM-MEG"), ("MEG", "teM", "teM-MEG"),
        ("biparental", "gbM", "gbM-biparental"), ("biparental", "UM", "UM-biparental")]

g = pd.read_csv(WORK / "08_rnaseq_allele" / "per_gene_maternal_fraction.tsv", sep="\t")
master = pd.read_csv(WORK / "02_gene_sets" / "master_table.csv")[
    ["b73_gene", "imp_class", "meth_class"]]
g = g.merge(master, left_on="gene", right_on="b73_gene", how="left")


def stats(d, a, b):
    """Paired per-gene delta in percentage points, with the signed-rank test."""
    d = d[d[a].notna() & d[b].notna()]
    if len(d) < 2:
        return None
    x, y = d[a].to_numpy(), d[b].to_numpy()
    delta = (y - x) * 100.0
    try:
        p = float(wilcoxon(x, y).pvalue)
    except ValueError:                      # all differences zero
        p = np.nan
    q1, q3 = np.percentile(delta, [25, 75])
    return dict(n=len(delta), median_wt=float(np.median(x)), median_mut=float(np.median(y)),
                median_delta_pts=float(np.median(delta)), iqr_lo_pts=float(q1),
                iqr_hi_pts=float(q3), p=p, n_up=int((delta > 0).sum()),
                n_down=int((delta < 0).sum()))


rows = []
for imp, meth, label in SETS:
    d = g[(g.imp_class == imp) & (g.meth_class == meth)]
    for side, a, b in (("maternal (bM - BM)", "BM_frac", "bM_frac"),
                       ("paternal control (Mb - MB)", "MB_frac", "Mb_frac")):
        s = stats(d, a, b)
        if s:
            rows.append(dict(set=label, side=side, n_total=len(d), **s))
out = pd.DataFrame(rows)
TAB.mkdir(parents=True, exist_ok=True)
out.to_csv(TAB / "xu_delta_by_methylation_class.tsv", sep="\t", index=False)
print(out.to_string(index=False))
