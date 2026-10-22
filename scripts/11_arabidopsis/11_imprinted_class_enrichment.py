#!/usr/bin/env python
"""Methylation classes of the Arabidopsis consensus PEGs and MEGs, and gbM enrichment in PEGs.

The consensus PEGs and MEGs are classified by the Bewick et al. 2016 Col-0 leaf classes of the
master table (UM, gbM, teM; genes without a class are not counted). gbM genes among classified
PEGs are compared with gbM genes among classified MEGs: fold difference of the two proportions,
odds ratio and two-sided Fisher exact test.

Input   $GI_RESULTS/tables/consensus_{pegs,megs}.csv   (03_consensus_sets.py)
        $GI_WORK/11_arabidopsis/06_master_table.csv    (05_master_table.py)
Output  $GI_RESULTS/tables/arabidopsis_imprinted_classes.csv   set, classified, UM, gbM, teM
        $GI_RESULTS/tables/arabidopsis_gbm_peg_vs_meg.csv      the enrichment test
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

master = pd.read_csv(WORK / "11_arabidopsis" / "06_master_table.csv")
cls = dict(zip(master.AGI, master.meth_class))
sets = {"PEG": set(pd.read_csv(TAB / "consensus_pegs.csv").gene),
        "MEG": set(pd.read_csv(TAB / "consensus_megs.csv").gene)}
assert not sets["PEG"] & sets["MEG"], "a gene is in both consensus sets"

rows = []
for name, genes in sets.items():
    c = pd.Series([cls.get(g) for g in genes]).dropna().value_counts()
    rows.append(dict(set=name, n=len(genes), classified=int(c.sum()), UM=int(c.get("UM", 0)),
                     gbM=int(c.get("gbM", 0)), teM=int(c.get("teM", 0))))
counts = pd.DataFrame(rows).set_index("set")

p, m = counts.loc["PEG"], counts.loc["MEG"]
table = [[p.gbM, p.classified - p.gbM], [m.gbM, m.classified - m.gbM]]
odds, pval = fisher_exact(table)
test = pd.DataFrame([dict(comparison="gbM, consensus PEGs vs consensus MEGs",
                          peg_gbM=int(p.gbM), peg_classified=int(p.classified),
                          meg_gbM=int(m.gbM), meg_classified=int(m.classified),
                          pct_peg=100 * p.gbM / p.classified, pct_meg=100 * m.gbM / m.classified,
                          fold=(p.gbM / p.classified) / (m.gbM / m.classified),
                          odds_ratio=odds, fisher_p=pval)])

TAB.mkdir(parents=True, exist_ok=True)
counts.reset_index().to_csv(TAB / "arabidopsis_imprinted_classes.csv", index=False)
test.to_csv(TAB / "arabidopsis_gbm_peg_vs_meg.csv", index=False)
print(counts.to_string())
print(test.T.to_string(header=False))
