#!/usr/bin/env python
"""Arabidopsis biparental genes by methylation class.

Biparental genes are every gene with a Bewick methylation class that is on none of the three
screens' PEG or MEG lists (listed_imprinted_genes.csv, which holds both consensus sets and
Pignatta's high-confidence calls too). A gene one study called imprinted has no place in a
biparental control, so this is stricter than the master table's impr_call (biparental unless a
Pignatta 2-of-3 call): 26,665 genes (20,771 UM, 4,851 gbM, 1,043 teM) against 27,093.

Input   $GI_WORK/11_arabidopsis/06_master_table.csv
        $GI_RESULTS/tables/listed_imprinted_genes.csv   (03_consensus_sets.py)
Output  $GI_RESULTS/tables/arabidopsis_biparental_class_counts.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

RES = WORK / "11_arabidopsis"

master = pd.read_csv(RES / "06_master_table.csv")
imprinted = set(pd.read_csv(TAB / "listed_imprinted_genes.csv").gene)

classified = master[master.meth_class.notna()]
bip = classified[~classified.AGI.isin(imprinted)]
loose = master[master.impr_call == "biparental"]

counts = bip.meth_class.value_counts()
rows = [dict(definition="strict (excludes every imprinted list; canonical)",
             total=len(bip), UM=int(counts.get("UM", 0)), gbM=int(counts.get("gbM", 0)),
             teM=int(counts.get("teM", 0)))]
lc = loose.meth_class.value_counts()
rows.append(dict(definition="loose (06_master_table impr_call == biparental)",
                 total=len(loose), UM=int(lc.get("UM", 0)), gbM=int(lc.get("gbM", 0)),
                 teM=int(lc.get("teM", 0))))
out = pd.DataFrame(rows)
TAB.mkdir(parents=True, exist_ok=True)
out.to_csv(TAB / "arabidopsis_biparental_class_counts.csv", index=False)
print(out.to_string(index=False))
print(f"\nimprinted genes excluded: {len(imprinted & set(master.AGI))} of the {len(imprinted)} "
      f"on the three screens' lists")
