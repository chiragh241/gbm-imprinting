#!/usr/bin/env python
"""Figure 2a-d source table: CG methylation metagenes of four gene sets in four W22 tissues.

For each tissue, gene set and 100 bp bin: the pooled methylation (methylated over total CG reads
summed over all genes in the set), which Figure 2a-d draws, and, per gene, a 500 bp rolling
ratio (5 bins, within each 4 kb half of the axis, at least 5 reads) summarised as median and IQR
with the number of genes reaching the bin.

Tissues: 15-DAP endosperm, 15-DAP embryo, developing second leaf, premeiotic tassel.
Gene sets: gbM-PEG, gbM-biparental, UM-PEG, UM-biparental (W22 intervals).

Input   $GI_WORK/04_emseq_methylation/metagene/<tissue>_metagene_cg_fixed_pergene.tsv
        (04_emseq_methylation/08_metagene_per_gene.sh)
Output  $GI_RESULTS/tables/fig2ad_source_data.csv   pooled, median, IQR and n per bin
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK as GI_WORK, TAB

WORK = GI_WORK / "04_emseq_methylation" / "metagene"
TAB.mkdir(parents=True, exist_ok=True)

NUP, NIN, NDOWN = 20, 20, 20
NBIN = NUP + 2 * NIN + NDOWN
HALF = NUP + NIN                 # bins 0-39 are one contiguous 4 Kb stretch, 40-79 the other
WIN = 5                          # 500 bp rolling window for the per-gene ratio
MIN_READS = 5

PANELS = [("PEG__gbM", "gbM-PEG"), ("biparental__gbM", "gbM-biparental"),
          ("PEG__UM", "UM-PEG"), ("biparental__UM", "UM-biparental")]
TISSUES = ["wt_15dap", "embryo", "leaf", "tassel"]

rows = []
for t in TISSUES:
    d = pd.read_csv(WORK / f"{t}_metagene_cg_fixed_pergene.tsv", sep="\t")
    for key, _ in PANELS:
        s = d[d.set == key]
        if s.empty:
            continue
        genes = s.gene.unique()
        gi = {g: i for i, g in enumerate(genes)}
        M = np.zeros((len(genes), NBIN)); T = np.zeros((len(genes), NBIN))
        M[s.gene.map(gi).to_numpy(), s.bin.to_numpy()] = s.mc.to_numpy()
        T[s.gene.map(gi).to_numpy(), s.bin.to_numpy()] = s.total.to_numpy()
        pooled = np.where(T.sum(0) > 0, M.sum(0) / np.maximum(T.sum(0), 1), np.nan)
        Mr = np.full_like(M, np.nan); Tr = np.full_like(T, np.nan)
        for a, b in ((0, HALF), (HALF, NBIN)):
            mm = pd.DataFrame(M[:, a:b]).T.rolling(WIN, center=True, min_periods=3).sum().T
            tt = pd.DataFrame(T[:, a:b]).T.rolling(WIN, center=True, min_periods=3).sum().T
            Mr[:, a:b] = mm.to_numpy(); Tr[:, a:b] = tt.to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            per_gene = np.where(Tr >= MIN_READS, Mr / np.maximum(Tr, 1), np.nan)
        with np.errstate(invalid="ignore"):
            med = np.nanmedian(per_gene, axis=0)
            q1 = np.nanpercentile(per_gene, 25, axis=0)
            q3 = np.nanpercentile(per_gene, 75, axis=0)
        ngene = np.isfinite(per_gene).sum(axis=0)
        for b in range(NBIN):
            rows.append(dict(tissue=t, set=key, bin=b, pooled=pooled[b], median=med[b],
                             q1=q1[b], q3=q3[b], n_genes=int(ngene[b])))
    print(f"{t} done")

st = pd.DataFrame(rows)
st.to_csv(TAB / "fig2ad_source_data.csv", index=False)
print(f"wrote {TAB / 'fig2ad_source_data.csv'} ({len(st)} rows)")
