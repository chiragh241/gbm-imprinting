#!/usr/bin/env python
"""Build the target-gene table: AGI coordinates (TAIR10) x imprinting call x
methylation class. Backbone for the DMR intersection + per-gene methylation test.

Inputs:
  - $GI_DATA/references/Arabidopsis_thaliana.TAIR10.57.gff3.gz (Ensembl Plants 57; 00_references)
  - $GI_WORK/11_arabidopsis/06_master_table.csv (impr_call, meth_class per AGI)
Output:
  - $GI_WORK/11_arabidopsis/gene_targets.csv
"""
import gzip, re, sys
from pathlib import Path
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, SHIPPED, WORK

GFF = DATA / "references" / "Arabidopsis_thaliana.TAIR10.57.gff3.gz"
MASTER = WORK / "11_arabidopsis" / "06_master_table.csv"
OUT = WORK / "11_arabidopsis" / "gene_targets.csv"

# --- parse gene coordinates from GFF ---
rows = []
gid_re = re.compile(r"gene_id=([^;]+)")
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        chrom = f[0]
        if chrom not in {"1", "2", "3", "4", "5"}:
            continue
        m = gid_re.search(f[8])
        if not m:
            continue
        rows.append((m.group(1), chrom, int(f[3]), int(f[4]), f[6]))
coords = pd.DataFrame(rows, columns=["AGI", "chr", "start", "end", "strand"])
print(f"GFF nuclear genes: {len(coords)}")

# --- merge with imprinting + methylation classes ---
master = pd.read_csv(MASTER)[["AGI", "impr_call", "meth_class", "impr_meth",
                              "n_crosses_support"]]
tab = coords.merge(master, on="AGI", how="inner")
tab["tss"] = tab.apply(lambda r: r.start if r.strand == "+" else r.end, axis=1)
tab.to_csv(OUT, index=False)
print(f"merged (coords + master): {len(tab)} genes  -> {OUT}")

print("\n=== imprinting call counts (with coords + meth class) ===")
print(tab.impr_call.value_counts())
print("\n=== PEG / MEG split by methylation class ===")
sub = tab[tab.impr_call.isin(["PEG", "MEG"])]
print(pd.crosstab(sub.impr_call, sub.meth_class, margins=True))
