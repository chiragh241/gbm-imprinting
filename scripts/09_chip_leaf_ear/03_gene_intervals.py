#!/usr/bin/env python
"""Gene intervals for Supplemental Fig. S1, in AGPv4 coordinates to match the GSE120304 BEDs.

v5 -> v4 via the MaizeGDB B73v5_to_B73v4.tsv cross-reference. Only genes with EXACTLY ONE v4 counterpart are
kept: a v5 gene mapping to several v4 models is usually a tandem array, and counting reads
across all of them would smear signal in a way that differs systematically between gene
classes, which is precisely the comparison this figure makes.

Chromosome naming: the ChIP BEDs use bare "1".."10"; the v4 GFF uses "Chr1".."Chr10". Output
uses the ChIP convention, and drops contigs / Mt / Pt.

Input   $GI_DATA/references/B73v5_to_B73v4.tsv, $GI_WORK/02_gene_sets/master_table.csv,
        $GI_WORK/09_chip_leaf_ear/B73v4_genes.bed   (02_annotation.sh)
Output  $GI_WORK/09_chip_leaf_ear/fig1e_gene_intervals.bed   chrom start end gene set
        $GI_WORK/09_chip_leaf_ear/fig1e_interval_counts.csv
"""
from pathlib import Path
import csv
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

D = WORK / "09_chip_leaf_ear"
MAP = DATA / "references" / "B73v5_to_B73v4.tsv"
MASTER = WORK / "02_gene_sets" / "master_table.csv"
V4BED = D / "B73v4_genes.bed"
OUT = D / "fig1e_gene_intervals.bed"

GROUPS = {("PEG", "gbM"): "gbM-PEG", ("PEG", "UM"): "UM-PEG",
          ("MEG", "UM"): "UM-MEG", ("MEG", "teM"): "teM-MEG",
          ("biparental", "gbM"): "gbM-biparental",
          ("biparental", "UM"): "UM-biparental"}

v5_to_v4 = {}
for ln in open(MAP):
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 2 and f[1]:
        ids = [x for x in f[1].split(",") if x]
        if len(ids) == 1:                 # unique only
            v5_to_v4[f[0]] = ids[0]

v4coord = {}
for ln in open(V4BED):
    c, s, e, g = ln.rstrip("\n").split("\t")
    if c.startswith("Chr") and c[3:].isdigit():
        v4coord[g] = (c[3:], int(s), int(e))   # "Chr1" -> "1"

rows, stats = [], {}
for r in csv.DictReader(open(MASTER)):
    key = (r["imp_class"], r["meth_class"])
    if key not in GROUPS:
        continue
    g = GROUPS[key]
    st = stats.setdefault(g, dict(total=0, unique_v4=0, on_chr=0))
    st["total"] += 1
    v4 = v5_to_v4.get(r["b73_gene"])
    if not v4:
        continue
    st["unique_v4"] += 1
    co = v4coord.get(v4)
    if not co:
        continue
    st["on_chr"] += 1
    rows.append((co[0], co[1], co[2], r["b73_gene"], g))

rows.sort(key=lambda x: (int(x[0]), x[1]))
with open(OUT, "w", newline="") as fh:
    for c, s, e, gene, g in rows:
        fh.write(f"{c}\t{s}\t{e}\t{gene}\t{g}\n")

df = pd.DataFrame([dict(set=k, **v) for k, v in stats.items()])
df["kept"] = df["set"].map(pd.Series(rows and pd.Series([r[4] for r in rows]).value_counts()))
df.to_csv(D / "fig1e_interval_counts.csv", index=False)
print(df.to_string(index=False))
print(f"\nwrote {OUT}  ({len(rows):,} intervals)")
print(f"median gene length: {int(pd.Series([r[2]-r[1] for r in rows]).median()):,} bp")
