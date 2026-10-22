#!/usr/bin/env python
"""Metagene bins with a compact gene body for Figure 4: 20 x 1 kb flanks, 10 body bins, B73 v5.

With 15 kb flanks, a scaled gene body of 60 of 100 bins lets a gene of a few kb fill most of
the axis, which misstates its length. Here the body is 10 bins, each a tenth of the gene, which
the figure draws at the width of one 1 kb flank bin.

Genes and groups from gene_meta_v5.tsv (04_gene_meta.py), 5' to 3' order. Bin index runs
0-49: 0-19 upstream, 20-29 body, 30-49 downstream. The flanks reach 20 kb; Figure 4 shows 15 kb,
the first distance at which the maternal and paternal curves of both PEG sets have met on both
sides.

Input   $GI_DATA/references/B73v5_genes.gff3, $GI_WORK/10_chip_endosperm_allele/gene_meta_v5.tsv
Output  $GI_WORK/10_chip_endosperm_allele/metagene_bins_v5_compact.bed
          chrom start end group gene scheme bin
"""
from pathlib import Path
import csv
import re

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

D = WORK / "10_chip_endosperm_allele"
GFF = DATA / "references" / "B73v5_genes.gff3"
NUP, NBODY, NDOWN, BINW = 20, 10, 20, 1000
SCHEME = "compact"

coords = {}
for ln in open(GFF):
    if ln.startswith("#"):
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    m = re.search(r"ID=([^;\t]+)", ln)
    if not m:
        continue
    c = f[0][3:] if f[0].startswith("chr") else f[0]
    coords[m.group(1)] = (c, int(f[3]) - 1, int(f[4]), f[6])

members = [(r["gene"], r["group"]) for r in
           csv.DictReader(open(D / "gene_meta_v5.tsv"), delimiter="\t")]
rows = []
for gene, grp in members:
    c, s, e, sd = coords[gene]
    L = e - s
    edges = [(s - (NUP - i) * BINW, s - (NUP - i - 1) * BINW) for i in range(NUP)]
    edges += [(s + L * i // NBODY, s + L * (i + 1) // NBODY) for i in range(NBODY)]
    edges += [(e + i * BINW, e + (i + 1) * BINW) for i in range(NDOWN)]
    if sd == "-":
        edges = [(min(a, b), max(a, b)) for a, b in reversed(edges)]
    for k, (bs, be) in enumerate(edges):
        bs = max(0, bs)
        if be > bs:
            rows.append((c, bs, be, grp, gene, SCHEME, k))

rows.sort(key=lambda r: (int(r[0]) if r[0].isdigit() else 99, r[1]))
with open(D / "metagene_bins_v5_compact.bed", "w", newline="\n") as fh:
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(f"wrote {len(rows):,} compact-scheme bins for {len(members):,} gene-group pairs")
