#!/usr/bin/env python
"""The genes and groups of the Figure 4 metagenes, with strand, length and nearest neighbour.

Gene groups are taken from the allele panel (03_allele_panel.py), so the metagene is built on
exactly the genes with enough allele-informative reads. Genes shorter than 60 bp are left out.

Input   $GI_DATA/references/B73v5_genes.gff3                     (00_references)
        $GI_WORK/10_chip_endosperm_allele/05_panel_source_data.csv  (03_allele_panel.py)
Output  $GI_WORK/10_chip_endosperm_allele/gene_meta_v5.tsv
          gene group strand length nb_dist nb_free_2k nb_free_5k
"""
from pathlib import Path
from collections import defaultdict
import csv
import re

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

D = WORK / "10_chip_endosperm_allele"
GFF = DATA / "references" / "B73v5_genes.gff3"
PANEL = D / "05_panel_source_data.csv"

MIN_LENGTH = 60      # bp; a guard only: the shortest panel gene is 609 bp

coords, strand, allgenes = {}, {}, defaultdict(list)
for ln in open(GFF):
    if ln.startswith("#"):
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "gene":
        continue
    m = re.search(r"ID=([^;\t]+)", ln)
    if not m:
        continue
    g = m.group(1)
    c = f[0][3:] if f[0].startswith("chr") else f[0]
    s, e = int(f[3]) - 1, int(f[4])
    coords[g] = (c, s, e)
    strand[g] = f[6]
    allgenes[c].append((s, e, g))
for c in allgenes:
    allgenes[c].sort()


def nearest(c, s, e, gid):
    best = None
    for os_, oe, og in allgenes.get(c, ()):
        if og == gid:
            continue
        d = s - oe if oe <= s else (os_ - e if os_ >= e else 0)
        best = d if best is None else min(best, d)
        if best == 0:
            break
    return best if best is not None else 10 ** 9


# Keyed on (gene, group), not gene -> group: the panel's groups overlap ("gbM-biparental (all)"
# contains the demethylated quartile), and a gene -> group dict would keep only the last group
# each gene appears in, silently dropping the quartile groups.
members = sorted({(r["gene"], r["group"]) for r in csv.DictReader(open(PANEL))})

meta_out = []
missing = 0
for gene, grp in members:
    if gene not in coords:
        missing += 1
        continue
    c, s, e = coords[gene]
    sd = strand[gene]
    L = e - s
    if L < MIN_LENGTH:
        continue
    nb = nearest(c, s, e, gene)
    meta_out.append((gene, grp, sd, L, nb, int(nb >= 2000), int(nb >= 5000)))

with open(D / "gene_meta_v5.tsv", "w", newline="\n") as fh:
    fh.write("gene\tgroup\tstrand\tlength\tnb_dist\tnb_free_2k\tnb_free_5k\n")
    for r in meta_out:
        fh.write("\t".join(map(str, r)) + "\n")

print(f"wrote {len(meta_out):,} gene-group rows ({missing} panel genes absent from the v5 GFF)")
byg = defaultdict(list)
for m in meta_out:
    byg[m[1]].append(m)
for grp in sorted(byg):
    v = byg[grp]
    f2 = sum(1 for m in v if m[5])
    f5 = sum(1 for m in v if m[6])
    print(f"  {grp:34s} n={len(v):5d}   flank clean at 2 kb: {f2:5d}   at 5 kb: {f5:5d}")
