#!/usr/bin/env python
"""Metagene bins for the Ricci 2019 leaf/ear H3K27me3 profiles, in AGPv4 coordinates.

Two bin schemes are emitted in one file so the read pass over 2.2 GB of BEDs happens once:

  fig2   20 upstream bins x 100 bp (2 kb) | 60 body bins, scaled | 20 downstream x 100 bp
         The scaled scheme of the earlier Figure 2a-d draft; kept because the profiles and
         domain-extent numbers are computed for both schemes.
  wide   20 upstream bins x 250 bp (5 kb) | 60 body bins, scaled | 20 downstream x 250 bp
         UM-PEG signal is still 6-7x enriched at -3 kb, so 2 kb of flank cannot show where
         the domain ends. This one can.

Bins are always ordered 5'->3', so plus and minus strand genes stack.

Genes also carry a neighbour-free flag. gbM-PEGs sit in gene-dense neighbourhoods and
UM-PEGs sit in gene deserts (37% vs 3% with another gene within 2 kb), so for 40% of
gbM-PEGs the "flank" is another gene's body. Any claim about flanks has to be checked on
the subset where the flank is really intergenic.

Supplemental Fig. S1 uses the wide scheme.

Input   $GI_DATA/references/B73v4.gff3.gz, B73v5_to_B73v4.tsv          (00_references)
        $GI_WORK/09_chip_leaf_ear/fig1e_gene_intervals.bed              (03_gene_intervals.py)
Output  $GI_WORK/09_chip_leaf_ear/metagene_bins_v4.bed   chrom start end set gene scheme bin
        $GI_WORK/09_chip_leaf_ear/gene_meta_v4.tsv       gene set strand length nb_dist
                                                         nb_free_2k nb_free_5k
"""
from pathlib import Path
import gzip
import re
from collections import defaultdict

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

D = WORK / "09_chip_leaf_ear"
GFF = DATA / "references" / "B73v4.gff3.gz"
INTERVALS = D / "fig1e_gene_intervals.bed"
MAP = DATA / "references" / "B73v5_to_B73v4.tsv"

# The six gene sets drawn in Supplemental Figure 1: both PEG sets, the two biparental
# comparison sets and the two MEG sets.
SETS = {"gbM-PEG", "UM-PEG", "gbM-biparental", "UM-biparental", "UM-MEG", "teM-MEG"}
SCHEMES = {"fig2": dict(nup=20, nbody=60, ndown=20, binw=100),
           "wide": dict(nup=20, nbody=60, ndown=20, binw=250)}

# ---- v4 annotation: strand, and every gene interval for the neighbour distance ----------
chrom_len, strand, allgenes = {}, {}, defaultdict(list)
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln.startswith("##sequence-region"):
            f = ln.split()
            if len(f) >= 4:
                chrom_len[f[1]] = int(f[3])
            continue
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        c = f[0][3:] if f[0].startswith("Chr") else f[0]
        # This GFF has literal tabs inside description=, so the attribute column is
        # split across fields 9..NF and gene_id= is often not in f[8] at all.
        # ID=gene:... always leads the attributes, so match it against the whole line.
        gid = re.search(r"ID=gene:([^;\t]+)", ln)
        if not gid:
            continue
        strand[gid.group(1)] = f[6]
        allgenes[c].append((int(f[3]) - 1, int(f[4]), gid.group(1)))
for c in allgenes:
    allgenes[c].sort()

v5_to_v4 = {}
for ln in open(MAP):
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 2 and f[1]:
        ids = [x for x in f[1].split(",") if x]
        if len(ids) == 1:
            v5_to_v4[f[0]] = ids[0]


def nearest(c, s, e, gid_v4):
    """Distance in bp to the nearest other annotated gene, 0 if one overlaps."""
    best = None
    for os_, oe, og in allgenes.get(c, ()):
        if og == gid_v4:
            continue
        d = s - oe if oe <= s else (os_ - e if os_ >= e else 0)
        best = d if best is None else min(best, d)
        if best == 0:
            break
    return best if best is not None else 10 ** 9


bins_out, meta_out = [], []
kept = defaultdict(int)
dropped_short = defaultdict(int)

for ln in open(INTERVALS):
    c, s, e, gene, st = ln.rstrip("\n").split("\t")
    if st not in SETS:
        continue
    s, e = int(s), int(e)
    v4 = v5_to_v4.get(gene, "")
    sd = strand.get(v4, "+")
    L = e - s
    if L < SCHEMES["fig2"]["nbody"]:
        dropped_short[st] += 1
        continue
    clen = chrom_len.get(c, 10 ** 12)
    nb = nearest(c, s, e, v4)
    kept[st] += 1
    meta_out.append((gene, st, sd, L, nb, int(nb >= 2000), int(nb >= 5000)))

    for scheme, p in SCHEMES.items():
        nup, nbody, ndown, bw = p["nup"], p["nbody"], p["ndown"], p["binw"]
        edges = []
        # upstream flank, then body, then downstream flank, in genomic order for + strand
        for i in range(nup):
            edges.append((s - (nup - i) * bw, s - (nup - i - 1) * bw))
        for i in range(nbody):
            edges.append((s + L * i // nbody, s + L * (i + 1) // nbody))
        for i in range(ndown):
            edges.append((e + i * bw, e + (i + 1) * bw))
        if sd == "-":
            edges = [(b, a) for a, b in reversed(edges)]      # keep 5'->3' bin order
            edges = [(min(a, b), max(a, b)) for a, b in edges]
        for k, (bs, be) in enumerate(edges):
            bs, be = max(0, bs), min(clen, be)
            if be > bs:
                bins_out.append((c, bs, be, st, gene, scheme, k))

bins_out.sort(key=lambda r: (int(r[0]) if r[0].isdigit() else 99, r[1]))
with open(D / "metagene_bins_v4.bed", "w", newline="\n") as fh:
    for r in bins_out:
        fh.write("\t".join(map(str, r)) + "\n")
with open(D / "gene_meta_v4.tsv", "w", newline="\n") as fh:
    fh.write("gene\tset\tstrand\tlength\tnb_dist\tnb_free_2k\tnb_free_5k\n")
    for r in meta_out:
        fh.write("\t".join(map(str, r)) + "\n")

print(f"wrote {len(bins_out):,} bin rows for {len(meta_out):,} genes")
for st in sorted(SETS):
    n = kept[st]
    f2 = sum(1 for m in meta_out if m[1] == st and m[5])
    f5 = sum(1 for m in meta_out if m[1] == st and m[6])
    print(f"  {st:16s} n={n:5d}   flank clean at 2 kb: {f2:5d} ({100*f2/max(n,1):3.0f}%)"
          f"   at 5 kb: {f5:5d} ({100*f5/max(n,1):3.0f}%)   dropped short: {dropped_short[st]}")
