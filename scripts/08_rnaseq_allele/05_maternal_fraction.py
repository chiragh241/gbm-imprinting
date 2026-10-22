#!/usr/bin/env python3
"""Per-gene maternal expression ratio in each of the four crosses of Xu et al. 2022.

Groups: BM (B73 wild-type mother x Mo17), bM (zmros1ab mother x Mo17), MB (Mo17 x B73 wild-type
father), Mb (Mo17 x zmros1ab father). For each gene and group, maternal and paternal fragment
counts are summed over the three replicates, and the maternal ratio (maternal / (maternal +
paternal)) is reported when that sum is at least 20; otherwise it is left blank.

Input   $GI_WORK/08_rnaseq_allele/counts/<run>.{B73,Mo17}.txt   (04_align_split_count.sh)
        data/rnaseq/samplesheet.tsv
Output  $GI_WORK/08_rnaseq_allele/per_gene_maternal_fraction.tsv
          gene, <group>_mat, <group>_pat, <group>_frac (4 decimals) for BM, bM, MB, Mb
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import SHIPPED, WORK

CNT = WORK / "08_rnaseq_allele" / "counts"
OUT = WORK / "08_rnaseq_allele" / "per_gene_maternal_fraction.tsv"
MIN_GROUP_READS = 20     # min maternal+paternal fragments in a group to score a gene

samples = []
with open(SHIPPED / "rnaseq" / "samplesheet.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        r["group"] = r["sample"].split("_")[0]      # BM/bM/MB/Mb
        samples.append(r)


def read_fc(path):
    """featureCounts output -> {geneid: count}."""
    d = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            c = line.rstrip("\n").split("\t")
            d[c[0]] = int(c[-1])
    return d


# gene -> group -> [maternal_sum, paternal_sum]
agg = defaultdict(lambda: defaultdict(lambda: [0, 0]))
groups_with_data = defaultdict(int)
for s in samples:
    p1 = CNT / f"{s['run']}.B73.txt"
    p2 = CNT / f"{s['run']}.Mo17.txt"
    if not (os.path.exists(p1) and os.path.exists(p2)):
        sys.exit(f"missing counts for {s['sample']} ({s['run']})")
    groups_with_data[s['group']] += 1
    b73 = read_fc(p1)
    mo17 = read_fc(p2)
    mat_is_b73 = s["maternal_genome"].startswith("B73")
    for g in set(b73) | set(mo17):
        b, m = b73.get(g, 0), mo17.get(g, 0)
        mat, pat = (b, m) if mat_is_b73 else (m, b)
        agg[g][s["group"]][0] += mat
        agg[g][s["group"]][1] += pat
print("samples per group:", dict(groups_with_data))


def frac(g, grp):
    mat, pat = agg[g][grp]
    tot = mat + pat
    return (mat / tot, tot) if tot >= MIN_GROUP_READS else (None, tot)


groups = ["BM", "bM", "MB", "Mb"]
with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene"] + [f"{g}_mat" for g in groups] + [f"{g}_pat" for g in groups]
               + [f"{g}_frac" for g in groups])
    for g in sorted(agg):
        row = [g]
        for grp in groups:
            row.append(agg[g][grp][0])
        for grp in groups:
            row.append(agg[g][grp][1])
        for grp in groups:
            fr, _ = frac(g, grp)
            row.append("" if fr is None else f"{fr:.4f}")
        w.writerow(row)
print(f"wrote {OUT} ({len(agg):,} genes)")
