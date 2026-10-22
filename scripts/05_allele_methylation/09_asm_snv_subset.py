#!/usr/bin/env python3
"""The variants allele-specific methylation is called at: genic variants plus a control sample.

cgmaptools asm costs about 0.09 s per variant per library and context; over all 5,158,268
variants that is about a month for the 18 jobs. So it runs on a subset:

  genic    variants inside the W22 gene intervals of all 16 gene sets, merged and extended by
           2 kb on each side (clipped to the sequence ends), merged again: 10,955 regions over
           128.7 Mb, 472,785 variants. asm links a cytosine to a variant within one read pair,
           so the 2 kb margin keeps every gene-body cytosine reachable from a kept variant.
  control  every 93rd of the remaining (non-genic) variants, in VCF order: 50,381 variants, a
           systematic genome-wide sample that is reproducible without a seed and covers every
           sequence (284 on "unmapped"). Their allele-specific methylation is the non-genic
           baseline (11_gene_summary.py).

Input   $GI_WORK/02_gene_sets/w22/*.bed                        (02_gene_sets/03_build_sets.py)
        $GI_WORK/05_allele_methylation/snv/endosperm_snv.vcf   (08_snv.sh)
        $GI_DATA/references/W22v2.fa.fai                       (04_emseq_methylation/02_bisulfite_reference.sh)
Output  $GI_WORK/05_allele_methylation/snv/
          asm_regions.bed        the genic regions
          .keep_genic.txt        genic variants, "chrom:pos", VCF order
          .keep_ctrl.txt         control variants, "chrom:pos", VCF order
          .keep.txt              both, sorted
          endosperm_snv.asm.vcf  the VCF restricted to these variants (523,166)
"""
import glob
import sys
from bisect import bisect_right
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

SNV = WORK / "05_allele_methylation" / "snv"
FLANK, EVERY = 2000, 93

size = {}
with open(DATA / "references" / "W22v2.fa.fai") as fh:
    for line in fh:
        f = line.split("\t")
        size[f[0]] = int(f[1])

intervals = []
for bed in sorted(glob.glob(str(WORK / "02_gene_sets" / "w22" / "*.bed"))):
    with open(bed) as fh:
        for line in fh:
            c, s, e = line.split("\t")[:3]
            intervals.append((c, int(s), int(e)))


def merge(ivs):
    """bedtools merge: overlapping or book-ended intervals, chromosomes in text order."""
    out = []
    for c, s, e in sorted(ivs):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


genes = merge(intervals)
regions = merge([(c, max(0, s - FLANK), min(size[c], e + FLANK)) for c, s, e in genes])
with open(SNV / "asm_regions.bed", "w", newline="\n") as fh:
    fh.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in regions)
print(f"{len(intervals):,} gene intervals -> {len(genes):,} merged -> {len(regions):,} regions "
      f"over {sum(e - s for _, s, e in regions) / 1e6:.1f} Mb")

starts = {}
for c, s, e in regions:
    starts.setdefault(c, ([], []))
    starts[c][0].append(s)
    starts[c][1].append(e)


def in_region(c, pos):
    """VCF POS is 1-based: the variant base is pos-1 in the 0-based half-open regions."""
    if c not in starts:
        return False
    ss, ee = starts[c]
    i = bisect_right(ss, pos - 1) - 1
    return i >= 0 and pos - 1 < ee[i]


header, rows, genic, pool = [], [], [], []
with open(SNV / "endosperm_snv.vcf") as fh:
    for line in fh:
        if line.startswith("#"):
            header.append(line)
            continue
        f = line.split("\t", 2)
        key = f"{f[0]}:{f[1]}"
        rows.append((key, line))
        (genic if in_region(f[0], int(f[1])) else pool).append(key)
ctrl = [k for i, k in enumerate(pool, 1) if i % EVERY == 0]
keep = sorted(set(genic) | set(ctrl))

for name, keys in ((".keep_genic.txt", genic), (".keep_ctrl.txt", ctrl), (".keep.txt", keep)):
    with open(SNV / name, "w", newline="\n") as fh:
        fh.writelines(k + "\n" for k in keys)
kset = set(keep)
with open(SNV / "endosperm_snv.asm.vcf", "w", newline="\n") as fh:
    fh.writelines(header)
    fh.writelines(line for k, line in rows if k in kset)
print(f"variants {len(rows):,}: genic {len(genic):,}, control {len(ctrl):,} "
      f"(of {len(pool):,} non-genic), kept {len(keep):,}")
