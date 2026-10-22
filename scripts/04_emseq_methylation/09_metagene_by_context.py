#!/usr/bin/env python
"""Per-gene, per-bin CG, CHG and CHH counts on the fixed-kb metagene axis. Pure Python.

Supplemental Figure 3 needs the CHG (a-d) and CHH (e-h) versions of Figure 2a-d. The CG
profile (08_metagene_per_gene.sh) streams the allc through bedtools; this does the same work in
one Python process with a flat memory profile, for all three contexts at once.

Memory is bounded by construction. Counts accumulate into two numpy arrays of length
n_bins * 3, one slot per (bin, context), so nothing grows with the number of cytosines read.
For 710,135 bins that is about 34 MB, whatever the size of the allc file.

Bin lookup is a bucket index rather than an interval tree. Bins are fixed 100 bp, so a
cytosine at position p can only fall in a bin whose start lies in bucket p//BUCKET or the one
below it. That makes assignment O(1) per cytosine with a small constant, and it does not
assume the allc file is sorted, which it is not: the first record is on chromosome 3.

Overlapping bins are all counted, matching `bedtools intersect -wa -wb`, which reports one row
per overlapping interval. Genes do overlap, so this matters.

allc columns: chrom, pos (1-based), strand, trinucleotide, mc, total, ...
Context from the trinucleotide: CG if the second base is G, else CHG if the third base is G,
else CHH. Same rule 08_metagene_per_gene.sh applies for CG.

Input   $GI_WORK/04_emseq_methylation/allc/allc_<tissue>_merged.tsv.gz   (04_merge_allc.sh)
        $GI_WORK/04_emseq_methylation/metagene/metagene_bins_fixed.bed   (07_metagene_bins.sh)
Output  $GI_WORK/04_emseq_methylation/metagene/<tissue>_metagene_bycontext_pergene.tsv
          context, set, gene, bin, mc, total
"""
import gzip
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK as GI_WORK

ALLC = GI_WORK / "04_emseq_methylation" / "allc"
WORK = GI_WORK / "04_emseq_methylation" / "metagene"
BINS = WORK / "metagene_bins_fixed.bed"

BUCKET = 1000          # bins are 100 bp, so a hit can only be in this bucket or the one below
CTX = {"CG": 0, "CHG": 1, "CHH": 2}
CTX_NAME = ["CG", "CHG", "CHH"]
TISSUES = sys.argv[1:] or ["wt_15dap", "embryo", "leaf", "tassel"]


def load_bins():
    chrom, start, end, label = [], [], [], []
    with BINS.open() as f:
        for line in f:
            c, s, e, st, g, b = line.rstrip("\n").split("\t")[:6]
            chrom.append(c)
            start.append(int(s))
            end.append(int(e))
            label.append((st, g, b))
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    index = {}
    for i, (c, s) in enumerate(zip(chrom, start)):
        index.setdefault((c, s // BUCKET), []).append(i)
    return start, end, label, index


def run(tissue, start, end, label, index):
    src = ALLC / f"allc_{tissue}_merged.tsv.gz"
    out = WORK / f"{tissue}_metagene_bycontext_pergene.tsv"
    if not src.exists():
        print(f"### {tissue}: no allc, skip")
        return
    if out.exists() and out.stat().st_size:
        print(f"### skip {tissue} (done)")
        return

    n = len(start)
    mc = np.zeros(n * 3, dtype=np.int64)
    tot = np.zeros(n * 3, dtype=np.int64)
    G = ord("G")
    t0, seen, hit = time.time(), 0, 0

    with gzip.open(src, "rb") as f:
        for line in f:
            seen += 1
            p = line.split(b"\t", 6)
            if len(p) < 6:
                continue
            c = p[0].decode()
            pos = int(p[1])
            k = pos // BUCKET
            rows = index.get((c, k))
            prev = index.get((c, k - 1))
            if rows is None and prev is None:
                continue
            tri = p[3]
            ctx = 0 if tri[1:2] == b"G" else (1 if tri[2:3] == b"G" else 2)
            m, t = int(p[4]), int(p[5])
            for group in (rows, prev):
                if not group:
                    continue
                for i in group:
                    # bed is half-open on 0-based coords: start <= pos-1 < end
                    if start[i] < pos <= end[i]:
                        j = i * 3 + ctx
                        mc[j] += m
                        tot[j] += t
                        hit += 1
            if seen % 50_000_000 == 0:
                print(f"    {seen:,} cytosines, {hit:,} bin hits, "
                      f"{time.time() - t0:.0f}s", flush=True)

    nz = np.nonzero(tot)[0]
    tmp = out.with_suffix(".tsv.part")
    with tmp.open("w", newline="\n") as fh:
        for j in nz:
            st, g, b = label[j // 3]
            fh.write(f"{CTX_NAME[j % 3]}\t{st}\t{g}\t{b}\t{mc[j]}\t{tot[j]}\n")
    tmp.replace(out)
    print(f"### {tissue}: {len(nz):,} context/set/gene/bin rows, "
          f"{seen:,} cytosines read, {time.time() - t0:.0f}s")


start, end, label, index = load_bins()
print(f"bins: {len(start):,} over {len({c for c, _ in index}):,} sequences\n")
for t in TISSUES:
    run(t, start, end, label, index)
