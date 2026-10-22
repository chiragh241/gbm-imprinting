"""Correct interval containment for gene-body assignment.

Replaces a `searchsorted(starts, pos) - 1` lookup that inspected only the LAST gene
starting before a position. Maize gene models overlap: measured on the 16 W22 class
BEDs, 275 of 12,842 genes (2.14%) overlap their predecessor and up to 3 can cover a
single base. Under the old lookup a cytosine inside an earlier, still-open gene failed
the `pos <= end` test against the wrong candidate and was DISCARDED rather than
assigned - silent data loss, not a misassignment.

Semantics: BED is 0-based half-open, asm/CGmap positions are 1-based, so a position p
is inside a gene when `start < p <= end`. A position inside several overlapping genes
is reported for EVERY one of them, which is what the interval query actually means and
what bedtools would return; callers that need one gene per site must say so explicitly.

The required backward scan depth is derived from the gene set rather than hardcoded, so
the module stays correct if the sets are rebuilt with different models.
"""
import numpy as np


def _lookback(starts, ends):
    """Max number of earlier genes that can still be open at a later gene's start."""
    maxend = np.maximum.accumulate(ends)
    worst = 0
    for i in range(1, len(starts)):
        j = i
        while j > 0 and maxend[j - 1] > starts[i]:
            j -= 1
        if i - j > worst:
            worst = i - j
    return worst


def build_index(genes, chrom_col="chrom", start_col="start", end_col="end"):
    """Per-chromosome sorted arrays plus the scan depth that chromosome needs.

    `genes` must have a unique positional index (call reset_index(drop=True) first);
    the returned gene ids are positions into it.
    """
    idx = {}
    for chrom, sub in genes.groupby(chrom_col, sort=False):
        sub = sub.sort_values(start_col)
        starts = sub[start_col].to_numpy()
        ends = sub[end_col].to_numpy()
        ids = sub.index.to_numpy()
        idx[str(chrom)] = (starts, ends, ids, _lookback(starts, ends))
    return idx


def assign(chrom, pos, index):
    """Every (query position, gene) containment pair.

    chrom : array of chromosome names, one per query
    pos   : array of 1-based positions
    index : output of build_index

    Returns (qi, gi): parallel arrays of query row numbers and gene ids. A query
    contained in n genes appears n times; a query in none does not appear.
    """
    chrom = np.asarray(chrom).astype(str)
    pos = np.asarray(pos)
    qi_out, gi_out = [], []

    for c in np.unique(chrom):
        if c not in index:
            continue
        starts, ends, ids, look = index[c]
        sel = np.flatnonzero(chrom == c)
        p = pos[sel]
        # last gene whose start is strictly less than p
        base = np.searchsorted(starts, p, side="left") - 1
        for k in range(look + 1):
            cand = base - k
            live = cand >= 0
            if not live.any():
                break
            hit = np.zeros(len(p), dtype=bool)
            hit[live] = p[live] <= ends[cand[live]]
            if hit.any():
                qi_out.append(sel[hit])
                gi_out.append(ids[cand[hit]])

    if not qi_out:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    return np.concatenate(qi_out), np.concatenate(gi_out)
