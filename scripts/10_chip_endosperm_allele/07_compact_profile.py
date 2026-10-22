#!/usr/bin/env python
"""Maternal and paternal H3K27me3 on the compact bins, profiles and per-gene body values.

Bins are 20 x 1 kb flanks and 10 length-scaled body bins (05_compact_bins.py); each read is
counted once, at its leftmost base (06_compact_counts.sh), so the summed body counts equal the
gene-level counts exactly.

Estimator: per allele, reads divided by allele-informative bp, normalised to that allele's
genome-wide density of assigned reads per informative bp (so 1 is the allele's genome average
and the 2:1 maternal dosage of endosperm divides out), the two reciprocal crosses averaged.
Profiles pool genes before the ratio, blank a bin below MIN_READS pooled reads and carry a 95%
bootstrap interval over genes. Per-gene values pool the 10 body bins into one value per gene,
for genes with at least MIN_PER_RUN assigned reads over the body in each cross.

Input   $GI_WORK/10_chip_endosperm_allele/counts_allele_compact/*.bincounts.tsv, *.total.txt
        $GI_WORK/10_chip_endosperm_allele/bin_informative_bp_compact.tsv, genome_informative_bp.txt
        $GI_WORK/10_chip_endosperm_allele/gene_meta_v5.tsv
Output  $GI_WORK/10_chip_endosperm_allele/matpat_profiles_compact.csv
        $GI_WORK/10_chip_endosperm_allele/matpat_per_gene_body_compact.csv
"""
from pathlib import Path
from collections import defaultdict
import csv
import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK

D = WORK / "10_chip_endosperm_allele"
CNT = D / "counts_allele_compact"
CROSS = {"SRR748861": ("genome1", "genome2"),    # B73 mother
         "SRR748862": ("genome2", "genome1")}    # Mo17 mother
GROUPS = ["gbM-PEG", "UM-PEG", "gbM-biparental demethylated",
          "gbM-biparental not demethylated", "UM-biparental", "gbM-biparental (all)"]
SCH = "compact"
NBIN, NUP, NBODY = 50, 20, 10
BODY = range(NUP, NUP + NBODY)
BOOT, MIN_READS, MIN_PER_RUN = 1000, 30, 10
RNG = np.random.default_rng(0)
GENOME_INF_BP = int(open(D / "genome_informative_bp.txt").read().split()[0])

# keyed on (group, gene): "gbM-biparental (all)" contains the demethylated quartile
meta = {(r["group"], r["gene"]) for r in
        csv.DictReader(open(D / "gene_meta_v5.tsv"), delimiter="\t")}
genes_of = defaultdict(list)
for grp, g in sorted(meta):
    genes_of[grp].append(g)
gidx = {(grp, g): i for grp in genes_of for i, g in enumerate(genes_of[grp])}

counts = {grp: {(run, side): np.zeros((len(genes_of[grp]), NBIN))
                for run in CROSS for side in ("mat", "pat")} for grp in genes_of}
infbp = {grp: np.zeros((len(genes_of[grp]), NBIN)) for grp in genes_of}
with open(D / "bin_informative_bp_compact.tsv") as fh:
    next(fh)
    for ln in fh:
        grp, gene, sch, b, bp = ln.rstrip("\n").split("\t")
        if (grp, gene) in gidx:
            infbp[grp][gidx[(grp, gene)], int(b)] += int(bp)
tot = {}
for run, (matg, patg) in CROSS.items():
    for gn in ("genome1", "genome2"):
        tot[(run, gn)] = int(open(CNT / f"{run}.{gn}.total.txt").read().split()[0])
    for gn, side in ((matg, "mat"), (patg, "pat")):
        with open(CNT / f"{run}.{gn}.bincounts.tsv") as fh:
            for ln in fh:
                grp, gene, sch, b, n = ln.rstrip("\n").split("\t")
                if (grp, gene) in gidx:
                    counts[grp][(run, side)][gidx[(grp, gene)], int(b)] += int(n)
side_gn = {"mat": 0, "pat": 1}
gdens = {(run, side): tot[(run, CROSS[run][side_gn[side]])] / GENOME_INF_BP
         for run in CROSS for side in ("mat", "pat")}


def curves(grp, rows):
    ibp = infbp[grp][rows].sum(axis=0)
    safe = np.where(ibp > 0, ibp, 1.0)
    per_run, depth = {"mat": [], "pat": []}, np.zeros(NBIN)
    for run in CROSS:
        m = counts[grp][(run, "mat")][rows].sum(axis=0)
        p = counts[grp][(run, "pat")][rows].sum(axis=0)
        depth += m + p
        for side, r in (("mat", m), ("pat", p)):
            per_run[side].append(np.where(ibp > 0, r / safe, np.nan) / gdens[(run, side)])
    return (np.nanmean(per_run["mat"], axis=0), np.nanmean(per_run["pat"], axis=0),
            depth, ibp)


# ---------------------------------------------------------------- profiles
rows_out = []
for grp in GROUPS:
    idx = np.arange(len(genes_of[grp]))
    mo, po, depth, ibp = curves(grp, idx)
    bs = {k: np.empty((BOOT, NBIN)) for k in ("mat", "pat")}
    for k in range(BOOT):
        bm, bp_, _, _ = curves(grp, RNG.choice(idx, len(idx), replace=True))
        bs["mat"][k], bs["pat"][k] = bm, bp_
    keep = depth >= MIN_READS
    ci = {k: np.nanpercentile(v, [2.5, 97.5], axis=0) for k, v in bs.items()}
    for side, obs, tag in (("mat", mo, "maternal"), ("pat", po, "paternal")):
        y = np.where(keep, obs, np.nan)
        for b in range(NBIN):
            rows_out.append(dict(group=grp, scheme=SCH, subset="all", bin=b, allele=tag,
                                 enrichment=y[b], lo=ci[side][0][b], hi=ci[side][1][b],
                                 reads=int(depth[b]), inf_bp=int(ibp[b]), n_genes=len(idx)))
pd.DataFrame(rows_out).to_csv(D / "matpat_profiles_compact.csv", index=False)
print(f"wrote {D / 'matpat_profiles_compact.csv'}")

# ---------------------------------------------------------------- per-gene body values
per_gene = []
for grp in genes_of:
    for g, i in ((g, gidx[(grp, g)]) for g in genes_of[grp]):
        bp = infbp[grp][i, list(BODY)].sum()
        if bp <= 0:
            continue
        vals, depth = {"mat": [], "pat": []}, {run: 0.0 for run in CROSS}
        for run in CROSS:
            for side in ("mat", "pat"):
                r = counts[grp][(run, side)][i, list(BODY)].sum()
                depth[run] += r
                vals[side].append((r / bp) / gdens[(run, side)])
        if min(depth.values()) < MIN_PER_RUN:
            continue
        m, p = float(np.mean(vals["mat"])), float(np.mean(vals["pat"]))
        per_gene.append(dict(gene=g, group=grp, maternal=m, paternal=p,
                             log2_ratio=float(np.log2(m / p)) if p > 0 else np.nan,
                             reads=int(sum(depth.values())), inf_bp=int(bp)))
pg = pd.DataFrame(per_gene)
pg.to_csv(D / "matpat_per_gene_body_compact.csv", index=False)
print(f"wrote {D / 'matpat_per_gene_body_compact.csv'} ({len(pg):,} rows)")
for grp in GROUPS:
    d = pg[pg.group == grp]
    diff = d.maternal - d.paternal
    print(f"  {grp:32s} n {len(d):5,} of {len(genes_of[grp]):5,}  mat - pat median {diff.median():+.2f}"
          f"  (IQR {diff.quantile(.25):+.2f} to {diff.quantile(.75):+.2f})"
          f"  mat > pat {int((diff > 0).sum())}/{len(d)}")
