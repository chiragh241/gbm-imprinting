#!/usr/bin/env python
"""Pooled H3K27me3 metagene profiles and per-gene matrices, B73 leaf and developing ear.

Two normalisations, each used where its denominator is actually populated:

  profile   per bin, pooled over the genes in a set:
                (sum chip reads / chip library total) / (sum input reads / input library total)
            Bin width cancels because ChIP and input share the bin. Pooling is what makes
            the input denominator usable: the input libraries are ~0.4x coverage, so a
            single gene's 100 bp bin holds well under one input read and a per-gene per-bin
            ratio would be noise. Confidence bands are bootstrapped over genes, 1000
            resamples, so both gene-to-gene spread and count noise sit in the band.

  per-gene  per bin, for the heatmap:
                (chip reads in bin / bin bp) / (that gene's input density over its whole
                 window) x (input total / chip total)
            Input is taken at window level, not bin level, for the same coverage reason,
            so each gene still carries its own mappability correction.

Replicates are pooled within tissue before any ratio is taken (leaf 2 ChIP + 2 input,
ear 2 ChIP + 1 input).

Input   $GI_WORK/09_chip_leaf_ear/bincounts/*.bincounts.tsv   (07_metagene_count.sh)
        $GI_WORK/09_chip_leaf_ear/counts/*.total.txt           (04_count.sh)
        $GI_WORK/09_chip_leaf_ear/gene_meta_v4.tsv, metagene_bins_v4.bed   (06_metagene_bins.py)
Output  $GI_WORK/09_chip_leaf_ear/profiles.csv        set tissue scheme subset bin ratio lo hi n_genes
        $GI_WORK/09_chip_leaf_ear/per_gene_bins.csv   gene set tissue scheme bin enrichment
        $GI_WORK/09_chip_leaf_ear/domain_extent.csv   flank extent of the enrichment
"""
from pathlib import Path
from collections import defaultdict
import csv
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK

D = WORK / "09_chip_leaf_ear"
CNT = D / "bincounts"
TOT = D / "counts"

TISSUE = {"leaf": dict(chip=["leaf_H3K27me3_rep1", "leaf_H3K27me3_rep2"],
                       inp=["leaf_input_rep1", "leaf_input_rep2"]),
          "ear": dict(chip=["ear_H3K27me3_rep1", "ear_H3K27me3_rep2"],
                      inp=["ear_input_rep1"])}
SETS = ["gbM-PEG", "UM-PEG", "gbM-biparental", "UM-biparental", "UM-MEG", "teM-MEG"]
NBIN = 100
BOOT = 1000
RNG = np.random.default_rng(0)

# ---- gene metadata ------------------------------------------------------------------
meta = {r["gene"]: r for r in
        csv.DictReader(open(D / "gene_meta_v4.tsv"), delimiter="\t")}
genes_of = defaultdict(list)
for g, r in meta.items():
    genes_of[r["set"]].append(g)
gidx = {}
for s in SETS:
    genes_of[s].sort()
    for i, g in enumerate(genes_of[s]):
        gidx[g] = i

# ---- bin widths, needed only for the per-gene density --------------------------------
binbp = defaultdict(float)
for ln in open(D / "metagene_bins_v4.bed"):
    c, bs, be, st, gene, scheme, b = ln.rstrip("\n").split("\t")
    binbp[(gene, scheme, int(b))] += int(be) - int(bs)


def libtotal(lib):
    return int(open(TOT / (lib + ".total.txt")).read().split()[0])


# ---- counts, library-size normalised to reads per million ----------------------------
mat = {s: {sch: {} for sch in ("fig2", "wide")} for s in SETS}
for s in SETS:
    for sch in ("fig2", "wide"):
        for t in TISSUE:
            for kind in ("chip", "inp"):
                mat[s][sch][(t, kind)] = np.zeros((len(genes_of[s]), NBIN))

for t, libs in TISSUE.items():
    for kind in ("chip", "inp"):
        for lib in libs[kind]:
            tot = libtotal(lib)
            with open(CNT / (lib + ".bincounts.tsv")) as fh:
                for ln in fh:
                    st, gene, sch, b, n = ln.rstrip("\n").split("\t")
                    if st not in mat or gene not in gidx:
                        continue
                    mat[st][sch][(t, kind)][gidx[gene], int(b)] += int(n) / tot * 1e6


def pooled(chip, inp, rows):
    c = chip[rows].sum(axis=0)
    i = inp[rows].sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(i > 0, c / i, np.nan)


# ---- pooled profile, bootstrapped over genes -----------------------------------------
prof_rows = []
for s in SETS:
    for sch in ("fig2", "wide"):
        key = "nb_free_2k" if sch == "fig2" else "nb_free_5k"
        for t in TISSUE:
            chip = mat[s][sch][(t, "chip")]
            inp = mat[s][sch][(t, "inp")]
            for subset in ("all", "flank_clean"):
                rows = np.array([gidx[g] for g in genes_of[s]
                                 if subset == "all" or meta[g][key] == "1"])
                if len(rows) < 5:
                    continue
                obs = pooled(chip, inp, rows)
                bs = np.empty((BOOT, NBIN))
                for k in range(BOOT):
                    bs[k] = pooled(chip, inp, RNG.choice(rows, len(rows), replace=True))
                with np.errstate(invalid="ignore"):
                    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
                for b in range(NBIN):
                    prof_rows.append(dict(set=s, tissue=t, scheme=sch, subset=subset,
                                          estimator="pooled",
                                          bin=b, ratio=obs[b], lo=lo[b], hi=hi[b],
                                          n_genes=len(rows)))

# ---- per-gene matrix: the heatmap, and a second profile estimator ---------------------
# The pooled ratio above weights every gene by its coverage, so wide body bins from long
# genes dominate the body while fixed-width flank bins weight all genes alike. That could
# manufacture a step at the TSS or TES on its own. The median over per-gene enrichments
# weights each gene once and cannot. If the two estimators agree, neither is the artefact.
pg_rows = []
for s in SETS:
    for sch in ("fig2", "wide"):
        widths = np.array([[binbp.get((g, sch, b), 0.0) for b in range(NBIN)]
                           for g in genes_of[s]])
        win_bp = widths.sum(axis=1)
        key = "nb_free_2k" if sch == "fig2" else "nb_free_5k"
        for t in TISSUE:
            chip = mat[s][sch][(t, "chip")]
            inp = mat[s][sch][(t, "inp")]
            win_inp = inp.sum(axis=1)
            dens_inp = np.where(win_bp > 0, win_inp / np.maximum(win_bp, 1), np.nan)
            with np.errstate(divide="ignore", invalid="ignore"):
                dens_chip = np.where(widths > 0, chip / np.maximum(widths, 1), np.nan)
                enr = dens_chip / dens_inp[:, None]
            ok = np.isfinite(dens_inp) & (dens_inp > 0)
            for gi, g in enumerate(genes_of[s]):
                if not ok[gi]:
                    continue
                for b in range(NBIN):
                    if np.isfinite(enr[gi, b]):
                        pg_rows.append(dict(gene=g, set=s, tissue=t, scheme=sch,
                                            bin=b, enrichment=enr[gi, b]))
            for subset in ("all", "flank_clean"):
                rows = np.array([gidx[g] for g in genes_of[s]
                                 if ok[gidx[g]] and
                                 (subset == "all" or meta[g][key] == "1")])
                if len(rows) < 5:
                    continue
                with np.errstate(invalid="ignore"):
                    obs = np.nanmedian(enr[rows], axis=0)
                    bs = np.empty((BOOT, NBIN))
                    for k in range(BOOT):
                        bs[k] = np.nanmedian(enr[RNG.choice(rows, len(rows), replace=True)],
                                             axis=0)
                    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
                for b in range(NBIN):
                    prof_rows.append(dict(set=s, tissue=t, scheme=sch, subset=subset,
                                          estimator="median_gene",
                                          bin=b, ratio=obs[b], lo=lo[b], hi=hi[b],
                                          n_genes=len(rows)))

# ---- domain extent, the numbers the flank claim rests on ------------------------------
byk = defaultdict(lambda: np.full(NBIN, np.nan))
for r in prof_rows:
    if r["estimator"] != "pooled":
        continue
    byk[(r["set"], r["tissue"], r["scheme"], r["subset"])][r["bin"]] = r["ratio"]

extent_rows = []
for (s, t, sch, subset), v in byk.items():
    body = float(np.nanmedian(v[20:80]))
    flank_bp = 100 if sch == "fig2" else 250
    # How far past the gene the signal stays above halfway between the body level and 1x.
    # Only defined for a set that is enriched over its body at all: for a set sitting at
    # 0.94x the half-max threshold is 0.97x, which any flank clears for free.
    half = 1.0 + (body - 1.0) / 2.0
    up = dn = ""
    if body >= 1.5:
        up = dn = 0
        for b in range(19, -1, -1):
            if np.isfinite(v[b]) and v[b] >= half:
                up = (20 - b) * flank_bp
            else:
                break
        for b in range(80, NBIN):
            if np.isfinite(v[b]) and v[b] >= half:
                dn = (b - 79) * flank_bp
            else:
                break
    extent_rows.append(dict(set=s, tissue=t, scheme=sch, subset=subset,
                            n_genes=int([r["n_genes"] for r in prof_rows
                                         if r["estimator"] == "pooled" and
                                         (r["set"], r["tissue"], r["scheme"],
                                          r["subset"]) == (s, t, sch, subset)][0]),
                            body_median=body,
                            upstream_edge=float(v[0]), downstream_edge=float(v[NBIN - 1]),
                            half_max_upstream_bp=up, half_max_downstream_bp=dn))

R = D
for name, rows in (("profiles.csv", prof_rows), ("per_gene_bins.csv", pg_rows),
                   ("domain_extent.csv", extent_rows)):
    with open(R / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("wrote {}  ({:,} rows)".format(R / name, len(rows)))

print("\ndomain extent, scheme=wide, subset=flank_clean")
for r in sorted(extent_rows, key=lambda r: (r["tissue"], r["set"])):
    if r["scheme"] == "wide" and r["subset"] == "flank_clean":
        reach = ("not enriched" if r["half_max_upstream_bp"] == "" else
                 "half-max reach {} bp up / {} bp down".format(
                     r["half_max_upstream_bp"], r["half_max_downstream_bp"]))
        print("  {:16s} {:5s} n={:4d}  body {:6.2f}x  -5kb edge {:5.2f}x  +5kb edge {:5.2f}x"
              "  {}".format(r["set"], r["tissue"], r["n_genes"], r["body_median"],
                            r["upstream_edge"], r["downstream_edge"], reach))
