"""Allele-specific endosperm CG methylation for imprinted genes, from the paternal
vs maternal Col allele (Kirkbride/Gehring ROS1 dataset, GSE280598).

Answers: does the leaf-based Bewick gbM/UM label hold on the biologically-relevant
PATERNAL endosperm allele? And are gbM-PEGs maternally demethylated at the allele level?

Design (reciprocal Col x C24 F1 WT endosperm, 7 DAP, 3 reps each):
  paternal Col allele = C24xCol F1 Col_spiked  (Col is the father)  GSM8601912/13/14
  maternal Col allele = ColxC24 F1 Col_spiked  (Col is the mother)  GSM8601909/10/11
Using the Col_spiked files keeps everything on TAIR10 coords (direct AGI mapping) and
compares the SAME Col sequence as paternal vs maternal (no ecotype confound).

Download (per sample, 3 contexts; we use CG):
  ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8601nnn/<GSM>/suppl/<GSM>_<cross>_<rep>_Col_spiked_CG.bed.gz
bed cols: chrom  start0  end1  unmeth  meth  ratio   (col5=methylated; confirmed genome mCG=0.31)

This script aggregates the beds (01_download.sh) to per-gene weighted mCG over gene bodies
(gene_targets.csv, 07_gene_targets.py), methylated and unmethylated reads summed over the three
replicates. Output:
  $GI_WORK/11_arabidopsis/10_paternal_allele_endosperm_mCG.csv   AGI, pat_mCG, pat_nCG, mat_mCG, mat_nCG
"""
import gzip, numpy as np, pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, SHIPPED, WORK

BED_DIR = DATA / "arabidopsis" / "gse280598"
GENES = WORK / "11_arabidopsis" / "gene_targets.csv"
OUT = WORK / "11_arabidopsis" / "10_paternal_allele_endosperm_mCG.csv"
PAT = ["pat_r1.CG.bed.gz", "pat_r2.CG.bed.gz", "pat_r3.CG.bed.gz"]   # C24xCol Col allele
MAT = ["mat_r1.CG.bed.gz", "mat_r2.CG.bed.gz", "mat_r3.CG.bed.gz"]   # ColxC24 Col allele

g = pd.read_csv(GENES)[["AGI", "chr", "start", "end", "strand"]]
g["chr"] = g["chr"].astype(str)
by_chr = {}
for c, sub in g.groupby("chr"):
    sub = sub.sort_values("start")
    by_chr[c] = (sub["start"].to_numpy(), sub["end"].to_numpy(), sub["AGI"].to_numpy())

def accumulate(files):
    acc = {c: (np.zeros(len(by_chr[c][2])), np.zeros(len(by_chr[c][2])),
               np.zeros(len(by_chr[c][2]))) for c in by_chr}
    for fn in files:
        with gzip.open(BED_DIR / fn, "rt") as fh:
            for chunk in pd.read_csv(fh, sep="\t", header=None,
                                     names=["chrom", "s0", "p", "um", "m", "r"],
                                     usecols=["chrom", "p", "um", "m"], chunksize=1_000_000):
                chunk["chrom"] = chunk["chrom"].astype(str).str.replace("Chr", "", regex=False)
                for c, cc in chunk.groupby("chrom"):
                    if c not in by_chr:
                        continue
                    starts, ends, _ = by_chr[c]
                    pos = cc["p"].to_numpy()
                    idx = np.searchsorted(starts, pos, side="right") - 1
                    ok = (idx >= 0) & (pos <= ends[np.clip(idx, 0, len(ends) - 1)])
                    gi = idx[ok]
                    mm = cc["m"].to_numpy()[ok]; cov = mm + cc["um"].to_numpy()[ok]
                    am, at, an = acc[c]
                    np.add.at(am, gi, mm); np.add.at(at, gi, cov); np.add.at(an, gi, (cov > 0).astype(float))
    rows = []
    for c in by_chr:
        _, _, agis = by_chr[c]; am, at, an = acc[c]
        rows += list(zip(agis, am, at, an))
    return pd.DataFrame(rows, columns=["AGI", "meth", "tot", "nsite"])

pat = accumulate(PAT); mat = accumulate(MAT)
df = pat.merge(mat, on="AGI", suffixes=("_pat", "_mat"))
df["pat_mCG"] = df.meth_pat / df.tot_pat.replace(0, np.nan)
df["mat_mCG"] = df.meth_mat / df.tot_mat.replace(0, np.nan)
df.rename(columns={"nsite_pat": "pat_nCG", "nsite_mat": "mat_nCG"})[
    ["AGI", "pat_mCG", "pat_nCG", "mat_mCG", "mat_nCG"]].to_csv(OUT, index=False)
print("wrote", OUT, len(df), "genes")
