#!/usr/bin/env python
"""Supplemental Dataset 1: one row per gene, everything measured about it.

Spine is master_table.csv (13,209 B73 v5 genes with their imprinting call, methylation
class, coordinates and W22 syntelogs). Adds two blocks that until now only existed
as per-set summaries:

  gene-body mCG per genotype   wt_15dap, embryo, mdr1_15dap, wt_mature, mdr1_mature,
                               dng102_mature, leaf, tassel
                               (from 04_emseq_methylation/06_region_hits.sh)
  allele-resolved endosperm    maternal and paternal mCG fractions, pooled over the six
                               reciprocal libraries (05_allele_methylation/11_gene_summary.py)

Coverage columns are kept alongside every measurement so a reader can see what each value
rests on rather than having to trust it. No coverage filter is applied here: this is the
dataset, not a figure, so thin genes are reported with their thinness visible.

A B73 gene with two W22 syntelogs has its cytosines from both summed into one value.

Input   $GI_WORK/02_gene_sets/master_table.csv
        $GI_WORK/04_emseq_methylation/violin/<genotype>_hits.tsv
        $GI_WORK/05_allele_methylation/gene_summary/12_per_gene_allele_methylation.csv
Output  $GI_RESULTS/tables/Supplemental_Dataset_1_per_gene.csv
"""
from pathlib import Path
import io
import subprocess
import sys
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK as GI_WORK, TAB

WORK = GI_WORK / "04_emseq_methylation" / "violin"
TAB.mkdir(parents=True, exist_ok=True)
OUT = TAB / "Supplemental_Dataset_1_per_gene.csv"

GENOTYPES = ["wt_15dap", "embryo", "mdr1_15dap", "wt_mature", "mdr1_mature", "dng102_mature",
             "leaf", "tassel"]

m = pd.read_csv(GI_WORK / "02_gene_sets" / "master_table.csv")
print(f"spine: {len(m):,} genes")

# ---- gene-body mCG per genotype -------------------------------------------------
# region ids in the hits files are the BED name column, "W22gene|B73gene"; the spine is
# keyed on the B73 id, so split and take the second half.
for g in GENOTYPES:
    f = WORK / f"{g}_hits.tsv"
    if not f.exists():
        print(f"  !! {g}: no hits file, column omitted")
        continue
    # Aggregate in awk: these files are ~2.4 GB / ~45M rows and pandas would need several GB
    # each. The B73 id is the second half of the "W22gene|B73gene" region name.
    cmd = ("awk -F'\\t' 'substr($3,2,1)==\"G\" && $1!=\"shuffled_control\" "
           "{split($2,a,\"|\"); k=a[2]; mc[k]+=$4; tot[k]+=$5; n[k]++} "
           "END{for(k in tot) print k\"\\t\"mc[k]\"\\t\"tot[k]\"\\t\"n[k]}' "
           f"'{f}'")
    out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"awk failed on {f}: {out.stderr[:400]}")
    agg = pd.read_csv(io.StringIO(out.stdout), sep="\t", header=None,
                      names=["b73_gene", f"mc_{g}", f"cov_{g}", f"nCG_{g}"])
    agg[f"mCG_{g}"] = agg[f"mc_{g}"] / agg[f"cov_{g}"]
    agg = agg.drop(columns=[f"mc_{g}"])
    m = m.merge(agg, on="b73_gene", how="left")
    print(f"  {g}: {agg[f'mCG_{g}'].notna().sum():,} genes")

# ---- allele-resolved endosperm --------------------------------------------------
ap = GI_WORK / "05_allele_methylation" / "gene_summary" / "12_per_gene_allele_methylation.csv"
if ap.exists():
    a = pd.read_csv(ap)
    a = a[a.context == "CG"]
    a = (a.groupby("gene")[["mat_met", "mat_un", "pat_met", "pat_un", "n_sites"]]
           .sum().reset_index())
    a["allele_mat_mCG"] = a.mat_met / (a.mat_met + a.mat_un)
    a["allele_pat_mCG"] = a.pat_met / (a.pat_met + a.pat_un)
    a["allele_mat_reads"] = a.mat_met + a.mat_un
    a["allele_pat_reads"] = a.pat_met + a.pat_un
    a = a.rename(columns={"gene": "w22_gene", "n_sites": "allele_n_sites"})[
        ["w22_gene", "allele_mat_mCG", "allele_pat_mCG",
         "allele_mat_reads", "allele_pat_reads", "allele_n_sites"]]
    # spine's w22_genes can list more than one W22 model per B73 gene; join on the first
    m["w22_gene"] = m.w22_genes.astype(str).str.split(";").str[0].str.split(",").str[0]
    m = m.merge(a, on="w22_gene", how="left")
    print(f"  allele-resolved: {m.allele_mat_mCG.notna().sum():,} genes")

m.to_csv(OUT, index=False)
print(f"\nwrote {OUT}  ({len(m):,} rows, {len(m.columns)} columns)")
print("\ncolumns:")
for c in m.columns:
    print(f"  {c}")

print("\n=== sanity: gbM-PEG vs gbM-biparental, endosperm vs embryo ===")
for imp in ("PEG", "biparental"):
    s = m[(m.imp_class == imp) & (m.meth_class == "gbM")]
    if "mCG_wt_15dap" in m.columns and "mCG_embryo" in m.columns:
        print(f"  {imp:<11} gbM  n {len(s):5d}  endosperm {s.mCG_wt_15dap.median():.4f}  "
              f"embryo {s.mCG_embryo.median():.4f}")
