#!/usr/bin/env python
"""Supplemental Dataset 2: rice, one row per gene, everything measured about it.

Rows are the MSU7 genes of the leaf methylation table, plus the persistent imprinted genes the
table lacks (MSU7 transposable element related, so no leaf data). Columns:

  gene, chrom, start, end, strand   MSU7 gene model
  tonosaki_imprinting               PEG or MEG when the gene is on Tonosaki et al.'s lists
                                    (Supplementary Data 3), with tonosaki_imprint_class, the
                                    stages it is imprinted at ("Persistent" is all three)
  msu7_is_te                        MSU7's transposable element flag (locus_brief_info)
  leaf_mCG, leaf_mCHG, leaf_mCHH    Nipponbare young leaf (SRR10763664), from the leaf table,
  leaf_CG_reads, leaf_CHG_reads     with the reads behind the CG and CHG values
  leaf_category                     UM, gbM, teM, ambiguous or "no leaf data" (06_gene_sets.py)
  gene_set                          the gene's set in the paper (e.g. gbM-PEG, UM-biparental):
                                    persistent PEGs and MEGs, and biparental genes on none of
                                    Tonosaki et al.'s lists, of a category; blank otherwise
  endo_mat_mCG, endo_pat_mCG,       5 DAP endosperm, maternal and paternal CG methylation as
  endo_pat_minus_mat,               means over the comparisons where both sides have at least 10
  endo_n_comparisons                covered CG sites (of four; 07_gene_matpat.py), for the genes
                                    of a set; Supplemental Fig. S6e plots pat - mat

No coverage filter is applied beyond those of the cited scripts; the read and comparison counts
are kept so a reader can see what each value rests on.

Input   data/rice/rice_leaf_methylation.csv (shipped)
        $GI_DATA/rice/supplementary/tonosaki2024_SupplementaryData3_imprint_classes.xlsx   (01)
        $GI_DATA/rice/ref/osa1_r7.locus_brief_info.txt.gz, osa1_r7.all_models.gff3.gz     (01)
        $GI_RESULTS/tables/rice_gene_sets.csv, rice_persistent_imprinted_genes.csv,
        rice_imprinting_counts.csv   (06_gene_sets.py)
        $GI_RESULTS/tables/rice_matpat_per_gene.csv   (07_gene_matpat.py)
Output  $GI_RESULTS/tables/Supplemental_Dataset_2_rice.csv
"""
import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, SHIPPED, TAB

SD3 = DATA / "rice" / "supplementary" / "tonosaki2024_SupplementaryData3_imprint_classes.xlsx"
BRIEF = DATA / "rice" / "ref" / "osa1_r7.locus_brief_info.txt.gz"
GFF = DATA / "rice" / "ref" / "osa1_r7.all_models.gff3.gz"
OUT = TAB / "Supplemental_Dataset_2_rice.csv"


def gene_coords():
    """MSU7 gene features, as 07_gene_matpat.py reads them."""
    rows = []
    with gzip.open(GFF, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) == 9 and f[2] == "gene":
                gid = f[8].split("ID=")[1].split(";")[0]
                rows.append((gid, f[0], int(f[3]), int(f[4]), f[6]))
    return pd.DataFrame(rows, columns=["gene", "chrom", "start", "end", "strand"])


leaf = pd.read_csv(SHIPPED / "rice" / "rice_leaf_methylation.csv")
per = pd.read_csv(TAB / "rice_persistent_imprinted_genes.csv")
sets = pd.read_csv(TAB / "rice_gene_sets.csv")
mp = pd.read_csv(TAB / "rice_matpat_per_gene.csv")

# ---- Tonosaki et al.'s lists, every stage ---------------------------------------
lists = []
for status, sheet in {"PEG": "PEGs", "MEG": "MEGs"}.items():
    t = pd.read_excel(SD3, sheet_name=sheet, header=None, skiprows=3)[[0, 1]]
    t.columns = ["gene", "tonosaki_imprint_class"]
    t = t[t.gene.astype(str).str.startswith("LOC_Os")]
    t["tonosaki_imprinting"] = status
    lists.append(t)
lists = pd.concat(lists, ignore_index=True)
lists["tonosaki_imprint_class"] = lists.tonosaki_imprint_class.astype(str).str.strip()
assert not lists.gene.duplicated().any(), "a gene on both of Tonosaki et al.'s lists"

# ---- one row per gene -------------------------------------------------------------
d = leaf[["gene"]].copy()
missing = per.loc[per.leaf_category == "no leaf data", ["gene"]]
d = pd.concat([d, missing], ignore_index=True)
assert not d.gene.duplicated().any()
d = d.merge(gene_coords(), on="gene", how="left")
d = d.merge(lists[["gene", "tonosaki_imprinting", "tonosaki_imprint_class"]], on="gene", how="left")
brief = pd.read_csv(BRIEF, sep="\t").groupby("locus").is_TE.first()
d["msu7_is_te"] = d.gene.map(brief)
lv = leaf.rename(columns={"CG_methy": "leaf_mCG", "CHG_methy": "leaf_mCHG", "CHH_methy": "leaf_mCHH",
                          "CG_totalC": "leaf_CG_reads", "CHG_totalC": "leaf_CHG_reads",
                          "epiallele": "leaf_category"})
d = d.merge(lv[["gene", "leaf_mCG", "leaf_mCHG", "leaf_mCHH", "leaf_CG_reads", "leaf_CHG_reads",
                "leaf_category"]], on="gene", how="left")
d["leaf_category"] = d.leaf_category.fillna("no leaf data")
d = d.merge(sets.rename(columns={"class": "gene_set"}), on="gene", how="left")
endo = mp[["gene", "mat_mCG", "pat_mCG", "pat_minus_mat", "n_comparisons"]].rename(columns={
    "mat_mCG": "endo_mat_mCG", "pat_mCG": "endo_pat_mCG", "pat_minus_mat": "endo_pat_minus_mat",
    "n_comparisons": "endo_n_comparisons"})
d = d.merge(endo, on="gene", how="left")

# ---- checks against the paper's tables ---------------------------------------------
# 06_gene_sets.py checked the leaf categories against the rules; these tie the rows to its tables
assert len(d) == len(leaf) + len(missing) and not d.gene.duplicated().any()
assert d.chrom.notna().all(), "a gene without an MSU7 gene model"
p = d[d.tonosaki_imprint_class == "Persistent"]
counts = pd.read_csv(TAB / "rice_imprinting_counts.csv").set_index("status")
for status in ("PEG", "MEG"):
    s = p[p.tonosaki_imprinting == status]
    c = counts.loc[status]
    assert len(s) == c.persistent, status
    assert (s.leaf_category == "no leaf data").sum() == c.te_without_leaf_data, status
    for cat in ("UM", "gbM", "teM"):
        assert (s.leaf_category == cat).sum() == c[cat], (status, cat)
pg = per.set_index("gene")
assert (p.set_index("gene").leaf_category.sort_index() == pg.leaf_category.sort_index()).all()
assert (d.dropna(subset=["gene_set"]).set_index("gene").gene_set.sort_index()
        == sets.set_index("gene")["class"].sort_index()).all()
assert d.endo_mat_mCG.notna().sum() == mp.mat_mCG.notna().sum()
bip = d[d.tonosaki_imprinting.isna() & d.leaf_category.isin(["UM", "gbM", "teM"])]
assert (bip.gene_set == bip.leaf_category + "-biparental").all()

for c in ("leaf_CG_reads", "leaf_CHG_reads", "endo_n_comparisons"):
    d[c] = d[c].astype("Int64")          # counts, blank where there is no value
d.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(d):,} rows, {len(d.columns)} columns)")
print(d.gene_set.value_counts().to_string())
print(f"with endosperm values: {d.endo_mat_mCG.notna().sum():,}")
