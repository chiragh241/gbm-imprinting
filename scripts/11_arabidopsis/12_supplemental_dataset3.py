#!/usr/bin/env python
"""Supplemental Dataset 3: arabidopsis, one row per gene, everything measured about it.

Rows are the genes with a Bewick et al. 2016 methylation class (06_master_table.csv), plus the
consensus imprinted genes without one. Columns:

  AGI, chr, start, end, strand      TAIR10 gene (07_gene_targets.py; blank if not there)
  bewick_class, meth_class          Bewick et al.'s class, and as used here (CHG and CHH
                                    methylated combined as teM)
  pignatta2014, gehring2011,        PEG or MEG when the gene is on that screen's list
  wolff2011                         (03_consensus_sets.py reads the same lists)
  consensus                         PEG or MEG when on that list in at least two of the three
  gene_set                          the gene's set in the paper (e.g. gbM-PEG, UM-biparental):
                                    consensus PEGs and MEGs of a class, and biparental genes on
                                    none of the three screens' lists, so with the three columns
                                    above blank (10_biparental_class_counts.py); blank otherwise
  endo_mat_mCG, endo_mat_nCG,       7 DAP endosperm, reciprocal Col-0 x C24 (GSE280598): CG
  endo_pat_mCG, endo_pat_nCG        methylation of the maternal and of the paternal Col-0 allele
                                    over the gene body, and covered CG sites (08); Supplemental
                                    Fig. S7e uses genes with at least 10 on each allele
  wt_mat_mCG, wt_pat_mCG,           Ibarra et al. 2012 (GSE38935) endosperm, wild-type Col x Ler
  dme_mat_mCG, dme_pat_mCG          and dme-2 Col-gl x Ler: maternal (Col) and paternal (Ler)
  and their _nCG                    CG methylation and covered CG sites, given when an allele
                                    has at least 5 sites (09); Supplemental Fig. S7f

No coverage filter is applied beyond those of the cited scripts; the site counts are kept so a
reader can see what each value rests on.

Input   $GI_WORK/11_arabidopsis/06_master_table.csv, gene_targets.csv, 01_imprinting_calls.csv,
        10_paternal_allele_endosperm_mCG.csv, 11_ibarra_per_gene.csv
        $GI_DATA/arabidopsis/supplementary/gehring2011_tableS6.xls, wolff2011_tableS7_PEG.xlsx,
        wolff2011_tableS4_MEG.xlsx   (01_download.sh)
        $GI_RESULTS/tables/consensus_{pegs,megs}.csv   (03)
        $GI_RESULTS/tables/arabidopsis_imprinted_classes.csv (11),
        arabidopsis_biparental_class_counts.csv (10), fig5cde_source.csv (12_figures/09)   checks
Output  $GI_RESULTS/tables/Supplemental_Dataset_3_arabidopsis.csv
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, TAB

RES = WORK / "11_arabidopsis"
RAW = DATA / "arabidopsis" / "supplementary"
OUT = TAB / "Supplemental_Dataset_3_arabidopsis.csv"
CLASSES = ("UM", "gbM", "teM")


def agi(s):
    return s.astype(str).str.strip().str.upper().str.split(".").str[0]


# ---- the three screens' lists, read as 03_consensus_sets.py reads them ----------------
pig = pd.read_csv(RES / "01_imprinting_calls.csv")
pig["gene"] = agi(pig["gene"])
geh = pd.read_excel(RAW / "gehring2011_tableS6.xls", header=1)
geh["gene"] = agi(geh["Gene ID"])
wp = pd.read_excel(RAW / "wolff2011_tableS7_PEG.xlsx", sheet_name="Tabelle1", header=3)
wm = pd.read_excel(RAW / "wolff2011_tableS4_MEG.xlsx", sheet_name=0, header=5)
STUDIES = {
    "pignatta2014": {"PEG": set(pig.loc[pig.call == "PEG", "gene"]),
                     "MEG": set(pig.loc[pig.call == "MEG", "gene"])},
    "gehring2011": {"PEG": set(geh.loc[geh["Favored Parent"] == "father", "gene"]),
                    "MEG": set(geh.loc[geh["Favored Parent"] == "mother", "gene"])},
    "wolff2011": {"PEG": {x for x in agi(wp["Locus"]) if x.startswith("AT")},
                  "MEG": {x for x in agi(wm["Locus"]) if x.startswith("AT")}},
}
cons = {c: set(pd.read_csv(TAB / f"consensus_{c.lower()}s.csv").gene) for c in ("PEG", "MEG")}
for c in ("PEG", "MEG"):
    n = pd.Series([g for s in STUDIES.values() for g in s[c]]).value_counts()
    assert set(n[n >= 2].index) == cons[c], f"consensus {c} lists do not reproduce"

# ---- one row per gene ------------------------------------------------------------------
master = pd.read_csv(RES / "06_master_table.csv")
d = master[["AGI", "raw_class", "meth_class"]].rename(columns={"raw_class": "bewick_class"})
extra = sorted((cons["PEG"] | cons["MEG"]) - set(d.AGI))
d = pd.concat([d, pd.DataFrame({"AGI": extra})], ignore_index=True)
assert not d.AGI.duplicated().any()
gt = pd.read_csv(RES / "gene_targets.csv")[["AGI", "chr", "start", "end", "strand"]]
d = d.merge(gt, on="AGI", how="left")
d = d[["AGI", "chr", "start", "end", "strand", "bewick_class", "meth_class"]]
for study, lists in STUDIES.items():
    d[study] = ["/".join(c for c in ("PEG", "MEG") if g in lists[c]) for g in d.AGI]
d["consensus"] = ["/".join(c for c in ("PEG", "MEG") if g in cons[c]) for g in d.AGI]
assert not d.consensus.str.contains("/").any(), "a gene in both consensus sets"

# gene sets: consensus imprinted genes of a class; biparental genes on none of the three
# screens' lists, the rule of 10_biparental_class_counts.py and 12_figures/09_fig5cde.py
listed = d[[*STUDIES]].ne("").any(axis=1)
cls = d.meth_class.where(d.meth_class.isin(CLASSES))
d["gene_set"] = ""
imp = (d.consensus != "") & cls.notna()
d.loc[imp, "gene_set"] = cls[imp] + "-" + d.consensus[imp]
bip = cls.notna() & ~listed
d.loc[bip, "gene_set"] = cls[bip] + "-biparental"

pa = pd.read_csv(RES / "10_paternal_allele_endosperm_mCG.csv").rename(columns={
    "mat_mCG": "endo_mat_mCG", "mat_nCG": "endo_mat_nCG",
    "pat_mCG": "endo_pat_mCG", "pat_nCG": "endo_pat_nCG"})
d = d.merge(pa[["AGI", "endo_mat_mCG", "endo_mat_nCG", "endo_pat_mCG", "endo_pat_nCG"]],
            on="AGI", how="left")
IB = {"matWT_col": "wt_mat", "patWT_ler": "wt_pat", "matDME_col": "dme_mat", "patDME_ler": "dme_pat"}
ib = pd.read_csv(RES / "11_ibarra_per_gene.csv")
ren = {}
for src, dst in IB.items():
    ren[src], ren[src + "_nsite"] = f"{dst}_mCG", f"{dst}_nCG"
ib = ib[["AGI", *ren]].rename(columns=ren)
d = d.merge(ib[["AGI", "wt_mat_mCG", "wt_mat_nCG", "wt_pat_mCG", "wt_pat_nCG",
                "dme_mat_mCG", "dme_mat_nCG", "dme_pat_mCG", "dme_pat_nCG"]], on="AGI", how="left")

# ---- checks against the paper's tables ----------------------------------------------------
assert not d.AGI.duplicated().any() and len(d) == len(master) + len(extra)
ic =pd.read_csv(TAB / "arabidopsis_imprinted_classes.csv").set_index("set")
for c in ("PEG", "MEG"):
    s = d[d.consensus == c]
    assert len(s) == ic.loc[c, "n"] and s.meth_class.isin(CLASSES).sum() == ic.loc[c, "classified"], c
    for k in CLASSES:
        assert (d.gene_set == f"{k}-{c}").sum() == ic.loc[c, k], (c, k)
bc = pd.read_csv(TAB / "arabidopsis_biparental_class_counts.csv").iloc[0]
for k in CLASSES:
    assert (d.gene_set == f"{k}-biparental").sum() == bc[k], k
src = pd.read_csv(TAB / "fig5cde_source.csv")


def n_used(gene_set, mat, pat, min_sites):
    s = d[(d.gene_set == gene_set) & (d[f"{mat}_nCG"] >= min_sites) & (d[f"{pat}_nCG"] >= min_sites)
          & d[f"{mat}_mCG"].notna() & d[f"{pat}_mCG"].notna()]
    return len(s)


for _, r in src.iterrows():
    if r.panel == "c":
        assert n_used(r.cls, "endo_mat", "endo_pat", 10) == r.n, (r.panel, r.cls)
    elif r.panel == "d":
        g = "wt" if r.genotype == "maternal WT" else "dme"
        assert n_used(r.cls, f"{g}_mat", f"{g}_pat", 5) == r.n, (r.panel, r.genotype, r.cls)

for c in ["chr", "start", "end"] + [c for c in d.columns if c.endswith("_nCG")]:
    d[c] = d[c].astype("Int64")          # positions and counts, blank where there is no value
d.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(d):,} rows, {len(d.columns)} columns)")
print(d.gene_set.replace("", "(none)").value_counts().to_string())
print(f"on any screen's list: {listed.sum():,}; of them in a biparental set: "
      f"{(listed & d.gene_set.str.endswith('-biparental')).sum():,}")
