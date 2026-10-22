#!/usr/bin/env python
"""Rice gene sets: methylation categories, imprinted and biparental genes, and the counts and
enrichment test of the rice Results paragraph.

CATEGORIES, from Nipponbare young leaf (SRR10763664, Zhao et al. 2020), summarised per MSU7 gene
by Yirui Sun (data/rice/rice_leaf_methylation.csv: methylated and unmethylated read counts over
the gene's CG, CHG and CHH sites). The rules are Zeng et al. 2023's, as for maize, without a
coverage minimum except for teM:
    UM   mCG <= 0.05 and mCHG <= 0.05
    gbM  mCG >= 0.2 and mCHG <= 0.05
    teM  mCG >= 0.4 and mCHG >= 0.4, with at least 10 reads over the gene's CHG sites
    anything else is ambiguous (no category).
The rules are applied here and must give the table's own epiallele column for every gene.

IMPRINTED, Tonosaki et al. 2024 Supplementary Data 3: the persistent PEGs and MEGs (imprinted at
3, 5 and 7 DAP; 245 and 95). The leaf table covers MSU7 genes not annotated as transposable
element related, so the persistent genes that MSU7 flags as TE related (locus_brief_info is_TE)
have no category: 7 PEGs and 20 MEGs, leaving 238 and 75. BIPARENTAL: the genes of a category on
none of Tonosaki's imprinted lists at any stage (every locus of Supplementary Data 3).

TEST, gbM among PEGs against MEGs: two-sided Fisher's exact test on the genes with a category
(gbM against UM and teM).

Input   data/rice/rice_leaf_methylation.csv (shipped)
        $GI_DATA/rice/supplementary/tonosaki2024_SupplementaryData3_imprint_classes.xlsx   (01)
        $GI_DATA/rice/ref/osa1_r7.locus_brief_info.txt.gz                                 (01)
Output  $GI_RESULTS/tables/rice_gene_sets.csv                 gene, class (e.g. gbM-PEG)
        $GI_RESULTS/tables/rice_persistent_imprinted_genes.csv every persistent gene: status,
                                                             MSU7 TE flag, leaf category
        $GI_RESULTS/tables/rice_imprinting_counts.csv         counts per status, the test
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, SHIPPED, TAB

SD3 = DATA / "rice" / "supplementary" / "tonosaki2024_SupplementaryData3_imprint_classes.xlsx"
BRIEF = DATA / "rice" / "ref" / "osa1_r7.locus_brief_info.txt.gz"
CATS = ["UM", "gbM", "teM"]

# ---------------------------------------------------------------- categories
leaf = pd.read_csv(SHIPPED / "rice" / "rice_leaf_methylation.csv")
cg, chg = leaf.CG_methy, leaf.CHG_methy
leaf["category"] = np.select(
    [(cg <= 0.05) & (chg <= 0.05),
     (cg >= 0.2) & (chg <= 0.05),
     (cg >= 0.4) & (chg >= 0.4) & (leaf.CHG_totalC >= 10)],
    CATS, "ambiguous")
bad = leaf[leaf.category != leaf.epiallele]
if len(bad):
    sys.exit(f"the rules do not give the table's categories for {len(bad)} genes:\n{bad.head()}")
print(f"leaf table: {len(leaf):,} genes; " + ", ".join(
    f"{c} {n:,}" for c, n in leaf.category.value_counts().items()))
category = dict(zip(leaf.gene, leaf.category))

# ---------------------------------------------------------------- imprinted genes
sheets = {"PEG": "PEGs", "MEG": "MEGs"}
imprinted, persistent = set(), []
for status, sheet in sheets.items():
    t = pd.read_excel(SD3, sheet_name=sheet, header=None, skiprows=3)[[0, 1]]
    t.columns = ["gene", "imprint_class"]
    t = t[t.gene.astype(str).str.startswith("LOC_Os")]
    imprinted |= set(t.gene)
    p = t[t.imprint_class.astype(str).str.strip() == "Persistent"]
    persistent += [(g, status) for g in p.gene]
brief = pd.read_csv(BRIEF, sep="\t").groupby("locus").is_TE.first()
per = pd.DataFrame(persistent, columns=["gene", "status"]).sort_values("gene", ignore_index=True)
per["msu7_is_te"] = per.gene.map(brief)
per["leaf_category"] = per.gene.map(category).fillna("no leaf data")
missing = per[per.leaf_category == "no leaf data"]
if (missing.msu7_is_te != "Y").any():
    sys.exit(f"persistent genes without leaf data that MSU7 does not flag as TE:\n{missing}")
print(f"imprinted at any stage {len(imprinted):,}; persistent "
      + ", ".join(f"{s} {n}" for s, n in per.status.value_counts().items())
      + f"; without leaf data (all MSU7 TE related) "
      + ", ".join(f"{s} {n}" for s, n in missing.status.value_counts().items()))

# ---------------------------------------------------------------- gene sets
sets = {}
for g, status, cat in zip(per.gene, per.status, per.leaf_category):
    if cat in CATS:
        sets[g] = f"{cat}-{status}"
for g, cat in category.items():
    if g not in imprinted and cat in CATS:
        sets[g] = f"{cat}-biparental"
gs = pd.DataFrame({"gene": list(sets), "class": list(sets.values())})
print(gs["class"].value_counts().to_string())

# ---------------------------------------------------------------- counts and the test
rows = []
for status in ("PEG", "MEG"):
    d = per[per.status == status]
    n_cat = {c: int((d.leaf_category == c).sum()) for c in CATS}
    rows.append(dict(status=status, persistent=len(d),
                     te_without_leaf_data=int((d.leaf_category == "no leaf data").sum()),
                     with_leaf_data=int((d.leaf_category != "no leaf data").sum()),
                     ambiguous=int((d.leaf_category == "ambiguous").sum()),
                     with_category=sum(n_cat.values()), **n_cat))
counts = pd.DataFrame(rows).set_index("status")
peg, meg = counts.loc["PEG"], counts.loc["MEG"]
table = [[peg.gbM, peg.with_category - peg.gbM], [meg.gbM, meg.with_category - meg.gbM]]
odds, p = fisher_exact(table)
counts["gbM_fraction"] = counts.gbM / counts.with_category
counts["gbM_fold_PEG_over_MEG"] = peg.gbM / peg.with_category / (meg.gbM / meg.with_category)
counts["fisher_odds_ratio"] = odds
counts["fisher_p"] = p
print(counts.to_string())

TAB.mkdir(parents=True, exist_ok=True)
gs.to_csv(TAB / "rice_gene_sets.csv", index=False)
per.to_csv(TAB / "rice_persistent_imprinted_genes.csv", index=False)
counts.to_csv(TAB / "rice_imprinting_counts.csv")
