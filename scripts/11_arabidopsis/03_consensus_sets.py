"""Build the canonical consensus imprinting sets.

A consensus imprinted gene is called by at least 2 of 3 published endosperm screens,
each contributing its own author list (no screen is re-called):
  - Pignatta 2014  full published list (all PEG/MEG rows of Fig 1 source data)
  - Gehring 2011   table S6 (Favored Parent = father/mother)
  - Wolff 2011     table S7 (PEG) / table S4 (MEG)

Input   $GI_WORK/11_arabidopsis/01_imprinting_calls.csv   (02_imprinting_calls.py)
        $GI_DATA/arabidopsis/supplementary/gehring2011_tableS6.xls, wolff2011_tableS7_PEG.xlsx,
        wolff2011_tableS4_MEG.xlsx   (01_download.sh)
Output  $GI_RESULTS/tables/consensus_pegs.csv   gene, n_studies, studies, call=PEG   (36)
        $GI_RESULTS/tables/consensus_megs.csv   gene, n_studies, studies, call=MEG   (70)
        (in Supplemental Dataset 3, 12_supplemental_dataset3.py)
        $GI_RESULTS/tables/listed_imprinted_genes.csv   gene, PEG, MEG (the studies listing
        it): every gene on any of the six lists (552), which every biparental set leaves out
"""
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, TAB

RAW = DATA / "arabidopsis" / "supplementary"
LEGACY_CALLS = WORK / "11_arabidopsis" / "01_imprinting_calls.csv"
SETS = TAB
SETS.mkdir(parents=True, exist_ok=True)


def agi(s):
    return s.astype(str).str.strip().str.upper().str.split(".").str[0]


# --- Pignatta 2014: full published lists (all crosses) ---
imp = pd.read_csv(LEGACY_CALLS)
imp["gene"] = agi(imp["gene"])
pig_peg = set(imp.loc[imp.call == "PEG", "gene"])
pig_meg = set(imp.loc[imp.call == "MEG", "gene"])

# --- Gehring 2011 ---
g = pd.read_excel(RAW / "gehring2011_tableS6.xls", header=1)
g["gene"] = agi(g["Gene ID"])
geh_peg = set(g.loc[g["Favored Parent"] == "father", "gene"])
geh_meg = set(g.loc[g["Favored Parent"] == "mother", "gene"])

# --- Wolff 2011 ---
wp = pd.read_excel(RAW / "wolff2011_tableS7_PEG.xlsx", sheet_name="Tabelle1", header=3)
wp["gene"] = agi(wp["Locus"])
wol_peg = {x for x in wp["gene"] if x.startswith("AT")}

wm = pd.read_excel(RAW / "wolff2011_tableS4_MEG.xlsx", sheet_name=0, header=5)
wm["gene"] = agi(wm["Locus"])
wol_meg = {x for x in wm["gene"] if x.startswith("AT")}


def consensus_table(sets_by_study, call):
    """sets_by_study: {study_name: set(genes)}. Return DataFrame of genes in >=2 studies."""
    support = defaultdict(list)
    for study, genes in sets_by_study.items():
        for gene in genes:
            support[gene].append(study)
    rows = [
        {"gene": gene, "n_studies": len(studies), "studies": ";".join(sorted(studies)), "call": call}
        for gene, studies in support.items()
        if len(studies) >= 2
    ]
    return pd.DataFrame(rows).sort_values(["n_studies", "gene"], ascending=[False, True]).reset_index(drop=True)


peg_studies = {"Pignatta2014": pig_peg, "Gehring2011": geh_peg, "Wolff2011": wol_peg}
meg_studies = {"Pignatta2014": pig_meg, "Gehring2011": geh_meg, "Wolff2011": wol_meg}

peg_tbl = consensus_table(peg_studies, "PEG")
meg_tbl = consensus_table(meg_studies, "MEG")

peg_tbl.to_csv(SETS / "consensus_pegs.csv", index=False)
meg_tbl.to_csv(SETS / "consensus_megs.csv", index=False)

# every gene on any screen's PEG or MEG list, which no biparental set may contain
listed = pd.DataFrame({"gene": sorted(set().union(*peg_studies.values(), *meg_studies.values()))})
for call, studies in (("PEG", peg_studies), ("MEG", meg_studies)):
    listed[call] = [";".join(s for s, genes in studies.items() if g in genes) for g in listed.gene]
listed.to_csv(SETS / "listed_imprinted_genes.csv", index=False)

print("Input list sizes:")
for lab, s in [("Pignatta PEG", pig_peg), ("Gehring PEG", geh_peg), ("Wolff PEG", wol_peg),
               ("Pignatta MEG", pig_meg), ("Gehring MEG", geh_meg), ("Wolff MEG", wol_meg)]:
    print(f"  {lab:14}: {len(s)}")
print(f"\nCONSENSUS PEG (>=2 of 3): {len(peg_tbl)}   (all-3: {(peg_tbl.n_studies==3).sum()})")
print(f"CONSENSUS MEG (>=2 of 3): {len(meg_tbl)}   (all-3: {(meg_tbl.n_studies==3).sum()})")
print(f"ON ANY LIST: {len(listed)}")
print(f"\nWrote {SETS/'consensus_pegs.csv'}, {SETS/'consensus_megs.csv'} and "
      f"{SETS/'listed_imprinted_genes.csv'}")
