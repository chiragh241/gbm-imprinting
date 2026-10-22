"""Pignatta et al. 2014 endosperm imprinting calls (eLife 3:e03198, Figure 1 source data 2).

Input: $GI_DATA/arabidopsis/supplementary/pignatta2014_endo.xlsx   (01_download.sh)
  - MEGs sheet: union of 286 MEG candidates with per-cross flags
    (Col-Ler MEG, Col-Cvi MEG, Ler-Cvi MEG) = "yes"/NaN
  - PEGs sheet: union of 104 PEG candidates with per-cross flags

Output ($GI_WORK/11_arabidopsis/):
  01_imprinting_calls.csv  per-gene call (PEG/MEG), n_crosses_support, support_tier
  01_high_conf_pegs.csv    PEGs called in at least 2 of the 3 crosses
  01_high_conf_megs.csv    MEGs called in at least 2 of the 3 crosses
"""
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

RAW = DATA / "arabidopsis" / "supplementary" / "pignatta2014_endo.xlsx"
OUT = WORK / "11_arabidopsis"
OUT.mkdir(parents=True, exist_ok=True)


def load_sheet(name, kind):
    df = pd.read_excel(RAW, sheet_name=name, header=2)
    cols = ["locus_name", f"Col-Ler {kind}", f"Col-Cvi {kind}", f"Ler-Cvi {kind}"]
    df = df[cols].copy()
    df.columns = ["gene", "ColLer", "ColCvi", "LerCvi"]
    for c in ["ColLer", "ColCvi", "LerCvi"]:
        df[c] = df[c].astype(str).str.lower().eq("yes").astype(int)
    df["n_crosses_support"] = df[["ColLer", "ColCvi", "LerCvi"]].sum(axis=1)
    df = df[df["gene"].astype(str).str.startswith("AT")].copy()
    df["call"] = kind
    return df


megs = load_sheet("MEGs", "MEG")
pegs = load_sheet("PEGs", "PEG")

print(f"Pignatta MEG candidates: {len(megs)}  (n_crosses dist: {megs.n_crosses_support.value_counts().to_dict()})")
print(f"Pignatta PEG candidates: {len(pegs)}  (n_crosses dist: {pegs.n_crosses_support.value_counts().to_dict()})")

both = pd.concat([megs, pegs], ignore_index=True)
both["support_tier"] = both["n_crosses_support"].map(
    {1: "single_cross", 2: "high_2of3", 3: "all_3"}
)

both.to_csv(OUT / "01_imprinting_calls.csv", index=False)
print(f"\nWrote {OUT/'01_imprinting_calls.csv'} ({len(both)} rows)")

high_pegs = pegs[pegs.n_crosses_support >= 2]
high_megs = megs[megs.n_crosses_support >= 2]
high_pegs.to_csv(OUT / "01_high_conf_pegs.csv", index=False)
high_megs.to_csv(OUT / "01_high_conf_megs.csv", index=False)
print(f"High-conf PEGs (>=2 crosses): {len(high_pegs)}")
print(f"High-conf MEGs (>=2 crosses): {len(high_megs)}")

dup = both[both.duplicated("gene", keep=False)]
if len(dup):
    print(f"\nWARN: {len(dup)} genes appear in both MEG and PEG sheets:")
    print(dup[["gene", "call", "n_crosses_support"]])
