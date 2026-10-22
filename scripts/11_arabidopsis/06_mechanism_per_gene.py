"""Every PEG on any of the three published lists, with its methylation class and DMR flags.

PEG set = union of Pignatta 2014, Gehring 2011 and Wolff 2011 (135 genes); consensus PEGs (in at
least 2 of the 3) are the 36 of 03_consensus_sets.py. Per gene:
  meth, mCG          Bewick 2016 class and leaf mCG
  demeth_DMR, dme_DMR Gehring 2011 endosperm-vs-embryo and dme-vs-WT DMRs (Table S6)

Input   $GI_WORK/11_arabidopsis/01_imprinting_calls.csv
        $GI_DATA/arabidopsis/supplementary/gehring2011_tableS6.xls, wolff2011_tableS7_PEG.xlsx,
        pnas.1604666113.sd01.xlsx   (01_download.sh)
Output  $GI_WORK/11_arabidopsis/09_mechanism_per_gene.csv
"""
import pandas as pd, numpy as np, warnings
from pathlib import Path
from scipy.stats import fisher_exact, mannwhitneyu
warnings.filterwarnings("ignore")

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

RAW = DATA / "arabidopsis" / "supplementary"
OUT = WORK / "11_arabidopsis"

def agi(s):
    return s.astype(str).str.strip().str.upper().str.split(".").str[0]

# ---------- PEG sets ----------
imp = pd.read_csv(OUT / "01_imprinting_calls.csv")
imp["gene"] = agi(imp["gene"])
pig_union = set(imp.loc[imp.call == "PEG", "gene"])

g = pd.read_excel(RAW / "gehring2011_tableS6.xls", header=1)
g["gene"] = agi(g["Gene ID"])
g["demeth_DMR"] = g["Col-gl embryo-endosp DMR"].notna() | g["Ler embryo-endosp DMR"].notna()
g["dme_DMR"] = g["Col-gl dme endosp-endosp DMR"].notna() | g["Ler dme endosp-endosp DMR"].notna()
geh_peg = set(g.loc[g["Favored Parent"] == "father", "gene"])

w = pd.read_excel(RAW / "wolff2011_tableS7_PEG.xlsx", sheet_name="Tabelle1", header=3)
w["gene"] = agi(w["Locus"])
wol_peg = {x for x in w["gene"] if x.startswith("AT")}

# consensus (>=2 studies) and union
from collections import Counter
cnt = Counter()
for s in (pig_union, geh_peg, wol_peg):
    cnt.update(s)
consensus = {x for x, n in cnt.items() if n >= 2}
union = pig_union | geh_peg | wol_peg

# ---------- Bewick methylation ----------
b = pd.read_excel(RAW / "pnas.1604666113.sd01.xlsx",
                  sheet_name="A. thaliana gbM genes - A. thal", header=1)
b["gene"] = agi(b["Gene"])
b["meth"] = b["Classification"].map({"UM": "UM", "gbM": "gbM", "mCHH": "teM", "mCHG": "teM"})
b["mCG"] = b["Percent CG Methylation"]
meth = b.set_index("gene")[["meth", "mCG"]]

# ---------- assemble per-gene PEG table ----------
gdmr = g.set_index("gene")[["demeth_DMR", "dme_DMR"]]
gdmr = gdmr[~gdmr.index.duplicated()]

def build(peg_set, label):
    df = pd.DataFrame({"gene": sorted(peg_set)})
    df = df.join(meth, on="gene")
    df = df.join(gdmr, on="gene")           # NaN => not in Gehring's annotated table
    df["peg_set"] = label
    return df

per_gene = build(union, "union")
per_gene.to_csv(OUT / "09_mechanism_per_gene.csv", index=False)
