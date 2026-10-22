"""Per-gene Arabidopsis master table: methylation class and the Pignatta call.

One row per gene in the Col-0 leaf methylation classification of Bewick et al. 2016 (PNAS
113:9111, Dataset S1), its CHG- and CHH-methylated classes combined as teM, with the Pignatta
2014 call at the 2-of-3-cross threshold (impr_call; biparental otherwise).

Input   $GI_DATA/arabidopsis/supplementary/pnas.1604666113.sd01.xlsx   (01_download.sh)
        $GI_WORK/11_arabidopsis/01_imprinting_calls.csv
Output  $GI_WORK/11_arabidopsis/06_master_table.csv
"""
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

RES = WORK / "11_arabidopsis"
BEW = DATA / "arabidopsis" / "supplementary" / "pnas.1604666113.sd01.xlsx"

imp = pd.read_csv(RES / "01_imprinting_calls.csv")
meth = pd.read_excel(BEW, sheet_name="A. thaliana gbM genes - A. thal", header=1)
meth["AGI"] = meth["Gene"].str.split(".").str[0]
meth = meth.rename(columns={"Classification": "raw_class"})[["AGI", "raw_class"]]
meth["meth_class"] = meth["raw_class"].map(
    {"UM": "UM", "gbM": "gbM", "mCHH": "teM", "mCHG": "teM"}
)

hc = imp[imp.n_crosses_support >= 2][["gene", "call", "n_crosses_support"]].rename(columns={"gene": "AGI"})
master = meth.merge(hc.rename(columns={"call": "impr_call"}), on="AGI", how="left")
master["impr_call"] = master["impr_call"].fillna("biparental")
master["impr_meth"] = master["impr_call"] + "_" + master["meth_class"]
master.to_csv(RES / "06_master_table.csv", index=False)
print(f"Master table: {master.shape}")
print(f"Imprinting call counts: {master.impr_call.value_counts().to_dict()}")
