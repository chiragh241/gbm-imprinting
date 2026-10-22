"""Allele-resolved endosperm CG methylation in dme and fie mutants (Ibarra 2012,
Science, GSE38935): do DME and FIE act on gbM-PEGs?

Design (reciprocal Col x Ler endosperm, single library each; both mutants MATERNAL):
  WT  Col x Ler  (GSM952441 TFH44 ): maternal = cg-col     , paternal = cg-ler
  dme Col x Ler  (GSM952440 TFH605): maternal = cg-dmecol  , paternal = cg-ler
  WT  Ler x Col  (GSM952443 TFH52 ): maternal = cg-ler.merged, paternal = cg-col.merged
  fie Ler x Col  (GSM952442 TFH113): maternal = cg-fieler  , paternal = cg-col

DME (labelled Claim2 in the summary): maternal COL allele, WT vs dme (same Col mother).
  DME normally demethylates the maternal genome, so dme should RETAIN/GAIN gene-body mCG.
  The paternal Ler allele is the internal negative control (father is WT in both).
FIE (labelled Claim3): maternal LER allele, WT vs fie (same Ler mother). Prediction: gbM-PEG
  gene-body mCG is unchanged (FIE = PRC2 = H3K27me3, not a DNA-methylation writer).
  gbM-biparental genes are the discriminating control.

The group tests here use the 30 gbM-PEGs among all 135 PEGs of 09_mechanism_per_gene.csv.
Supplemental Figure S7f (dme; the fie panel is not in the paper), and the numbers quoted
with them, use the consensus PEG-36 instead and are
computed from 11_ibarra_per_gene.csv by 12_figures/09_fig5cde.py.

Per-cytosine GFF cols: chr . CG start end frac . . c=<meth>;t=<unmeth>;n=<reads>
Weighted gene-body mCG = sum(c) / sum(c+t) over CG sites in [start,end]; a gene needs
>= MIN_SITES covered CG on an allele to get a call (else NaN = "no data").
Coords: gene_targets.csv (TAIR10; identical assembly coords to the TAIR8 GFF).

Input   $GI_DATA/arabidopsis/gse38935/*.gff.gz   (01_download.sh)
        $GI_WORK/11_arabidopsis/gene_targets.csv, 09_mechanism_per_gene.csv, 01_imprinting_calls.csv
Outputs ($GI_WORK/11_arabidopsis/):
  11_ibarra_per_gene.csv       per-gene allele mCG for every sample/allele
  11_dme_gbmpeg_maternal.csv   the 12 gbM-PEGs with a Gehring 2011 DMR annotation: maternal
                               WT vs dme delta
  11_summary.csv               group-level medians + paired tests
  11_offset_corrected.csv      per-gene mutant deltas against the gbM-biparental background
"""
import gzip, re, warnings
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import wilcoxon, mannwhitneyu
warnings.filterwarnings("ignore")

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, SHIPPED, WORK

MDIR = DATA / "arabidopsis" / "gse38935"
OUT = WORK / "11_arabidopsis"
MIN_SITES = 5                      # covered CG sites in gene body to make a call
ATTR = re.compile(r"c=(\d+);t=(\d+)")

# ---------- gene coordinates + classes ----------
gt = pd.read_csv(OUT / "gene_targets.csv")
gt["chr"] = gt["chr"].astype(str)
gt = gt[gt["chr"].isin(list("12345"))]
coords = gt.set_index("AGI")[["chr", "start", "end"]]

pg = pd.read_csv(OUT / "09_mechanism_per_gene.csv")
imp = pd.read_csv(OUT / "01_imprinting_calls.csv")
imp["gene"] = imp["gene"].astype(str).str.upper().str.split(".").str[0]
# every gene we (or gene_targets) ever call imprinted -> keep out of the biparental control
imprinted = set(pg.gene) | set(imp.loc[imp.call.isin(["PEG", "MEG"]), "gene"]) \
            | set(gt.loc[gt.impr_call != "biparental", "AGI"])
SETS = {
    "gbM_PEG": set(pg[pg.meth == "gbM"].gene),
    "UM_PEG":  set(pg[pg.meth == "UM"].gene),
    # gene_targets' narrow call labels 24/30 of our union gbM-PEGs "biparental_gbM";
    # subtract every imprinted gene so the control is genuinely biparental gbM.
    "gbM_biparental": set(gt.loc[gt.impr_meth == "biparental_gbM", "AGI"]) - imprinted,
}
assert len(SETS["gbM_PEG"]) == 30, len(SETS["gbM_PEG"])

# per-chromosome sorted intervals for fast assignment
by_chr = {}
for c, sub in gt.sort_values("start").groupby("chr"):
    by_chr[c] = (sub["start"].to_numpy(), sub["end"].to_numpy(), sub["AGI"].to_numpy())

FILES = {
    "matWT_col":  "GSM952441_TFH44.all.cg-col.gff.gz",
    "patWT_ler":  "GSM952441_TFH44.all.cg-ler.gff.gz",
    "matDME_col": "GSM952440_TFH605.all.cg-dmecol.gff.gz",
    "patDME_ler": "GSM952440_TFH605.all.cg-ler.gff.gz",
    "matWT_ler":  "GSM952443_TFH52.all.cg-ler.merged.gff.gz",
    "patWT_col":  "GSM952443_TFH52.all.cg-col.merged.gff.gz",
    "matFIE_ler": "GSM952442_TFH113.all.cg-fieler.gff.gz",
    "patFIE_col": "GSM952442_TFH113.all.cg-col.gff.gz",
}

def accumulate(fname):
    """Return DataFrame AGI -> meth, cov, nsite for one allele GFF."""
    acc = {c: (np.zeros(len(by_chr[c][2])), np.zeros(len(by_chr[c][2])),
               np.zeros(len(by_chr[c][2]))) for c in by_chr}
    with gzip.open(MDIR / fname, "rt") as fh:
        for chunk in pd.read_csv(fh, sep="\t", header=None,
                                 names=["chrom", "s", "ctx", "pos", "e", "f", "st", "fr", "attr"],
                                 usecols=["chrom", "pos", "attr"], chunksize=1_000_000):
            chunk["chrom"] = chunk["chrom"].str.replace("chr", "", regex=False)
            ct = chunk["attr"].str.extract(ATTR).astype(float)
            chunk["c"] = ct[0]; chunk["t"] = ct[1]
            for c, cc in chunk.groupby("chrom"):
                if c not in by_chr:
                    continue
                starts, ends, _ = by_chr[c]
                pos = cc["pos"].to_numpy()
                idx = np.searchsorted(starts, pos, side="right") - 1
                ok = (idx >= 0) & (pos <= ends[np.clip(idx, 0, len(ends) - 1)])
                gi = idx[ok]
                m = cc["c"].to_numpy()[ok]; cov = m + cc["t"].to_numpy()[ok]
                am, at, an = acc[c]
                np.add.at(am, gi, m); np.add.at(at, gi, cov)
                np.add.at(an, gi, (cov > 0).astype(float))
    rows = []
    for c in by_chr:
        _, _, agis = by_chr[c]; am, at, an = acc[c]
        rows += list(zip(agis, am, at, an))
    d = pd.DataFrame(rows, columns=["AGI", "meth", "cov", "nsite"]).set_index("AGI")
    return d

print("accumulating 8 allele files ...")
data = {}
for k, fn in FILES.items():
    d = accumulate(fn)
    data[k] = d
    print(f"  {k:12} {fn:42} total CG reads={int(d['cov'].sum()):>10,}")

# ---------- assemble per-gene mCG (NaN if < MIN_SITES) ----------
def mcg(key):
    d = data[key]
    v = (d["meth"] / d["cov"].replace(0, np.nan)).where(d["nsite"] >= MIN_SITES)
    return v.rename(key)

per = pd.concat([mcg(k) for k in FILES], axis=1)
for k in FILES:
    per[k + "_nsite"] = data[k]["nsite"].reindex(per.index)
# label set membership
per["set"] = "other"
for name, s in SETS.items():
    per.loc[per.index.isin(s), "set"] = name
per.to_csv(OUT / "11_ibarra_per_gene.csv")

# ---------- helpers ----------
def paired(df, a, b):
    sub = df[[a, b]].dropna()
    if len(sub) < 3:
        return dict(n=len(sub), median_a=np.nan, median_b=np.nan, median_delta=np.nan,
                    n_up=np.nan, p=np.nan)
    delta = sub[b] - sub[a]
    try:
        p = wilcoxon(sub[a], sub[b]).pvalue
    except ValueError:
        p = np.nan
    return dict(n=len(sub), median_a=round(sub[a].median(), 4), median_b=round(sub[b].median(), 4),
                median_delta=round(delta.median(), 4), n_up=int((delta > 0).sum()), p=p)

rows = []
def add(label, group, a, b):
    df = per[per["set"] == group] if group in SETS else per
    r = paired(df, a, b)
    rows.append({"test": label, "group": group, "allele_WT": a, "allele_MUT": b, **r})

# --- Claim 2: DME, maternal Col allele (WT vs dme) ---
for grp in ["gbM_PEG", "UM_PEG", "gbM_biparental"]:
    add("Claim2_dme_MATERNAL_col", grp, "matWT_col", "matDME_col")
# paternal Ler negative control (father WT in both)
for grp in ["gbM_PEG", "UM_PEG", "gbM_biparental"]:
    add("Claim2_dme_PATERNAL_ler_ctrl", grp, "patWT_ler", "patDME_ler")

# --- Claim 3: FIE, maternal Ler allele (WT vs fie) ---
for grp in ["gbM_PEG", "UM_PEG", "gbM_biparental"]:
    add("Claim3_fie_MATERNAL_ler", grp, "matWT_ler", "matFIE_ler")
for grp in ["gbM_PEG", "UM_PEG", "gbM_biparental"]:
    add("Claim3_fie_PATERNAL_col_ctrl", grp, "patWT_col", "patFIE_col")

summary = pd.DataFrame(rows)
summary.to_csv(OUT / "11_summary.csv", index=False)

# ---------- offset-corrected: per-gene delta, gbM-PEG vs biparental background ----------
# single-library methylomes carry a global offset (esp. fie); comparing the gbM-PEG delta
# distribution to the genome-wide biparental delta controls for it.
def delta_vs_background(a, b, label):
    d = (per[b] - per[a])
    out = {}
    bg = d[per["set"] == "gbM_biparental"].dropna()
    for grp in ["gbM_PEG", "UM_PEG"]:
        gd = d[per["set"] == grp].dropna()
        p = mannwhitneyu(gd, bg, alternative="two-sided").pvalue if len(gd) >= 3 else np.nan
        out[grp] = dict(n=len(gd), median_delta=round(gd.median(), 4),
                        bg_n=len(bg), bg_median_delta=round(bg.median(), 4),
                        excess_over_bg=round(gd.median() - bg.median(), 4),
                        p_vs_bg=p)
    return pd.DataFrame(out).T.assign(contrast=label).reset_index(names="group")

offset = pd.concat([
    delta_vs_background("matWT_col", "matDME_col", "dme_maternal_col"),
    delta_vs_background("patWT_ler", "patDME_ler", "dme_paternal_ler_ctrl"),
    delta_vs_background("matWT_ler", "matFIE_ler", "fie_maternal_ler"),
    delta_vs_background("patWT_col", "patFIE_col", "fie_paternal_col_ctrl"),
], ignore_index=True)
offset.to_csv(OUT / "11_offset_corrected.csv", index=False)

# ---------- Claim 3: does the maternal<paternal asymmetry at gbM-PEGs survive in fie? ----------
def asymmetry(mat, pat, tag):
    sub = per[per["set"] == "gbM_PEG"][[mat, pat]].dropna()
    delta = sub[pat] - sub[mat]           # paternal minus maternal (demethylation gap)
    p = wilcoxon(sub[pat], sub[mat]).pvalue if len(sub) >= 3 else np.nan
    return dict(tag=tag, n=len(sub), maternal=round(sub[mat].median(), 4),
                paternal=round(sub[pat].median(), 4),
                pat_minus_mat=round(delta.median(), 4), p=p)

asym = pd.DataFrame([
    asymmetry("matWT_ler", "patWT_col", "gbM_PEG_WT (Ler x Col)"),
    asymmetry("matFIE_ler", "patFIE_col", "gbM_PEG_fie (Ler x Col)"),
])

# ---------- the 12 Gehring-annotated gbM-PEGs: quantitative maternal dme delta ----------
geh12 = pg[(pg.meth == "gbM") & (pg.demeth_DMR.notna())][["gene", "demeth_DMR", "dme_DMR"]].copy()
geh12 = geh12.set_index("gene")
tab = per.loc[per.index.isin(geh12.index),
              ["matWT_col", "matDME_col", "matWT_col_nsite", "matDME_col_nsite",
               "patWT_ler", "patDME_ler"]].copy()
tab = tab.join(geh12)
tab["mat_delta_dme"] = tab["matDME_col"] - tab["matWT_col"]
tab["pat_delta_dme"] = tab["patDME_ler"] - tab["patWT_ler"]
def call(r):
    if pd.isna(r["mat_delta_dme"]):
        return "no_data"
    return "gain" if r["mat_delta_dme"] > 0.02 else ("loss" if r["mat_delta_dme"] < -0.02 else "flat")
tab["dme_maternal_call"] = tab.apply(call, axis=1)
tab.sort_values("mat_delta_dme", ascending=False).to_csv(OUT / "11_dme_gbmpeg_maternal.csv")

# ---------- report ----------
pd.set_option("display.width", 200); pd.set_option("display.max_columns", 20)
print("\n===== SUMMARY (paired WT vs mutant, maternal allele) =====")
print(summary.to_string(index=False))
print("\n===== 12 Gehring gbM-PEGs: maternal Col allele, WT vs dme =====")
show = tab[["demeth_DMR", "dme_DMR", "matWT_col", "matDME_col", "mat_delta_dme",
            "matDME_col_nsite", "dme_maternal_call", "pat_delta_dme"]]
print(show.round(3).to_string())
print("\nmaternal dme call counts:", tab["dme_maternal_call"].value_counts().to_dict())
print("of the 7 previously-unflagged (dme_DMR False):")
sub7 = tab[tab["dme_DMR"] == False]
print("  calls:", sub7["dme_maternal_call"].value_counts().to_dict())

print("\n===== OFFSET-CORRECTED: per-gene delta vs biparental background =====")
print(offset[["contrast", "group", "n", "median_delta", "bg_median_delta",
              "excess_over_bg", "p_vs_bg"]].to_string(index=False))
print("\n===== Claim 3: gbM-PEG maternal<paternal asymmetry, WT vs fie =====")
print(asym.to_string(index=False))
