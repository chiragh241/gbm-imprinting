#!/usr/bin/env python
"""Per-gene endosperm-minus-tissue CG methylation deltas, the numbers quoted with Figure 2a-d.

For gbM-PEGs and gbM-biparental genes, the whole-gene CG methylation of each B73 gene in 15 DAP
endosperm is paired with the same gene in embryo, leaf and tassel, and the median per-gene
difference (endosperm minus the other tissue) is reported with its interquartile range and a
paired two-sided Wilcoxon signed-rank test. The shuffled control regions are paired the same
way, embryo against endosperm.

Genes and control regions need at least MIN_SITES covered CG sites and MIN_COV reads in both
tissues of a pair, the filter of 02_mutant_deltas.py.

The genes' values come from Supplemental Dataset 1; the control regions' from the per-region
hits files, CG cytosines (second base of the trinucleotide G) summed per region.

Input   $GI_RESULTS/tables/Supplemental_Dataset_1_per_gene.csv   (01_supplemental_dataset1.py)
        $GI_WORK/04_emseq_methylation/violin/{embryo,wt_15dap}_hits.tsv
                                                            (04_emseq_methylation/06_region_hits.sh)
Output  $GI_RESULTS/tables/tissue_deltas.csv
          comparison, gene_set, n, tissue and endosperm medians (%), delta and IQR (points), p
"""
import io
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB
import lib.stats_lib as S

V = WORK / "04_emseq_methylation" / "violin"
MIN_COV, MIN_SITES = 10, 5
SETS = {"gbM-PEG": ("PEG", "gbM"), "gbM-biparental": ("biparental", "gbM")}


def awk(program, path):
    """Run an awk aggregation over one hits file (~2.4 GB, ~45M rows) and parse its output."""
    out = subprocess.run(["bash", "-c", f"awk -F'\\t' '{program}' '{path}'"],
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"awk failed on {path}: {out.stderr[:400]}")
    return pd.read_csv(io.StringIO(out.stdout), sep="\t", header=None,
                       names=["id", "mc", "tot", "nsite"])


CTRL_AWK = ('$1=="shuffled_control" && substr($3,2,1)=="G" '
            '{mc[$2]+=$4; tot[$2]+=$5; n[$2]++} '
            'END{for(k in tot) print k"\\t"mc[k]"\\t"tot[k]"\\t"n[k]}')

m = pd.read_csv(TAB / "Supplemental_Dataset_1_per_gene.csv")

rows = []


def emit(comparison, name, s):
    rows.append(dict(comparison=comparison, gene_set=name, n=s["n"],
                     tissue_median_pct=round(s["median_a"] * 100, 1),
                     endosperm_median_pct=round(s["median_b"] * 100, 1),
                     delta_points=round(s["median_delta"] * 100, 1),
                     iqr_lo=round(s["iqr_lo"] * 100, 1), iqr_hi=round(s["iqr_hi"] * 100, 1),
                     wilcoxon_p=s["p_wilcoxon"]))
    print(S.line(f"{name}, {comparison}", s))


for other in ("embryo", "leaf", "tassel"):
    for name, (imp, meth) in SETS.items():
        d = m[(m.imp_class == imp) & (m.meth_class == meth)]
        ok = (d[f"mCG_{other}"].notna() & d["mCG_wt_15dap"].notna()
              & (d[f"cov_{other}"] >= MIN_COV) & (d["cov_wt_15dap"] >= MIN_COV)
              & (d[f"nCG_{other}"] >= MIN_SITES) & (d["nCG_wt_15dap"] >= MIN_SITES))
        emit(f"{other} to 15 DAP endosperm", name,
             S.paired(d.loc[ok, f"mCG_{other}"].to_numpy(), d.loc[ok, "mCG_wt_15dap"].to_numpy()))

ctrl = {}
for g in ("embryo", "wt_15dap"):
    t = awk(CTRL_AWK, V / f"{g}_hits.tsv")
    t = t[(t.tot >= MIN_COV) & (t.nsite >= MIN_SITES)]
    ctrl[g] = (t.mc / t.tot).set_axis(t.id).rename(g)
j = pd.concat([ctrl["embryo"], ctrl["wt_15dap"]], axis=1, join="inner")
emit("embryo to 15 DAP endosperm", "shuffled control",
     S.paired(j.embryo.to_numpy(), j.wt_15dap.to_numpy()))

TAB.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(TAB / "tissue_deltas.csv", index=False)
print(f"wrote {TAB / 'tissue_deltas.csv'}")
