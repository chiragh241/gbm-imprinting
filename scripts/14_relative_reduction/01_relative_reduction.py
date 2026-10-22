#!/usr/bin/env python
"""Relative maternal reduction in gbM-biparental CG methylation in maize, rice and arabidopsis:
1 - mean(maternal mCG) / mean(paternal mCG) over the genes of the set, the "on average" figure of
the Introduction and the three Results paragraphs. The per-gene median and mean of
1 - maternal/paternal are kept as columns for comparison; the per-gene mean is unstable where
paternal mCG is near zero (a maize gene scores below -3,000%), the ratio of means is not.

The genes are those of each paragraph's delta, with n and the median delta checked against the
tables the figures are drawn from:
  maize        Fig. 2e gbM-biparentals (12_figures/03_fig2.py's data section, which itself
               checks 05_allele_methylation's 12_class_summary.csv)
  rice         Supplemental Fig. S6e (rice_matpat_per_gene.csv, against rice_matpat_summary.csv)
  arabidopsis  Supplemental Fig. S7 (12_figures/09_fig5cde.py's genes, against fig5cde_source.csv)

Output  $GI_RESULTS/tables/relative_reduction.csv
"""
import contextlib
import importlib.util
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import TAB

FIGS = Path(__file__).resolve().parents[1] / "12_figures"
SET = "gbM-biparental"
genes = {}

# maize: the genes of Figure 2e, from the figure script's own data section
src = FIGS / "03_fig2.py"
head = src.read_text(encoding="utf-8").split(
    "# ---------------------------------------------------------------- a to d")[0]
mod = {"__file__": str(src), "__name__": "_fig2_data"}
with contextlib.redirect_stdout(io.StringIO()), matplotlib.rc_context():
    exec(compile(head, str(src), "exec"), mod)
d = mod["pooled"]
d = d[d["set"] == "biparental__gbM"]
genes["maize"] = (d.mat.to_numpy(), d.pat.to_numpy())

# rice: the genes of Supplemental Fig. S6e (at least 10 covered CG on both sides)
g = pd.read_csv(TAB / "rice_matpat_per_gene.csv")
g = g[(g.n_comparisons > 0) & (g["class"] == SET)]
s = pd.read_csv(TAB / "rice_matpat_summary.csv")
s = s[s.min_sites == 10].set_index("cls").loc[SET]
assert len(g) == s.n, (len(g), s.n)
assert abs(np.median(g.pat_minus_mat) * 100 - s.median_pat_minus_mat_points) < 5e-3
genes["rice"] = (g.mat_mCG.to_numpy(), g.pat_mCG.to_numpy())

# arabidopsis: the genes of panel c of 09_fig5cde.py (Supplemental Fig. S7)
spec = importlib.util.spec_from_file_location("_fig5", FIGS / "09_fig5cde.py")
f5 = importlib.util.module_from_spec(spec)
with matplotlib.rc_context():
    spec.loader.exec_module(f5)
a = pd.read_csv(f5.RES / "10_paternal_allele_endosperm_mCG.csv")
a["cls"] = a.AGI.map(f5.gene_classes())
a = a[(a.pat_nCG >= f5.MIN_SITES_A) & (a.mat_nCG >= f5.MIN_SITES_A)
      & a.pat_mCG.notna() & a.mat_mCG.notna() & (a.cls == SET)]
s = pd.read_csv(TAB / "fig5cde_source.csv")
s = s[s.panel == "c"].set_index("cls").loc[SET]
assert len(a) == s.n, (len(a), s.n)
assert abs(np.median(a.pat_mCG - a.mat_mCG) * 100 - s.median_pat_minus_mat_points) < 1e-9
genes["arabidopsis"] = (a.mat_mCG.to_numpy(), a.pat_mCG.to_numpy())

rows = []
for species, (mat, pat) in genes.items():
    with np.errstate(divide="ignore", invalid="ignore"):
        raw = 1 - mat / pat
    # the per-gene median counts a gene with both values zero as no difference and one with
    # paternal zero under maternal methylation as an unbounded gain, which ranks last; the
    # per-gene mean is over the genes with paternal mCG above zero
    r = np.where((mat == 0) & (pat == 0), 0.0, raw)
    rows.append(dict(species=species, gene_set=SET, n=len(mat),
                     mean_maternal_mCG=mat.mean(), mean_paternal_mCG=pat.mean(),
                     reduction_pct=100 * (1 - mat.mean() / pat.mean()),
                     per_gene_median_pct=100 * np.median(r),
                     per_gene_mean_pct=100 * raw[pat > 0].mean(),
                     n_paternal_zero=int((pat == 0).sum())))

out = pd.DataFrame(rows)
TAB.mkdir(parents=True, exist_ok=True)
out.to_csv(TAB / "relative_reduction.csv", index=False)
print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
print("\n1 - mean maternal / mean paternal: " + ", ".join(
    f"{x.species} {x.reduction_pct:.0f}%" for x in out.itertuples()))
