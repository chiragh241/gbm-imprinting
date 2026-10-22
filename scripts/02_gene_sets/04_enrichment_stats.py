#!/usr/bin/env python
"""Fisher tests and confidence intervals for the methylation-category enrichments.

A test statistic for the gbM enrichment in PEGs over MEGs, and for the teM enrichment
quoted in the abstract.

Denominators are CATEGORIZED genes only, gbM + UM + teM, excluding "undefined". That is
the basis the manuscript text uses: it is what makes the MEG denominator 17 and the
biparental denominator 6378 in B73. Using all genes instead would give 64.3% -> 41.2%
for the PEG gbM fraction, so the two bases are not interchangeable.

Counts come from counts_grid.csv (03_build_sets.py), not retyped.

Odds-ratio CIs use the Haldane-Anscombe correction (+0.5 to every cell), needed because
the MEG gbM cell is 1. The Fisher p-values are exact and uncorrected.

Input   $GI_WORK/02_gene_sets/counts_grid.csv
Output  $GI_RESULTS/tables/enrichment_stats.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

OUT = TAB
OUT.mkdir(parents=True, exist_ok=True)

CATEGORIES = ["gbM", "UM", "teM"]      # "undefined" deliberately excluded
Z = 1.96


def wilson(k, n):
    """Wilson score interval. Behaves sensibly at k=1, which normal-approx does not."""
    p = k / n
    d = 1 + Z * Z / n
    centre = (p + Z * Z / (2 * n)) / d
    half = Z * np.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return centre - half, centre + half


def or_ci(a, b, c, d):
    """Haldane-Anscombe corrected odds-ratio CI."""
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    lor = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return np.exp(lor - Z * se), np.exp(lor + Z * se)


grid = pd.read_csv(WORK / "02_gene_sets" / "counts_grid.csv", index_col=[0, 1])
grid.index.names = ["coord_system", "imprint_class"]

TESTS = [
    # (methylation category, group A, group B, what it is for)
    ("gbM", "PEG", "MEG",        "gbM enrichment in PEGs over MEGs (the claim in the text)"),
    ("gbM", "PEG", "biparental", "gbM enrichment in PEGs over background (not claimed)"),
    ("teM", "MEG", "biparental", "teM enrichment in MEGs over background (abstract)"),
    ("teM", "MEG", "PEG",        "teM enrichment in MEGs over PEGs (Results)"),
]

rows = []
for coord in ["B73_genes", "W22_pairs"]:
    sub = grid.loc[coord]
    n_cat = sub[CATEGORIES].sum(axis=1)          # categorized total per imprint class
    for cat, ga, gb in [(t[0], t[1], t[2]) for t in TESTS]:
        note = next(t[3] for t in TESTS if t[:3] == (cat, ga, gb))
        a, n1 = int(sub.loc[ga, cat]), int(n_cat[ga])
        c, n2 = int(sub.loc[gb, cat]), int(n_cat[gb])
        b, d = n1 - a, n2 - c
        orv, p = fisher_exact([[a, b], [c, d]])
        lo_or, hi_or = or_ci(a, b, c, d)
        lo_a, hi_a = wilson(a, n1)
        lo_b, hi_b = wilson(c, n2)
        rows.append({
            "coord_system": "B73 v5" if coord == "B73_genes" else "W22 v2",
            "methylation_category": cat,
            "group_a": ga, "a_count": a, "a_n_categorized": n1,
            "a_fraction": round(a / n1, 4),
            "a_ci95_low": round(lo_a, 4), "a_ci95_high": round(hi_a, 4),
            "group_b": gb, "b_count": c, "b_n_categorized": n2,
            "b_fraction": round(c / n2, 4),
            "b_ci95_low": round(lo_b, 4), "b_ci95_high": round(hi_b, 4),
            "fold_enrichment": round((a / n1) / (c / n2), 2),
            "odds_ratio": round(orv, 2),
            "or_ci95_low": round(lo_or, 2), "or_ci95_high": round(hi_or, 2),
            "fisher_p": f"{p:.3g}",
            "significant_at_0.05": "yes" if p < 0.05 else "no",
            "note": note,
        })

out = pd.DataFrame(rows)
out.to_csv(OUT / "enrichment_stats.csv", index=False)

pd.set_option("display.width", 200)
for coord in out.coord_system.unique():
    print(f"\n=== {coord}  (categorized genes only: gbM + UM + teM) ===")
    for _, r in out[out.coord_system == coord].iterrows():
        print(f"  {r.methylation_category:>3}  {r.group_a} {r.a_count}/{r.a_n_categorized} "
              f"({r.a_fraction:.1%})  vs  {r.group_b} {r.b_count}/{r.b_n_categorized} "
              f"({r.b_fraction:.1%})")
        print(f"       {r.fold_enrichment}-fold   OR {r.odds_ratio} "
              f"(95% CI {r.or_ci95_low} - {r.or_ci95_high})   Fisher p = {r.fisher_p}"
              f"   {'SIGNIFICANT' if r['significant_at_0.05'] == 'yes' else 'not significant'}")

print(f"\nwrote {OUT / 'enrichment_stats.csv'}")
