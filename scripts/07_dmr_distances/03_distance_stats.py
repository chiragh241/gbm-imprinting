#!/usr/bin/env python
"""Distances from genes to their nearest DMR: per-gene table, per-class summary and tests.

Classes are Figure 2's gene sets plus MEGs and undefined-PEGs: gbM-PEG, gbM-biparental, UM-PEG,
UM-biparental, MEG, undefined-PEG. Each PEG class is compared with the biparental genes of its
own methylation category, because gbM and UM genes differ in their distance to DMRs regardless
of imprinting.

  Holm family (6 tests, corrected together): for gbM-PEG vs gbM-biparental, UM-PEG vs
  UM-biparental and gbM-PEG vs UM-PEG, absolute distance by two-sided Mann-Whitney U and the
  fraction of genes overlapping a DMR by two-sided Fisher's exact test.
  Within 2 kb (uncorrected, separate questions in the text): the fraction of genes whose nearest
  DMR is within 2,000 bp (overlap included), by two-sided Fisher's exact test, for UM-PEG vs
  UM-biparental, gbM-PEG vs gbM-biparental, gbM-PEG vs MEG and UM-PEG vs MEG.

Each gene is placed in one of five zones by where its nearest DMR lies: overlapping the gene,
within 2 kb upstream, within 2 kb downstream, more than 2 kb upstream, more than 2 kb downstream.

Input   $GI_WORK/07_dmr_distances/closest.tsv, gene_category.tsv   (01, 02)
        $GI_WORK/02_gene_sets/master_table.csv
Output  $GI_RESULTS/tables/fig_dmr_distance_source_data.csv   per-gene distance and zone
        $GI_RESULTS/tables/fig_dmr_distance_summary.csv       per-class medians and percentages
        $GI_RESULTS/tables/fig_dmr_distance_stats.csv         the 6 Holm-family tests
        $GI_RESULTS/tables/dmr_within_2kb_tests.csv           the within-2 kb tests
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, TAB

W = WORK / "07_dmr_distances"
TAB.mkdir(parents=True, exist_ok=True)
PROMOTER = 2000
MAIN = ["gbM-PEG", "gbM-biparental", "UM-PEG", "UM-biparental"]
EXTRA = ["MEG", "undefined-PEG"]
CONTRASTS = [("gbM-PEG", "gbM-biparental"), ("UM-PEG", "UM-biparental"), ("gbM-PEG", "UM-PEG")]
WITHIN_2KB = [("UM-PEG", "UM-biparental"), ("gbM-PEG", "gbM-biparental"),
              ("gbM-PEG", "MEG"), ("UM-PEG", "MEG")]


def load():
    cl = pd.read_csv(W / "closest.tsv", sep="\t", header=None,
                     names=["chrom", "g_start", "g_end", "b73_gene", "score", "strand",
                            "d_chr", "d_start", "d_end", "dist"])
    cat = pd.read_csv(W / "gene_category.tsv", sep="\t")
    d = cl.merge(cat.rename(columns={"category": "cat"}), on="b73_gene", how="left")
    d["dist"] = d["dist"].astype(float)
    m = pd.read_csv(WORK / "02_gene_sets" / "master_table.csv",
                    usecols=["b73_gene", "imp_class", "meth_class"])
    d = d.merge(m, on="b73_gene", how="left")

    def cls(r):
        ic, mc = r["imp_class"], r["meth_class"]
        if ic == "MEG":
            return "MEG"
        if ic == "PEG":
            return {"gbM": "gbM-PEG", "UM": "UM-PEG", "undefined": "undefined-PEG"}.get(mc)
        if ic == "biparental":
            return {"gbM": "gbM-biparental", "UM": "UM-biparental"}.get(mc)
        return None

    d["cls"] = d.apply(cls, axis=1)
    d = d.dropna(subset=["cls"]).copy()
    d["abs_dist"] = d["dist"].abs()
    d["zone"] = np.select(
        [d["dist"] == 0,
         (d["dist"] < 0) & (d["abs_dist"] <= PROMOTER),
         d["dist"] < 0,
         (d["dist"] > 0) & (d["abs_dist"] <= PROMOTER)],
        ["gene body (overlap)", "promoter (<=2 Kb up)", "upstream >2 Kb",
         "3' (<=2 Kb down)"],
        default="downstream >2 Kb")
    return d


def holm(pvals):
    """Holm-Bonferroni adjusted p-values, in the order given."""
    order = np.argsort(pvals)
    adj = np.empty(len(pvals))
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(pvals) - rank) * pvals[i])
        adj[i] = min(1.0, running)
    return adj


def summarise(d, classes):
    rows = []
    for c in classes:
        s = d[d["cls"] == c]
        rows.append(dict(
            gene_class=c, n=len(s),
            median_signed_bp=float(s["dist"].median()),
            median_abs_bp=float(s["abs_dist"].median()),
            q1_abs_bp=float(s["abs_dist"].quantile(0.25)),
            q3_abs_bp=float(s["abs_dist"].quantile(0.75)),
            pct_body_overlap=100 * float((s["dist"] == 0).mean()),
            pct_promoter=100 * float((s["zone"] == "promoter (<=2 Kb up)").mean()),
            pct_within_2kb=100 * float((s["abs_dist"] <= PROMOTER).mean())))
    return pd.DataFrame(rows)


def tests(d):
    rows = []
    for a, b in CONTRASTS:
        x, y = d.loc[d["cls"] == a], d.loc[d["cls"] == b]
        u, pu = mannwhitneyu(x["abs_dist"], y["abs_dist"], alternative="two-sided")
        rows.append(dict(test="Mann-Whitney, |distance|", group_a=a, group_b=b,
                         stat_a=float(x["abs_dist"].median()),
                         stat_b=float(y["abs_dist"].median()),
                         n_a=len(x), n_b=len(y), p_raw=float(pu)))
        tab = [[int((x["dist"] == 0).sum()), int((x["dist"] != 0).sum())],
               [int((y["dist"] == 0).sum()), int((y["dist"] != 0).sum())]]
        orr, pf = fisher_exact(tab)
        rows.append(dict(test="Fisher, DMR overlaps gene body", group_a=a, group_b=b,
                         stat_a=100 * float((x["dist"] == 0).mean()),
                         stat_b=100 * float((y["dist"] == 0).mean()),
                         n_a=len(x), n_b=len(y), p_raw=float(pf), odds_ratio=float(orr)))
    t = pd.DataFrame(rows)
    t["p_holm"] = holm(t["p_raw"].to_numpy())
    return t


def within_2kb(d):
    rows = []
    for a, b in WITHIN_2KB:
        x = d.loc[d["cls"] == a, "abs_dist"] <= PROMOTER
        y = d.loc[d["cls"] == b, "abs_dist"] <= PROMOTER
        orr, p = fisher_exact([[int(x.sum()), int((~x).sum())], [int(y.sum()), int((~y).sum())]])
        fa, fb = 100 * float(x.mean()), 100 * float(y.mean())
        rows.append(dict(group_a=a, group_b=b, n_a=len(x), n_b=len(y),
                         within_a=int(x.sum()), within_b=int(y.sum()),
                         pct_a=fa, pct_b=fb, fold=fa / fb, odds_ratio=float(orr), fisher_p=p))
    return pd.DataFrame(rows)


def main():
    d = load()
    all_keys = MAIN + EXTRA
    stats = tests(d)
    src = d.loc[d["cls"].isin(all_keys),
                ["b73_gene", "cls", "chrom", "g_start", "g_end", "strand",
                 "d_start", "d_end", "dist", "abs_dist", "zone"]]
    src = src.rename(columns={"cls": "gene_class"}).sort_values(["gene_class", "b73_gene"])
    src.to_csv(TAB / "fig_dmr_distance_source_data.csv", index=False)

    summ = summarise(d, all_keys).round({"median_signed_bp": 0, "median_abs_bp": 0,
                                         "q1_abs_bp": 0, "q3_abs_bp": 0,
                                         "pct_body_overlap": 2, "pct_promoter": 2,
                                         "pct_within_2kb": 2})
    summ.to_csv(TAB / "fig_dmr_distance_summary.csv", index=False)
    stats.to_csv(TAB / "fig_dmr_distance_stats.csv", index=False)
    w2 = within_2kb(d)
    w2.to_csv(TAB / "dmr_within_2kb_tests.csv", index=False)

    print("\n" + summ.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\n" + stats.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\n" + w2.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
