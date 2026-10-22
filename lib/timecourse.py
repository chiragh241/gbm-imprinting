"""Maternal expression ratios of the B73/W22 time course, averaged over both cross directions.

One direction alone mixes parent of origin with the B73 versus W22 allele; the two directions'
per-gene ratios correlate at only r = 0.04 to 0.26 by stage, and the mean of the two cancels
the allele effect. Used by Figure 1a and Supplemental Figure S5.

The two directions are keyed on different gene IDs in B73_W22_maternalpreference_with_avg_rpm.csv:
  genome == "B73"   B73 x W22, B73 mother, feature = B73 gene, maternal allele = A
  genome == "W22"   W22 x B73, W22 mother, feature = W22 gene, maternal allele = B
Genes are paired through b73v5_to_w22.tsv (02_gene_sets/02_b73v5_to_w22.py), keeping strict
one-to-one pairs only: a ratio cannot be shared out over tandem copies. At each stage a gene
needs rpm_min in EACH direction, the calling pipeline's per-record floor.
"""
import numpy as np
import pandas as pd

from lib.paths import DATA, WORK

STAGES = (11, 14, 17, 21)
MATERNAL = {"B73": "A", "W22": "B"}


def _one_direction(raw, genome):
    """Long table of one direction: gene, stage, rpm, maternal ratio."""
    d = raw[raw.genome == genome]
    assert d.feature.is_unique, f"{genome} rows repeat a gene"
    out = []
    for st in STAGES:
        a, b = d[f"A_mean{st}"].to_numpy(), d[f"B_mean{st}"].to_numpy()
        tot = a + b
        mat = a if MATERNAL[genome] == "A" else b
        ratio = np.where(tot > 0, mat / np.where(tot > 0, tot, 1.0), np.nan)
        pref = d[f"pref_{st}"].to_numpy()
        ok = (tot > 0) & ~np.isnan(pref)
        diff = float(np.abs(ratio[ok] - pref[ok]).max())
        assert diff < 1e-9, f"{genome} rows, {st} DAP: maternal ratio off the file by {diff:.1e}"
        out.append(pd.DataFrame({"gene": d.feature.to_numpy(), "stage": st,
                                 "rpm": tot, "ratio": ratio}))
    return pd.concat(out, ignore_index=True)


def averaged_ratios(rpm_min=5):
    """Per gene and stage, both directions' maternal ratios and their mean.

    Returns (table, n_pairs): one row per B73 gene and stage passing rpm_min in both
    directions, and the number of one-to-one B73/W22 pairs present in both directions.
    """
    raw = pd.read_csv(DATA / "maize_tables" / "B73_W22_maternalpreference_with_avg_rpm.csv")
    bm, wm = _one_direction(raw, "B73"), _one_direction(raw, "W22")
    m = pd.read_csv(WORK / "02_gene_sets" / "b73v5_to_w22.tsv", sep="\t").dropna()
    m["w22"] = m.w22.str.split(",")
    m = m.explode("w22").drop_duplicates()
    m = m[m.b73_v5.isin(bm.gene) & m.w22.isin(wm.gene)]
    m = m[m.groupby("b73_v5").w22.transform("nunique").eq(1)
          & m.groupby("w22").b73_v5.transform("nunique").eq(1)]
    pairs = m.rename(columns={"b73_v5": "b73_gene", "w22": "w22_gene"})
    t = (pairs
         .merge(bm.rename(columns={"gene": "b73_gene", "rpm": "rpm_b73_mother",
                                   "ratio": "ratio_b73_mother"}), on="b73_gene")
         .merge(wm.rename(columns={"gene": "w22_gene", "rpm": "rpm_w22_mother",
                                   "ratio": "ratio_w22_mother"}), on=["w22_gene", "stage"]))
    t = t[(t.rpm_b73_mother >= rpm_min) & (t.rpm_w22_mother >= rpm_min)].copy()
    t["mat_ratio"] = (t.ratio_b73_mother + t.ratio_w22_mother) / 2
    return t.sort_values(["stage", "b73_gene"]).reset_index(drop=True), len(pairs)
