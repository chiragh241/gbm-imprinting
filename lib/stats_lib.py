"""Shared summary + test helpers for the manuscript number fills.

One place for the conventions, so every number in the Results section is built the same way:

  * every contrast is PAIRED WITHIN GENE (or within control region/window) and summarised as
    the median of the per-gene differences, with the interquartile range of those differences;
  * the paired test is a two-sided Wilcoxon signed-rank on the per-gene differences;
  * a comparison of two sets of differences is a two-sided Mann-Whitney U, reported with
    Cliff's delta and the probability of superiority, because the draft's claim is about
    OVERLAP ("the two distributions overlapped strongly"), which a p-value cannot express;
  * a comparison of two proportions is a two-sided Fisher exact test.
"""
import numpy as np
from scipy.stats import wilcoxon, mannwhitneyu, fisher_exact


def q(v):
    v = np.asarray(v, float)
    return np.percentile(v, [25, 50, 75])


def paired(a, b):
    """Summarise per-element differences b - a. Returns a dict, everything on the 0-1 scale."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = b - a
    q1, med, q3 = q(d)
    try:
        p = float(wilcoxon(a, b).pvalue)
    except ValueError:          # all differences zero
        p = float("nan")
    return dict(n=len(d), median_a=float(np.median(a)), median_b=float(np.median(b)),
                median_delta=float(med), iqr_lo=float(q1), iqr_hi=float(q3),
                mean_delta=float(np.mean(d)), p_wilcoxon=p,
                n_neg=int(np.sum(d < 0)), n_pos=int(np.sum(d > 0)), deltas=d)


def compare(d1, d2):
    """Two-sided MWU on two sets of differences, with overlap-aware effect sizes.

    prob_superiority is P(|d1| > |d2|) + 0.5 P(=) after orienting both to magnitude, i.e. the
    AUC. 0.5 = the distributions are interchangeable; 1.0 = no overlap at all.
    """
    d1, d2 = np.asarray(d1, float), np.asarray(d2, float)
    u = mannwhitneyu(d1, d2, alternative="two-sided")
    auc = float(u.statistic) / (len(d1) * len(d2))
    return dict(n1=len(d1), n2=len(d2), p_mwu=float(u.pvalue),
                prob_superiority=auc, cliffs_delta=2 * auc - 1)


def prop(hits, n):
    return dict(hits=int(hits), n=int(n), pct=100.0 * hits / n if n else float("nan"))


def prop_test(h1, n1, h2, n2):
    odds, p = fisher_exact([[h1, n1 - h1], [h2, n2 - h2]])
    return dict(odds_ratio=float(odds), p_fisher=float(p))


def pct(x, nd=1):
    return f"{100 * x:.{nd}f}%"


def line(label, s, unit="pp"):
    """One human-readable line for a paired summary."""
    return (f"{label:<34} n={s['n']:>6,}  {s['median_a']*100:5.1f}% -> {s['median_b']*100:5.1f}%"
            f"   delta {s['median_delta']*100:+6.1f} {unit}"
            f"  (IQR {s['iqr_lo']*100:+.1f} to {s['iqr_hi']*100:+.1f})"
            f"  p={s['p_wilcoxon']:.3g}")
