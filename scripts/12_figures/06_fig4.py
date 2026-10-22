#!/usr/bin/env python
"""Figure 4: allele-resolved H3K27me3 in endosperm, as ten images.

a-e  metagene profiles of maternal and paternal H3K27me3 enrichment, gbM-PEG, strongly
     demethylated gbM-biparental, all gbM-biparental, UM-PEG and UM-biparental. 15 kb flanks in
     1 kb bins; the gene body as 10 length-scaled bins, each drawn at the width of one flank
     bin, so the body is not to scale. Thick lines are a centred 5-bin rolling mean run within
     each flank and within the body, never across the TSS or polyA; thin lines unsmoothed.
     Bands are 95% bootstrap intervals over genes. Maternal red, paternal blue, no key. Linear axis 0 to 15 on all five; values above 15 are clipped.
     Each is as narrow as its TSS and polyA labels allow, and only a carries the y axis title
     and numbers, so the five fit one 6.3 in row.
f-j  per-gene maternal minus paternal body enrichment, one violin per set with its box, on
     -10 to 15 shared by all five. Above each: genes with maternal > paternal and the
     Holm-adjusted paired Wilcoxon p.

The 14 tests (a paired Wilcoxon within each of the five sets; for gbM-PEG against the
demethylated and against all gbM-biparentals, and UM-PEG against UM-biparentals, Mann-Whitney U
on the log2 maternal/paternal ratio, the maternal value and the paternal value) are
Holm-corrected together.

Input   $GI_WORK/10_chip_endosperm_allele/matpat_profiles_compact.csv,
        matpat_per_gene_body_compact.csv   (10_chip_endosperm_allele/07_compact_profile.py)
Output  $GI_RESULTS/figures/fig4i_{a..j}_*.{png,pdf}
        $GI_RESULTS/tables/fig4_pergene_tests.csv   the 14 tests, raw and Holm
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK, FIG, TAB

# ---------------------------------------------------------------- house style
FULL_WIDTH = 6.3
HALF_WIDTH = FULL_WIDTH / 2
BASE = 8
LETTER = 10
if "Arial" not in {f.name for f in fm.fontManager.ttflist}:
    raise SystemExit("Arial is not installed; matplotlib would silently fall back to "
                     "DejaVu Sans and the figure would not match the published one")
matplotlib.rcParams.update({
    "font.family": "Arial", "font.size": BASE,
    "axes.labelsize": BASE, "axes.titlesize": BASE,
    "xtick.labelsize": BASE, "ytick.labelsize": BASE, "legend.fontsize": BASE,
    "text.color": "black", "axes.labelcolor": "black",
    "xtick.color": "black", "ytick.color": "black",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

RES = WORK / "10_chip_endosperm_allele"

SCHEME, SUBSET = "wide", "all"
SMOOTH = 5
MAT, PAT = "#c0392b", "#2c7fb8"

PANELS = [("gbM-PEG", "gbM-PEG", "fig4h_1_gbm_peg"),
          ("gbM-biparental demethylated", "gbM-biparental (demethylated)",
           "fig4h_2_gbm_bip_demeth"),
          ("gbM-biparental (all)", "gbM-biparental (all)", "fig4h_3_gbm_bip_all"),
          ("UM-PEG", "UM-PEG", "fig4h_4_um_peg"),
          ("UM-biparental", "UM-biparental", "fig4h_5_um_biparental")]
CONTRASTS = [("gbM-PEG", "gbM-biparental demethylated"),
             ("gbM-PEG", "gbM-biparental (all)"),
             ("UM-PEG", "UM-biparental")]

def fmt_p(p, digits=0):
    """Plain decimals where a reader can hold the number in their head, scientific below 0.01.

    A fixed "%.0e" turned a Holm-saturated p into "1e+00", which reads as a typo rather than as
    "not significant at all".
    """
    return f"{p:.2f}" if p >= 0.01 else f"{p:.{digits}e}".replace("e-0", "e-")


def roll(y, k):
    """Centred rolling mean that keeps NaN bins NaN and does not run off the ends."""
    return pd.Series(y).rolling(k, center=True, min_periods=1).mean().to_numpy()


def holm(p):
    """Holm step-down adjusted p-values, in the input order."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    adj, running = np.empty(len(p)), 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(p) - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def build_tests(gene):
    """Five paired within-group tests plus three contrasts x three metrics; Holm over all 14.

    The family is the fourteen tests behind this figure, whichever way the panels are
    arranged, so re-sizing the images changes no p-value.
    """
    by = {g: d for g, d in gene.groupby("group")}
    rows, pvals = [], []
    for grp, _, _ in PANELS:
        g = by[grp]
        m, p_ = g["maternal"].to_numpy(), g["paternal"].to_numpy()
        lr = g["log2_ratio"].to_numpy(dtype=float)
        pw = float(wilcoxon(m, p_, alternative="two-sided").pvalue)
        rows.append(dict(test="Wilcoxon, paired within gene", group_a=grp, group_b="",
                         metric="maternal vs paternal", n_a=len(g), n_b=0,
                         median_a=float(np.median(m)), median_b=float(np.median(p_)),
                         median_log2_ratio=float(np.nanmedian(lr)),
                         n_a_maternal_higher=int((m > p_).sum()), p_raw=pw))
        pvals.append(pw)
    for metric in ("log2_ratio", "maternal", "paternal"):
        for a, b in CONTRASTS:
            x = by[a][metric].to_numpy(dtype=float)
            y = by[b][metric].to_numpy(dtype=float)
            x, y = x[~np.isnan(x)], y[~np.isnan(y)]
            pm = float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
            rows.append(dict(test="Mann-Whitney, PEG vs matched control", group_a=a, group_b=b,
                             metric=metric, n_a=len(x), n_b=len(y),
                             median_a=float(np.median(x)), median_b=float(np.median(y)),
                             median_log2_ratio=np.nan, n_a_maternal_higher="", p_raw=pm))
            pvals.append(pm)
    t = pd.DataFrame(rows)
    t["p_holm"] = holm(pvals)
    return t



# Linear y axis from 0. ONE CAP, 15, ON ALL FIVE PROFILES, so the panels are comparable by
# eye; the top of UM-PEG's maternal band is clipped, and how much is reported by
# profile_image() on every run.
YLIM = (0, 15)
TICKS = {(0, 15): [0, 5, 10, 15]}
# 15 kb flanks and the compact gene body (10_chip_endosperm_allele/05_compact_bins.py)
FLANK_KB, NUP, NBODY, NBIN = 15, 20, 10, 50
FIRST_BIN = NUP - FLANK_KB
DIFF_YLIM, DIFF_TICKS = (-10, 15), [-10, -5, 0, 5, 10, 15]


def ylim_for(grp):
    return YLIM


def segment_roll(v):
    """The centred rolling mean, run within each flank and within the 10 body bins.

    Never across the TSS or polyA, where it would mix body and flank values. A 5-bin window is
    5 kb in a flank but half the gene in the compact body, so the body's shape is flattened:
    the thin lines keep every per-bin value."""
    out = np.array(v, dtype=float)
    for seg in (slice(0, NUP), slice(NUP, NUP + NBODY), slice(NUP + NBODY, None)):
        out[seg] = roll(out[seg], SMOOTH)
    return out

# Abbreviated set names: "gbM-biparental (demethylated)" is 1.52 in against a 0.94 in violin.
SHORT = {"gbM-PEG": "gbM-PEG",
         "gbM-biparental (demethylated)": "gbM-bip demeth",
         "gbM-biparental (all)": "gbM-bip all",
         "UM-PEG": "UM-PEG",
         "UM-biparental": "UM-bip"}

LEFT, VIOL_W = 0.50, 0.94
# The five profiles go in one row: a carries the y title and numbers, b to e only the tick
# marks. At 8 pt the TSS and polyA labels, a quarter of the axis apart, touch below about
# 1.0 in of axis; 1.09 in keeps them about one space (0.03 in) apart, checked in save().
PROF_LEFT_A, PROF_LEFT, PROF_W, PROF_RPAD = 0.38, 0.05, 1.09, 0.02
MIN_LABEL_GAP = 0.03
ROW_W = PROF_LEFT_A + 4 * PROF_LEFT + 5 * (PROF_W + PROF_RPAD)
assert ROW_W <= FULL_WIDTH, f"five profiles need {ROW_W:.2f} in, the page is {FULL_WIDTH}"
# TITLE_H holds the tallest title, which is the violins' three lines. The profiles use two, so
# they leave a gap above it; that is deliberate, because it keeps all ten images exactly the
# same height with their axes at the same y, so they tile.
AX_H, TITLE_H, XTICK_H, PAD = 1.35, 0.50, 0.32, 0.05
VIOL_W_TOTAL = LEFT + VIOL_W + 0.04
IMG_H = PAD + XTICK_H + AX_H + TITLE_H


def save(fig, name):
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    w, h = fig.get_size_inches()
    for ax in fig.axes:
        arts = (list(ax.get_xticklabels()) + list(ax.get_yticklabels())
                + [ax.xaxis.label, ax.yaxis.label])
        for a in arts:
            if not a.get_text():
                continue
            bb = a.get_window_extent(r)
            if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                    or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
                raise SystemExit(
                    f"!!! {name}: {a.get_text()[:40]!r} runs off the {w:.2f} x {h:.2f} in "
                    f"sheet, x {bb.x0 / fig.dpi:.2f}-{bb.x1 / fig.dpi:.2f} "
                    f"y {bb.y0 / fig.dpi:.2f}-{bb.y1 / fig.dpi:.2f}")
        # Neighbouring x labels must stay apart: the profiles are as narrow as TSS and polyA
        # allow.
        boxes = sorted((a.get_window_extent(r) for a in ax.get_xticklabels() if a.get_text()),
                       key=lambda b: b.x0)
        for b0, b1 in zip(boxes, boxes[1:]):
            gap = (b1.x0 - b0.x1) / fig.dpi
            if gap < MIN_LABEL_GAP:
                raise SystemExit(f"!!! {name}: x labels {gap:.3f} in apart, under "
                                 f"{MIN_LABEL_GAP} in")
            print(f"   {name}: x labels {gap:.3f} in apart")
    # Figure-level text counts too. These titles are fig.text rather than ax.set_title, and
    # the first version of this check walked only the axes, so a title 0.36 in off the sheet
    # passed silently.
    for a in fig.texts:
        bb = a.get_window_extent(r)
        if (bb.x0 < -0.5 or bb.x1 / fig.dpi > w + 1e-6
                or bb.y0 < -0.5 or bb.y1 / fig.dpi > h + 1e-6):
            raise SystemExit(
                f"!!! {name}: title {a.get_text().splitlines()[0]!r} runs off the "
                f"{w:.2f} x {h:.2f} in sheet, x {bb.x0 / fig.dpi:.2f}-{bb.x1 / fig.dpi:.2f} "
                f"y {bb.y0 / fig.dpi:.2f}-{bb.y1 / fig.dpi:.2f}")
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", dpi=600)
    print(f"wrote {name}.png / .pdf at {w:.2f} x {h:.2f} in")
    plt.close(fig)


def linear_axis(ax, ylim):
    """Plain linear ticks, with the null at 1.0 still drawn by the caller."""
    ax.set_ylim(*ylim)
    ax.set_yticks(TICKS[ylim])


def profile_image(grp, label, prof, wil, stem, first):
    left = PROF_LEFT_A if first else PROF_LEFT
    fig = plt.figure(figsize=(left + PROF_W + PROF_RPAD, IMG_H))
    W, H = fig.get_size_inches()
    ax = fig.add_axes([left / W, (PAD + XTICK_H) / H, PROF_W / W, AX_H / H])
    p = prof[prof["group"] == grp]
    n_prof = int(p["n_genes"].iloc[0])
    for allele, col in (("maternal", MAT), ("paternal", PAT)):
        s = p[p["allele"] == allele].sort_values("bin")
        # smoothed over the full 20 kb, then trimmed, so the outermost kept bins are not
        # averaged over a truncated window
        y, lo, hi = s["enrichment"].to_numpy(), s["lo"].to_numpy(), s["hi"].to_numpy()
        ys, los, his = (segment_roll(v) for v in (y, lo, hi))
        keep = slice(FIRST_BIN, NBIN - FIRST_BIN)
        y, ys, los, his = y[keep], ys[keep], los[keep], his[keep]
        x = np.arange(len(y))
        ax.fill_between(x, los, his, color=col, alpha=0.16, lw=0)
        ax.plot(x, y, color=col, lw=0.6, alpha=0.35)
        ax.plot(x, ys, color=col, lw=1.6)
        over = [int((v > YLIM[1]).sum()) for v in (ys, y, his)]
        if any(over):
            print(f"   {label} {allele}: above {YLIM[1]} in {over[0]} thick-line, {over[1]} "
                  f"thin-line and {over[2]} band bins of {len(y)} (band peak {np.nanmax(his):.1f})")
    linear_axis(ax, ylim_for(grp))
    ax.set_xlim(0, len(x) - 1)
    ax.axhline(1.0, color="0.5", ls=":", lw=1)
    # on the bin edges: the body is bins 15 to 24 of the 40 drawn
    tss = NUP - FIRST_BIN - 0.5
    polya = tss + NBODY
    for b in (tss, polya):
        ax.axvline(b, color="0.45", ls="--", lw=1)
    ax.set_xticks([tss, polya])
    ax.set_xticklabels(["TSS", "polyA"])
    if first:
        ax.set_ylabel("H3K27me3 enrichment")
    else:
        ax.tick_params(axis="y", labelleft=False)     # read off a's axis in the row
    ax.grid(alpha=0.2, lw=0.6)
    fig.text(left / W, (PAD + XTICK_H + AX_H + 0.05) / H, f"{label}\nn = {n_prof:,}",
             va="bottom", ha="left", linespacing=1.3)
    save(fig, stem)


def violin_image(grp, label, gene, wil, stem):
    fig = plt.figure(figsize=(VIOL_W_TOTAL, IMG_H))
    W, H = fig.get_size_inches()
    ax = fig.add_axes([LEFT / W, (PAD + XTICK_H) / H, VIOL_W / W, AX_H / H])
    # One violin of per-gene maternal minus paternal, with a box inside it. Positive is a
    # maternal excess.
    g = gene[gene["group"] == grp]
    d = (g["maternal"] - g["paternal"]).to_numpy()
    parts = ax.violinplot([d], positions=[0], widths=0.78, showextrema=False)
    parts["bodies"][0].set_facecolor(MAT); parts["bodies"][0].set_alpha(0.35)
    parts["bodies"][0].set_edgecolor(MAT)
    bx = ax.boxplot([d], positions=[0], widths=0.2, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.4))
    bx["boxes"][0].set_facecolor(MAT); bx["boxes"][0].set_alpha(0.85)
    ax.axhline(0, color="0.5", ls=":", lw=1)
    ax.set_ylim(*DIFF_YLIM)
    ax.set_yticks(DIFF_TICKS)
    ax.set_yticklabels([f"{t:g}".replace("-", "−") for t in DIFF_TICKS])
    ax.set_xlim(-0.6, 0.6)
    ax.set_xticks([])
    ax.set_xlabel("per gene", labelpad=2)
    ax.set_ylabel("mat − pat enrichment")
    ax.grid(axis="y", alpha=0.2, lw=0.6)
    r = wil.loc[grp]
    # Three lines, not two: over 0.94 in of panel, "15/16 mat > pat, p = 3e-4" is 1.3 in and
    # ran off the sheet. The check in save() should have caught that and did not, because it
    # walked only the axes artists and these titles are fig.text.
    fig.text(LEFT / W, (PAD + XTICK_H + AX_H + 0.05) / H,
             f"{label}\n{int(r['n_a_maternal_higher'])}/{int(r['n_a'])} mat > pat\n"
             f"p = {fmt_p(r['p_holm'])}",
             va="bottom", ha="left", linespacing=1.3)
    save(fig, stem)


def main():
    prof = pd.read_csv(RES / "matpat_profiles_compact.csv")
    prof = prof[(prof["scheme"] == "compact") & (prof["subset"] == SUBSET)]
    gene = pd.read_csv(RES / "matpat_per_gene_body_compact.csv")
    keep = [g for g, _, _ in PANELS]
    missing = sorted(set(keep) - set(gene["group"]))
    if missing:
        raise SystemExit(f"missing group(s) {missing}: rerun 10_chip_endosperm_allele")
    gene = gene[gene["group"].isin(keep)]
    prof = prof[prof["group"].isin(keep)]
    tests = build_tests(gene)
    TAB.mkdir(parents=True, exist_ok=True)
    tests.to_csv(TAB / "fig4_pergene_tests.csv", index=False)
    wil = tests[tests["test"].str.startswith("Wilcoxon")].set_index("group_a")
    # Clipping is reported per set, because the cap is per set and a single total would hide
    # which panel is clipped.
    print(f"per-gene mat - pat outside {DIFF_YLIM[0]} to {DIFF_YLIM[1]}:")
    for grp, label, _ in PANELS:
        sub = gene[gene["group"] == grp]
        d = (sub["maternal"] - sub["paternal"]).dropna()
        out = d[(d < DIFF_YLIM[0]) | (d > DIFF_YLIM[1])]
        print(f"   {SHORT[label]:16s} {len(out):>3} of {len(d):>5} off the scale"
              + (f"  ({', '.join(f'{v:+.1f}' for v in out)})" if len(out) else "")
              + f"   median {d.median():+.2f}")
    FIG.mkdir(parents=True, exist_ok=True)
    for i, (grp, label, _stem) in enumerate(PANELS):
        key = _stem.replace("fig4h_", "").lstrip("12345_")
        short = SHORT[label]
        profile_image(grp, short, prof, wil, f"fig4i_{'abcde'[i]}_{key}_profile", i == 0)
        violin_image(grp, short, gene, wil, f"fig4i_{'fghij'[i]}_{key}_pergene")
    print("\nten images, no panel letters. a-e are the metagenes and f-j the per-gene "
          "violins, matching the (Fig. 4a-e) and (Fig. 4f-j) citations in the text.")
    print(f"the five profiles share {YLIM[0]} to {YLIM[1]} with {FLANK_KB} kb flanks and make a "
          f"row {ROW_W:.2f} in wide (a {PROF_LEFT_A + PROF_W + PROF_RPAD:.2f} in, b to e "
          f"{PROF_LEFT + PROF_W + PROF_RPAD:.2f} in each); the five per-gene panels share "
          f"{DIFF_YLIM[0]} to {DIFF_YLIM[1]}.")


if __name__ == "__main__":
    main()
