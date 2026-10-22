#!/usr/bin/env python
"""Maternal and paternal mCG metagene from the B73/W22 hybrid endosperm, parallel to Figure 2a-d
(a supplement no longer in the paper).

SAME DATA AS 2e. PRJNA759188, reciprocal B73 x W22 and W22 x B73 MATURE endosperm, three
libraries each way, mapped to W22 v2 with BS-Seeker2 and split by allele with cgmaptools asm
(05_allele_methylation). Only cytosines linked by reads to a B73/W22 SNV can be given an allele,
so coverage is a fraction of what 2a-d draws on and short sets are gappy.

SAME AXIS AS 2a-d. The fixed 100 bp bins of 04_emseq_methylation/07_metagene_bins.sh: 2 kb
upstream, 2 kb inward from the TSS, 2 kb inward from the TES, 2 kb downstream, 80 bins, 5' to 3'.

MATERNAL FOLLOWS THE CROSS. The reference is W22, so Allele1 (= Ref) is W22 and Allele2 is B73,
as in 05_allele_methylation/11_gene_summary.py:
    bxw  B73 x W22   mother B73   maternal = Allele2
    wxb  W22 x B73   mother W22   maternal = Allele1

ESTIMATOR. Per direction, counts are summed over replicate libraries and over the genes in a
set, and each bin is sum(mC) / sum(C) per allele. The two directions are then averaged with
EQUAL weight, which is what cancels strain and reference bias (bxw has about three times the
linked cytosines of wxb, so summing would let it dominate). 2e's per-gene values sum counts
over all six libraries instead; the two agree in sign and are not meant to agree to the digit.

WHY A 500 bp WINDOW. At 100 bp a gbM-PEG bin has a median of 9 linked cytosines in the W22 x
B73 direction and a UM-PEG bin 2.5, so the raw 100 bp trace is mostly gaps. Each point is
therefore the pooled count over a centred five-bin (500 bp) window, still plotted on 2a-d's
100 bp axis, never pooling across the mid-gene break. A point needs MIN_SITES linked cytosines
in EACH direction or it is left as a gap; UM-PEG, n = 34, stays largely gaps.
Shading is a 95% bootstrap interval over genes (1,000 resamples, multinomial weights).

Input   $GI_WORK/04_emseq_methylation/metagene/metagene_bins_fixed.bed   (04_emseq_methylation/07)
        $GI_WORK/05_allele_methylation/asm/<library>.CG.asm.tsv           (05_allele_methylation/10)
Output  $GI_RESULTS/figures/figS_matpat_metagene.{png,pdf}
        $GI_RESULTS/tables/figS_matpat_metagene_source_data.csv
"""
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D

BASE = 8
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
    "xtick.major.size": 2.0, "ytick.major.size": 2.0,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK, FIG, TAB
BINS = WORK / "04_emseq_methylation" / "metagene" / "metagene_bins_fixed.bed"
ASM = WORK / "05_allele_methylation" / "asm"
LIBS = {"bxw": ["bxw_endo_rep1", "bxw_endo_rep2", "bxw_endo_rep3"],
        "wxb": ["wxb_endo_rep1", "wxb_endo_rep2", "wxb_endo_rep3"]}
MATERNAL_ALLELE = {"bxw": "a2", "wxb": "a1"}

SETS = [("PEG__gbM", "gbM-PEG"), ("biparental__gbM", "gbM-bip"),
        ("PEG__UM", "UM-PEG"), ("biparental__UM", "UM-bip")]
NUP, NIN, NDOWN = 20, 20, 20
NBIN = NUP + 2 * NIN + NDOWN
MID = NUP + NIN
MIN_SITES, BOOT, WIN = 20, 1000, 5
MAT, PAT = "#c0392b", "#2c7fb8"
RNG = np.random.default_rng(0)
FIELDS = ("mat_met", "mat_tot", "pat_met", "pat_tot", "sites")


def load_asm(lib, direction):
    """One library's linked CGs as maternal / paternal counts, one row per cytosine."""
    df = pd.read_csv(ASM / f"{lib}.CG.asm.tsv", sep="\t", dtype={"Chr": str})
    n0 = len(df)
    bad = (df.Allele1 != df.Ref)
    df = df[~bad]
    counts = {}
    for a in ("Allele1", "Allele2"):
        sp = df[f"{a}_linked_C"].astype(str).str.split("-", n=1, expand=True)
        counts[a] = (pd.to_numeric(sp[0], errors="coerce"), pd.to_numeric(sp[1], errors="coerce"))
    df = df.assign(a1_met=counts["Allele1"][0], a1_un=counts["Allele1"][1],
                   a2_met=counts["Allele2"][0], a2_un=counts["Allele2"][1]).dropna(
        subset=["a1_met", "a1_un", "a2_met", "a2_un"])
    # asm writes one row per (SNV, linked C); a C linked to three SNVs would count three times
    df["_cov"] = df.a1_met + df.a1_un + df.a2_met + df.a2_un
    # pandas' default sort is not stable, so which row survives a coverage tie can differ by
    # platform: the paper's values are from WSL, and a Windows run differs by up to 0.0006. Kept
    # as run, so the output matches the paper's.
    df = df.sort_values("_cov", ascending=False).drop_duplicates(["Chr", "C_Pos"])
    m, p = (("a2", "a1") if MATERNAL_ALLELE[direction] == "a2" else ("a1", "a2"))
    out = pd.DataFrame({
        "chrom": df.Chr.str.replace("^chr", "", regex=True),
        "pos": df.C_Pos.astype(np.int64),
        "mat_met": df[f"{m}_met"], "mat_tot": df[f"{m}_met"] + df[f"{m}_un"],
        "pat_met": df[f"{p}_met"], "pat_tot": df[f"{p}_met"] + df[f"{p}_un"]})
    print(f"  {lib}: {n0:,} asm rows, {int(bad.sum()):,} with Allele1 != Ref dropped, "
          f"{len(out):,} unique linked CGs")
    return out


def bin_sums(cyt, bins):
    """Sum each count over the cytosines inside every bin. asm positions are 1-based and the
    BED is 0-based half-open, so a C at P sits in [start, end) when start < P <= end."""
    out = np.zeros((len(bins), len(FIELDS)))
    for chrom, b in bins.groupby("chrom"):
        c = cyt[cyt.chrom == chrom].sort_values("pos")
        if c.empty:
            continue
        p = c.pos.to_numpy()
        cs = np.vstack([np.zeros(len(FIELDS)),
                        np.cumsum(np.column_stack([c[f].to_numpy(float) if f != "sites"
                                                   else np.ones(len(c)) for f in FIELDS]),
                                  axis=0)])
        lo = np.searchsorted(p, b.start.to_numpy() + 1, "left")
        hi = np.searchsorted(p, b.end.to_numpy(), "right")
        out[b.index.to_numpy()] = cs[hi] - cs[lo]
    return out


def window(a):
    """Centred WIN-bin sums along the last axis, inside each side of the mid-gene break."""
    out = np.zeros_like(a)
    half = WIN // 2
    for s0, s1 in ((0, MID), (MID, NBIN)):
        seg = a[..., s0:s1]
        cs = np.concatenate([np.zeros(seg.shape[:-1] + (1,)), np.cumsum(seg, axis=-1)], axis=-1)
        idx = np.arange(s1 - s0)
        lo, hi = np.clip(idx - half, 0, s1 - s0), np.clip(idx + half + 1, 0, s1 - s0)
        out[..., s0:s1] = cs[..., hi] - cs[..., lo]
    return out


def main():
    bins = pd.read_csv(BINS, sep="\t", header=None,
                       names=["chrom", "start", "end", "set", "name", "bin"],
                       dtype={"chrom": str})
    bins = bins[bins.set.isin([s for s, _ in SETS])].reset_index(drop=True)
    bins["gene"] = bins.name.str.split("|").str[0]

    per_dir = {}
    for direction, libs in LIBS.items():
        tot = np.zeros((len(bins), len(FIELDS)))
        for lib in libs:
            tot += bin_sums(load_asm(lib, direction), bins)
        per_dir[direction] = tot

    rows, prof = [], {}
    for key, label in SETS:
        b = bins[bins.set == key]
        genes = sorted(b.gene.unique())
        gi = {g: i for i, g in enumerate(genes)}
        arr = {}
        for direction in LIBS:
            a = np.zeros((len(FIELDS), len(genes), NBIN))
            v = per_dir[direction][b.index.to_numpy()]
            for k in range(len(FIELDS)):
                np.add.at(a[k], (b.gene.map(gi).to_numpy(), b.bin.to_numpy()), v[:, k])
            arr[direction] = window(a)
        # point estimate, and the same through bootstrap weights over genes
        W = RNG.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=BOOT)
        est = {}
        for allele in ("mat", "pat"):
            pt, bt, ok = [], [], np.ones(NBIN, bool)
            for direction in LIBS:
                a = arr[direction]
                met, tot = a[FIELDS.index(f"{allele}_met")], a[FIELDS.index(f"{allele}_tot")]
                sites = a[FIELDS.index("sites")].sum(0)
                ok &= sites >= MIN_SITES
                with np.errstate(invalid="ignore", divide="ignore"):
                    pt.append(met.sum(0) / tot.sum(0))
                    bt.append((W @ met) / (W @ tot))
            mean_pt = np.where(ok, (pt[0] + pt[1]) / 2, np.nan)
            mean_bt = (bt[0] + bt[1]) / 2
            with warnings.catch_warnings():
                # a gap point has no finite resample; it is masked below either way
                warnings.simplefilter("ignore", RuntimeWarning)
                lo, hi = np.nanpercentile(mean_bt, [2.5, 97.5], axis=0)
            est[allele] = (mean_pt, np.where(ok, lo, np.nan), np.where(ok, hi, np.nan))
        s_both = arr["bxw"][FIELDS.index("sites")] + arr["wxb"][FIELDS.index("sites")]
        n_genes = (s_both > 0).sum(0)
        n_linked = int((s_both.sum(1) > 0).sum())
        for k in range(NBIN):
            rows.append({"set": label, "bin": k, "n_genes_with_sites_window": int(n_genes[k]),
                         "sites_bxw_window": int(arr["bxw"][FIELDS.index("sites")][:, k].sum()),
                         "sites_wxb_window": int(arr["wxb"][FIELDS.index("sites")][:, k].sum()),
                         "mat": est["mat"][0][k], "mat_lo": est["mat"][1][k],
                         "mat_hi": est["mat"][2][k], "pat": est["pat"][0][k],
                         "pat_lo": est["pat"][1][k], "pat_hi": est["pat"][2][k]})
        prof[label] = (est, n_linked)
        body = slice(NUP, NUP + 2 * NIN)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            bm, bp = np.nanmean(est["mat"][0][body]), np.nanmean(est["pat"][0][body])
        print(f"{label:8s} genes {len(genes):5,}, {n_linked:5,} with a linked CG  "
              f"body mat {bm:.3f} pat {bp:.3f}  points drawn "
              f"{int(np.isfinite(est['mat'][0]).sum())}/{NBIN}")
    src = pd.DataFrame(rows)
    src.to_csv(TAB / "figS_matpat_metagene_source_data.csv", index=False)
    draw(prof)


def draw(prof):
    # Same geometry as Figure 2a-d (03_fig2.py).
    W, L_LAB, R_PAD, PAD = 1.62, 0.34, 0.03, 0.03
    TITLE, AX, GAP, LEG = 0.14, 0.60, 0.05, 0.30
    h = PAD + 4 * (TITLE + AX) + 3 * GAP + LEG
    fig = plt.figure(figsize=(W, h))
    y = PAD
    for k, (label, (est, n)) in enumerate(prof.items()):
        y += TITLE
        ax = fig.add_axes([L_LAB / W, (h - y - AX) / h, (W - L_LAB - R_PAD) / W, AX / h])
        x = np.arange(NBIN)
        for allele, colour in (("mat", MAT), ("pat", PAT)):
            m, lo, hi = est[allele]
            for a, b in ((0, MID), (MID, NBIN)):
                ax.fill_between(x[a:b], lo[a:b], hi[a:b], color=colour, alpha=0.18, lw=0)
                ax.plot(x[a:b], m[a:b], color=colour, lw=0.9)
        ax.axvline(NUP - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(NBIN - NDOWN - 0.5, color="0.55", ls="--", lw=0.6)
        ax.axvline(MID - 0.5, color="0.75", ls=":", lw=0.6)
        ax.set_xticks([])
        ax.set_xlim(0, NBIN - 1)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_yticklabels(["0", "0.5", "1"])
        ax.set_ylabel("mCG", labelpad=1)
        ax.set_title(f"{label}, n = {n:,}", pad=1.5, loc="left")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        y += AX + (GAP if k < 3 else 0)
    fig.legend(handles=[Line2D([], [], color=MAT, lw=1.2, label="mat"),
                        Line2D([], [], color=PAT, lw=1.2, label="pat")],
               ncol=2, frameon=False, loc="upper left",
               bbox_to_anchor=(0.02, (h - y - 0.04) / h),
               handletextpad=0.35, columnspacing=0.8, handlelength=1.2)
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"figS_matpat_metagene.{ext}", dpi=600)
    print(f"wrote figS_matpat_metagene.png / .pdf at {W} x {h:.2f} in")
    plt.close(fig)


if __name__ == "__main__":
    main()
