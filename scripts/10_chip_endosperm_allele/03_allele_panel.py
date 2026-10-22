#!/usr/bin/env python3
"""Per-gene maternal H3K27me3 fraction in endosperm, and the gene groups of Figure 4.

Genes need at least 20 allele-informative reads in each of the two reciprocal libraries
(51 gbM-PEGs, 16 UM-PEGs, 1,830 gbM-biparentals, 676 UM-biparentals). Strongly demethylated
gbM-biparentals are the lowest quartile of maternal minus paternal CG methylation, computed per
cross direction from the allele-resolved EM-seq when both alleles have at least 10 reads and
averaged over the directions available; "not demethylated" is the top quartile, and "gbM-
biparental (all)" contains both.

Reciprocal design: SRR748861 has a B73 mother, SRR748862 a Mo17 mother. Averaging the two
per-gene maternal fractions cancels reference and strain bias; their difference measures it.
Endosperm is 2 maternal : 1 paternal, so the null is ~0.667, not 0.5. The null used here is
the OBSERVED genome-wide maternal fraction, which absorbs dosage and residual reference bias
together rather than assuming either.

Input   $GI_WORK/10_chip_endosperm_allele/counts/<run>.{B73,Mo17}.txt   (02_align_split_count.sh)
        $GI_WORK/02_gene_sets/master_table.csv
        $GI_WORK/05_allele_methylation/gene_summary/12_per_gene_allele_methylation.csv
Output  $GI_WORK/10_chip_endosperm_allele/05_panel_source_data.csv   gene, group, fractions, mCG delta
"""
import csv, os, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import WORK
import lib.stats_lib as S

W = str(WORK / "10_chip_endosperm_allele")
OUT = W
TAG = ""
MASTER = str(WORK / "02_gene_sets" / "master_table.csv")
METH = str(WORK / "05_allele_methylation" / "gene_summary" / "12_per_gene_allele_methylation.csv")
MINREADS = 20   # allele-informative reads per gene per library
RUNS = {"SRR748861": "B73", "SRR748862": "Mo17"}   # value = the MOTHER


def read_fc(path):
    """featureCounts output -> {gene: count}. Col 0 = GeneID, last col = count."""
    d = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#") or ln.startswith("Geneid") or ln.startswith("GeneID"):
                continue
            f = ln.rstrip("\n").split("\t")
            d[f[0]] = int(f[-1])
    return d


# ---------------------------------------------------------------- H3K27me3
counts = {}
for run in RUNS:
    for lab in ("B73", "Mo17"):
        p = os.path.join(W, "counts", f"{run}.{lab}.txt")
        if not os.path.exists(p):
            sys.exit(f"missing {p} -- run 02_align_split_count.sh first")
        counts[(run, lab)] = read_fc(p)

genes = sorted(counts[("SRR748861", "B73")])   # sorted, so the output rows come in a fixed order

# empirical null: genome-wide maternal fraction, per library
null = {}
for run, mother in RUNS.items():
    father = "Mo17" if mother == "B73" else "B73"
    m = sum(counts[(run, mother)].values())
    p = sum(counts[(run, father)].values())
    null[run] = m / (m + p)
    print(f"{run}: mother {mother:<4} genome-wide maternal fraction {null[run]:.4f} "
          f"({m:,} / {m+p:,})   [2m:1p dosage expectation 0.667]")
NULL = float(np.mean(list(null.values())))
print(f"reciprocal-mean null = {NULL:.4f}\n")

k27 = {}
for g in genes:
    fs, ok = [], True
    for run, mother in RUNS.items():
        father = "Mo17" if mother == "B73" else "B73"
        m, p = counts[(run, mother)].get(g, 0), counts[(run, father)].get(g, 0)
        if m + p < MINREADS:
            ok = False
            break
        fs.append(m / (m + p))
    if ok:
        k27[g] = dict(f_bar=float(np.mean(fs)), f_861=fs[0], f_862=fs[1],
                      excess=float(np.mean(fs)) - NULL, strain_gap=abs(fs[0] - fs[1]))
print(f"genes with >= {MINREADS} allele-informative reads in BOTH libraries: {len(k27):,}\n")

# ---------------------------------------------------------------- classes + mCG delta
cls, w22map = {}, defaultdict(list)
for r in csv.DictReader(open(MASTER)):
    cls[r["b73_gene"]] = (r["imp_class"], r["meth_class"])
    for w in r["w22_genes"].replace(";", ",").split(","):
        if w.strip():
            w22map[w.strip()].append(r["b73_gene"])

mc = defaultdict(lambda: defaultdict(lambda: [0, 0, 0, 0]))   # gene -> dir -> counts
for r in csv.DictReader(open(METH)):
    if r["context"] != "CG":
        continue
    a = mc[r["gene"]][r["direction"]]
    for i, k in enumerate(("mat_met", "mat_un", "pat_met", "pat_un")):
        a[i] += int(r[k])

mcg_delta = {}
for w22g, dirs in mc.items():
    per = []
    for d, (mm, mu, pm, pu) in dirs.items():
        if mm + mu < 10 or pm + pu < 10:
            continue
        per.append(mm / (mm + mu) - pm / (pm + pu))
    if not per:
        continue
    for b73g in w22map.get(w22g, []):
        mcg_delta[b73g] = float(np.mean(per))
print(f"genes with a maternal-minus-paternal mCG delta: {len(mcg_delta):,}\n")

# ---------------------------------------------------------------- assemble groups
gbm_bip = [g for g in k27 if cls.get(g) == ("biparental", "gbM") and g in mcg_delta]
d = np.array([mcg_delta[g] for g in gbm_bip])
lo, hi = np.percentile(d, 25), np.percentile(d, 75)     # lo = most negative = demethylated

groups = {
    "gbM-PEG":                    [g for g in k27 if cls.get(g) == ("PEG", "gbM")],
    "UM-PEG":                     [g for g in k27 if cls.get(g) == ("PEG", "UM")],
    "gbM-biparental demethylated":    [g for g in gbm_bip if mcg_delta[g] <= lo],
    "gbM-biparental not demethylated": [g for g in gbm_bip if mcg_delta[g] >= hi],
    # The complete gbM-biparental set, which Figure 4 compares with the strongly demethylated
    # set. NB this group CONTAINS the demethylated quartile, so the two are not disjoint.
    "gbM-biparental (all)":       list(gbm_bip),
    "UM-biparental":              [g for g in k27 if cls.get(g) == ("biparental", "UM")],
}

print(f"gbM-biparental mCG delta quartiles: demethylated <= {lo:+.3f}, "
      f"not demethylated >= {hi:+.3f}\n")
print(f"{'group':<34}{'n':>6}{'med excess':>12}{'IQR':>22}{'p vs null':>12}")
print("-" * 88)
res = {}
for name, gl in groups.items():
    if not gl:
        print(f"{name:<34}{0:>6}   (empty)")
        continue
    e = np.array([k27[g]["excess"] for g in gl])
    s = S.paired(np.zeros_like(e), e)
    res[name] = e
    iqr = f"{s['iqr_lo']:+.4f} to {s['iqr_hi']:+.4f}"
    print(f"{name:<34}{len(e):>6}{s['median_delta']:>+12.4f}{iqr:>22}"
          f"{s['p_wilcoxon']:>12.3g}")

print("\n=== the test: gbM-PEG vs gbM-biparental demethylated ===")
if "gbM-PEG" in res and "gbM-biparental demethylated" in res:
    c = S.compare(res["gbM-PEG"], res["gbM-biparental demethylated"])
    print(f"  MWU p={c['p_mwu']:.3g}  probability of superiority {c['prob_superiority']:.3f}")

# ---------------------------------------------------------------- source data
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, f"05_panel_source_data{TAG}.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene", "group", "imp_class", "meth_class", "mat_frac_B73xMo17",
                "mat_frac_Mo17xB73", "mat_frac_mean", "maternal_excess",
                "strain_gap", "mCG_delta"])
    for name, gl in groups.items():
        for g in gl:
            k = k27[g]
            w.writerow([g, name, *cls.get(g, ("", "")), f"{k['f_861']:.4f}",
                        f"{k['f_862']:.4f}", f"{k['f_bar']:.4f}", f"{k['excess']:.4f}",
                        f"{k['strain_gap']:.4f}",
                        f"{mcg_delta[g]:.4f}" if g in mcg_delta else ""])
print(f"\nwrote {OUT}/05_panel_source_data.csv")
