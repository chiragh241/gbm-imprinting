#!/usr/bin/env python3
"""Per-gene maternal and paternal methylation from the allele-specific tables, per gene set.

Input   $GI_WORK/05_allele_methylation/asm/<library>.<context>.asm.tsv   (10_asm.sh)
        $GI_WORK/05_allele_methylation/libs.tsv                           (02_merge_libraries.sh)
        $GI_WORK/02_gene_sets/w22/*.bed                                    (02_gene_sets)
        $GI_WORK/05_allele_methylation/snv/.keep_ctrl.txt                 (09_asm_snv_subset.py)
Output  $GI_WORK/05_allele_methylation/gene_summary/
          12_per_gene_allele_methylation.csv   per gene, library and context: maternal and
                                               paternal methylated and unmethylated reads
          12_class_summary.csv                 per gene set: median maternal, paternal, delta
          12_control_baseline.csv              the non-genic control variants, pooled
Genes need at least 3 linked cytosines and 10 reads on each allele for the class summary.

WHAT THIS ANSWERS: for each gene, is the MATERNAL allele less methylated than the PATERNAL one, and
does that differ by imprinting class? The PEG body-demethylation model predicts maternal-specific
loss of gene-body mCG at gbM PEGs, with biparental gbM genes flat.

MATERNAL/PATERNAL COMES FROM THE CROSS, NOT THE GENOTYPE. The reference is W22v2, so cgmaptools'
Allele1 (= REF) is the W22 allele and Allele2 (= ALT) is the B73 allele. libs.tsv records the cross
as mother x father, so:
    B73xW22  -> mother is B73 -> maternal = Allele2, paternal = Allele1
    W22xB73  -> mother is W22 -> maternal = Allele1, paternal = Allele2
That Allele1 == Ref identity is an assumption about how cgmaptools orders alleles for a 0/1 GT, so
it is CHECKED per row rather than trusted, and violations are counted and reported.

WHY THE RECIPROCAL SPLIT IS REPORTED SEPARATELY: pooling both directions is what cancels strain
effects, but pooling also HIDES them. A real parent-of-origin effect has the same sign in both
directions; a B73-vs-W22 strain effect flips sign when the cross flips. The per-direction deltas are
therefore always emitted alongside the pooled one - if they disagree in sign, the pooled number is
not measuring imprinting and should not be reported as such.

DE-DUPLICATION IS NOT OPTIONAL: asm emits one row per (SNP, linked C) pair, so a cytosine linked to
three nearby htSNPs appears three times with the same reads behind it. Summing naively would triple
that cytosine's weight. Rows are collapsed to one per (Chr, C_Pos) keeping the best-covered, and the
number collapsed is reported - if it is a large fraction, SNP density is doing more of the weighting
than coverage is.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import gene_intervals
from lib.paths import WORK

ASM_COLS = ["Chr", "SNP_Pos", "Ref", "Allele1", "Allele2", "C_Pos",
            "Allele1_linked_C", "Allele2_linked_C",
            "Allele1_linked_C_met", "Allele2_linked_C_met", "pvalue", "fdr", "ASM"]


def load_gene_sets(bed_dir):
    """Read the 16 imprint x methylation BEDs into one gene table."""
    rows = []
    for path in sorted(glob.glob(os.path.join(bed_dir, "*.bed"))):
        imprint, meth = os.path.basename(path)[:-4].split("__")
        bed = pd.read_csv(path, sep="\t", header=None,
                          names=["chrom", "start", "end", "name", "score", "strand"])
        if bed.empty:
            continue
        bed["imprint_class"] = imprint
        bed["meth_class"] = meth
        # name is "W22gene|B73gene"; the W22 id is the one matching these coordinates
        bed["gene"] = bed["name"].str.split("|").str[0]
        rows.append(bed[["gene", "chrom", "start", "end", "imprint_class", "meth_class"]])
    if not rows:
        sys.exit(f"no BED files found in {bed_dir}")
    genes = pd.concat(rows, ignore_index=True)
    dup = genes["gene"].duplicated().sum()
    if dup:
        print(f"  note: {dup} gene ids appear in more than one class file; keeping first "
              f"(all carry identical coordinates, so only the class label is arbitrary)")
        genes = genes.drop_duplicates(subset="gene", keep="first")
    # gene_intervals returns index labels that are used positionally downstream, so the
    # index must be a clean 0..n-1 range after the drop above leaves gaps in it.
    return genes.reset_index(drop=True)


def parse_counts(series):
    """'8-2' -> (methylated, unmethylated) as int arrays."""
    parts = series.str.split("-", n=1, expand=True)
    met = pd.to_numeric(parts[0], errors="coerce")
    unmet = pd.to_numeric(parts[1], errors="coerce")
    return met, unmet


def assign_genes(df, genes, index=None):
    """Map each C_Pos to EVERY gene whose body contains it.

    BED is 0-based half-open, asm C_Pos is 1-based, so position p is inside when start < p <= end.

    Delegates to gene_intervals, which scans back far enough to catch earlier
    overlapping genes. The previous inline version inspected only the last gene
    starting before the position, so a cytosine inside an earlier still-open gene
    failed the end test against the wrong candidate and was dropped outright. Maize
    models overlap enough for that to matter: 275 of 12,842 genes (2.14%) overlap
    their predecessor. A cytosine inside two overlapping bodies is now counted for
    both, which is what the containment query means.
    """
    if index is None:
        index = gene_intervals.build_index(genes)
    qi, gidx = gene_intervals.assign(df["Chr"].to_numpy(), df["C_Pos"].to_numpy(), index)
    if len(qi) == 0:
        return pd.DataFrame(columns=list(df.columns) + ["gene"])
    hit = df.iloc[qi].copy()
    hit["gene"] = genes["gene"].to_numpy()[gidx]
    return hit.reset_index(drop=True)


def load_library(path, direction, genes, ctrl=None, index=None):
    """One asm .tsv -> per-gene maternal/paternal read counts.

    `ctrl` is the set of "chr:pos" htSNPs that are NON-GENIC controls (09_asm_snv_subset.py's
    .keep_ctrl.txt). They are the randomised baseline and are deliberately in the
    10_asm.sh input VCF, but every row here is assigned to a gene and non-hits are
    dropped, so without this they leave no trace in any output of this script. Their
    weighted counts ride along in `stats` because this function is the only place
    the asm file is read.
    """
    df = pd.read_csv(path, sep="\t", dtype={"Chr": str})
    missing = [c for c in ASM_COLS if c not in df.columns]
    if missing:
        # ASM is only present in some cgmaptools builds; everything before it is required
        if [c for c in missing if c != "ASM"]:
            sys.exit(f"{path}: unexpected columns {list(df.columns)}")
    n_raw = len(df)

    # allele identity check: Allele1 should be the reference (= W22) allele
    bad = int((df["Allele1"] != df["Ref"]).sum())

    a1_met, a1_un = parse_counts(df["Allele1_linked_C"])
    a2_met, a2_un = parse_counts(df["Allele2_linked_C"])
    df = df.assign(a1_met=a1_met, a1_un=a1_un, a2_met=a2_met, a2_un=a2_un)
    df = df.dropna(subset=["a1_met", "a1_un", "a2_met", "a2_un"])

    # collapse duplicate cytosines (same C linked to several htSNPs), keeping best covered.
    # pandas' default sort is not stable, so which row survives a coverage tie can differ by
    # platform: the paper's values are from WSL. Kept as run, so the output matches the paper's.
    df["_cov"] = df["a1_met"] + df["a1_un"] + df["a2_met"] + df["a2_un"]
    before = len(df)
    df = (df.sort_values("_cov", ascending=False)
            .drop_duplicates(subset=["Chr", "C_Pos"], keep="first"))
    collapsed = before - len(df)

    # Allele1 = W22 (reference), Allele2 = B73
    if direction == "B73xW22":       # mother is B73
        mat_met, mat_un, pat_met, pat_un = df.a2_met, df.a2_un, df.a1_met, df.a1_un
    elif direction == "W22xB73":     # mother is W22
        mat_met, mat_un, pat_met, pat_un = df.a1_met, df.a1_un, df.a2_met, df.a2_un
    else:
        sys.exit(f"unrecognised cross direction {direction!r} (expected B73xW22 or W22xB73)")
    df = df.assign(mat_met=mat_met, mat_un=mat_un, pat_met=pat_met, pat_un=pat_un)

    # non-genic control baseline, taken before the gene filter discards these rows
    cstat = dict(ctrl_sites=0, ctrl_mat_met=0, ctrl_mat_un=0,
                 ctrl_pat_met=0, ctrl_pat_un=0)
    if ctrl:
        isc = (df.Chr.astype(str) + ":" + df.SNP_Pos.astype(str)).isin(ctrl)
        c = df[isc]
        cstat = dict(ctrl_sites=len(c),
                     ctrl_mat_met=int(c.mat_met.sum()), ctrl_mat_un=int(c.mat_un.sum()),
                     ctrl_pat_met=int(c.pat_met.sum()), ctrl_pat_un=int(c.pat_un.sum()))

    hit = assign_genes(df, genes, index)
    per_gene = (hit.groupby("gene")[["mat_met", "mat_un", "pat_met", "pat_un"]].sum()
                   .join(hit.groupby("gene").size().rename("n_sites")))
    # a cytosine inside overlapping bodies is counted for each, so hit rows can exceed
    # the distinct cytosines that landed in any gene
    multi = len(hit) - hit.drop_duplicates(subset=["Chr", "C_Pos"]).shape[0]
    stats = dict(rows=n_raw, collapsed=collapsed, allele1_ne_ref=bad,
                 sites_in_genes=len(hit), multi_gene=multi,
                 genes=len(per_gene), **cstat)
    return per_gene.reset_index(), stats


def summarise(per_gene, min_sites, min_reads):
    """Per-gene fractions with the coverage filter applied."""
    d = per_gene.copy()
    d["mat_n"] = d.mat_met + d.mat_un
    d["pat_n"] = d.pat_met + d.pat_un
    d = d[(d.n_sites >= min_sites) & (d.mat_n >= min_reads) & (d.pat_n >= min_reads)]
    d["mat_frac"] = d.mat_met / d.mat_n
    d["pat_frac"] = d.pat_met / d.pat_n
    d["delta"] = d.mat_frac - d.pat_frac      # negative = maternal allele less methylated
    return d


def class_stats(d, label):
    """Paired maternal-vs-paternal test across the genes of one class."""
    n = len(d)
    if n == 0:
        return dict(group=label, n_genes=0, med_mat=np.nan, med_pat=np.nan,
                    med_delta=np.nan, pvalue=np.nan)
    p = np.nan
    if n >= 6 and np.any(d.delta != 0):
        try:
            p = wilcoxon(d.mat_frac, d.pat_frac, zero_method="wilcox").pvalue
        except ValueError:
            p = np.nan
    return dict(group=label, n_genes=n,
                med_mat=d.mat_frac.median(), med_pat=d.pat_frac.median(),
                med_delta=d.delta.median(), pvalue=p)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--asm-dir", default=str(WORK / "05_allele_methylation" / "asm"))
    ap.add_argument("--libs", default=str(WORK / "05_allele_methylation" / "libs.tsv"))
    ap.add_argument("--bed-dir", default=str(WORK / "02_gene_sets" / "w22"))
    ap.add_argument("--out-dir", default=str(WORK / "05_allele_methylation" / "gene_summary"))
    ap.add_argument("--contexts", default="CG,CHG,CHH")
    ap.add_argument("--min-sites", type=int, default=3,
                    help="min informative cytosines per gene [3]")
    ap.add_argument("--min-reads", type=int, default=10,
                    help="min reads per allele per gene [10]")
    ap.add_argument("--ctrl-list",
                    default=str(WORK / "05_allele_methylation" / "snv" / ".keep_ctrl.txt"),
                    help="09_asm_snv_subset.py's .keep_ctrl.txt of non-genic control htSNPs; "
                         "without it the randomised baseline is not reported")
    args = ap.parse_args()

    ctrl = None
    if args.ctrl_list and os.path.exists(args.ctrl_list):
        with open(args.ctrl_list) as fh:
            ctrl = {ln.strip() for ln in fh if ln.strip()}
        print(f"control htSNPs: {len(ctrl)} (non-genic randomised baseline)")

    os.makedirs(args.out_dir, exist_ok=True)
    genes = load_gene_sets(args.bed_dir)
    print(f"gene sets: {len(genes)} genes across "
          f"{genes.imprint_class.nunique()} imprint x {genes.meth_class.nunique()} methylation classes")
    # built once; rebuilding it per library x context would repeat the overlap scan 18 times
    gindex = gene_intervals.build_index(genes)
    print(f"  interval index: backward scan depth per chromosome "
          f"{sorted({v[3] for v in gindex.values()})}")

    libs = pd.read_csv(args.libs, sep="\t", header=None,
                       names=["lib", "tissue", "direction", "r1", "r2"])
    libs = libs[libs.tissue == "endosperm"]

    all_gene_rows, summary_rows, ctrl_rows = [], [], []
    for ctx in args.contexts.split(","):
        pooled = []
        for _, lib in libs.iterrows():
            path = os.path.join(args.asm_dir, f"{lib.lib}.{ctx}.asm.tsv")
            if not os.path.exists(path):
                print(f"  {ctx} {lib.lib}: MISSING {os.path.basename(path)} - skipped")
                continue
            per_gene, st = load_library(path, lib.direction, genes, ctrl, gindex)
            if ctrl:
                ctrl_rows.append(dict(context=ctx, lib=lib.lib,
                                      direction=lib.direction, **{
                    k: st[k] for k in ("ctrl_sites", "ctrl_mat_met", "ctrl_mat_un",
                                       "ctrl_pat_met", "ctrl_pat_un")}))
            print(f"  {ctx} {lib.lib} ({lib.direction}): {st['rows']} rows, "
                  f"{st['collapsed']} dup C collapsed, {st['sites_in_genes']} in genes "
                  f"({st['multi_gene']} in overlapping bodies), {st['genes']} genes"
                  + (f", WARNING {st['allele1_ne_ref']} rows Allele1!=Ref" if st["allele1_ne_ref"] else ""))
            per_gene["lib"] = lib.lib
            per_gene["direction"] = lib.direction
            per_gene["context"] = ctx
            pooled.append(per_gene)

        if not pooled:
            print(f"  {ctx}: no libraries - skipping")
            continue
        pooled = pd.concat(pooled, ignore_index=True)
        all_gene_rows.append(pooled)

        # pooled over all six libraries, and each direction on its own
        for scope, sub in [("pooled", pooled),
                           ("B73xW22", pooled[pooled.direction == "B73xW22"]),
                           ("W22xB73", pooled[pooled.direction == "W22xB73"])]:
            agg = sub.groupby("gene")[["mat_met", "mat_un", "pat_met", "pat_un", "n_sites"]].sum()
            d = summarise(agg.reset_index(), args.min_sites, args.min_reads)
            d = d.merge(genes[["gene", "imprint_class", "meth_class"]], on="gene", how="left")
            for (imp, meth), grp in d.groupby(["imprint_class", "meth_class"]):
                row = class_stats(grp, f"{imp}__{meth}")
                row.update(context=ctx, scope=scope, imprint_class=imp, meth_class=meth)
                summary_rows.append(row)

    if not all_gene_rows:
        sys.exit("no asm results found - run 10_asm.sh first")

    per_gene_out = os.path.join(args.out_dir, "12_per_gene_allele_methylation.csv")
    pd.concat(all_gene_rows, ignore_index=True).to_csv(per_gene_out, index=False)

    summary = pd.DataFrame(summary_rows)[
        ["context", "scope", "imprint_class", "meth_class", "n_genes",
         "med_mat", "med_pat", "med_delta", "pvalue"]]
    summary_out = os.path.join(args.out_dir, "12_class_summary.csv")
    summary.to_csv(summary_out, index=False)

    print(f"\nwrote {per_gene_out}\nwrote {summary_out}")

    # Randomised baseline. Reported as COVERAGE-WEIGHTED fractions, not as medians of
    # per-gene fractions like the class summary, because control htSNPs have no gene to
    # aggregate within. Kept in its own file so the two estimators are never mixed in
    # one table.
    if ctrl_rows:
        c = pd.DataFrame(ctrl_rows)
        agg = c.groupby("context")[["ctrl_sites", "ctrl_mat_met", "ctrl_mat_un",
                                    "ctrl_pat_met", "ctrl_pat_un"]].sum().reset_index()
        agg["mat_frac"] = agg.ctrl_mat_met / (agg.ctrl_mat_met + agg.ctrl_mat_un)
        agg["pat_frac"] = agg.ctrl_pat_met / (agg.ctrl_pat_met + agg.ctrl_pat_un)
        agg["delta"] = agg.mat_frac - agg.pat_frac
        ctrl_out = os.path.join(args.out_dir, "12_control_baseline.csv")
        agg.to_csv(ctrl_out, index=False, float_format="%.5f")
        print(f"wrote {ctrl_out}")
        print("\n=== non-genic control baseline (weighted; NOT medians) ===")
        print(agg[["context", "ctrl_sites", "mat_frac", "pat_frac", "delta"]]
              .to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    else:
        print("\n*** WARNING: no --ctrl-list given, so the non-genic randomised "
              "baseline\n*** is NOT in any "
              "output. Every row\n*** here is gene-assigned; control htSNPs are "
              "silently dropped.")

    # headline: the comparison this whole arm was built to make
    print("\n=== gbM genes, CG context, pooled (delta = maternal - paternal) ===")
    head = summary[(summary.context == "CG") & (summary.scope == "pooled")
                   & (summary.meth_class == "gbM")]
    if head.empty:
        print("  no gbM/CG rows survived the coverage filter")
    else:
        print(head.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        print("\n  reciprocal check (a parent-of-origin effect keeps its sign in both directions):")
        recip = summary[(summary.context == "CG") & (summary.meth_class == "gbM")
                        & (summary.scope != "pooled")]
        print(recip[["imprint_class", "scope", "n_genes", "med_delta"]].to_string(
            index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
