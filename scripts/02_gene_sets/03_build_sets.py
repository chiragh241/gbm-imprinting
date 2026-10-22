#!/usr/bin/env python3
"""
Every evaluable B73 gene on two axes, the 16 gene sets and their BED files in B73 and W22.

  imprinting axis (01_imprinting_calls, final_call):
      MEG | PEG | biparental | undefined      (undefined = unknown + insufficient_data)
  methylation axis (01_categorize_methylation.py, methylation_category):
      gbM | UM | teM | undefined              (undefined = gene not categorised)

  -> 4 x 4 = 16 mutually exclusive sets over the 13,209 evaluable genes.

mCG, mCHG, cCG and cCHG are Zeng's per-gene CDS values, given for every gene in Zeng's table
whatever its class, so the table shows why an undefined gene is undefined.

Input   $GI_WORK/01_imprinting_calls/b73_imprinting_calls_all.csv
        $GI_WORK/02_gene_sets/genes_all_categorized.csv
        $GI_DATA/maize_tables/B73.all.csv                      (00_references)
        $GI_WORK/02_gene_sets/b73v5_to_w22.tsv                 (02_b73v5_to_w22.py)
        $GI_DATA/references/B73v5_genes.gff3, W22v2.gff3.gz     (00_references)
Output  $GI_WORK/02_gene_sets/master_table.csv   one row per gene: both classes, B73 coordinates,
                                                 W22 syntelogs, mCG, mCHG, cCG and cCHG
        $GI_WORK/02_gene_sets/b73/<imp>__<meth>.bed   B73 v5 coordinates, name = Zm00001eb id
        $GI_WORK/02_gene_sets/w22/<imp>__<meth>.bed   W22 v2 coordinates, name = Zm00004b|Zm00001eb,
                                                 one line per W22 syntelog
        $GI_WORK/02_gene_sets/counts_grid.csv    set sizes in both genomes
BED = 6 columns, 0-based half-open (start = GFF start - 1), score 0, strand.
"""
import re
from pathlib import Path
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

ROOT = WORK / "02_gene_sets"
REF = DATA / "references"
B73DIR = ROOT / "b73"; W22DIR = ROOT / "w22"
for d in (B73DIR, W22DIR): d.mkdir(parents=True, exist_ok=True)

IMP_ORDER  = ["MEG", "PEG", "biparental", "undefined"]
METH_ORDER = ["gbM", "UM", "teM", "undefined"]

# ---------------------------------------------------------------- imprinting axis
imp = pd.read_csv(WORK/"01_imprinting_calls"/"b73_imprinting_calls_all.csv")[["b73_gene", "final_call"]]
imp["imp_class"] = imp["final_call"].replace(
    {"unknown": "undefined", "insufficient_data": "undefined"})
imp = imp[["b73_gene", "imp_class"]]

# ---------------------------------------------------------------- methylation axis
cat = pd.read_csv(ROOT/"genes_all_categorized.csv", usecols=["gene", "methylation_category"])
assert set(cat["methylation_category"]) <= {"gbM", "UM", "teM"}
meth_of = dict(zip(cat["gene"], cat["methylation_category"]))
zeng = (pd.read_csv(DATA/"maize_tables"/"B73.all.csv", usecols=["gene", "mCG", "mCHG", "cCG", "cCHG"])
        .drop_duplicates("gene").set_index("gene"))

# ---------------------------------------------------------------- B73 v5 coords
def parse_gff(path, opener=open):
    coords = {}
    with opener(path) as fh:
        for line in fh:
            if isinstance(line, bytes): line = line.decode()
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene": continue
            gid = re.search(r"ID=([^;]+)", f[8])
            if gid: coords[gid.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return coords
b73_coord = parse_gff(REF/"B73v5_genes.gff3")

import gzip
w22_coord = parse_gff(REF/"W22v2.gff3.gz",
                      opener=lambda p: gzip.open(p, "rt"))

# ---------------------------------------------------------------- B73 -> W22 map
mp = pd.read_csv(ROOT/"b73v5_to_w22.tsv", sep="\t")
w22_of = {b: [x for x in str(w).split(",") if x.startswith("Zm00004b")]
          for b, w in zip(mp["b73_v5"], mp["w22"]) if isinstance(w, str)}

# ---------------------------------------------------------------- master table
imp["meth_class"] = imp["b73_gene"].map(meth_of).fillna("undefined")
imp["b73_chr"]    = imp["b73_gene"].map(lambda g: b73_coord.get(g, (None,))[0])
imp["b73_start0"] = imp["b73_gene"].map(lambda g: b73_coord[g][1]-1 if g in b73_coord else None)
imp["b73_end"]    = imp["b73_gene"].map(lambda g: b73_coord[g][2] if g in b73_coord else None)
imp["b73_strand"] = imp["b73_gene"].map(lambda g: b73_coord.get(g, (None,None,None,None))[3])
imp["w22_genes"]  = imp["b73_gene"].map(lambda g: ",".join(w22_of.get(g, [])))
for c in ("mCG", "mCHG", "cCG", "cCHG"):
    imp[c] = imp["b73_gene"].map(zeng[c])
imp.to_csv(ROOT/"master_table.csv", index=False)

n_no_b73 = imp["b73_chr"].isna().sum()
print(f"imprinting-analyzed genes: {len(imp)}")
print(f"  without B73 v5 GFF coords (dropped from B73 beds): {n_no_b73}")
print(f"  without any W22 ortholog  (dropped from W22 beds): {(imp['w22_genes']=='').sum()}")

# ---------------------------------------------------------------- write 16 + 16 BEDs
def bed_line(chrom, start0, end, name, strand):
    return f"{chrom}\t{start0}\t{end}\t{name}\t0\t{strand}\n"

b73_counts, w22_counts = {}, {}
for ic in IMP_ORDER:
    for mc in METH_ORDER:
        sub = imp[(imp.imp_class == ic) & (imp.meth_class == mc)]
        stem = f"{ic}__{mc}.bed"
        # B73
        nb = 0
        with open(B73DIR/stem, "w", newline="\n") as fh:
            for _, r in sub.iterrows():
                if r.b73_gene in b73_coord:
                    c, s, e, st = b73_coord[r.b73_gene]
                    fh.write(bed_line(c, s-1, e, r.b73_gene, st)); nb += 1
        # W22 (one line per ortholog with coords)
        nw = 0
        with open(W22DIR/stem, "w", newline="\n") as fh:
            for _, r in sub.iterrows():
                for w in w22_of.get(r.b73_gene, []):
                    if w in w22_coord:
                        c, s, e, st = w22_coord[w]
                        fh.write(bed_line(c, s-1, e, f"{w}|{r.b73_gene}", st)); nw += 1
        b73_counts[(ic, mc)] = nb; w22_counts[(ic, mc)] = nw

# ---------------------------------------------------------------- counts grid
def grid(counts, tag):
    g = pd.DataFrame(index=IMP_ORDER, columns=METH_ORDER, dtype=int)
    for (ic, mc), n in counts.items(): g.loc[ic, mc] = n
    g["All"] = g.sum(axis=1); g.loc["All"] = g.sum(axis=0)
    print(f"\n=== {tag} BED line counts ===\n{g}")
    return g
gb = grid(b73_counts, "B73 (genes)"); gw = grid(w22_counts, "W22 (gene-ortholog pairs)")
out = pd.concat({"B73_genes": gb, "W22_pairs": gw}, axis=0)
out.to_csv(ROOT/"counts_grid.csv")
print("\nwrote 16 B73 + 16 W22 BEDs, master_table.csv, counts_grid.csv ->", ROOT)
