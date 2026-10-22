#!/usr/bin/env python3
"""B73 v5 to W22 v2 syntelogs, from the MaizeGDB 2020 pan-gene table.

Each row of the pan-gene table is one pan-gene: every gene, in every genome, assigned to it.
Every B73 v5 gene (Zm00001eb) in a row is paired with all W22 v2 genes (Zm00004b) of that row,
sorted and comma-joined; a B73 gene with no W22 gene in its row keeps an empty w22 field. Rows
are written in table order, each B73 gene once.

Input   $GI_DATA/references/MaizeGDB_maize_pangene_2020_08.tsv.gz
        (00_references/01_download_references.sh)
Output  $GI_WORK/02_gene_sets/b73v5_to_w22.tsv   b73_v5, w22   (36,810 genes)
"""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from lib.paths import DATA, WORK

OUT = WORK / "02_gene_sets"
OUT.mkdir(parents=True, exist_ok=True)

out, seen = [], set()
with gzip.open(DATA / "references" / "MaizeGDB_maize_pangene_2020_08.tsv.gz", "rt") as fh:
    for line in fh:
        row = line.rstrip("\n").split("\t")
        w22 = sorted(g for g in row if g.startswith("Zm00004b"))
        for g in row:
            if g.startswith("Zm00001eb") and g not in seen:
                seen.add(g)
                out.append(f"{g}\t{','.join(w22)}")

with open(OUT / "b73v5_to_w22.tsv", "w", newline="\n") as fh:
    fh.write("b73_v5\tw22\n" + "\n".join(out) + "\n")
n_w22 = sum(1 for r in out if not r.endswith("\t"))
print(f"wrote {len(out):,} B73 v5 genes, {n_w22:,} with at least one W22 syntelog")
