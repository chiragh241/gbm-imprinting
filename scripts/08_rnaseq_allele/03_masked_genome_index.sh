#!/usr/bin/env bash
# B73 v5 N-masked at the B73/Mo17 variants and its HISAT2 index, so reads of either allele align
# without reference bias; and exon (SAF) and gene (BED) annotations from the B73 v5 GFF3.
#
# Input   $GI_DATA/references/B73v5.fa, B73v5.gff3   (00_references)
#         $GI_WORK/08_rnaseq_allele/snps/snps.bed      (02_build_snps.sh)
# Output  $GI_WORK/08_rnaseq_allele/genome/B73v5.Nmasked.fa
#         $GI_WORK/08_rnaseq_allele/index/B73v5_Nmasked.*.ht2
#         $GI_WORK/08_rnaseq_allele/annot/exons.saf, genes.bed
# The masked genome and index are also used by the ChIP-seq allele module.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

REF=$GI_DATA/references/B73v5.fa
GFF=$GI_DATA/references/B73v5.gff3
W=$GI_WORK/08_rnaseq_allele
GEN=$W/genome; IDX=$W/index; ANN=$W/annot
mkdir -p "$GEN" "$IDX" "$ANN"

[ -s "$REF.fai" ] || samtools faidx "$REF"
NMASK=$GEN/B73v5.Nmasked.fa
# The masked genome is written under a temporary name and the index marked when built, so an
# interrupted run never leaves a truncated genome or a partial index for the next run to skip.
if [ ! -s "$NMASK" ]; then
    bedtools maskfasta -fi "$REF" -bed "$W/snps/snps.bed" -fo "$NMASK.part"
    mv -f "$NMASK.part" "$NMASK"
    samtools faidx "$NMASK"
fi

python3 - "$GFF" "$ANN" <<'PY'
import re, sys
gff, ann = sys.argv[1], sys.argv[2]
exon_rows = []; genes = {}
for line in open(gff):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9: continue
    ftype = f[2]; attrs = f[8]
    if ftype == "exon":
        m = re.search(r'Parent=([^;]+)', attrs)
        if not m: continue
        parent = m.group(1).split(",")[0]
        gene = re.sub(r'^gene:', '', re.sub(r'_T\d+$', '', parent))
        exon_rows.append((gene, f[0], f[3], f[4], f[6]))
    elif ftype == "gene":
        m = re.search(r'ID=([^;]+)', attrs)
        if m:
            gid = re.sub(r'^gene:', '', m.group(1))
            genes[gid] = (f[0], int(f[3]), int(f[4]), f[6])
with open(f"{ann}/exons.saf", "w") as o:
    o.write("GeneID\tChr\tStart\tEnd\tStrand\n")
    for g, c, s, e, st in exon_rows: o.write(f"{g}\t{c}\t{s}\t{e}\t{st}\n")
with open(f"{ann}/genes.bed", "w") as o:
    for gid, (c, s, e, st) in genes.items(): o.write(f"{c}\t{s-1}\t{e}\t{gid}\t.\t{st}\n")
print(f"exon rows: {len(exon_rows)}  genes: {len(genes)}")
PY

[ -f "$IDX/.build_done" ] || { hisat2-build -p "$GI_THREADS" "$NMASK" "$IDX/B73v5_Nmasked"; touch "$IDX/.build_done"; }
ls -lh "$IDX"
