#!/usr/bin/env bash
# B73/Mo17 single-nucleotide variants in B73 v5 coordinates, from whole-genome alignment,
# independent of any published variant set.
#
#   Mo17 CAU-2.0, split into 1 Mb pieces to bound memory -> minimap2 2.28 -cx asm5 --cs to B73 v5
#   -> paftools.js call -L 5000 -> biallelic SNVs (4,636,099)
#
# Input   $GI_DATA/references/B73v5.fa, Mo17_CAU2.fa.gz   (00_references)
# Output  $GI_WORK/08_rnaseq_allele/snps/
#           snps.bed            0-based, used to N-mask the genome and to find allele-informative reads
#           snps.snpsplit.txt   SNPsplit format: id, chrom, pos, 1, B73/Mo17 (genome1 = B73)
#           snps.tsv            chrom, pos, B73 base, Mo17 base
#           stats.txt           counts per sequence
# Also used by the ChIP-seq allele module (10_chip_endosperm_allele).
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

REF=$GI_DATA/references/B73v5.fa
W=$GI_WORK/08_rnaseq_allele
SNPD=$W/snps
mkdir -p "$SNPD" "$W/ref"
MO17=$W/ref/Mo17_CAU2.fa
MO17CHUNK=$W/ref/Mo17_CAU2.1Mb.fa

# Each step that is skipped when its output exists writes under a temporary name and renames on
# success, so an interrupted run never leaves a truncated file for the next run to skip.
[ -s "$MO17" ] || { zcat "$GI_DATA/references/Mo17_CAU2.fa.gz" > "$MO17.part"; mv -f "$MO17.part" "$MO17"; }

# 1 Mb pieces: whole 200 Mb chromosomes as single queries need more than 24 GB in minimap2.
if [ ! -s "$MO17CHUNK" ]; then
python3 - "$MO17" "$MO17CHUNK.part" <<'PY'
import sys
inp, out = sys.argv[1], sys.argv[2]
CH = 1000000
def rd(fh):
    n = None; s = []
    for l in fh:
        if l[0] == '>':
            if n: yield n, ''.join(s)
            n = l[1:].split()[0]; s = []
        else:
            s.append(l.strip())
    if n: yield n, ''.join(s)
with open(inp) as fi, open(out, 'w') as fo:
    for n, seq in rd(fi):
        for i in range(0, len(seq), CH):
            fo.write(">%s_%d\n%s\n" % (n, i, seq[i:i + CH]))
PY
mv -f "$MO17CHUNK.part" "$MO17CHUNK"
fi

PAF=$SNPD/mo17_vs_b73v5.srt.paf
if [ ! -s "$PAF" ]; then
    minimap2 -cx asm5 --cs -t 6 -I 4G "$REF" "$MO17CHUNK" 2> "$SNPD/minimap2.log" \
        | sort -k6,6 -k8,8n -S 1G > "$PAF.part"
    mv -f "$PAF.part" "$PAF"
fi
VAR=$SNPD/mo17_b73v5.var.txt
[ -s "$VAR" ] || { paftools.js call -L 5000 -f "$REF" "$PAF" > "$VAR.part" 2> "$SNPD/paftools.log"
                   mv -f "$VAR.part" "$VAR"; }

python3 - "$VAR" "$SNPD" <<'PY'
import sys
var, outd = sys.argv[1], sys.argv[2]
ACGT = set("ACGT")
bed = open(f"{outd}/snps.bed", "w"); ss = open(f"{outd}/snps.snpsplit.txt", "w")
tsv = open(f"{outd}/snps.tsv", "w")
n = 0; per = {}
for line in open(var):
    if line.startswith("#"): continue          # paftools -f emits VCF
    f = line.rstrip("\n").split("\t")
    if len(f) < 5: continue
    chrom = f[0]; ref = f[3].upper(); alt = f[4].upper()
    if len(ref) != 1 or len(alt) != 1 or ref not in ACGT or alt not in ACGT or ref == alt: continue
    pos1 = int(f[1]); start = pos1 - 1; end = pos1
    n += 1; per[chrom] = per.get(chrom, 0) + 1
    bed.write(f"{chrom}\t{start}\t{end}\n")
    ss.write(f"{n}\t{chrom}\t{pos1}\t1\t{ref}/{alt}\n")   # Ref = B73 (genome1), SNP = Mo17 (genome2)
    tsv.write(f"{chrom}\t{pos1}\t{ref}\t{alt}\n")
bed.close(); ss.close(); tsv.close()
with open(f"{outd}/stats.txt", "w") as st:
    st.write(f"total_SNVs\t{n}\n")
    for c in sorted(per): st.write(f"{c}\t{per[c]}\n")
print(f"biallelic SNVs: {n}")
PY
