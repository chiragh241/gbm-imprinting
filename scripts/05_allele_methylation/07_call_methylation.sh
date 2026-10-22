#!/usr/bin/env bash
# Per-cytosine ATCGmap files of the pooled endosperm BAM, one per sequence, for SNV calling.
#
#   bs_seeker2-call_methylation.py --sorted, run on each sequence's slice of the pooled BAM
#
# The caller must carry the patch in envs/patches (applied by envs/setup_envs.sh): unpatched,
# its pysam pileup drops every read without the proper-pair flag, which BS-Seeker2 never sets,
# and it writes empty ATCGmaps while exiting 0. The check below refuses to run without it.
#
# Input   $GI_WORK/05_allele_methylation/combined/endosperm_combined.bam   (06_pool_endosperm.sh)
#         $GI_WORK/05_allele_methylation/bsseeker_index/W22v2.fa_bowtie2/   (03_bsseeker2_index.sh)
# Output  $GI_WORK/05_allele_methylation/atcgmap/<sequence>.ATCGmap.gz
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_bsseeker2

W=$GI_WORK/05_allele_methylation
IN=$W/combined/endosperm_combined.bam
OUT=$W/atcgmap
TMP=$W/tmp
DB=$W/bsseeker_index/W22v2.fa_bowtie2
MAXJ=${GI_CALL_JOBS:-2}      # each caller holds one chromosome's reference in memory
mkdir -p "$OUT" "$TMP"

CALLER=$(command -v bs_seeker2-call_methylation.py)
grep -q 'ignore_orphans=False' "$CALLER" \
    || { echo "!!! $CALLER is unpatched; run envs/setup_envs.sh"; exit 1; }

call_one() {
    local chr=$1 cbam=$TMP/${1}.bam
    samtools view -b -o "$cbam" "$IN" "$chr" && samtools index "$cbam" \
      && bs_seeker2-call_methylation.py --sorted -i "$cbam" -d "$DB" -o "$OUT/$chr" \
             > "$OUT/${chr}.call.log" 2>&1 \
      && [ "$(stat -c %s "$OUT/${chr}.ATCGmap.gz")" -ge 100000 ] \
      && touch "$OUT/${chr}.done" \
      || echo "!!! $chr failed (see $OUT/${chr}.call.log)"
    rm -f "$cbam" "$cbam.bai"
}

# largest sequence first, so the longest job starts immediately
mapfile -t CHRS < <(samtools idxstats "$IN" | awk '$1!="*" && $3>0 {print $2"\t"$1}' \
                    | sort -rn | cut -f2)
for chr in "${CHRS[@]}"; do
    [ -f "$OUT/${chr}.done" ] && continue
    while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
    call_one "$chr" &
done
wait
n=0; for chr in "${CHRS[@]}"; do [ -f "$OUT/${chr}.done" ] && n=$((n + 1)); done
echo "$n of ${#CHRS[@]} sequences called"
[ "$n" -eq "${#CHRS[@]}" ]
