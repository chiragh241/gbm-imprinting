#!/usr/bin/env bash
# BS-Seeker2 alignment of each endosperm library to W22 v2, in pieces.
#
#   bs_seeker2-align.py --aligner=bowtie2 --bt2-p 4 -f bam   (BS-Seeker2 2.1.7, bowtie2 2.5.4)
#
# BS-Seeker2 keeps every internal chunk's temporary files until a library finishes, which for a
# deep library approaches a terabyte. Each library is therefore split into pieces of 4,000,000
# read pairs, each piece is aligned as its own job (its temporary files are removed when it ends),
# and the piece BAMs are concatenated in order. Reads are aligned independently, so the result
# does not depend on the split. Two pieces run at a time (GI_MAP_JOBS to change).
#
# Input   $GI_WORK/05_allele_methylation/libs.tsv                  (02_merge_libraries.sh)
#         $GI_WORK/05_allele_methylation/bsseeker_index/           (03_bsseeker2_index.sh)
# Output  $GI_WORK/05_allele_methylation/map/<library>.bam         (read order, unsorted)
#
# About 1.5 days per library on 16 cores; the allele-specific tables this leads to are deposited
# on Zenodo (see README).
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_bsseeker2

W=$GI_WORK/05_allele_methylation
GENOME=$GI_DATA/references/W22v2.fa
IDX=$W/bsseeker_index
OUT=$W/map
SPLIT=$W/split
JTMP=$W/tmp
CHUNK=16000000            # lines per piece = 4,000,000 read pairs
MAXJ=${GI_MAP_JOBS:-2}
mkdir -p "$OUT" "$SPLIT" "$JTMP"

FAILED=0
while IFS=$'\t' read -r lib tissue direction r1 r2; do
    [ -f "$OUT/${lib}.align_done" ] && { echo "### skip $lib"; continue; }
    sdir=$SPLIT/$lib
    if [ ! -f "$sdir/.split_done" ]; then
        rm -rf "$sdir"; mkdir -p "$sdir"
        # identical -l on both mates keeps the pairs in step; zcat's own status is checked
        zcat "$r1" | split -l "$CHUNK" -d -a 4 --additional-suffix=.fq \
            --filter='gzip -1 > $FILE.gz' - "$sdir/R1_"
        rc1=${PIPESTATUS[0]}
        zcat "$r2" | split -l "$CHUNK" -d -a 4 --additional-suffix=.fq \
            --filter='gzip -1 > $FILE.gz' - "$sdir/R2_"
        rc2=${PIPESTATUS[0]}
        n1=$(ls "$sdir"/R1_*.fq.gz | wc -l); n2=$(ls "$sdir"/R2_*.fq.gz | wc -l)
        if [ "$rc1" -ne 0 ] || [ "$rc2" -ne 0 ] || [ "$n1" -ne "$n2" ]; then
            echo "!!! $lib: split failed"; rm -rf "$sdir"; FAILED=1; continue
        fi
        touch "$sdir/.split_done"
    fi
    echo "### $lib: $(ls "$sdir"/R1_*.fq.gz | wc -l) pieces  $(date)"
    for f1 in "$sdir"/R1_*.fq.gz; do
        idx=$(basename "$f1" .fq.gz); idx=${idx#R1_}
        [ -f "$OUT/${lib}.${idx}.done" ] && continue
        while [ "$(jobs -rp | wc -l)" -ge "$MAXJ" ]; do wait -n; done
        (
            t=$JTMP/${lib}_${idx}; rm -rf "$t"; mkdir -p "$t"
            if bs_seeker2-align.py -1 "$f1" -2 "$sdir/R2_${idx}.fq.gz" -g "$GENOME" -d "$IDX" \
                   --aligner=bowtie2 --bt2-p 4 -f bam -o "$OUT/${lib}.${idx}.bam" --temp_dir "$t" \
                   > "$OUT/${lib}.${idx}.log" 2>&1 \
               && samtools quickcheck "$OUT/${lib}.${idx}.bam"; then
                touch "$OUT/${lib}.${idx}.done"
            else
                echo "!!! $lib piece $idx failed (see $OUT/${lib}.${idx}.log)"
            fi
            rm -rf "$t"
        ) &
    done
    wait
    n=$(ls "$sdir"/R1_*.fq.gz | wc -l)
    ndone=$(ls "$OUT/${lib}".[0-9]*.done 2>/dev/null | wc -l)
    if [ "$ndone" -ne "$n" ]; then echo "!!! $lib: $ndone of $n pieces done"; FAILED=1; continue; fi
    pieces=("$OUT/${lib}".[0-9]*.bam)    # glob order is piece order: indices are zero-padded
    samtools cat -o "$OUT/${lib}.bam" "${pieces[@]}" && samtools quickcheck "$OUT/${lib}.bam" \
        && touch "$OUT/${lib}.align_done" \
        && rm -f "$OUT/${lib}".[0-9]*.bam "$OUT/${lib}".[0-9]*.done && rm -rf "$sdir"
done < "$W/libs.tsv"
exit $FAILED
