#!/usr/bin/env bash
# Allele-assigned H3K27me3 reads and allele-informative bp in each compact metagene bin.
#
# Each read is counted once, in the bin that holds its leftmost aligned base: the position the
# informative windows are defined on ([variant - 35, variant], the starts from which a 36 bp read
# covers the variant), so reads and informative bp count the same thing. Informative bp per bin is
# the bin's overlap with the merged informative windows.
#
# Genome-wide denominators, written here too:
#   <run>.genome{1,2}.total.txt   allele-assigned reads in each split BAM
#   genome_informative_bp.txt     total length of the informative windows
#
# Input   $GI_WORK/10_chip_endosperm_allele/split/<run>.inf.genome{1,2}.bam   (02_align_split_count.sh)
#         $GI_WORK/10_chip_endosperm_allele/annot/informative_windows.bed     (02_align_split_count.sh)
#         $GI_WORK/10_chip_endosperm_allele/metagene_bins_v5_compact.bed      (05_compact_bins.py)
# Output  $GI_WORK/10_chip_endosperm_allele/counts_allele_compact/<run>.genome{1,2}.bincounts.tsv
#           group gene scheme bin reads
#         $GI_WORK/10_chip_endosperm_allele/bin_informative_bp_compact.tsv
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads

W=$GI_WORK/10_chip_endosperm_allele
SPLIT=$W/split
INF=$W/annot/informative_windows.bed
BINS=$W/metagene_bins_v5_compact.bed
OUT=$W/counts_allele_compact
mkdir -p "$OUT"

awk 'BEGIN{FS="\t"} {s+=$3-$2} END {print s}' "$INF" > "$W/genome_informative_bp.txt"

for run in SRR748861 SRR748862; do
    for gn in 1 2; do
        bam=$SPLIT/${run}.inf.genome${gn}.bam
        o=$OUT/${run}.genome${gn}.bincounts.tsv
        samtools view -c "$bam" > "$OUT/${run}.genome${gn}.total.txt"
        [ -s "$o" ] && continue
        # leftmost aligned base as a 1 bp interval; FS="\t": group names contain spaces
        bedtools bamtobed -i "$bam" \
            | awk 'BEGIN{OFS="\t"} {print $1, $2, $2 + 1}' \
            | bedtools intersect -a stdin -b "$BINS" -wa -wb \
            | awk 'BEGIN{FS=OFS="\t"} {n[$7"\t"$8"\t"$9"\t"$10]++}
                   END{for(k in n) print k,n[k]}' > "$o.part"
        mv -f "$o.part" "$o"
        echo "$run genome$gn: $(cat "$OUT/${run}.genome${gn}.total.txt") assigned reads"
    done
done

IBP=$W/bin_informative_bp_compact.tsv
printf 'group\tgene\tscheme\tbin\tinf_bp\n' > "$IBP.part"
bedtools intersect -a "$BINS" -b "$INF" -wao \
    | awk 'BEGIN{FS=OFS="\t"} {n[$4 OFS $5 OFS $6 OFS $7] += $11}
           END {for (k in n) print k, n[k]}' \
    | sort -k1,1 -k2,2 -k3,3 -k4,4n >> "$IBP.part"
mv -f "$IBP.part" "$IBP"
echo "informative-bp rows: $(( $(wc -l < "$IBP") - 1 )) (bins: $(wc -l < "$BINS"))"
