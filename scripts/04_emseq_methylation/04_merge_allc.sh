#!/usr/bin/env bash
# Sum the three replicates of each genotype and tissue per cytosine (methylpy merge-allc).
#
# Input   $GI_WORK/04_emseq_methylation/allc/allc_<tag>.tsv.gz   (03_map.sh)
# Output  $GI_WORK/04_emseq_methylation/allc/allc_<sample>_merged.tsv.gz, eight samples:
#           wt_15dap, mdr1_15dap, wt_mature, mdr1_mature, dng102_mature, embryo, leaf, tassel
#         These eight files are deposited on Zenodo (see README): placing them here lets every
#         later step run without mapping any reads.
set -u
source "$(dirname "$0")/../../config.sh"
gi_env gi_reads
OUT=$GI_WORK/04_emseq_methylation/allc

merge() {  # sample rep_tags...
    local samp=$1; shift
    local files=""
    if [ -s "$OUT/allc_${samp}_merged.tsv.gz" ]; then echo "### $samp already merged"; return 0; fi
    for r in "$@"; do
        local fp="$OUT/allc_${r}.tsv.gz"
        # a complete gzip, not just a non-empty file: a library still being written is truncated
        gzip -t "$fp" 2>/dev/null || { echo "!!! $samp: $r missing or incomplete"; return 1; }
        files="$files $fp"
    done
    echo "### merge $samp $(date '+%H:%M')"
    methylpy merge-allc --allc-files $files \
        --output-file "$OUT/allc_${samp}_merged.tsv.gz" \
        --num-procs 4 --compress-output True
}

rc=0
merge dng102_mature dng102_mature_rep1 dng102_mature_rep2 dng102_mature_rep3 || rc=1
merge mdr1_mature   mdr1_mature_rep1   mdr1_mature_rep2   mdr1_mature_rep3   || rc=1
merge wt_mature     wt_mature_rep1     wt_mature_rep2     wt_mature_rep3     || rc=1
merge wt_15dap      wt_15dap_rep1      wt_15dap_rep2      wt_15dap_rep3      || rc=1
merge mdr1_15dap    mdr1_15dap_rep1    mdr1_15dap_rep2    mdr1_15dap_rep4    || rc=1
merge tassel        tassel_rep1        tassel_rep2        tassel_rep3        || rc=1
merge embryo        embryo_rep1        embryo_rep2        embryo_rep3        || rc=1
merge leaf          leaf_rep1          leaf_rep2          leaf_rep3          || rc=1
ls -la "$OUT"/allc_*_merged.tsv.gz
exit $rc
