#!/usr/bin/env bash
# Published per-gene maize tables, used as deposited (no reads are mapped for these).
#
#   B73.all.csv      Zeng et al. 2023: per-gene CDS methylation (mCG, mCHG and covered sites) and
#                    expression in ten B73 tissues (TPM, from the NAM RNA-seq of Hufford et al. 2021)
#                    https://github.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize
#                    Data/B73.all.csv
#                    -> methylation categories (02_gene_sets)
#   b73.all.6_38.txt Sun et al. 2025: the same table merged with an endosperm time course; its ten
#                    tissue columns are those of B73.all.csv
#                    https://github.com/yiruiS/endo-teM-genes  01gene_sets_identification/01.3_b73.all.6_38.txt
#                    -> tissue specificity (Fig. 1b, 1c)
#   B73_W22_maternalpreference_with_avg_rpm.csv
#                    Higgins et al. 2025 (G3): allele-resolved RPM in reciprocal B73 x W22 endosperm at
#                    11, 14, 17 and 21 DAP, as redistributed (byte-identical) at
#                    https://github.com/yiruiS/endo-teM-genes
#                    01gene_sets_identification/04.2_B73_W22_maternalpreference_with_avg_rpm.csv
#                    -> Fig. 1a and Supplemental Fig. S5
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
T=$GI_DATA/maize_tables
mkdir -p "$T"

ZENG=https://raw.githubusercontent.com/dawelab/Natural-methylation-epialleles-correlate-with-gene-expression-in-maize/main
SUN=https://raw.githubusercontent.com/yiruiS/endo-teM-genes/main/01gene_sets_identification

gi_fetch "$ZENG/Data/B73.all.csv" "$T/B73.all.csv"
gi_fetch "$SUN/01.3_b73.all.6_38.txt" "$T/b73.all.6_38.txt"
gi_fetch "$SUN/04.2_B73_W22_maternalpreference_with_avg_rpm.csv" \
    "$T/B73_W22_maternalpreference_with_avg_rpm.csv"
wc -l "$T"/*
