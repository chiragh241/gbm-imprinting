#!/usr/bin/env Rscript
# Imprinting calls: maternally (MEG) and paternally (PEG) expressed B73 genes from allele-resolved
# endosperm RNA-seq of reciprocal B73 x NAM-founder hybrids at 14 DAP.
#
# Input   data/maize/B73_NAM_syntelog_maternalpref.txt.gz  (ships with the repository)
#         Allele-specific expression of Higgins et al. 2025 (G3), mat_pref_NAM.txt.gz at
#         https://github.com/kmhiggins/Imprinting_timeseries_2024 (444,878 rows). This table holds
#         237,338 of those rows, allele values unchanged, with two columns added:
#           feature_2    the B73 v5 syntelog(s) of the row's gene: for rows aligned to the NAM
#                        founder (type A), every Zm00001eb gene in that gene's row of the MaizeGDB
#                        2020 pan-gene table, in table order; for B73-aligned rows (type B), the
#                        gene itself
#           cont_genome  genome and contrast, "<genome> _ <contrast>"
#         The kept rows are those with A_mean + B_mean >= 2 RPM, type A rows also needing a B73
#         syntelog; that rule reproduces all 237,338 rows but would add 848 more that the table
#         does not contain, so the table is shipped as it was used rather than rebuilt.
#
# Output  $GI_WORK/01_imprinting_calls/
#           b73_imprinting_calls_all.csv        one row per evaluable gene with its call
#           imprinted_3fold_genes.csv            one row per gene x hybrid x alignment, MEG/PEG
#           imprinted_3fold_genes_collapsed.csv  one row per MEG/PEG, hybrids spread wide
#
# Classification rules:
#   maternal_pref = maternal_RPM / (maternal_RPM + paternal_RPM)
#     type == "A"  -> maternal = A_mean, paternal = B_mean
#     type == "B"  -> maternal = B_mean, paternal = A_mean
#   MEG  if maternal_pref >= 6/7  (~0.857)  3-fold above 2/3 expectation
#   PEG  if maternal_pref <= 2/5  (0.4)     3-fold below 2/3 expectation
#   else biparental
#
# Fold deviation vs expected 2:1 maternal:paternal:
#   ratio               = maternal_rpm / paternal_rpm  (expected = 2)
#   maternal_fold_excess = ratio / 2                    (MEG when >=3)
#   paternal_fold_excess = 2 / ratio                    (PEG when >=3)
#   fold_vs_expected     = max(ratio/2, 2/ratio)        (>=3 == 3-fold deviation)
#   signed_fold          = +ratio/2 if maternal skew, -2/ratio if paternal skew

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
repo     <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."))
run_dir  <- Sys.getenv("GI_RUN", file.path(Sys.getenv("HOME"), "gbm_imprinting_run"))
shipped  <- Sys.getenv("GI_SHIPPED", file.path(repo, "data"))
work     <- Sys.getenv("GI_WORK", file.path(run_dir, "work"))
in_path  <- file.path(shipped, "maize", "B73_NAM_syntelog_maternalpref.txt.gz")
out_dir  <- file.path(work, "01_imprinting_calls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

MEG_CUT   <- 6 / 7
PEG_CUT   <- 2 / 5
RPM_MIN   <- 5
MIN_HYBRIDS <- 2  # at least two hybrids; the support tier (2, 3, 4 or more) is in the output
EXPR_OUTLIER_FOLD <- 3  # Step 5b: drop rows whose rpm_sum is >3x below the gene's mean rpm_sum for that type
DROP_CONTRASTS <- c("OB_MN")

classify <- function(mp) {
  dplyr::case_when(
    mp >= MEG_CUT ~ "MEG",
    mp <= PEG_CUT ~ "PEG",
    TRUE          ~ "biparental"
  )
}

message("Reading ", in_path)
raw <- read_tsv(
  in_path,
  col_types = cols(
    feature             = col_character(),
    genome              = col_character(),
    contrast            = col_character(),
    type                = col_character(),
    A_mean              = col_double(),
    B_mean              = col_double(),
    maternal_preference = col_double(),
    expression          = col_character(),
    feature_2           = col_character(),
    cont_genome         = col_character()
  )
)
message("  rows read: ", nrow(raw))

tidy <- raw %>%
  filter(!contrast %in% DROP_CONTRASTS) %>%
  separate_rows(feature_2, sep = " ") %>%
  mutate(feature_2 = str_trim(feature_2)) %>%
  filter(str_starts(feature_2, "Zm00001eb")) %>%
  mutate(
    maternal_rpm = if_else(type == "A", A_mean, B_mean),
    paternal_rpm = if_else(type == "A", B_mean, A_mean),
    rpm_sum      = maternal_rpm + paternal_rpm,
    mat_pref     = maternal_rpm / rpm_sum,
    mat_pat_ratio         = maternal_rpm / paternal_rpm,
    maternal_fold_excess  = mat_pat_ratio / 2,
    paternal_fold_excess  = 2 / mat_pat_ratio,
    fold_vs_expected      = pmax(maternal_fold_excess, paternal_fold_excess),
    signed_fold           = if_else(mat_pat_ratio >= 2,
                                    maternal_fold_excess,
                                    -paternal_fold_excess)
  ) %>%
  filter(rpm_sum >= RPM_MIN, is.finite(mat_pref))

message("  rows after rpm_sum filter: ", nrow(tidy))

# Step 5b: per (gene, alignment type) drop measurements where total expression
# (rpm_sum) is more than EXPR_OUTLIER_FOLD below that gene's mean rpm_sum. Suppresses
# hybrids where gene is poorly expressed in that NAM background (both alleles low,
# noisy mat_pref looks biparental). Does NOT punish zero-allele rows, which are
# legitimate strong MEG/PEG evidence.
tidy <- tidy %>%
  group_by(feature_2, type) %>%
  mutate(
    mean_rpm_sum      = mean(rpm_sum),
    low_expr_outlier  = rpm_sum < mean_rpm_sum / EXPR_OUTLIER_FOLD
  ) %>%
  ungroup() %>%
  filter(!low_expr_outlier) %>%
  select(-mean_rpm_sum, -low_expr_outlier) %>%
  mutate(class = classify(mat_pref))

message("  rows after expression-outlier filter: ", nrow(tidy))

# sanity check: recomputed maternal_pref should match file's maternal_preference
check <- tidy %>%
  mutate(diff = abs(mat_pref - maternal_preference)) %>%
  summarise(max_diff = max(diff, na.rm = TRUE),
            median_diff = median(diff, na.rm = TRUE))
message(sprintf("  mat_pref recompute vs file: max_diff=%.3g  median_diff=%.3g",
                check$max_diff, check$median_diff))

# Per-row detail: one row per (gene, hybrid, type alignment).
row_detail <- tidy %>%
  transmute(
    b73_gene      = feature_2,
    hybrid        = contrast,
    alignment_type = type,
    alignment_ref  = if_else(type == "A", "NAM-reference", "B73-reference"),
    nam_stock      = if_else(type == "A", genome, cont_genome),
    maternal_feature_id = feature,
    A_mean,
    B_mean,
    maternal_rpm,
    paternal_rpm,
    rpm_sum,
    mat_pref,
    mat_pat_ratio,
    fold_vs_expected,
    signed_fold,
    class
  )

# Step 6: per (B73 gene, hybrid), require both type A and type B rows and same class.
gene_hybrid <- tidy %>%
  group_by(feature_2, contrast) %>%
  summarise(
    n_types         = n_distinct(type),
    n_classes       = n_distinct(class),
    first_class     = first(class),
    mat_pref_typeA  = mat_pref[type == "A"][1],
    mat_pref_typeB  = mat_pref[type == "B"][1],
    mat_rpm_typeA   = maternal_rpm[type == "A"][1],
    pat_rpm_typeA   = paternal_rpm[type == "A"][1],
    mat_rpm_typeB   = maternal_rpm[type == "B"][1],
    pat_rpm_typeB   = paternal_rpm[type == "B"][1],
    fold_typeA      = fold_vs_expected[type == "A"][1],
    fold_typeB      = fold_vs_expected[type == "B"][1],
    signed_fold_typeA = signed_fold[type == "A"][1],
    signed_fold_typeB = signed_fold[type == "B"][1],
    min_mat_pref    = min(mat_pref),
    max_mat_pref    = max(mat_pref),
    min_fold        = min(fold_vs_expected),
    .groups     = "drop"
  ) %>%
  filter(n_types == 2L, n_classes == 1L) %>%
  transmute(
    b73_gene = feature_2,
    hybrid   = contrast,
    class    = first_class,
    mat_pref_typeA, mat_pref_typeB,
    mat_rpm_typeA, pat_rpm_typeA,
    mat_rpm_typeB, pat_rpm_typeB,
    fold_typeA, fold_typeB,
    signed_fold_typeA, signed_fold_typeB,
    min_mat_pref, max_mat_pref, min_fold
  )

message("  (gene, hybrid) pairs surviving type A/B agreement: ", nrow(gene_hybrid))

# Steps 7-8: per B73 gene, require >= MIN_HYBRIDS and unanimous class.
# Summaries use the *minimum* fold magnitude across all hybrid x type measurements
# used for the gene's call (i.e. the weakest supporting point). Also report mean + max.
gene_summary <- gene_hybrid %>%
  group_by(b73_gene) %>%
  summarise(
    n_hybrids       = n(),
    hybrids         = paste(sort(unique(hybrid)), collapse = ","),
    unique_classes  = n_distinct(class),
    consensus_class = if (n_distinct(class) == 1L) first(class) else NA_character_,
    min_mat_pref    = min(min_mat_pref),
    max_mat_pref    = max(max_mat_pref),
    mean_mat_pref   = mean(c(mat_pref_typeA, mat_pref_typeB)),
    min_fold        = min(min_fold),
    max_fold        = max(fold_typeA, fold_typeB),
    mean_fold       = mean(c(fold_typeA, fold_typeB)),
    .groups         = "drop"
  ) %>%
  mutate(
    final_call = dplyr::case_when(
      n_hybrids < MIN_HYBRIDS          ~ "insufficient_data",
      is.na(consensus_class)           ~ "unknown",
      consensus_class == "biparental"  ~ "biparental",
      TRUE                             ~ consensus_class
    ),
    # Support tier: high = >=4 hybrids, medium = 3, low = 2
    # Low-tier calls have higher false-positive risk and should be flagged downstream.
    support_tier = dplyr::case_when(
      n_hybrids >= 4 ~ "high",
      n_hybrids == 3 ~ "medium",
      n_hybrids == 2 ~ "low",
      TRUE           ~ NA_character_
    )
  )

counts <- gene_summary %>% count(final_call, name = "n_genes") %>% arrange(desc(n_genes))
message("Gene call distribution:")
print(counts)

# Build single consolidated output: one row per (gene, hybrid, alignment) for the
# MEGs and PEGs. Adds gene-level summary columns and a
# plain-English reason for the call.
fmt_num <- function(x) formatC(x, format = "fg", digits = 3)

imprinted <- gene_summary %>%
  filter(final_call %in% c("MEG", "PEG")) %>%
  mutate(
    gene_reason = sprintf(
      "%s: unanimous across %d hybrids (%s); weakest supporting measurement %s-fold, strongest %s-fold, mean %s-fold. Both NAM-reference and B73-reference alignments agree for every hybrid used. Threshold: mat_pref %s.",
      final_call,
      n_hybrids,
      hybrids,
      fmt_num(min_fold),
      fmt_num(max_fold),
      fmt_num(mean_fold),
      ifelse(final_call == "MEG", ">= 6/7", "<= 2/5")
    )
  ) %>%
  arrange(final_call, desc(min_fold))

imprinted_rows <- row_detail %>%
  # keep only (gene, hybrid) pairs that survived Step 7 (type A/B agreement)
  semi_join(gene_hybrid, by = c("b73_gene", "hybrid")) %>%
  inner_join(
    imprinted %>%
      select(b73_gene, final_call, support_tier, n_hybrids, hybrids_used = hybrids,
             gene_min_fold = min_fold, gene_max_fold = max_fold,
             gene_mean_fold = mean_fold, gene_reason),
    by = "b73_gene"
  ) %>%
  transmute(
    b73_gene,
    final_call,
    support_tier,
    n_hybrids,
    hybrids_used,
    gene_min_fold,
    gene_max_fold,
    gene_mean_fold,
    hybrid,
    alignment_type,
    alignment_ref,
    nam_stock,
    maternal_feature_id,
    A_mean,
    B_mean,
    maternal_rpm,
    paternal_rpm,
    rpm_sum,
    mat_pref,
    mat_pat_ratio,
    fold_vs_expected,
    signed_fold,
    row_class = class,
    gene_reason
  ) %>%
  arrange(final_call, desc(gene_min_fold), b73_gene, hybrid, alignment_type)

single_out <- file.path(out_dir, "imprinted_3fold_genes.csv")
write_csv(imprinted_rows, single_out)

# Collapsed version: one row per gene with per-hybrid columns spread wide.
# Uses gene_hybrid (one row per gene x hybrid with type A/B values already side-by-side),
# pivots hybrid into column names. Fixed 8 hybrids: BB, BK, BO, BW, CB, DB, EB, FB.
collapsed <- imprinted %>%
  transmute(
    b73_gene, final_call, support_tier, n_hybrids, hybrids_used = hybrids,
    min_mat_pref, max_mat_pref, mean_mat_pref,
    gene_min_fold = min_fold, gene_max_fold = max_fold, gene_mean_fold = mean_fold,
    gene_reason
  ) %>%
  left_join(
    gene_hybrid %>%
      semi_join(imprinted, by = "b73_gene") %>%
      select(b73_gene, hybrid,
             mat_pref_typeA, mat_pref_typeB,
             mat_rpm_typeA, pat_rpm_typeA,
             mat_rpm_typeB, pat_rpm_typeB,
             fold_typeA, fold_typeB) %>%
      pivot_wider(
        names_from  = hybrid,
        values_from = c(mat_pref_typeA, mat_pref_typeB,
                        mat_rpm_typeA, pat_rpm_typeA,
                        mat_rpm_typeB, pat_rpm_typeB,
                        fold_typeA, fold_typeB),
        names_glue  = "{hybrid}_{.value}"
      ),
    by = "b73_gene"
  ) %>%
  arrange(final_call, desc(gene_min_fold))

collapsed_out <- file.path(out_dir, "imprinted_3fold_genes_collapsed.csv")
write_csv(collapsed, collapsed_out)

# keep the raw all-genes call sheet as a secondary reference
all_out <- file.path(out_dir, "b73_imprinting_calls_all.csv")
write_csv(gene_summary, all_out)

message("Wrote:")
message("  ", single_out,    "  (", nrow(imprinted_rows), " rows, ",
        dplyr::n_distinct(imprinted_rows$b73_gene), " genes)")
message("  ", collapsed_out, "  (", nrow(collapsed),      " rows, ",
        ncol(collapsed), " cols)")
message("  ", all_out)
