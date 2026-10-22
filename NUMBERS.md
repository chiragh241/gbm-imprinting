# Where each number in the text comes from

Every quantity the Results quote, the script that computes it and the file it is written to.
Paths are under `$GI_RESULTS` (`tables/`, `figures/`) or `$GI_WORK`; scripts are under `scripts/`.
"points" are percentage points. Where a figure plots a difference with the opposite sign to the
text (paternal minus maternal against maternal minus paternal, or wild type minus mutant against
mutant minus wild type), the file holds the figure's sign.

## Identification of imprinted genes in maize

| Quantity | Script | Output |
|---|---|---|
| MEGs, PEGs and biparental genes, and genes with an imprinting call | `01_imprinting_calls/01_identify_imprinted.R` | `$GI_WORK/01_imprinting_calls/b73_imprinting_calls_all.csv` (`final_call`); counts printed in the log |

## Gene body methylation and expression of PEGs

| Quantity | Script | Output |
|---|---|---|
| PEGs, MEGs and biparental genes in each methylation category | `02_gene_sets/03_build_sets.py` | `$GI_WORK/02_gene_sets/counts_grid.csv` |
| teM and gbM enrichments: percentages, odds ratios, Fisher p | `02_gene_sets/04_enrichment_stats.py` | `tables/enrichment_stats.csv` |
| Median Tau of gbM-PEGs and UM-PEGs, Mann-Whitney p | `12_figures/01_fig1.py` | `tables/fig1_tau_source_data.csv`; `figures/fig1_legend_numbers.txt` |
| Median endosperm specificity of gbM-PEGs and UM-PEGs, Mann-Whitney p | `03_expression_specificity/01_endosperm_specificity.py` | `tables/endosperm_specificity_summary.csv` |
| Leaf and ear H3K27me3: fraction of UM-PEGs and gbM-PEGs enriched twofold; Mann-Whitney p, UM-PEGs against gbM-PEGs, on log2 enrichment | `09_chip_leaf_ear/05_gene_enrichment.py` | `tables/figS1_gene_enrichment.csv` |

## Demethylation is not a consequence of maternal repression

| Quantity | Script | Output |
|---|---|---|
| Median CG methylation of gbM-PEGs in embryo and endosperm; per-gene deltas and IQRs, endosperm against embryo, leaf and tassel, for gbM-PEGs and gbM-biparentals; the shuffled control, embryo against endosperm | `06_per_gene_methylation/05_tissue_deltas.py` | `tables/tissue_deltas.csv` |
| Allele-resolved paternal and maternal delta, IQR, n and Wilcoxon p for gbM-PEGs and gbM-biparentals | `12_figures/03_fig2.py` | `figures/fig2_legend_numbers.txt` (paternal minus maternal) |
| Relative maternal reduction of gbM-biparentals, 1 - mean maternal / mean paternal mCG (also the Introduction's figure) | `14_relative_reduction/01_relative_reduction.py` | `tables/relative_reduction.csv` (`reduction_pct`, maize) |

## Cis-regulatory element demethylation

| Quantity | Script | Output |
|---|---|---|
| Fraction of MEGs, UM-PEGs and gbM-PEGs with a DMR within 2 kb; fold enrichments and Fisher p against biparentals and MEGs | `07_dmr_distances/03_distance_stats.py` | `tables/dmr_within_2kb_tests.csv` |

## DNGs mediate maternal demethylation and repression

| Quantity | Script | Output |
|---|---|---|
| *mdr1* (15 DAP and mature) and *dng102* (mature) minus wild type: per-gene delta and IQR for gbM-PEGs and gbM-biparentals | `06_per_gene_methylation/02_mutant_deltas.py` | `tables/mutant_deltas_per_gene.csv` |
| The same for the shuffled control regions | `06_per_gene_methylation/04_control_deltas_summary.py` | `tables/control_deltas.csv` |
| *mdr1 dng103*: per-gene change in maternal expression ratio (median, IQR, Wilcoxon p) for each gene class, mutant as mother and as father | `08_rnaseq_allele/06_delta_by_class.py` | `tables/xu_delta_by_methylation_class.tsv`; also `figures/fig3_legend_numbers.txt` |
| Maternal expression ratio of gbM, UM and teM biparental genes at four stages | `12_figures/13_figS4.py` | `tables/figS_matpref_source.csv` |

## Gene body demethylation and H3K27me3

| Quantity | Script | Output |
|---|---|---|
| Per-gene maternal minus paternal H3K27me3 enrichment (median, IQR) for each gene group | `10_chip_endosperm_allele/07_compact_profile.py` | printed in the log; per-gene values in `$GI_WORK/10_chip_endosperm_allele/matpat_per_gene_body_compact.csv` |
| Wilcoxon and Mann-Whitney p, raw and Holm-corrected | `12_figures/06_fig4.py` | `tables/fig4_pergene_tests.csv` |
| The CG delta threshold defining strongly demethylated gbM-biparentals | `10_chip_endosperm_allele/03_allele_panel.py` | printed in the log (lower quartile of maternal minus paternal mCG) |

## Rice

| Quantity | Script | Output |
|---|---|---|
| Persistent PEGs and MEGs (245, 95), those without leaf methylation (MSU7 transposable element related), those with a category and in each (UM, gbM, teM); gbM enrichment in PEGs over MEGs (fold, Fisher p) | `13_rice/06_gene_sets.py` | `tables/rice_imprinting_counts.csv`; per gene `tables/rice_persistent_imprinted_genes.csv` |
| Biparental genes in each methylation category | `13_rice/06_gene_sets.py` | `tables/rice_gene_sets.csv` |
| Per-gene paternal minus maternal CG methylation (n, median, IQR, Wilcoxon p) for each set; Mann-Whitney p between sets | `13_rice/07_gene_matpat.py` | `tables/rice_matpat_summary.csv`, `tables/rice_matpat_tests.csv`; also `figures/rice_matpat_legend_numbers.txt` |
| Relative maternal reduction of gbM-biparentals | `14_relative_reduction/01_relative_reduction.py` | `tables/relative_reduction.csv` (`reduction_pct`, rice) |
| Kasalath and 93-11 SNPs (Methods) | `13_rice/03_call_snps.sh` | `$GI_WORK/13_rice/snps/snp_counts.tsv` |
| Bisulfite conversion from chloroplast CHH (Methods) | `13_rice/09_qc_table.py` | `tables/rice_qc_vs_tonosaki.csv` (`conversion_chh`) |

## Arabidopsis

| Quantity | Script | Output |
|---|---|---|
| Relative maternal reduction of gbM-biparentals | `14_relative_reduction/01_relative_reduction.py` | `tables/relative_reduction.csv` (`reduction_pct`, arabidopsis) |
| Consensus PEGs and MEGs | `11_arabidopsis/03_consensus_sets.py` | `tables/consensus_pegs.csv`, `tables/consensus_megs.csv` |
| Consensus PEGs and MEGs in each methylation category; gbM enrichment in PEGs over MEGs (fold, odds ratio, Fisher p) | `11_arabidopsis/11_imprinted_class_enrichment.py` | `tables/arabidopsis_imprinted_classes.csv`, `tables/arabidopsis_gbm_peg_vs_meg.csv` |
| Biparental genes in each methylation category | `11_arabidopsis/10_biparental_class_counts.py` | `tables/arabidopsis_biparental_class_counts.csv` |
| Paternal minus maternal CG methylation in wild-type, *dme* and *fie* endosperm (median, IQR, p) | `12_figures/09_fig5cde.py` | `figures/fig5_legend_numbers.txt`; `tables/fig5cde_source.csv` |
