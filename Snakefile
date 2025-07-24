#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

import pandas as pd


##### load config file #####
configfile: "config.yaml"


##### Define total chromosome numbers #####
CHROM = [f"chr{i}" for i in range(1, config["chromosome_number"] + 1)]


##### set wildcard constraints #####
wildcard_constraints:
    chrom="chr\d+",
    info_cutoff="\d+.\d+",
    maf="\d+.\d+",
    maf_cutoff="\d+.\d+",
    info="\d+.\d+",
    canid_subset="\w+",
    cov_cutoff="\d+.\d+",
    hom_win_het="\d+",

# Bams used for the validation
samples_con_df = pd.read_table(
    "sample_lists/concordance_bams_published.tsv", dtype=str, delimiter="\t"
).set_index("Original_ID", drop=False)
SAMPLE_CON = list(samples_con_df["Original_ID"])
BAM_CON = list(samples_con_df["bam_path"])

#Bams which will be imputed (might need to be updated from the online spreadsheet)
samples_df = pd.read_table(
    "sample_lists/bams_published_imputation_metadata_cutoff.tsv", dtype=str,delimiter="\t",
).set_index("Sample", drop=False)
SAMPLE = list(samples_df["Sample"])
BAM = list(samples_df["bam_path"])

# Define coverage values to downsample target bams to
COVERAGE_VAL = ["2", "1", "0.5", "0.2", "0.1", "0.05"]

# Define INFO score cutoffs for filtering imputed sites
INFO_CUTOFF = ["0.0", "0.8", "0.9", "0.95"]

# datasets to use
CANID_SUBSET = ["dogs", "wolves", "dogwolf"]


##### Rules to include #####
include: "rules/genetic_map.smk"
include: "rules/fetch_data.smk"
include: "rules/ref_panel.smk"
include: "rules/ref_panel_concordance_only_dogs.smk"
include: "rules/GLIMPSE_concordance.smk"
include: "rules/GLIMPSE_concordance_only_dogs.smk"
include: "rules/GLIMPSE_concordance_HC.smk"
include: "rules/GLIMPSE_concordance_transversions.smk"
include: "rules/GLIMPSE_concordance_estimate_sites_per_MAF_bin.smk"
include: "rules/filter_MAF_INFO_concordance.smk"
include: "rules/validation_filtering_pseudohaploid.smk"
include: "rules/smartpca_concordance_pseudohaploid_HC_genotyped.smk"
include: "rules/smartpca_concordance_pseudohaploid_no_filters.smk"
include: "rules/plink_ROH_validation.smk"
include: "rules/plink_ROH_concordance.smk"
include: "rules/GLIMPSE_impute.smk"
include: "rules/filter_MAF_INFO_imputation.smk"
include: "rules/smartpca_imputed.smk"
include: "rules/plink_ROH_phased.smk"
include: "rules/GO_analysis.smk"


rule all:
    input:
        #Full datasets results
        expand(
            [
            "output/GLIMPSE_imputation/plots/PCA/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_P1P2.png",
            "output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_only.png",
            "output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_only.png",
            "output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs.png",
            "output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs_map.png",
            "output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_all_ROHs.png",
            "output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_ancient_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_main.png",
            "output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/cnv_window_freq_density_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs.png",
            "output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/ROH_prevelance_cnv_freq_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs.png",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_GO_terms.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_candidate_genes.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_GO_terms_noDLA.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_candidate_genes_noDLA.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/cnv_window_freq_density_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves.png",
            "output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/ROH_prevelance_cnv_freq_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves.png",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_GO_terms.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_candidate_genes.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_GO_terms_noDLA.txt",
            "output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_candidate_genes_noDLA.txt",
            ],
            info=config["info_cutoff"],
            maf_cutoff=config["maf_cutoff"],
            canid_subset=CANID_SUBSET,
            hom_win_het=config["roh"]["homozyg_window_het"],
            chrom=CHROM
        ),
        #Benchmarking plots:
        expand(
            ["output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_allchrom_filtered-main.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_allchrom_filtered-{sample_con}.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_allchrom_filtered_dogs_full.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_allchrom_filtered_wolves_full.png",
            "output/GLIMPSE_concordance_only_dogs/plots/glimpse_concordance/concordance_allchrom_0.5x_1x_filtered-all_sites_only_dogs.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_dogs_full_transversions.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_wolves_full_transversions.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_dogs_full_all_sites_transversions_comparison.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_wolves_full_all_sites_transversions_comparison.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_0.5x_1x_filtered-all_sites_transversions.png",
            "output/GLIMPSE_concordance/plots/glimpse_concordance_MAF_bins_reference_panel/phased_annotated.NGDG_CGG32_allchrom_INFO_all.png",
            "output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample_con}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.png",
            "output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample_con}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_accuracy_HC.png",
            "output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample_con}_allchrom_corr_{canid_subset}.png",
            "output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample_con}_allchrom_corr_{canid_subset}_accuracy_HC.png",
            "output/GLIMPSE_concordance/plots/ROH_concordance/{sample_con}_allchrom_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_band_accuracy.png",
            "output/GLIMPSE_concordance/plots/ROH_concordance/{sample_con}_allchrom_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_band_accuracy.png",
            "output/GLIMPSE_concordance/plots/ROH_concordance_ROHan/{sample_con}_allchrom_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy_imputed_ROHan.png",
            ],
            sample_con=SAMPLE_CON,
            coverage_val=COVERAGE_VAL,
            info=config["info_cutoff"],
            maf=config["maf_cutoff"],
            hom_win_het=config["roh"]["homozyg_window_het"],
            canid_subset=CANID_SUBSET,
        ),
       #Publication figures
        expand(
            ["output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_allchrom_filtered-main.pdf",
            "output/GLIMPSE_concordance/plots/ROH_concordance/NGDG_allchrom_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_band_accuracy.pdf",
            "output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs_map.pdf",
            "output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_ancient_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_main.pdf", 
            ],
            info=config["info_cutoff"], 
            maf=config["maf_cutoff"], 
            hom_win_het=config["roh"]["homozyg_window_het"],
            maf_cutoff=config["maf_cutoff"],
            canid_subset=CANID_SUBSET,
        )