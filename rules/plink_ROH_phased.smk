#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

#######################################
#  Plink ROH for full imputed dataset #
#######################################

global CHROM, CANID_SUBSET

# define output files of make plink
DOCS = ["bed", "bim", "fam"]


rule transversions_phased_subset:
    """
    Only take transversions from recalibrated phased data 
    Doing this for dogs and wolves seperately
    """
    input:
        #merged_phased_vcf_maf_info = "output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info}.vcf.gz",
        merged_phased_vcf_maf_recalibrated_info="output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz",
        imputed_canid_subset="sample_lists/names_imputed_{canid_subset}.tsv",
    output:
        #tranversion_sites = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.tsv.gz",
        #phased_transversions = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.vcf.gz"
        tranversion_sites="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.tsv.gz",
        phased_transversions="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.vcf.gz",
    log:
        #"output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.log"
        "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.log",
    threads: 4
    shell:
        """
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.merged_phased_vcf_maf_recalibrated_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.merged_phased_vcf_maf_recalibrated_info} \
        -S {input.imputed_canid_subset} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.phased_transversions} 2> {log}

        bcftools index -f {output.phased_transversions}
        """


rule make_plink_transversions_phased:
    """
    Prepare file format for plink
    """
    input:
        #phased_transversions = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.vcf.gz"
        phased_transversions="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.vcf.gz",
    output:
        #expand("output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.{doc}", doc=DOCS, allow_missing=True)
        expand(
            "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        #prefix = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}"
        prefix="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}",
    log:
        #"output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.log"
        "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.log",
    shell:
        """
        plink \
        --vcf {input.phased_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule roh_transversions_phased:
    """
    Run ROH estimation with plink
    """
    input:
        #bim = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.bim"
        bim="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}.bim",
    output:
        #phased_roh = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        #phased_roh_sum = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
        phased_roh="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_roh_sum="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    params:
        #prefix_in = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}",
        #prefix_out = "output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}"
        prefix_in="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_{canid_subset}",
        prefix_out="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}",
    log:
        #"output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.log"
        "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.log",
    shell:
        """
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out} 2> {log}
        """


rule merge_chrom_ROH_phased_transversions:
    """
    TODO add block header
    """
    input:
        #phased_roh = expand("output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom", chrom=CHROM, allow_missing=True),
        #phased_roh_sum = expand('output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
        phased_roh=expand(
            "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom",
            chrom=CHROM,
            allow_missing=True,
        ),
        phased_roh_sum=expand(
            "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        #phased_roh_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        #phased_roh_sum_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        phased_roh_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_roh_sum_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    shell:
        """
        awk 'FNR>1 || NR==1' {input.phased_roh} > {output.phased_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.phased_roh_sum} > {output.phased_roh_sum_allchrom}
        """


####################################################
#### Run ROHs for transversions and transitions ####
####################################################


rule all_sites_phased_subset:
    """
    Take all sites from recalibrated phased data
    Doing this for dogs and wolves seperately
    """
    input:
        #merged_phased_vcf_maf_info = 'output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info}.vcf.gz',
        merged_phased_vcf_maf_recalibrated_info="output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz",
        imputed_canid_subset="sample_lists/names_imputed_{canid_subset}.tsv",
    output:
        #phased_all_sites_subset = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz'
        phased_all_sites_subset="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.vcf.gz",
    log:
        #'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.log'
        "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.log",
    threads: 4
    shell:
        """
        bcftools view {input.merged_phased_vcf_maf_recalibrated_info} \
        -S {input.imputed_canid_subset} \
        --threads {threads} \
        -Oz -o {output.phased_all_sites_subset} 2> {log}

        bcftools index -f {output.phased_all_sites_subset}
        """


rule make_plink_all_sites_phased:
    """
    Prepare file format for plink
    """
    input:
        #phased_all_sites_subset = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz'
        phased_all_sites_subset="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.vcf.gz",
    output:
        #expand('output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
        expand(
            "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        #prefix = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}'
        prefix="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}",
    log:
        #'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.log'
        "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.log",
    shell:
        """
        plink \
        --vcf {input.phased_all_sites_subset} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule roh_all_sites_phased:
    """
    Run ROH estimation with plink
    """
    input:
        #bim = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.bim'
        bim="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.bim",
    output:
        #phased_roh = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        #phased_roh_ind = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        #phased_roh_sum = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary'
        phased_roh="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_roh_ind="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
        phased_roh_sum="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    params:
        #prefix_in = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}',
        #prefix_out = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}'
        prefix_in="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}",
        prefix_out="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}",
    log:
        #'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
        "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log",
    shell:
        """
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out} 2> {log}
        """


rule merge_chrom_ROH_phased:
    """
    TODO add block header
    """
    input:
        #phased_roh = expand('output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
        #phased_roh_sum = expand('output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
        phased_roh=expand(
            "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
            chrom=CHROM,
            allow_missing=True,
        ),
        phased_roh_sum=expand(
            "output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        #phased_roh_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        #phased_roh_sum_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        phased_roh_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_roh_sum_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    shell:
        """
        awk 'FNR>1 || NR==1' {input.phased_roh} > {output.phased_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.phased_roh_sum} > {output.phased_roh_sum_allchrom}
        """


####################################################
#### Run ROHs for transversions reference panel ####
####################################################


rule transversions_ref_panel_subset:
    """
    Only take transversions from ref panel
    """
    input:
        #ref_sample_snp_filltags_filter = 'output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        ref_sample_snp_filltags_filter_maf_vcf="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.phased.vcf.gz",
        modern_canid_subset="sample_lists/ref_panel_filt_{canid_subset}.tsv",
    output:
        tranversion_sites="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.tsv.gz",
        ref_panel_transversions="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.vcf.gz",
    log:
        "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.log",
    threads: 4
    shell:
        """
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.ref_sample_snp_filltags_filter_maf_vcf} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.ref_sample_snp_filltags_filter_maf_vcf} \
        -S {input.modern_canid_subset} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.ref_panel_transversions} 2> {log}

        bcftools index -f {output.ref_panel_transversions}
        """


rule make_plink_transversions_ref_panel:
    """
    Prepare file format for plink
    """
    input:
        ref_panel_transversions="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.vcf.gz",
    output:
        expand(
            "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}",
    log:
        "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.log",
    shell:
        """
        plink \
        --vcf {input.ref_panel_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule roh_transversions_ref_panel:
    """
    Run ROH estimation with plink
    """
    input:
        bim="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.bim",
    output:
        phased_roh="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_roh_sum="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
        phased_roh_ind="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
    params:
        prefix_in="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}",
        prefix_out="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}",
    log:
        "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom.log",
    shell:
        """
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out} 2> {log}
        """


rule merge_chrom_ROH_transversions_ref_panel:
    """
    TODO add block header
    """
    input:
        modern_roh=expand(
            "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom",
    shell:
        """
        awk 'FNR>1 || NR==1' {input.modern_roh} > {output.modern_roh_allchrom}
        """


################################################
#### Run ROHs for all sites reference panel ####
################################################


rule all_sites_ref_panel_subset:
    """
    Take all sites from ref panel
    """
    input:
        #ref_sample_snp_filltags_filter = 'output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        ref_sample_snp_filltags_filter_maf_vcf="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.phased.vcf.gz",
        modern_canid_subset="sample_lists/ref_panel_filt_{canid_subset}.tsv",
    output:
        ref_panel_all_sites_subset="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.vcf.gz",
    log:
        "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.log",
    threads: 4
    shell:
        """
        bcftools view {input.ref_sample_snp_filltags_filter_maf_vcf} \
        -S {input.modern_canid_subset} \
        --threads {threads} \
        -Oz -o {output.ref_panel_all_sites_subset} 2> {log}

        bcftools index -f {output.ref_panel_all_sites_subset}
        """


rule make_plink_all_sites_ref_panel:
    """
    Prepare file format for plink
    """
    input:
        ref_panel_all_sites_subset="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.vcf.gz",
    output:
        expand(
            "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}",
    log:
        "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.log",
    shell:
        """
        plink \
        --vcf {input.ref_panel_all_sites_subset} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule roh_all_sites_ref_panel:
    """
    Run ROH estimation with plink
    """
    input:
        bim="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.bim",
    output:
        phased_roh="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_roh_sum="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
        phased_roh_ind="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
    params:
        prefix_in="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}",
        prefix_out="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}",
    log:
        "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log",
    shell:
        """
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out} 2> {log}
        """


rule merge_chrom_ROH_ref_panel:
    """
    TODO add block header
    """
    input:
        modern_roh=expand(
            "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
            chrom=CHROM,
            allow_missing=True,
        ),
        modern_roh_sum=expand(
            "output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        modern_roh_allchrom_sum="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    shell:
        """
        awk 'FNR>1 || NR==1' {input.modern_roh} > {output.modern_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.modern_roh_sum} > {output.modern_roh_allchrom_sum}
        """


###########################################################################################################################################
#### Estimating ROHS for merged modern and ancient (since we want the .hom.summary file which is estimated based on the input file samples)
#### This is for the heatmap for all dogs and all wolves (seperately)
#### Doing this only for all sites
###########################################################################################################################################


rule merge_phased_modern_all_sites:
    """
    Merge phased recalibrated with selected samples from reference panel
    """
    input:
        #phased_all_sites_subset = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz',
        phased_all_sites_subset="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.vcf.gz",
        ref_panel_all_sites_subset="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.vcf.gz",
    output:
        #phased_modern_all_sites_subset = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz',
        phased_modern_all_sites_subset="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.vcf.gz",
    log:
        #'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz.log'
        "output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.vcf.gz.log",
    threads: 4
    shell:
        """
        bcftools merge \
        {input.phased_all_sites_subset} {input.ref_panel_all_sites_subset} \
        --threads {threads} \
        -Oz -o {output.phased_modern_all_sites_subset} 2> {log}

        bcftools index -f {output.phased_modern_all_sites_subset}
        """


rule make_plink_all_sites_phased_ref_panel:
    """
    Convert to plink
    """
    input:
        #phased_modern_all_sites_subset = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz',
        phased_modern_all_sites_subset="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.vcf.gz",
    output:
        #expand('output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
        expand(
            "output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        #prefix = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}'
        prefix="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}",
    log:
        #'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.plink.log'
        "output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.plink.log",
    shell:
        """
        plink \
        --vcf {input.phased_modern_all_sites_subset} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule roh_all_sites_phased_ref_panel:
    """
    Run ROH estimation with plink
    """
    input:
        #bim = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.bim'
        bim="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}.bim",
    output:
        #phased_modern_roh = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        #phased_modern_roh_ind = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        #phased_modern_roh_summary = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary'
        phased_modern_roh="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_modern_roh_ind="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
        phased_modern_roh_summary="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    params:
        #prefix_in = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}',
        #prefix_out = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}'
        prefix_in="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_{canid_subset}",
        prefix_out="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}",
    log:
        #'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
        "output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log",
    shell:
        """
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out} 2> {log}
        """


rule merge_chrom_ROH_phased_ref_panel:
    """
    Merge chromosomes
    """
    input:
        #phased_modern_roh = expand('output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
        #phased_modern_roh_sum = expand('output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
        phased_modern_roh=expand(
            "output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
            chrom=CHROM,
            allow_missing=True,
        ),
        phased_modern_roh_sum=expand(
            "output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        #phased_modern_roh_allchrom = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        #phased_modern_roh_sum_allchrom_sum = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary'
        phased_modern_roh_allchrom="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom",
        phased_modern_roh_sum_allchrom_sum="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
    shell:
        """
        awk 'FNR>1 || NR==1' {input.phased_modern_roh} > {output.phased_modern_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.phased_modern_roh_sum} > {output.phased_modern_roh_sum_allchrom_sum}
        """


##################################################################################
#### Steps for plotting all ROH for imputed and selected modern for all sites ####
##################################################################################


rule estimate_chrom_size_per_chr:
    """
    Estimate chromosome sizes for plotting x axis 
    """
    input:
        ref_fasta_chr="output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}.fasta",
    output:
        ref_fasta_size="output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}_size.genome",
    shell:
        """
        faidx {input.ref_fasta_chr} -i chromsizes > {output.ref_fasta_size}
        """


rule estimate_allchrom_size:
    """
    Estimate all autosome size 
    """
    input:
        ref_fasta_size=expand(
            "output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}_size.genome",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        ref_fasta_allchr_size="output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome",
    shell:
        """
        cat {input.ref_fasta_size} | awk '{{Total=Total+$2}} END{{print "allchrom " Total}}' > {output.ref_fasta_allchr_size}
        """


rule plot_transversions_ROH_all_chr:
    """
    Using only the dogwolf files for these plots (contains all imputed samples)
    """
    input:
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        #phased_roh_allchrom_tranversions_merged = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        phased_roh_sum_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_panel_roh_allchrom_transversions_merged = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_dogwolf.hom",
    output:
        #plot_dogs = 'output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_only.png',
        #plot_wolves = 'output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_only.png'
        plot_dogs="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_only.png",
        plot_wolves="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_only.png",
    params:
        sites="transversions",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_bands_all_chr.R"


rule plot_all_sites_ROH_all_chr:
    """
    Using only the dogwolf files for these plots (contains all imputed samples)
    """
    input:
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        #phased_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        phased_roh_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_panel_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom",
    output:
        plot_dogs="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_only.png",
        plot_wolves="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_only.png",
    params:
        sites="all sites",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_bands_all_chr.R"


rule plot_all_sites_ROH_count_length_dogs:
    """
    TODO add block header
    """
    input:
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        #phased_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        phased_roh_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_panel_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_fasta_size = 'output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size="output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome",
    output:
        plot_dogs_length_count="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs.png",
        plot_dogs_length_count_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs-labelled.png",
        plot_dogs_coeff="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs.png",
        plot_dogs_coeff_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs-labelled.png",
        plot_dogs_length_count_long="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs.png",
        plot_dogs_length_count_long_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs-labelled.png",
        plot_dogs_coeff_long="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs.png",
        plot_dogs_coeff_long_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs-labelled.png",
        plot_dogs_length_count_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs.png",
        plot_dogs_length_count_short_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs-labelled.png",
        plot_dogs_coeff_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs.png",
        plot_dogs_coeff_short_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs-labelled.png",
        froh_test="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_tests.tsv",
        roh_results_all="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_all_roh_results.tsv",
        plot_dogs_coeff_map_main="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs_map.png",
        plot_dogs_coeff_long_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs.png",
        plot_dogs_coeff_boxplot="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs_boxplot.png",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_count_length_coeff_dogs.R"


rule plot_transversions_ROH_count_length_dogs:
    """
    TODO add block header
    """
    input:
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        #phased_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        phased_roh_sum_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_panel_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_fasta_size = 'output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size="output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome",
    output:
        plot_dogs_length_count="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs.png",
        plot_dogs_length_count_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs-labelled.png",
        plot_dogs_coeff="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs.png",
        plot_dogs_coeff_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs-labelled.png",
        plot_dogs_length_count_long="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs.png",
        plot_dogs_length_count_long_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs-labelled.png",
        plot_dogs_coeff_long="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs.png",
        plot_dogs_coeff_long_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs-labelled.png",
        plot_dogs_length_count_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs.png",
        plot_dogs_length_count_short_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs-labelled.png",
        plot_dogs_coeff_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs.png",
        plot_dogs_coeff_short_laebelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs-labelled.png",
        froh_test_dogs="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_tests.tsv",
        roh_results_all="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_all_roh_results.tsv",
        plot_dogs_coeff_map_main="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs_map.png",
        plot_dogs_coeff_long_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs.png",
        plot_dogs_coeff_boxplot="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs_boxplot.png",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_count_length_coeff_dogs.R"


rule plot_all_sites_ROH_count_length_wolves:
    """
    TODO add block header
    """
    input:
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        #phased_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        phased_roh_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_panel_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_fasta_size = 'output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size="output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome",
    output:
        plot_wolves_length_count="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_all_ROHs.png",
        plot_wolves_length_count_long="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_long_ROHs.png",
        plot_wolves_length_count_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_short_ROHs.png",
        plot_wolves_coeff="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs.png",
        plot_wolves_coeff_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs-labelled.png",
        #plot_wolves_coeff_boxplot_modern = 'output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_boxplot.png',
        plot_wolves_coeff_boxplot_modern_ancient="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_ancient_boxplot.png",
        plot_pleistocene_wolves_coeff="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene.png",
        plot_pleistocene_wolves_coeff_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene-labelled.png",
        roh_results_all="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_all_roh_results.tsv",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_count_length_coeff_wolves.R"


rule plot_transversions_ROH_count_length_wolves:
    """
    TODO add block header
    """
    input:
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        #phased_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        phased_roh_sum_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_panel_roh_allchrom_all_sites_merged = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_dogwolf.hom",
        #ref_fasta_size = 'output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size="output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome",
    output:
        plot_wolves_length_count="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_count_length_all_ROHs.png",
        plot_wolves_length_count_long="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_count_length_long_ROHs.png",
        plot_wolves_length_count_short="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_count_length_short_ROHs.png",
        plot_wolves_coeff="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs.png",
        plot_wolves_coeff_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs-labelled.png",
        #plot_wolves_coeff_boxplot_modern = 'output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_boxplot.png',
        plot_wolves_coeff_boxplot_modern_ancient="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_ancient_boxplot.png",
        plot_pleistocene_wolves_coeff="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene.png",
        plot_pleistocene_wolves_coeff_labelled="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene-labelled.png",
        roh_results_all="output/GLIMPSE_imputation/plots/ROHs/merged_phased_annotated_ref_panel.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_all_roh_results.tsv",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_count_length_coeff_wolves.R"


##################################################################################
#### Steps for plotting all ROH for imputed and selected modern for all sites ####
##################################################################################


rule estimate_ROH_windows:
    """
    Get file with 500KB windows, to estimate coverage after. These will be the same, since the genome size is the same for all samples. 
    I'll further estimate coverages using either only the ancient wolf bams or only the ancient dog bams, or both dogs and wolves (seperate window coverage masks)
    """
    input:
        ref_fasta_size="output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}_size.genome",
    output:
        ROH_500_kb_windows="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/CanFam31_{chrom}_500_kb_windows.txt",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_500_kb_window.R"


rule bamlist_subset:
    """
    Create bam list for either dogs or wolves or dogs/wolves
    """
    input:
        bam_meta="sample_lists/bams_published_imputation_metadata_cutoff.tsv",
        imputed_canid_subset="sample_lists/names_imputed_{canid_subset}.tsv",
    output:
        bam_list_subset="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/bam_list_{canid_subset}_window_cov.txt",
    shell:
        """
        awk -F '\t' 'NR==FNR{{a[$1]; next}} FNR==1 || $1 in a' {input.imputed_canid_subset} {input.bam_meta} | awk -F'\t' 'FNR>1{{print $3}}' > {output.bam_list_subset}
        """


rule window_coverage_estimate:
    """
    Esimate window coverage for imputed bams
    """
    input:
        bam_list_subset="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/bam_list_{canid_subset}_window_cov.txt",
        ROH_500_kb_windows="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/CanFam31_{chrom}_500_kb_windows.txt",
    output:
        windows_cov_subset="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_{chrom}_windows_cov_500kb.txt",
    shell:
        """
        while read -r bam_line
        do
        echo $bam_line > temp_{wildcards.chrom}_{wildcards.canid_subset}.txt 
        while read -r line
        do
        samtools coverage -r $line --bam-list temp_{wildcards.chrom}_{wildcards.canid_subset}.txt    | grep -v '#' >> {output.windows_cov_subset}
        done < {input.ROH_500_kb_windows}
        done < {input.bam_list_subset}

        rm temp_{wildcards.chrom}_{wildcards.canid_subset}.txt 
        """


rule merge_chrom_window_coverage_estimate:
    """
    Merge all chromosomes 
    """
    input:
        windows_cov_subset=expand(
            "output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_{chrom}_windows_cov_500kb.txt",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        windows_cov_subset_allchrom="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_allchrom_windows_cov_500kb.txt",
    shell:
        """
        cat {input.windows_cov_subset} >> {output.windows_cov_subset_allchrom}
        """


rule plot_ROH_window_prevelance_depth:
    """
    Plotting all dogs, all wolves, and dogs/wolves
    Getting overlapping windows
    """
    input:
        #phased_roh_sum_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        #phased_roh_sum_ind = 'output/GLIMPSE_imputation/ROH_phased/merged_phased.chr1_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        phased_roh_sum_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
        phased_roh_allchrom="output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.chr1_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
        modern_roh_allchrom_sum="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
        modern_roh_ind="output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_chr1_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
        #phased_modern_roh_sum_allchrom_sum = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        #phased_modern_roh_sum_ind = 'output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.chr1_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        phased_modern_roh_sum_allchrom_sum="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary",
        phased_modern_roh_allchrom="output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_annotated_modern.chr1_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv",
        windows_cov_subset_allchrom="output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_allchrom_windows_cov_500kb.txt",
    output:
        windows_prev_depth_plot="output/GLIMPSE_imputation/plots/ROH_islands_deserts/windows_prevelance_cov_500kb_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        heatmap_imputed_modern="output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_modern_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        heatmap_density_imputed_modern="output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_modern_heatmap_density_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        heatmap_imputed="output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        heatmap_density_imputed="output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_heatmap_density_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        heatmap_modern="output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        heatmap_density_modern="output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_heatmap_density_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png",
        ROH_windows_modern_imputed="output/GLIMPSE_imputation/ROH_islands_deserts/modern_imputed_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.txt",
        ROH_windows_imputed="output/GLIMPSE_imputation/ROH_islands_deserts/imputed_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.txt",
        ROH_windows_modern="output/GLIMPSE_imputation/ROH_islands_deserts/modern_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.txt",
        imputed_windows_bed_islands="output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed",
        imputed_windows_bed_deserts="output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed",
        modern_windows_bed_islands="output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed",
        modern_windows_bed_deserts="output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed",
        imputed_modern_windows_bed_islands="output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed",
        imputed_modern_windows_bed_deserts="output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed",
        top_go_terms="output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts_GO_terms.txt",
        main_figure="output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_ancient_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_main.png",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/ROH_deserts_islands.R"
