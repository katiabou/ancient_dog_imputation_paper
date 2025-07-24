#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

########################################################
#  Plink ROH to test target imputed against validation #
########################################################

# define output files of make plink
DOCS = ["bed", "bim", "fam"]


rule transversions_validation:
    """
    Only take transversions from filtered validation data
    """
    input:
        validation_sample_filt_allelic="output/GLIMPSE_concordance/validation_bams/{sample_con}_{chrom}_validation_filt_qual_dp_ab.bcf",
    output:
        tranversion_sites="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_transversions.tsv.gz",
        validation_transversions="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions.vcf.gz",
    log:
        "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions.vcf.gz.log",
    shell:
        """
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.validation_sample_filt_allelic} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.validation_sample_filt_allelic} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.validation_transversions} 2> {log}

        bcftools index -f {output.validation_transversions}
        """


rule make_plink_transversions_validation:
    """
    Prepare file format for plink
    """
    input:
        validation_transversions="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions.vcf.gz",
    output:
        expand(
            "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_plink.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_plink",
    log:
        "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_plink.log",
    shell:
        """
        plink \
        --vcf {input.validation_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule roh_transversions_validation:
    """
    Run ROH estimation with plink
    """
    input:
        bim="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_plink.bim",
    output:
        validation_roh="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink.hom",
    params:
        prefix_in="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_plink",
        prefix_out="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink",
    log:
        "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_transversions_hom_win_het_{hom_win_het}_plink.hom.log",
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


rule merge_roh_transversions_validation:
    """
    Prepare plink ROH output for validation sample transversions
    """
    input:
        validation_roh="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink.hom",
    output:
        temp1_roh_validation=temp(
            "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink-temp1.hom"
        ),
        temp_roh_validation="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink-temp.hom",
    shell:
        """
        scp {input.validation_roh} {output.temp1_roh_validation}
        grep -qxF '{wildcards.sample_con}' {output.temp1_roh_validation} || printf "{wildcards.sample_con}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output.temp1_roh_validation}
        awk -F'\t' '{{print $0 "\t" "HC_genotyped"}}' {output.temp1_roh_validation} > {output.temp_roh_validation}
        sed -i '1s/HC_genotyped/cov/' {output.temp_roh_validation}
        """


####################################################
#### Run ROHs for transversions and transitions ####
####################################################


rule make_plink_all_sites_validation:
    """
    Prepare file format for plink
    """
    input:
        validation_sample_filt_allelic="output/GLIMPSE_concordance/validation_bams/{sample_con}_{chrom}_validation_filt_qual_dp_ab.bcf",
    output:
        expand(
            "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_plink.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_plink",
    log:
        "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_plink.log",
    shell:
        """
        plink \
        --bcf {input.validation_sample_filt_allelic} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 \
        --allow-extra-chr 2> {log}
        """


rule roh_all_sites_validation:
    """
    Run ROH estimation with plink
    """
    input:
        bim="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_plink.bim",
    output:
        validation_roh="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink.hom",
    params:
        prefix_in="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_plink",
        prefix_out="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink",
    log:
        "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink.hom.log",
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


rule merge_roh_all_sites_validation:
    """
    Prepare plink ROH output for validation sample all sites
    """
    input:
        validation_roh="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink.hom",
    output:
        temp1_roh_validation=temp(
            "output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink-temp1.hom"
        ),
        temp_roh_validation="output/GLIMPSE_concordance/ROH_validation/{sample_con}_{chrom}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink-temp.hom",
    shell:
        """
        scp {input.validation_roh} {output.temp1_roh_validation}
        grep -qxF '{wildcards.sample_con}' {output.temp1_roh_validation} || printf "{wildcards.sample_con}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output.temp1_roh_validation}
        awk -F'\t' '{{print $0 "\t" "HC_genotyped"}}' {output.temp1_roh_validation} > {output.temp_roh_validation}
        sed -i '1s/HC_genotyped/cov/' {output.temp_roh_validation}
        """
