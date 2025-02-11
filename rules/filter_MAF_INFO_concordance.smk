#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

###############################################################################################################
#  Filter imputed downsampled and HC datasets based on MAF of reference panel, INFO score and low cov samples #
###############################################################################################################

# define output files of make plink
DOCS = ["bed", "bim", "fam"]

### Downsampled imputed ###


rule MAF_sites_ref_pan:
    """
    Extract sites using a MAF filter from the reference panel 
    """
    input:
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
    output:
        ref_concordance_sample_excl_filltags_filter_maf_vcf="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz",
    params:
        maf=config["maf_cutoff"],
    log:
        "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.bcf.log",
    shell:
        """
        bcftools view \
        -q {params.maf}:minor \
        {input.ref_concordance_sample_excl_filltags_filter} \
        -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf} 2> {log}
        
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv}

        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv}
        """


rule maf_INFO_sites_concordance_phased:
    """
    Extract MAF and INFO filtered sites from phased VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz",
        phased_vcf_annotate="output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_{coverage_val}x.vcf.gz",
    output:
        phased_maf_info="output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz",
    params:
        info=config["info_cutoff"],
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz.log",
    threads: 4
    shell:
        """
        bcftools view {input.phased_vcf_annotate} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.phased_maf_info} 2> {log}

        bcftools index -f {output.phased_maf_info}
        """


rule maf_INFO_sites_concordance_HC_phased:
    """
    Extract MAF sites from HC phased VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz",
        phased_vcf_annotate="output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}.vcf.gz",
    output:
        phased_maf_info="output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_INFO_{info}_MAF_{maf}.vcf.gz",
    params:
        info=config["info_cutoff"],
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_INFO_{info}_MAF_{maf}.vcf.gz.log",
    threads: 4
    shell:
        """
        bcftools view {input.phased_vcf_annotate} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.phased_maf_info} 2> {log}

        bcftools index -f {output.phased_maf_info}
        """
