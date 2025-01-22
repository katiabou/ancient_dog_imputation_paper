#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

######################################################
# Prepare reference panel for concordance imputation #
######################################################

global CHROM


rule vcf_coverage:
    """
    Estimate depth of coverage of reference panel samples to filter out low coverage ones
    """
    input:
        ref=config["reference_panel"],
    output:
        ref_depth="output/reference_panel/reference_panel_depth.idepth",
    log:
        "output/reference_panel/reference_panel_depth.log",
    resources:
        mem_mb=50 * 1024,
    benchmark:
        "benchmarks/reference_panel/reference_panel_depth.tsv"
    shell:
        """
        vcftools \
        --gzvcf {input.ref} \
        --depth \
        --out {wildcards.path}/output/reference_panel/reference_panel_depth 2> {log}
        """


rule make_vcf_sample_list:
    """
    Make list of samples to remove from reference panel
    """
    input:
        ref_depth="output/reference_panel/reference_panel_depth.idepth",
        reseq_names="sample_lists/reseq_samples.tsv",
        boxer_names="sample_lists/boxer_samples.tsv",
    output:
        sample_list_exclude="output/reference_panel/remove_samples.txt",
    params:
        depth_cutoff=config["depth_cutoff"],
    shell:
        """
        awk -v cutoff="{params.depth_cutoff}" 'NR > 1 && $3 < cutoff {{print $1}}' {input.ref_depth} > temp.txt
        cat temp.txt {input.reseq_names} {input.boxer_names} | awk '!seen[$0]++' > {output.sample_list_exclude}
        rm temp.txt
        """


rule remove_sample_indels_multiallelic_snps:
    """
    Remove samples, trim-alt alleles created, normalise indels, only keep biallelic snps
    """
    input:
        ref=config["reference_panel"],
        sample_list_exclude="output/reference_panel/remove_samples.txt",
    output:
        ref_sample_snp=temp(
            "output/reference_panel/ref-panel_{chrom}_sample-snp.vcf.gz"
        ),
    log:
        "output/reference_panel/ref-panel_{chrom}_sample-snp.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/ref-panel_{chrom}_sample-snp.tsv"
    shell:
        """
        (
        bcftools view \
        -r {wildcards.chrom} \
        -S ^{input.sample_list_exclude} \
        --trim-alt-alleles \
        {input.ref} -Ou | \
        bcftools norm -a -Ou | \
        bcftools view -m 2 -M 2 -v snps  \
        -Oz -o {output.ref_sample_snp}
        ) 2> {log}
        """


rule fill_tags:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING, which is the fraction of missing genotypes)
    """
    input:
        ref_sample_snp="output/reference_panel/ref-panel_{chrom}_sample-snp.vcf.gz",
    output:
        ref_sample_snp_filltags=temp(
            "output/reference_panel/ref-panel_{chrom}_sample-snp_filltags.vcf.gz"
        ),
    log:
        "output/reference_panel/ref-panel_{chrom}_sample-snp_filltags.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/ref-panel_{chrom}_sample-snp_filltags.tsv"
    shell:
        """
        (
        bcftools +fill-tags {input.ref_sample_snp} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags} \
        -- -t all,F_MISSING
        ) 2> {log} 
        """


# rule get_F_missing_hist:
#     input:
#         ref_sample_snp_filltags="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags.vcf.gz",
#     output:
#         ref_panel_f_missing_hist=temp(
#             "output/reference_panel/ref-panel_{chrom}_f_missing.txt"
#         ),
#     # conda:
#     # "../envs/environment.yaml"
#     shell:
#         """
#        bcftools query -f '%INFO/F_MISSING\n' {input.ref_sample_snp_filltags} > {output.ref_panel_f_missing_hist}
#        """


# rule merge_F_missing:
#     input:
#         ref_panel_f_missing_hist=expand(
#             "output/reference_panel/ref-panel_{chrom}_f_missing.txt", chrom=CHROM
#         ),
#     output:
#         ref_panel_f_missing_hist_merged="output/reference_panel/ref-panel_allchrom_f_missing.txt",
#     shell:
#         """
#        cat {input.ref_panel_f_missing_hist} > {output.ref_panel_f_missing_hist_merged}
#        """


rule filter_sites:
    """
    Filter for only PASS sites and missingness (F_MISSING)
    """
    input:
        ref_sample_snp_filltags="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags.vcf.gz",
    output:
        ref_sample_snp_filltags_filter="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz",
        ref_sample_snp_filltags_filter_tbi="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz.tbi",
    params:
        f_missing=config["F_MISSING"],
    log:
        "output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.tsv"
    shell:
        """
        bcftools view -i 'FILTER=="PASS" & F_MISSING<{params.f_missing}' {input.ref_sample_snp_filltags} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags_filter} 2> {log}

        bcftools index --tbi {output.ref_sample_snp_filltags_filter}
        """


# rule phase_modern_data:
#     """
#     Phase filtered reference panel with shapeit
#     """
#     input:
#         ref_sample_snp_filltags_filter="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz",
#         gen_map="data/gen_map/{chrom}_average_canFam3.1_modified.tsv",
#     output:
#         ref_panel_phased="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
#         ref_panel_phased_tbi="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.tbi",
#     log:
#         log_shapeit="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.log",
#     resources:
#         mem_mb=720 * 1024,
#     threads: 96
#     benchmark:
#         "benchmarks/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.tsv"
#     shell:
#         """
#         /projects/racimolab/people/qcj125/programmes/shapeit4/bin/shapeit4.2 \
#         --input {input.ref_sample_snp_filltags_filter} \
#         --map {input.gen_map} \
#         --region {wildcards.chrom} \
#         --output {output.ref_panel_phased} \
#         --thread {threads} \
#         --sequencing &> {log.log_shapeit}
#
#         bcftools index --tbi {output.ref_panel_phased}
#         """


rule phase_modern_data_shapeit5:
    """
    Phase filtered reference panel with shapeit5
    """
    input:
        ref_sample_snp_filltags_filter="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz",
        gen_map="data/gen_map/{chrom}_average_canFam3.1_modified.tsv",
    output:
        ref_panel_phased="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        ref_panel_phased_tbi="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.tbi",
    log:
        log_shapeit="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.log",
    resources:
        mem_mb=720 * 1024,
    threads: 96
    benchmark:
        "benchmarks/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.tsv"
    shell:
        """
        SHAPEIT5_phase_common \
        --input {input.ref_sample_snp_filltags_filter} \
        --map {input.gen_map} \
        --region {wildcards.chrom} \
        --output {output.ref_panel_phased} \
        --thread {threads} &> {log.log_shapeit}

        bcftools index --tbi {output.ref_panel_phased}
        """
