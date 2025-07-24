#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

###############################
# Running GLIMPSE concordance #
###############################

global CHROM, COVERAGE_VAL, INFO_CUTOFF, SAMPLE_CON


rule transversions_validation_concordance:
    """
    Only take transversions from filtered validation data
    """
    input:
        validation_sample_filt_allelic="output/GLIMPSE_concordance/validation_bams/{sample_con}_{chrom}_validation_filt_qual_dp_ab.bcf",
    output:
        tranversion_sites_allelic="output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_{chrom}_validation_filt_transversions.tsv.gz",
        validation_transversions_allelic="output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_{chrom}_validation_filt_transversions.bcf",
    log:
        "output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_{chrom}_validation_filt_transversions.log",
    shell:
        """
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.validation_sample_filt_allelic} | bgzip -c > {output.tranversion_sites_allelic}

        tabix -s1 -b2 -e2 {output.tranversion_sites_allelic}

        bcftools view {input.validation_sample_filt_allelic} \
        --regions-file {output.tranversion_sites_allelic} \
        -Ob -o {output.validation_transversions_allelic} 2> {log}

        bcftools index -f {output.validation_transversions_allelic}
        """


rule prepare_merged_chr_list_trans_validation:
    """ 
    Prepare list to merge chromosomes for validation
    """
    input:
        validation_transversions_allelic=expand(
            "output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_{chrom}_validation_filt_transversions.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list_validation="output/GLIMPSE_concordance/validation_bams_transversions/chr_list.{sample_con}_validation_filt_transversions.txt",
    shell:
        """
        ls -v {input.validation_transversions_allelic} >> {output.chr_list_validation}
        """


rule merge_chr_trans_validation:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        chr_list_validation="output/GLIMPSE_concordance/validation_bams_transversions/chr_list.{sample_con}_validation_filt_transversions.txt",
    output:
        validation_sample_filt_allelic_allchrom="output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_allchrom_validation_filt_transversions.bcf",
        validation_sample_filt_allelic_allchrom_csi="output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_allchrom_validation_filt_transversions.bcf.csi",
    log:
        "output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_allchrom_validation_filt_transversions.bcf.log",
    threads: 4
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list_validation} \
        -Ob -o {output.validation_sample_filt_allelic_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.validation_sample_filt_allelic_allchrom}
        """


rule filter_transversions_imputed:
    """
    Filter for transversions only
    """
    input:
        phased_vcf_info="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf",
    output:
        tranversion_sites="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}_transversions.tsv.gz",
        phased_transversions="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf.log",
    threads: 4
    shell:
        """
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.phased_vcf_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.phased_vcf_info} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Ob -o {output.phased_transversions} 2> {log}

        bcftools index -f {output.phased_transversions}
        """


#######################################
#                                     #
#  Run GLIMPSE concordance all chrom  #
#                                     #
#######################################


rule prepare_merged_chr_list_trans_concordance:
    """ 
    Prepare list to merge chromosomes for reference, imputed and validation
    """
    input:
        phased_transversions=expand(
            "output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/chr_list.{sample_con}_{coverage_val}x-INFO_{info_cutoff}_transversions.txt",
    shell:
        """
        ls -v {input.phased_transversions} >> {output.chr_list}
        """


rule merge_chr_trans_concordance:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/chr_list.{sample_con}_{coverage_val}x-INFO_{info_cutoff}_transversions.txt",
    output:
        phased_info_allchrom = 'output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf',
        phased_info_allchrom_csi = 'output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf.csi'
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered_transversions/merged_ligated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.log",
    threads: 4
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Ob -o {output.phased_info_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.phased_info_allchrom}
        """


rule get_ID_for_targets_allchrom_trans:
    """
    Get sample name for concordance
    """
    input:
        phased_info_allchrom = 'output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf',
    output:
        sm_samples="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/sm_phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.txt",
    shell:
        """
        bcftools query -l {input.phased_info_allchrom} > {output.sm_samples}
        """


rule prepare_concordance_lst_info_score_filtered_allchrom_trans:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter_allchrom="output/GLIMPSE_concordance/reference_panel/allchrom_ref_panel_filltags_filter.phased.bcf",
        validation_sample_filt_allelic_allchrom="output/GLIMPSE_concordance/validation_bams_transversions/{sample_con}_allchrom_validation_filt_transversions.bcf",
        phased_info_allchrom = 'output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf',
    output:
        concordance_lst_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.lst",
    shell:
        """
        echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chr31,chr32,chr33,chr34,chr35,chr36,chr37,chr38" {input.ref_concordance_sample_excl_filltags_filter_allchrom} {input.validation_sample_filt_allelic_allchrom} {input.phased_info_allchrom} > {output.concordance_lst_info_score_filtered}
        """


rule GLIMPSE_concordance_info_score_filtered_allchrom_trans:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.lst",
        sm_samples="output/GLIMPSE_concordance/GLIMPSE_phased_INFO_filtered_transversions/sm_phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_transversions.txt",
    output:
        concordance_output_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.rsquare.grp.txt.gz",
        concordance_output_discordance_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.error.spl.txt.gz",
    params:
        prefix="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions",
    log:
        "output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.log",
    threads: 4
    shell:
        """
        GLIMPSE_concordance \
        --input {input.concordance_lst_info_score_filtered} \
        --minDP 8 \
        --output {params.prefix} \
        --minPROB 0.9 \
        --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
        --sample {input.sm_samples} \
        --af-tag AF \
        --thread {threads} 2> {log}
        """


rule prepare_concordance_output_filt_allchrom_trans:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.error.spl.txt.gz",
    output:
        concordance_output_discordance_filtered_temp=temp(
            "output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/temp_concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.txt"
        ),
        concordance_output_discordance_filtered_prep="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.txt",
    shell:
        """
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        """


rule merge_concordance_output_filt_allchrom_trans:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.txt",
            coverage_val=COVERAGE_VAL,
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        concordance_output_discordance_filtered_per_sample="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_filtered_transversions.txt",
    shell:
        """
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        """


rule plot_discordance_filt_allchrom_trans:
    """
    Plot genotype discordances
    """
    input:
        concordance_output_discordance_filtered_per_sample=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_filtered_transversions.txt",
            sample_con=SAMPLE_CON,
            allow_missing=True,
        ),
        concordance_metadata="sample_lists/concordance_bams_published.tsv",
    params:
        path_script="scripts",
        discordance_phased=lambda wildcards, input: ",".join(
            input.concordance_output_discordance_filtered_per_sample
        ),
    output:
        discordance_dogs_full="output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_dogs_full_transversions.png",
        discordance_wolves_full="output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_wolves_full_transversions.png",
    shell:
        """
        Rscript {params.path_script}/genotype_discordance.R \
        {params.discordance_phased} \
        {input.concordance_metadata} \
        {output.discordance_dogs_full} \
        {output.discordance_wolves_full}
        """


rule plot_discordance_filt_allsites_transversions_comparison_allchrom:
    """
    Plot genotype discordances of all sites against transversions
    """
    input:
        concordance_output_discordance_filtered_per_sample=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample_con}_allchrom_filtered.txt",
            sample_con=SAMPLE_CON,
            allow_missing=True,
        ),
        concordance_output_discordance_filtered_per_sample_trans=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_filtered_transversions.txt",
            sample_con=SAMPLE_CON,
            allow_missing=True,
        ),
        concordance_metadata="sample_lists/concordance_bams_published.tsv",
    params:
        path_script="scripts",
        discordance_phased=lambda wildcards, input: ",".join(
            input.concordance_output_discordance_filtered_per_sample
        ),
        discordance_phased_trans=lambda wildcards, input: ",".join(
            input.concordance_output_discordance_filtered_per_sample_trans
        ),
    output:
        discordance_dogs_full="output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_dogs_full_all_sites_transversions_comparison.png",
        discordance_wolves_full="output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_filtered_wolves_full_all_sites_transversions_comparison.png",
    shell:
        """
        Rscript {params.path_script}/genotype_discordance_allsites_transversions.R \
        {params.discordance_phased} \
        {params.discordance_phased_trans} \
        {input.concordance_metadata} \
        {output.discordance_dogs_full} \
        {output.discordance_wolves_full}
        """


rule prepare_transversion_files:
    """
    Prepare files for comparison plot
    """
    input:
        concordance_output_info_score_filtered_trans="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.rsquare.grp.txt.gz",
    output:
        concordance_output_info_score_filtered_trans_mod="output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.rsquare-mod.grp.txt.gz",
    shell:
        """
        zcat {input.concordance_output_info_score_filtered_trans} | awk -v FS=' ' -v OFS=' ' '{{$7={wildcards.coverage_val}}} {{$8={wildcards.info_cutoff}}} {{$9="{wildcards.sample_con}"}} 1' > {output.concordance_output_info_score_filtered_trans_mod}
        """


rule plot_concordance_filt_allsites_transversions_0_5x_1x:
    """
    Plot accuracy of full ref panel vs only dog ref panel
    """
    input:
        concordance_output_info_score_filtered_mod=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
            sample_con=SAMPLE_CON,
            info_cutoff=INFO_CUTOFF,
            coverage_val=["0.5", "1"],
            allow_missing=True,
        ),
        concordance_output_info_score_filtered_trans_mod=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.rsquare-mod.grp.txt.gz",
            sample_con=SAMPLE_CON,
            info_cutoff=INFO_CUTOFF,
            coverage_val=["0.5", "1"],
            allow_missing=True,
        ),
        concordance_metadata="sample_lists/concordance_bams_published.tsv",
    output:
        concordance_all_sites_trans="output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_0.5x_1x_filtered-all_sites_transversions.png",
    params:
        files_all_sites=lambda wildcards, input: ",".join(
            input.concordance_output_info_score_filtered_mod
        ),
        files_trans=lambda wildcards, input: ",".join(
            input.concordance_output_info_score_filtered_trans_mod
        ),
    shell:
        """
        Rscript scripts/rsquare_accuracy_all_sites_transversions.R \
        {params.files_all_sites} \
        {params.files_trans} \
        {input.concordance_metadata} \
        {output.concordance_all_sites_trans}
        """
