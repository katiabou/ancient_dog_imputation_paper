#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

##############################################
# Running GLIMPSE concordance for HC imputed #
##############################################

global CHROM


rule compute_GLs_HC_samples_concordance:
    """
    Compute GLs of HC target bams 
    """
    input:
        target_bams_chr="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}.bam",
        ref_panel_sites_vcf="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.tsv.gz",
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
    output:
        GL_vcf_HC_target_bams="output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}.vcf.gz",
        GL_vcf_HC_target_bams_csi="output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}.vcf.gz.csi",
    log:
        "output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}.log",
    threads: 8
    benchmark:
        "benchmarks/GLs_target_bams/{sample}_{chrom}.tsv"
    shell:
        """
        bcftools mpileup -f {input.ref_fasta_chr} -I -E -a 'FORMAT/DP' -T {input.ref_panel_sites_vcf} -r {wildcards.chrom} {input.target_bams_chr} -Ou | \
        bcftools call -Aim -C alleles -T {input.ref_panel_sites_tsv} -Oz -o {output.GL_vcf_HC_target_bams} --threads {threads} 2> {log}
        
        bcftools index -f {output.GL_vcf_HC_target_bams}
        """


rule impute_HC_concordance:
    """
    Impute each sample seperately
    """
    input:
        GL_vcf_HC_target_bams="output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}.vcf.gz",
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
        chunks="output/GLIMPSE_concordance/chunks/{chrom}_chunks.txt",
        gen_map="data/gen_map/{chrom}_average_canFam3.1_modified.tsv",
    output:
        imputed="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}.00.bcf",
        imputed_csi="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}.00.bcf.csi",
    params:
        prefix="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}",
    threads: 2
    log:
        "output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}.log",
    benchmark:
        "benchmarks/GLIMPSE_imputed/{sample}_{chrom}.tsv"
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT={params.prefix}.${{ID}}.bcf
            GLIMPSE_phase \
            --input {input.GL_vcf_HC_target_bams} --reference {input.ref_concordance_sample_excl_filltags_filter} --map {input.gen_map} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} --output ${{OUT}} \
            --thread {threads}
            bcftools index -f ${{OUT}}
        done < {input.chunks} 2> {log}
        """


rule ligate_HC_list_concordance:
    """
    Create list of imputed output files for each chunk to merge later
    """
    input:
        chunks="output/GLIMPSE_concordance/chunks/{chrom}_chunks.txt",
        imputed="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}.00.bcf",
    output:
        ligated_list="output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom}.txt",
    params:
        prefix="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}",
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        """


rule ligate_HC_concordance:
    """
    Merge all imputed chunks
    """
    input:
        ligated_list="output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom}.txt",
    output:
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.bcf",
        ligated_bcf_csi="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.log",
    threads: 8
    benchmark:
        "benchmarks/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.tsv"
    shell:
        """
        GLIMPSE_ligate \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.ligated_bcf}
        """


rule phase_HC_concordance:
    """
    Phase!!!
    """
    input:
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.bcf",
    output:
        phased_bcf="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.bcf",
        phased_bcf_csi="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.log",
    threads: 8
    benchmark:
        "benchmarks/GLIMPSE_phased/phased.{sample}_{chrom}.tsv"
    shell:
        """
        GLIMPSE_sample \
        --input {input.ligated_bcf} --solve --output {output.phased_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.phased_bcf}
        """


rule filter_info_score_HC:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.bcf",
    output:
        info_imputed_info="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}-INFO_{info_cutoff}.bcf",
        info_imputed_info_csi="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}-INFO_{info_cutoff}.bcf.csi",
    params:
        info_val="{info_cutoff}",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}-INFO_{info_cutoff}.log",
    threads: 8
    shell:
        """
        bcftools view {input.ligated_bcf} \
        --include 'INFO/INFO >= {params.info_val}' \
        --threads {threads} \
        -Ob -o {output.info_imputed_info} 2> {log}

        bcftools index -f {output.info_imputed_info}
        """


#######################################
#                                     #
#  Run GLIMPSE concordance all chrom  #
#                                     #
#######################################


rule prepare_merged_chr_list_HC_concordance:
    """ 
    Prepare list to merge chromosomes for reference, imputed and validation
    """
    input:
        info_imputed_info=expand(
            "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}-INFO_{info_cutoff}.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/chr_list.{sample}-INFO_{info_cutoff}.txt",
    shell:
        """
        ls -v {input.info_imputed_info} >> {output.chr_list}
        """


rule merge_chr_HC_concordance:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/chr_list.{sample}-INFO_{info_cutoff}.txt",
    output:
        info_imputed_info_allchrom="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom-INFO_{info_cutoff}.bcf",
        info_imputed_info_allchrom_csi="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom-INFO_{info_cutoff}.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom-INFO_{info_cutoff}.log",
    threads: 8
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Ob -o {output.info_imputed_info_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.info_imputed_info_allchrom}
        """
