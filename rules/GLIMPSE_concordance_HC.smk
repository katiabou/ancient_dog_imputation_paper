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
    threads: 4
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
        "benchmarks/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}.tsv"
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
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.tsv"
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
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.tsv"
    shell:
        """
        GLIMPSE_sample \
        --input {input.ligated_bcf} \
        --solve \
        --output {output.phased_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.phased_bcf}
        """


rule annotate_fields_concordance_HC:
    """
    SOS: This step is essential to retain the fields from the ligated samples in the phased files (not automatically transferred)
    GP is needed to be annotatedm since it's used downstream for recalibrating the INFO score
    """
    input:
        phased_bcf="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.bcf",
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}.bcf",
    output:
        phased_vcf_annotate="output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}.vcf.gz",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}.vcf.gz.log",
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}.tsv"
    shell:
        """
        bcftools annotate \
        --annotations {input.ligated_bcf} \
        --columns FORMAT/DS,FORMAT/GP,FORMAT/HS \
        --output-type z \
        --output {output.phased_vcf_annotate} {input.phased_bcf} 2> {log}
        
        bcftools index -f {output.phased_vcf_annotate}
        """