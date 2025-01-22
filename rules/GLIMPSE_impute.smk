#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

##############################
# Impute and phase bam files #
##############################

global bams_df

### Make sure bamfile list is correct !!!!


rule extract_chrom_ref_fast:
    """
    Extract chromosome from fasta reference file
    """
    input:
        ref_fasta=config["ref_fasta_file"],
    output:
        ref_fasta_chr="output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}.fasta",
        ref_fasta_chr_fai="output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}.fasta.fai",
    shell:
        """
        samtools faidx {input.ref_fasta} {wildcards.chrom} > {output.ref_fasta_chr}
        samtools faidx {output.ref_fasta_chr}
        """


rule extract_var_pos:
    """
    Extract sites from reference panel to use for GLs estimation from bams afterwards
    """
    input:
        ref_panel_phased="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
    output:
        ref_panel_sites_vcf="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.tsv.gz",
    log:
        "output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.log",
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.tsv"
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps \
        {input.ref_panel_phased} \
        --threads {threads} \
        -Oz -o {output.ref_panel_sites_vcf} 2> {log}

        bcftools index -f {output.ref_panel_sites_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_panel_sites_vcf} | \
        bgzip -c > {output.ref_panel_sites_tsv}

        tabix -s1 -b2 -e2 {output.ref_panel_sites_tsv}
        """


rule compute_GLs_imputed_samples:
    """
    Compute GLs of bams (do not need to extract exact chrom from bams)
    """
    input:
        imputation_bams=lambda wildcards: bams_df.loc[wildcards.bam_imputation, "Bam"],
        ref_panel_sites_vcf="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.tsv.gz",
        ref_fasta_chr="output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}.fasta",
    output:
        GL_imputed_bams="output/GLIMPSE_imputation/GLs_imputed_bams/{bam_imputation}_{chrom}.vcf.gz",
    log:
        "output/GLIMPSE_imputation/GLs_imputed_bams/{bam_imputation}_{chrom}.log",
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_imputation/GLs_imputed_bams/{bam_imputation}_{chrom}.tsv"
    shell:
        """
        bcftools mpileup -f {input.ref_fasta_chr} -I -E -a 'FORMAT/DP' -T {input.ref_panel_sites_vcf} -r {wildcards.chrom} {input.imputation_bams} -Ou | \
        bcftools call -Aim -C alleles -T {input.ref_panel_sites_tsv} -Oz -o {output.GL_imputed_bams} --threads {threads} 2> {log}
        
        bcftools index -f {output.GL_imputed_bams}
        """


rule chunk_spliting:
    """
    Split chromosome into chunks for imputation
    """
    input:
        ref_panel_sites_vcf="output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_sites.vcf.gz",
    output:
        chunks="output/GLIMPSE_imputation/chunks/{chrom}_chunks.txt",
    log:
        "output/GLIMPSE_imputation/chunks/{chrom}_chunks.log",
    params:
        window_size=config["window_size"],
        buffer_size=config["buffer_size"],
    shell:
        """
        GLIMPSE_chunk \
        --input {input.ref_panel_sites_vcf} \
        --region {wildcards.chrom} \
        --window-size {params.window_size} --buffer-size {params.buffer_size} \
        --output {output.chunks} 2> {log}
        """


rule imput_phase:
    """
    Impute all samples individually
    """
    input:
        GL_imputed_bams="output/GLIMPSE_imputation/GLs_imputed_bams/{bam_imputation}_{chrom}.vcf.gz",
        #ref_sample_snp_filltags_filter = "output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz",
        ref_panel_phased="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        chunks="output/GLIMPSE_imputation/chunks/{chrom}_chunks.txt",
        gen_map="data/gen_map/{chrom}_average_canFam3.1_modified.tsv",
    output:
        imputed=expand(
            "output/GLIMPSE_imputation/GLIMPSE_imputed/{bam_imputation}_imputed.{chrom}.00.bcf",
            allow_missing=True,
        ),
    params:
        prefix="output/GLIMPSE_imputation/GLIMPSE_imputed/{bam_imputation}_imputed.{chrom}",
    threads: 2
    log:
        "output/GLIMPSE_imputation/GLIMPSE_imputed/{bam_imputation}_imputed.{chrom}.log",
    benchmark:
        "benchmarks/GLIMPSE_imputation/GLIMPSE_imputed/{bam_imputation}_imputed.{chrom}.tsv"
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT={params.prefix}.${{ID}}.bcf
            GLIMPSE_phase \
            --input {input.GL_imputed_bams} --reference {input.ref_panel_phased} --map {input.gen_map} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} --output ${{OUT}} \
            --thread {threads}
            bcftools index -f ${{OUT}}
        done < {input.chunks} 2> {log}
        """


rule ligate_list:
    """
    Create list of imputed output files for each chunk to merge later
    """
    input:
        chunks="output/GLIMPSE_imputation/chunks/{chrom}_chunks.txt",
        imputed="output/GLIMPSE_imputation/GLIMPSE_imputed/{bam_imputation}_imputed.{chrom}.00.bcf",
    output:
        ligated_list="output/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated_list.{chrom}.txt",
    params:
        prefix="output/GLIMPSE_imputation/GLIMPSE_imputed/{bam_imputation}_imputed.{chrom}",
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        """


rule ligate:
    """
    Merge all imputed chunks
    """
    input:
        ligated_list="output/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated_list.{chrom}.txt",
    output:
        ligated_bcf="output/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated.{chrom}.bcf",
    log:
        "output/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated.{chrom}.log",
    threads: 8
    benchmark:
        "benchmarks/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated.{chrom}.tsv"
    shell:
        """
        GLIMPSE_ligate \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.ligated_bcf}
        """


rule sample_haplotype:
    """
    Phase!!!
    """
    input:
        ligated_bcf="output/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated.{chrom}.bcf",
    output:
        phased_bcf="output/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased.{chrom}.bcf",
    log:
        "output/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased.{chrom}.log",
    threads: 8
    benchmark:
        "benchmarks/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased.{chrom}.tsv"
    shell:
        """
        GLIMPSE_sample \
        --input {input.ligated_bcf} \
        --solve \
        --output {output.phased_bcf} \
        --thread {threads} 2> {log}
        
        bcftools index -f {output.phased_bcf}
        """


rule annotate_fields:
    """
    SOS: This step is essential to retain the fields from the ligated samples in the phased files (not automatically transferred)
    GP is needed to be annotatedm since it's used downstream for recalibrating the INFO score
    """
    input:
        phased_bcf="output/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased.{chrom}.bcf",
        ligated_bcf="output/GLIMPSE_imputation/GLIMPSE_ligated/{bam_imputation}_ligated.{chrom}.bcf",
    output:
        phased_vcf_annotate="output/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased_annotated.{chrom}.vcf.gz",
    log:
        "output/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased_annotated.{chrom}.bcf.log",
    benchmark:
        "benchmarks/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased_annotated.{chrom}.bcf.tsv"
    shell:
        """
        bcftools annotate \
        --annotations {input.ligated_bcf} \
        --columns FORMAT/DS,FORMAT/GP,FORMAT/HS \
        --output-type z \
        --output {output.phased_vcf_annotate} {input.phased_bcf} 2> {log}
        
        bcftools index --tbi {output.phased_vcf_annotate}
        """
