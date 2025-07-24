#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

##############################################
# Running imputation for GLIMPSE concordance #
##############################################

global CHROM, COVERAGE_VAL, INFO_CUTOFF, SAMPLE_CON


rule prepare_ref_panel_only_dogs:
    """
    Removes overlapping target and reference panel samples from the reference panel
    """
    input:
        ref_sample_snp_filltags_filter="output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        ref_val_sample_file="output/GLIMPSE_concordance/reference_panel/ref_val_sample.txt",
    output:
        ref_concordance_sample_excl=temp(
            "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel.phased.vcf.gz"
        ),
    log:
        "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/{chrom}_ref_panel.phased.tsv"
    shell:
        """
        bcftools view \
        -r {wildcards.chrom} \
        -S ^{input.ref_val_sample_file} \
        --trim-alt-alleles \
        {input.ref_sample_snp_filltags_filter} -Ou | \
        bcftools norm -a -Ou | \
        bcftools view -m 2 -M 2 -v snps  \
        -Oz -o {output.ref_concordance_sample_excl}
        """


rule fill_tags_ref_sample_excl_only_dogs:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING)
    """
    input:
        ref_concordance_sample_excl="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel.phased.vcf.gz",
    output:
        ref_concordance_sample_excl_filltags=temp(
            "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags.phased.vcf.gz"
        ),
    log:
        "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/{chrom}_ref_panel_filltags.tsv"
    shell:
        """
        bcftools +fill-tags {input.ref_concordance_sample_excl} \
        -Oz -o {output.ref_concordance_sample_excl_filltags} \
        --threads {threads} \
        -- -t all,F_MISSING 2> {log}
        """


rule filter_sites_ref_sample_excl_only_dogs:
    """
    Filter for missingness (F_MISSING) again, since we removed individuals
    """
    input:
        ref_concordance_sample_excl_filltags="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags.phased.vcf.gz",
    output:
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.bcf",
        ref_concordance_sample_excl_filltags_filter_csi="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.bcf.csi",
    params:
        f_missing=config["F_MISSING"],
    log:
        "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.tsv"
    shell:
        """
        bcftools view -i 'F_MISSING<{params.f_missing}' \
        {input.ref_concordance_sample_excl_filltags} \
        --threads {threads} \
        -Ob -o {output.ref_concordance_sample_excl_filltags_filter} 2> {log}

        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter}
        """


rule prepare_merged_chr_list_ref_only_dogs:
    """ 
    Prepare list to merge chromosomes for reference panel
    """
    input:
        ref_concordance_sample_excl_filltags_filter=expand(
            "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list_ref="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/chr_list.ref_panel_filltags_filter.phased.txt",
    shell:
        """
        ls -v {input.ref_concordance_sample_excl_filltags_filter} >> {output.chr_list_ref}
        """


rule merge_chr_ref_only_dogs:
    """
    Merge chromosomes of filtered reference panel
    """
    input:
        chr_list_ref="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/chr_list.ref_panel_filltags_filter.phased.txt",
    output:
        ref_concordance_sample_excl_filltags_filter_allchrom="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/allchrom_ref_panel_filltags_filter.phased.bcf",
        ref_concordance_sample_excl_filltags_filter_allchrom_csi="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/allchrom_ref_panel_filltags_filter.phased.bcf.csi",
    log:
        "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/allchrom_ref_panel_filltags_filter.phased.bcf.log",
    threads: 4
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list_ref} \
        -Ob -o {output.ref_concordance_sample_excl_filltags_filter_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_allchrom}
        """


rule extract_var_pos_concordance_only_dogs:
    """
    Extract sites from reference panel to use for GLs estimation from bams afterwards
    """
    input:
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.bcf",
    output:
        ref_panel_sites_vcf="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.tsv.gz",
    log:
        "output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.tsv"
    shell:
        """
        bcftools view -G -m 2 -M 2 -v snps \
        {input.ref_concordance_sample_excl_filltags_filter} \
        --threads {threads} \
        -Oz -o {output.ref_panel_sites_vcf} 2> {log}

        bcftools index -f {output.ref_panel_sites_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_panel_sites_vcf} | \
        bgzip -c > {output.ref_panel_sites_tsv}
        
        tabix -s1 -b2 -e2 {output.ref_panel_sites_tsv}
        """


rule compute_GLs_downsampled_samples_concordance_only_dogs:
    """
    Compute GLs of target bams 
    """
    input:
        downsampled_bam="output/GLIMPSE_concordance/target_bams/{sample_con}_{chrom}_{coverage_val}x.bam",
        ref_panel_sites_vcf="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.tsv.gz",
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
    output:
        GL_vcf_target_bams="output/GLIMPSE_concordance_only_dogs/GLs_target_bams/{sample_con}_{chrom}_{coverage_val}x.vcf.gz",
        GL_vcf_target_bams_csi="output/GLIMPSE_concordance_only_dogs/GLs_target_bams/{sample_con}_{chrom}_{coverage_val}x.vcf.gz.csi",
    log:
        "output/GLIMPSE_concordance_only_dogs/GLs_target_bams/{sample_con}_{chrom}_{coverage_val}x.log",
    threads: 4
    benchmark:
        "benchmarks/GLs_target_bams/{sample_con}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        bcftools mpileup -f {input.ref_fasta_chr} -I -E -a 'FORMAT/DP' -T {input.ref_panel_sites_vcf} -r {wildcards.chrom} {input.downsampled_bam} -Ou | \
        bcftools call -Aim -C alleles -T {input.ref_panel_sites_tsv} -Oz -o {output.GL_vcf_target_bams} --threads {threads} 2> {log}
        
        bcftools index -f {output.GL_vcf_target_bams}
        """


rule chunk_spliting_concordance_only_dogs:
    """
    Split chromosome into chunks for imputation
    """
    input:
        ref_panel_sites_vcf="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_sites.phased.vcf.gz",
    output:
        chunks="output/GLIMPSE_concordance_only_dogs/chunks/{chrom}_chunks.txt",
    params:
        window_size=config["window_size"],
        buffer_size=config["buffer_size"],
    log:
        "output/GLIMPSE_concordance_only_dogs/chunks/{chrom}_chunks.log",
    benchmark:
        "benchmarks/chunks/{chrom}_chunks.tsv"
    shell:
        """
        GLIMPSE_chunk \
        --input {input.ref_panel_sites_vcf} \
        --region {wildcards.chrom} \
        --window-size {params.window_size} --buffer-size {params.buffer_size} \
        --thread 5 \
        --output {output.chunks} 2> {log}
        """


rule impute_concordance_only_dogs:
    """
    Impute each sample seperately
    """
    input:
        GL_vcf_target_bams="output/GLIMPSE_concordance_only_dogs/GLs_target_bams/{sample_con}_{chrom}_{coverage_val}x.vcf.gz",
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom}_ref_panel_filltags_filter.phased.bcf",
        chunks="output/GLIMPSE_concordance_only_dogs/chunks/{chrom}_chunks.txt",
        gen_map="data/gen_map/{chrom}_average_canFam3.1_modified.txt",
    output:
        imputed="output/GLIMPSE_concordance_only_dogs/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x.00.bcf",
        imputed_csi="output/GLIMPSE_concordance_only_dogs/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x.00.bcf.csi",
    params:
        prefix="output/GLIMPSE_concordance_only_dogs/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x",
    threads: 2
    log:
        "output/GLIMPSE_concordance_only_dogs/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x.log",
    benchmark:
        "benchmarks/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT={params.prefix}.${{ID}}.bcf
            GLIMPSE_phase \
            --input {input.GL_vcf_target_bams} --reference {input.ref_concordance_sample_excl_filltags_filter} --map {input.gen_map} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} --output ${{OUT}} \
            --thread {threads}
            bcftools index -f ${{OUT}}
        done < {input.chunks} 2> {log}
        """


rule ligate_list_concordance_only_dogs:
    """
    Create list of imputed output files for each chunk to merge later
    """
    input:
        chunks="output/GLIMPSE_concordance_only_dogs/chunks/{chrom}_chunks.txt",
        imputed="output/GLIMPSE_concordance_only_dogs/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x.00.bcf",
    output:
        ligated_list="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/ligated_list_{sample_con}_{chrom}_{coverage_val}x.txt",
    params:
        prefix="output/GLIMPSE_concordance_only_dogs/GLIMPSE_imputed/{sample_con}_{chrom}_{coverage_val}x",
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        """


rule ligate_concordance_only_dogs:
    """
    Merge all imputed chunks
    """
    input:
        ligated_list="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/ligated_list_{sample_con}_{chrom}_{coverage_val}x.txt",
    output:
        ligated_bcf="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample_con}_{chrom}_{coverage_val}x.bcf",
        ligated_bcf_csi="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample_con}_{chrom}_{coverage_val}x.bcf.csi",
    log:
        "output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample_con}_{chrom}_{coverage_val}x.log",
    threads: 10
    benchmark:
        "benchmarks/GLIMPSE_ligated/merged_ligated.{sample_con}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        GLIMPSE_ligate \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.ligated_bcf}
        """


rule phase_concordance_only_dogs:
    """
    Phase!!!
    """
    input:
        ligated_bcf="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample_con}_{chrom}_{coverage_val}x.bcf",
    output:
        phased_bcf="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased.{sample_con}_{chrom}_{coverage_val}x.bcf",
        phased_bcf_csi="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased.{sample_con}_{chrom}_{coverage_val}x.bcf.csi",
    log:
        "output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased.{sample_con}_{chrom}_{coverage_val}x.log",
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_phased/phased.{sample_con}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        GLIMPSE_sample \
        --input {input.ligated_bcf} --solve --output {output.phased_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.phased_bcf}
        """

rule annotate_fields_concordance_only_dogs:
    """
    SOS: This step is essential to retain the fields from the ligated samples in the phased files (not automatically transferred)
    GP is needed to be annotatedm since it's used downstream for recalibrating the INFO score
    """
    input:
        phased_bcf = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased.{sample_con}_{chrom}_{coverage_val}x.bcf',
        ligated_bcf = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample_con}_{chrom}_{coverage_val}x.bcf'
    output:
        phased_vcf_annotate = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased_annotated.{sample_con}_{chrom}_{coverage_val}x.vcf.gz',
    log:
        'output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased_annotated.{sample_con}_{chrom}_{coverage_val}x.vcf.gz.log'
    benchmark:
        'benchmarks/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased_annotated.{sample_con}_{chrom}_{coverage_val}x.tsv'
    shell:
        '''
        bcftools annotate \
        --annotations {input.ligated_bcf} \
        --columns FORMAT/DS,FORMAT/GP,FORMAT/HS \
        --output-type z \
        --output {output.phased_vcf_annotate} {input.phased_bcf} 2> {log}
        
        bcftools index --tbi {output.phased_vcf_annotate}
        '''


rule filter_info_score_only_dogs:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        phased_vcf_annotate = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased/phased_annotated.{sample_con}_{chrom}_{coverage_val}x.vcf.gz',
    output:
        phased_vcf_info = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf',
        phased_vcf_csi = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf.csi'
    params:
        info_val="{info_cutoff}",
    log:
        "output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.log",
    threads: 4
    shell:
        """
        bcftools view {input.phased_vcf_annotate} \
        --include 'INFO/INFO >= {params.info_val}' \
        --threads {threads} \
        -Ob -o {output.phased_vcf_info} 2> {log}

        bcftools index -f {output.phased_vcf_info}
        """


#######################################
#                                     #
#  Run GLIMPSE concordance all chrom  #
#                                     #
#######################################


rule prepare_merged_chr_list_concordance_only_dogs:
    """ 
    Prepare list to merge chromosomes for reference, imputed and validation
    """
    input:
        phased_vcf_info=expand(
            "output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/chr_list.{sample_con}_{coverage_val}x-INFO_{info_cutoff}.txt",
    shell:
        """
        ls -v {input.phased_vcf_info} >> {output.chr_list}
        """


rule merge_chr_concordance_only_dogs:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        chr_list="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/chr_list.{sample_con}_{coverage_val}x-INFO_{info_cutoff}.txt",
    output:
        phased_vcf_info_allchrom="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf",
        phased_vcf_info_allchrom_csi="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf.csi",
    log:
        "output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.log",
    threads: 4
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Ob -o {output.phased_vcf_info_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.phased_vcf_info_allchrom}
        """


rule get_ID_for_targets_allchrom_only_dogs:
    """
    Get sample ID for targets
    """
    input:
        phased_vcf_info_allchrom="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf",
    output:
        sm_samples="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/sm_phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.txt",
    shell:
        """
        bcftools query -l {input.phased_vcf_info_allchrom} > {output.sm_samples}
        """


rule prepare_concordance_lst_info_score_filtered_allchrom_only_dogs:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter_allchrom="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/allchrom_ref_panel_filltags_filter.phased.bcf",
        validation_sample_filt_allelic_allchrom="output/GLIMPSE_concordance/validation_bams/{sample_con}_allchrom_validation_filt_qual_dp_ab.bcf",
        phased_vcf_info_allchrom="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf",
    output:
        concordance_lst_info_score_filtered="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.lst",
    shell:
        """
        echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chr31,chr32,chr33,chr34,chr35,chr36,chr37,chr38" {input.ref_concordance_sample_excl_filltags_filter_allchrom} {input.validation_sample_filt_allelic_allchrom} {input.phased_vcf_info_allchrom} > {output.concordance_lst_info_score_filtered}
        """


rule GLIMPSE_concordance_info_score_filtered_allchrom_only_dogs:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.lst",
        sm_samples="output/GLIMPSE_concordance_only_dogs/GLIMPSE_phased_INFO_filtered/sm_phased_annotated.{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}.txt",
    output:
        concordance_output_info_score_filtered="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz",
        concordance_output_discordance_filtered="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz",
    params:
        prefix="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered",
    log:
        "output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.log",
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


rule prepare_concordance_output_filt_allchrom_only_dogs:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz",
    output:
        concordance_output_discordance_filtered_temp=temp(
            "output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/temp_concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.txt"
        ),
        concordance_output_discordance_filtered_prep="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.txt",
    shell:
        """
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        """


rule merge_concordance_output_filt_allchrom_only_dogs:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep=expand(
            "output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.txt",
            coverage_val=COVERAGE_VAL,
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        concordance_output_discordance_filtered_per_sample="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_filtered.txt",
    shell:
        """
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        """


rule prepare_only_dogs_files:
    """
    Prepare files for comparison plot
    """
    input:
        concordance_output_info_score_filtered_only_dogs="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz",
    output:
        concordance_output_info_score_filtered_only_dogs_mod="output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
    shell:
        """
        zcat {input.concordance_output_info_score_filtered_only_dogs} | awk -v FS=' ' -v OFS=' ' '{{$7={wildcards.coverage_val}}} {{$8={wildcards.info_cutoff}}} {{$9="{wildcards.sample_con}"}} 1' > {output.concordance_output_info_score_filtered_only_dogs_mod}
        """


rule plot_concordance_filt_allsites_only_dogs_0_5x_1x:
    """
    Plot concordance for 0.5x and 1x comparing the two panels 
    """
    input:
        concordance_output_info_score_filtered_mod=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
            sample_con=SAMPLE_CON,
            info_cutoff=INFO_CUTOFF,
            coverage_val=["0.5", "1"],
            allow_missing=True,
        ),
        concordance_output_info_score_filtered_only_dogs_mod=expand(
            "output/GLIMPSE_concordance_only_dogs/concordance_INFO_filtered/concordance_{sample_con}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
            sample_con=SAMPLE_CON,
            info_cutoff=INFO_CUTOFF,
            coverage_val=["0.5", "1"],
            allow_missing=True,
        ),
        concordance_metadata="sample_lists/concordance_bams_published.tsv",
    output:
        concordance_all_sites_only_dogs="output/GLIMPSE_concordance_only_dogs/plots/glimpse_concordance/concordance_allchrom_0.5x_1x_filtered-all_sites_only_dogs.png",
    params:
        files_all=lambda wildcards, input: ",".join(
            input.concordance_output_info_score_filtered_mod
        ),
        files_dogs=lambda wildcards, input: ",".join(
            input.concordance_output_info_score_filtered_only_dogs_mod
        ),
    shell:
        """
        Rscript scripts/rsquare_accuracy_all_sites_only_dogs.R \
        {params.files_all} \
        {params.files_dogs} \
        {input.concordance_metadata} \
        {output.concordance_all_sites_only_dogs}
        """
