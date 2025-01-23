#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

##############################################
# Running imputation for GLIMPSE concordance #
##############################################

global samples_df, CHROM, INFO_CUTOFF, COVERAGE_VAL, SAMPLE


rule extract_chrom_ref_fast_concordance:
    """
    Extract chromosome from fasta reference file
    """
    input:
        ref_fasta=config["ref_fasta_file"],
    output:
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
        ref_fasta_chr_fai="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta.fai",
    shell:
        """
        samtools faidx {input.ref_fasta} {wildcards.chrom} > {output.ref_fasta_chr}
        samtools faidx {output.ref_fasta_chr}
        """


rule extract_chr_target_bams:
    """
    Extract chromosomes from target bam files 
    """
    input:
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
        target_bams=lambda wildcards: samples_df.loc[wildcards.sample, "bam_path"],
    output:
        target_bams_chr="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}.bam",
        target_bams_chr_bai="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}.bam.bai",
    shell:
        """
        samtools view -T {input.ref_fasta_chr} \
        -bo {output.target_bams_chr} \
        {input.target_bams} \
        {wildcards.chrom}

        samtools index {output.target_bams_chr}
        """


rule estimate_coverage_fraction:
    """
    Estimate coverage fraction to be used for downsampling the target bams to lower coverages (specified in Snakefile)
    Based on here: https://davemcg.github.io/post/easy-bam-downsampling/
    """
    input:
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
        target_bams_chr="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}.bam",
    output:
        seed_frac=temp(
            "output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}.txt"
        ),
    params:
        coverage_val="{coverage_val}",
    shell:
        """
        LENGTH=$(cut -f1,2 {input.ref_fasta_chr}.fai | awk -F ' ' '{{print $2}}')
        COV=$(samtools depth -a {input.target_bams_chr}  |  awk -v var="$LENGTH" '{{sum+=$3}} END {{ print sum/var}}')
        NUM_READS=$(samtools view -c {input.target_bams_chr})
        coverage_val={params.coverage_val}
        reads_down=$(echo "$coverage_val $NUM_READS $COV" | awk '{{print ($1 * $2)/$3}}')
        frac=$( samtools idxstats {input.target_bams_chr} | cut -f3 | awk -v var="$reads_down" 'BEGIN {{total=0}} {{total += $1}} END {{frac=var/total; if (frac > 1) {{print 1}} else {{print frac}}}}' )
        seed_frac=$(echo "$frac 1" | awk '{{print $1 + $2}}')
        echo $seed_frac > {output.seed_frac}
        """


rule downsample_target_bam:
    """
    Downsample target bams to lower coverages
    """
    input:
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
        target_bams_chr="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}.bam",
        seed_frac="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}.txt",
    output:
        downsampled_bam="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}x.bam",
        downsampled_bam_bai="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}x.bam.bai",
    log:
        "output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}.log",
    benchmark:
        "benchmarks/target_bams/{sample}_{chrom}_{coverage_val}.tsv"
    shell:
        """
        s=$(cat {input.seed_frac}) 

        samtools view -T {input.ref_fasta_chr} \
        -s $s \
        -bo {output.downsampled_bam} \
        {input.target_bams_chr} 2> {log}

        samtools index {output.downsampled_bam}
        """


rule prepare_list_validation_reference_removal:
    """
    Create a list of target samples which are also in the reference VCF but under another name (have to provide file with two columns for each sample)
    Only needed if samples overlap in target bams and VCF panel
    """
    input:
        ref_val_samples="sample_lists/target_names_reference_remove_published.tsv",
    output:
        ref_val_sample_file="output/GLIMPSE_concordance/reference_panel/ref_val_sample.txt",
    shell:
        """
        cat {input.ref_val_samples} | awk -F '\t' '{{print $2}}' >> {output.ref_val_sample_file}
        """


rule prepare_ref_panel:
    """
    Remove samples, trim-alt alleles created, normalise indels, only keep biallelic snps
    """
    input:
        ref_panel_phased="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        ref_val_sample_file="output/GLIMPSE_concordance/reference_panel/ref_val_sample.txt",
    output:
        ref_concordance_sample_excl=temp(
            "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel.phased.vcf.gz"
        ),
    log:
        "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/{chrom}_ref_panel.tsv"
    shell:
        """
        (
         bcftools view \
         -r {wildcards.chrom} \
         -S ^{input.ref_val_sample_file} \
         --trim-alt-alleles \
         {input.ref_panel_phased} -Ou | \
         bcftools norm -a -Ou | \
         bcftools view -m 2 -M 2 -v snps  \
         -Oz -o {output.ref_concordance_sample_excl}
        ) 2> {log}
        """


rule fill_tags_ref_sample_excl:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING, which is the fraction of missing genotypes)
    """
    input:
        ref_concordance_sample_excl="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel.phased.vcf.gz",
    output:
        ref_concordance_sample_excl_filltags=temp(
            "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags.phased.vcf.gz"
        ),
    log:
        "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/{chrom}_ref_panel_filltags.phased.tsv"
    shell:
        """
        bcftools +fill-tags {input.ref_concordance_sample_excl} \
        -Oz -o {output.ref_concordance_sample_excl_filltags} \
        --threads {threads} \
        -- -t all,F_MISSING 2> {log}
        """


rule filter_sites_ref_sample_excl:
    """
    Filter for missingness (F_MISSING) again, since we removed individuals
    """
    input:
        ref_concordance_sample_excl_filltags="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags.phased.vcf.gz",
    output:
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
        ref_concordance_sample_excl_filltags_filter_csi="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf.csi",
    params:
        f_missing=config["F_MISSING"],
    log:
        "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/{chrom}_ref_panel_filltags_filter.phased.tsv"
    shell:
        """
        bcftools view -i 'F_MISSING<{params.f_missing}' \
        {input.ref_concordance_sample_excl_filltags} \
        --threads {threads} \
        -Ob -o {output.ref_concordance_sample_excl_filltags_filter} 2> {log}
    
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter}
        """


rule prepare_merged_chr_list_ref:
    """ 
    Prepare list to merge chromosomes for reference panel
    """
    input:
        ref_concordance_sample_excl_filltags_filter=expand(
            "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list_ref="output/GLIMPSE_concordance/reference_panel/chr_list.ref_panel_filltags_filter.phased.txt",
    shell:
        """
        ls -v {input.ref_concordance_sample_excl_filltags_filter} >> {output.chr_list_ref}
        """


rule merge_chr_ref:
    """
    Merge chromosomes of filtered reference panel
    """
    input:
        chr_list_ref="output/GLIMPSE_concordance/reference_panel/chr_list.ref_panel_filltags_filter.phased.txt",
    output:
        ref_concordance_sample_excl_filltags_filter_allchrom="output/GLIMPSE_concordance/reference_panel/allchrom_ref_panel_filltags_filter.phased.bcf",
        ref_concordance_sample_excl_filltags_filter_allchrom_csi="output/GLIMPSE_concordance/reference_panel/allchrom_ref_panel_filltags_filter.phased.bcf.csi",
    log:
        "output/GLIMPSE_concordance/reference_panel/allchrom_ref_panel_filltags_filter.phased.bcf.log",
    threads: 4
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list_ref} \
        -Ob -o {output.ref_concordance_sample_excl_filltags_filter_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_allchrom}
        """


rule extract_var_pos_concordance:
    """
    Extract sites from reference panel to use for GLs estimation from bams afterwards
    """
    input:
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
    output:
        ref_panel_sites_vcf="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.tsv.gz",
    log:
        "output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel/{chrom}_ref_panel_sites.phased.tsv"
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


rule compute_GLs_downsampled_samples_concordance:
    """
    Compute GLs of target bams 
    """
    input:
        downsampled_bam="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}x.bam",
        ref_panel_sites_vcf="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.tsv.gz",
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
    output:
        GL_vcf_target_bams="output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}_{coverage_val}x.vcf.gz",
        GL_vcf_target_bams_csi="output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}_{coverage_val}x.vcf.gz.csi",
    log:
        "output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}_{coverage_val}x.log",
    threads: 4
    benchmark:
        "benchmarks/GLs_target_bams/{sample}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        bcftools mpileup -f {input.ref_fasta_chr} -I -E -a 'FORMAT/DP' -T {input.ref_panel_sites_vcf} -r {wildcards.chrom} {input.downsampled_bam} -Ou | \
        bcftools call -Aim -C alleles -T {input.ref_panel_sites_tsv} -Oz -o {output.GL_vcf_target_bams} --threads {threads} 2> {log}
        
        bcftools index -f {output.GL_vcf_target_bams}
        """


rule chunk_spliting_concordance:
    """
    Split chromosome into chunks for imputation
    """
    input:
        ref_panel_sites_vcf="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.vcf.gz",
    output:
        chunks="output/GLIMPSE_concordance/chunks/{chrom}_chunks.txt",
    params:
        window_size=config["window_size"],
        buffer_size=config["buffer_size"],
    log:
        "output/GLIMPSE_concordance/chunks/{chrom}_chunks.log",
    threads: 4
    benchmark:
        "benchmarks/chunks/{chrom}_chunks.tsv"
    shell:
        """
        GLIMPSE_chunk \
        --input {input.ref_panel_sites_vcf} \
        --region {wildcards.chrom} \
        --window-size {params.window_size} --buffer-size {params.buffer_size} \
        --thread {threads} \
        --output {output.chunks} 2> {log}
        """


rule impute_concordance:
    """
    Impute each sample seperately
    """
    input:
        GL_vcf_target_bams="output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom}_{coverage_val}x.vcf.gz",
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
        chunks="output/GLIMPSE_concordance/chunks/{chrom}_chunks.txt",
        gen_map="data/gen_map/{chrom}_average_canFam3.1_modified.tsv",
    output:
        imputed="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x.00.bcf",
        imputed_csi="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x.00.bcf.csi",
    params:
        prefix="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x",
    threads: 2
    log:
        "output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x.log",
    benchmark:
        "benchmarks/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x.tsv"
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


rule ligate_list_concordance:
    """
    Create list of imputed output files for each chunk to merge later
    """
    input:
        chunks="output/GLIMPSE_concordance/chunks/{chrom}_chunks.txt",
        imputed="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x.00.bcf",
    output:
        ligated_list="output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom}_{coverage_val}x.txt",
    params:
        prefix="output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom}_{coverage_val}x",
    shell:
        """
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        """


rule ligate_concordance:
    """
    Merge all imputed chunks
    """
    input:
        ligated_list="output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom}_{coverage_val}x.txt",
    output:
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf",
        ligated_bcf_csi="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.log",
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        GLIMPSE_ligate \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.ligated_bcf}
        """


rule phase_concordance:
    """
    Phase!!!
    """
    input:
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf",
    output:
        phased_bcf="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.bcf",
        phased_bcf_csi="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.log",
    threads: 4
    benchmark:
        "benchmarks/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        GLIMPSE_sample \
        --input {input.ligated_bcf} \
        --solve \
        --output {output.phased_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.phased_bcf}
        """


rule annotate_fields_concordance:
    """
    SOS: This step is essential to retain the fields from the ligated samples in the phased files (not automatically transferred)
    GP is needed to be annotatedm since it's used downstream for recalibrating the INFO score
    """
    input:
        phased_bcf="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.bcf",
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf",
    output:
        phased_vcf_annotate="output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_{coverage_val}x.{chrom}.vcf.gz",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_{coverage_val}x.{chrom}.vcf.gz.log",
    benchmark:
        "benchmarks/GLIMPSE_concordance/GLIMPSE_phased/phased_annotated.{sample}_{chrom}_{coverage_val}x.tsv"
    shell:
        """
        bcftools annotate \
        --annotations {input.ligated_bcf} \
        --columns FORMAT/DS,FORMAT/GP,FORMAT/HS \
        --output-type z \
        --output {output.phased_vcf_annotate} {input.phased_bcf} 2> {log}
        
        bcftools index --tbi {output.phased_vcf_annotate}
        """


#####################################################
#                                                   #
#   Getting validation genotypes (the ground truth) #
#                                                   #
#####################################################


rule estimate_coverage_target_sample:
    """
    TODO add block header
    """
    input:
        target_bams=lambda wildcards: samples_df.loc[wildcards.sample, "bam_path"],
    output:
        cov_depth_cutoff="output/GLIMPSE_concordance/target_bams/{sample}_cov_depth_cutoff.txt",
    log:
        "output/GLIMPSE_concordance/target_bams/{sample}_genome_coverage.txt.log",
    shell:
        """
        (
        #estimate the length of the regions within the bam file 
        LENGTH_GENOME=$(samtools view -H {input.target_bams} \
        | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{{sum+=$1}} END {{print sum}}')

        #estimate sample coverage
        COV_GENOME=$(samtools depth -a {input.target_bams}  \
        | awk -v var="$LENGTH_GENOME" '{{sum+=$3}} END {{ print sum/var}}')

        #avoid "invalid arithmetic operator" error when trying to add variables in Bash script
        C=$(echo $COV_GENOME/3 | bc)
        D=$(echo $COV_GENOME*2 | bc)

        #find max of cov/3 and 8 and add twice avg coverage to the same file
        echo $(($C>8 ? $C : 8)) > {output.cov_depth_cutoff}
        echo $D >> {output.cov_depth_cutoff}
        ) 2> {log}
        """


rule prepare_validation_samples_filt:
    """
    Take target bams (initial coverage) and call genotypes on the same site as the filtered reference panel, then filter sites based on Sousa da Mota 2022
    """
    input:
        target_bams_chr="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}.bam",
        ref_panel_sites_vcf="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.vcf.gz",
        ref_panel_sites_tsv="output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_sites.phased.tsv.gz",
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
        cov_depth_cutoff="output/GLIMPSE_concordance/target_bams/{sample}_cov_depth_cutoff.txt",
    output:
        validation_sample_filt="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp.bcf",
        validation_sample_filt_csi="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp.bcf.csi",
    log:
        "output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp.log",
    threads: 10
    benchmark:
        "benchmarks/validation_bams/{sample}_{chrom}_validation_filt_qual_dp.tsv"
    shell:
        """
        min=$(sed -n '1p' {input.cov_depth_cutoff})
        max=$(sed -n '2p' {input.cov_depth_cutoff})
        
        bcftools mpileup -f {input.ref_fasta_chr} \
        -I -E \
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
        -T {input.ref_panel_sites_vcf} \
        -r {wildcards.chrom} \
        -q 30 -Q 30 -C 50 \
        --threads {threads} \
        {input.target_bams_chr} -Ou | \
        bcftools call -Aim -C alleles \
        -T {input.ref_panel_sites_tsv} --threads {threads} -Ou | \
        bcftools filter -i "%QUAL>=30 && FORMAT/DP<$max && FORMAT/DP>$min" \
        -Ob -o {output.validation_sample_filt} 2> {log}
        
        bcftools index -f {output.validation_sample_filt}
        """


rule prepare_validation_samples_filt_allelic:
    """
    Filter sites for allelic imbalance
    """
    input:
        validation_sample_filt="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp.bcf",
        ref_fasta_chr="output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom}.fasta",
    output:
        validation_sample_filt_allelic="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp_ab.bcf",
        validation_sample_filt_csi_allelic="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp_ab.bcf.csi",
    log:
        "output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp_ab.log",
    threads: 10
    shell:
        """
        bcftools view --exclude 'GT="het" && ((INFO/AD[1] / INFO/DP < 0.15) || (INFO/AD[1] / INFO/DP > 0.85))' \
        {input.validation_sample_filt} \
        --threads {threads} \
        -Ob -o {output.validation_sample_filt_allelic} 2> {log}
        
        bcftools index -f {output.validation_sample_filt_allelic}
        """


rule prepare_merged_chr_list_validation:
    """ 
    Prepare list to merge chromosomes for validation
    """
    input:
        validation_sample_filt_allelic=expand(
            "output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp_ab.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list_validation="output/GLIMPSE_concordance/validation_bams/chr_list.{sample}_validation_filt_qual_dp_ab.txt",
    shell:
        """
        ls -v {input.validation_sample_filt_allelic} >> {output.chr_list_validation}
        """


rule merge_chr_validation:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        chr_list_validation="output/GLIMPSE_concordance/validation_bams/chr_list.{sample}_validation_filt_qual_dp_ab.txt",
    output:
        validation_sample_filt_allelic_allchrom="output/GLIMPSE_concordance/validation_bams/{sample}_allchrom_validation_filt_qual_dp_ab.bcf",
        validation_sample_filt_allelic_allchrom_csi="output/GLIMPSE_concordance/validation_bams/{sample}_allchrom_validation_filt_qual_dp_ab.bcf.csi",
    log:
        "output/GLIMPSE_concordance/validation_bams/{sample}_allchrom_validation_filt_qual_dp_ab.bcf.log",
    threads: 8
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list_validation} \
        -Ob -o {output.validation_sample_filt_allelic_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.validation_sample_filt_allelic_allchrom}
        """


#############################
#                           #
#  Run GLIMPSE concordance  #
#                           #
#############################


rule filter_info_score:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        ligated_bcf="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf",
    output:
        info_imputed_info="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf",
        info_imputed_info_csi="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf.csi",
    params:
        info_val="{info_cutoff}",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.log",
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


rule prepare_merged_chr_list_concordance:
    """ 
    Prepare list to merge chromosomes for reference, imputed and validation
    """
    input:
        info_imputed_info=expand(
            "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom}_{coverage_val}x-INFO_{info_cutoff}.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/chr_list.{sample}_{coverage_val}x-INFO_{info_cutoff}.txt",
    shell:
        """
        ls -v {input.info_imputed_info} >> {output.chr_list}
        """


rule merge_chr_concordance:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/chr_list.{sample}_{coverage_val}x-INFO_{info_cutoff}.txt",
    output:
        info_imputed_info_allchrom="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf",
        info_imputed_info_allchrom_csi="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.log",
    threads: 8
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Ob -o {output.info_imputed_info_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.info_imputed_info_allchrom}
        """


rule get_ID_for_targets_allchrom:
    """
    TODO add block header
    """
    input:
        info_imputed_info_allchrom="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf",
    output:
        sm_samples="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/sm_merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.txt",
    shell:
        """
        bcftools query -l {input.info_imputed_info_allchrom} > {output.sm_samples}
        """


rule prepare_concordance_lst_info_score_filtered_allchrom:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter_allchrom="output/GLIMPSE_concordance/reference_panel/allchrom_ref_panel_filltags_filter.phased.bcf",
        validation_sample_filt_allelic_allchrom="output/GLIMPSE_concordance/validation_bams/{sample}_allchrom_validation_filt_qual_dp_ab.bcf",
        info_imputed_info_allchrom="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.bcf",
    output:
        concordance_lst_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.lst",
    shell:
        """
        echo "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chr31,chr32,chr33,chr34,chr35,chr36,chr37,chr38" {input.ref_concordance_sample_excl_filltags_filter_allchrom} {input.validation_sample_filt_allelic_allchrom} {input.info_imputed_info_allchrom} > {output.concordance_lst_info_score_filtered}
        """


rule GLIMPSE_concordance_info_score_filtered_allchrom:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.lst",
        sm_samples="output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/sm_merged_ligated.{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}.txt",
    output:
        concordance_output_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz",
        concordance_output_discordance_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz",
    params:
        prefix="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered",
    log:
        "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.log",
    threads: 8
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


rule plot_rsquare_accuracy_filtered_allchrom:
    """
    Plot accuracy 
    """
    input:
        concordance_output_info_score_1="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_0.8_filtered.rsquare.grp.txt.gz",
        concordance_output_info_score_2="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_0.9_filtered.rsquare.grp.txt.gz",
        concordance_output_info_score_3="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_0.95_filtered.rsquare.grp.txt.gz",
        concordance_output_info_score_4="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_0.0_filtered.rsquare.grp.txt.gz",
    output:
        plot="output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_{sample}_allchrom_{coverage_val}x_filtered.png",
    params:
        chr="all autosomes",
        name="{sample}",
        cov="{coverage_val}",
    script:
        # TODO replace with `shell: Rscript script/*.R --arg1 etc
        "../scripts/rsquare_accuracy.R"


rule plot_rsquare_accuracy_filtered_allchrom2_prep:
    """
    Prepare files for main plot
    """
    input:
        concordance_output_info_score_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz",
    output:
        concordance_output_info_score_filtered_mod="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
    shell:
        """
        zcat {input.concordance_output_info_score_filtered} | awk -v FS=' ' -v OFS=' ' '{{$7={wildcards.coverage_val}}} {{$8={wildcards.info_cutoff}}} {{$9="{wildcards.sample}"}} 1' > {output.concordance_output_info_score_filtered_mod}
        """


rule plot_rsquare_accuracy_filtered_allchrom2:
    """
    Plot accuracy 
    """
    input:
        concordance_output=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
            sample=["NGDG", "PortauChoix", "TRF.05.05", "CGG33"],
            coverage_val=["0.5", "1", "2"],
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        plot="output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_allchrom_filtered-main.png",
    params:
        files=lambda wildcards, input: ",".join(input.concordance_output),
    shell:
        """
        Rscript scripts/rsquare_accuracy-main.R \
        {params.files} \
        {output.plot}
        """


rule plot_accuray_per_sample:
    """
    TODO add block header
    """
    input:
        concordance_output=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare-mod.grp.txt.gz",
            coverage_val=COVERAGE_VAL,
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        plot="output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_allchrom_filtered-{sample}.png",
    params:
        files=lambda wildcards, input: ",".join(input.concordance_output),
        name="{sample}",
        info_sample=lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
    shell:
        """
        Rscript scripts/rsquare_accuracy_per_sample.R \
        {params.files} \
        {params.name} \
        {params.info_sample} \
        {output.plot}
        """


rule prepare_concordance_output_filt_allchrom:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz",
    output:
        concordance_output_discordance_filtered_temp=temp(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/temp_concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.txt"
        ),
        concordance_output_discordance_filtered_prep="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.txt",
    shell:
        """
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        """


rule merge_concordance_output_filt_allchrom:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered.txt",
            coverage_val=COVERAGE_VAL,
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        concordance_output_discordance_filtered_per_sample="output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_filtered.txt",
    shell:
        """
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        """


rule plot_discordance_filt_allchrom:
    """
    Plot genotype discordances
    """
    input:
        concordance_output_discordance_filtered_per_sample=expand(
            "output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_allchrom_filtered.txt",
            sample=SAMPLE,
            allow_missing=True,
        ),
        concordance_metadata="sample_lists/concordance_bams_published.tsv",
    params:
        path_script="scripts",
        discordance_phased=lambda wildcards, input: ",".join(
            input.concordance_output_discordance_filtered_per_sample
        ),
    output:
        discordance_dogs_full="output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_allchrom_filtered_dogs_full.png",
        discordance_dogs="output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_allchrom_filtered_dogs.png",
        discordance_wolves_full="output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_allchrom_filtered_wolves_full.png",
        discordance_wolves="output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_allchrom_filtered_wolves.png",
    shell:
        """
        Rscript {params.path_script}/genotype_discordance.R \
        {params.discordance_phased} \
        {input.concordance_metadata} \
        {output.discordance_dogs_full} \
        {output.discordance_dogs} \
        {output.discordance_wolves_full} \
        {output.discordance_wolves}
        """
