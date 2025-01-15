############################
# Filter validation sample #
############################

rule estimate_coverage_target_sample:
    input:
        target_bams = lambda wildcards: samples_df.loc[wildcards.sample, "bam_path"]
    output:
        cov_depth_cutoff = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_cov_depth_cutoff.txt'
    log:
        '{path}/output/GLIMPSE_concordance/target_bams/{sample}_genome_coverage.txt.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
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
        '''

rule prepare_validation_samples_filt: 
    """
    Take target bams (initial coverage) and call genotypes on the same site as the filtered reference panel, then filter sites based on Sousa da Mota 2022
    """
    input:
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam',
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz',
        ref_panel_sites_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.tsv.gz',
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta',
        cov_depth_cutoff = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_cov_depth_cutoff.txt'
    output:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.bcf',
        validation_sample_filt_csi = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.tsv'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        min=$(sed -n '1p' {input.cov_depth_cutoff})
        max=$(sed -n '2p' {input.cov_depth_cutoff})
        
        bcftools mpileup -f {input.ref_fasta_chr} \
        -I -E \
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
        -T {input.ref_panel_sites_vcf} \
        -r {wildcards.chrom_con} \
        -q 30 -Q 30 -C 50 \
        --threads {threads} \
        {input.target_bams_chr} -Ou | \
        bcftools call -Aim -C alleles \
        -T {input.ref_panel_sites_tsv} --threads {threads} -Ou | \
        bcftools filter -i "%QUAL>=30 && FORMAT/DP<$max && FORMAT/DP>$min" \
        -Ob -o {output.validation_sample_filt} 2> {log}
        
        bcftools index -f {output.validation_sample_filt}
        '''

rule plot_ad_dp_prep:
    input:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.bcf'
    output:
        validation_sample_filt_het = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.bcf',
        validation_sample_filt_het_sites = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.txt'
    shell:
        '''
        bcftools view --include 'GT="het"' {input.validation_sample_filt} -Ob -o {output.validation_sample_filt_het}
        bcftools query -f '%CHROM %POS %INFO/DP %INFO/AD{{1}}\n' {output.validation_sample_filt_het} > {output.validation_sample_filt_het_sites}
        '''

rule plot_ad_dp:
    input:
        validation_sample_filt_het_sites = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.txt'
    output:
        validation_sample_filt_het_plot = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.png'
    script:
        "../scripts/allelic_depth_read_depth_plot.R"
        

rule prepare_validation_samples_filt_allelic: 
    """
    Filter sites for allelic imbalance
    """
    input:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.bcf',
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta' 
    output:
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        validation_sample_filt_csi_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view --exclude 'GT="het" && ((INFO/AD[1] / INFO/DP < 0.15) || (INFO/AD[1] / INFO/DP > 0.85))' \
        {input.validation_sample_filt} \
        --threads {threads} \
        -Ob -o {output.validation_sample_filt_allelic} 2> {log}
        
        bcftools index -f {output.validation_sample_filt_allelic}
        '''


rule filter_info_score:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
    output:
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf',
        info_imputed_info_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf.csi'
    params:
        info_val = '{info_cutoff}'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.log'
    threads: 8
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view {input.ligated_bcf} \
        --include 'INFO/INFO >= {params.info_val}' \
        --threads {threads} \
        -Ob -o {output.info_imputed_info} 2> {log}

        bcftools index -f {output.info_imputed_info}
        '''

rule get_ID_for_targets:
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
    output:
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}_{coverage_val}x.txt'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools query -l {input.ligated_bcf} > {output.sm_samples}
        '''

rule prepare_concordance_lst_info_score_filtered:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf',
    output:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.lst'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter} {input.validation_sample_filt_allelic} {input.info_imputed_info} > {output.concordance_lst_info_score_filtered}
        '''

rule GLIMPSE_concordance_info_score_filtered:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}_{coverage_val}x.txt'
    output:
        concordance_output_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz',
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered'
    log:
        '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.log'
    threads: 8
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        {glimpse_concordance} \
        --input {input.concordance_lst_info_score_filtered} \
        --minDP 8 \
        --output {params.prefix} \
        --minPROB 0.9 \
        --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
        --sample {input.sm_samples} \
        --af-tag AF \
        --thread {threads} 2> {log}
        '''

rule plot_rsquare_accuracy_filtered:
    """
    Plot accuracy 
    """
    input:
        concordance_output_info_score_1 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.8_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_2 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.9_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_3 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.95_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_4 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.0_filtered.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x_filtered.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = '{coverage_val}'
    script:
        "../scripts/rsquare_accuracy.R"


rule prepare_concordance_output_filt:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    output:
        concordance_output_discordance_filtered_temp = temp('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/temp_concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.txt'),
        concordance_output_discordance_filtered_prep = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.txt'
    shell:
        '''
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        '''

rule merge_concordance_output_filt:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.txt', coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        concordance_output_discordance_filtered_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filtered.txt'
    shell:
        '''
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        '''

rule plot_discordance_filt:
    """
    Plot genotype discordances
    """
    input:
        concordance_output_discordance_filtered_per_sample = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filtered.txt', sample=SAMPLE, allow_missing=True),
        concordance_metadata = config['bam_targets']
    params:
        path_script = '{path}/scripts',
        discordance_phased=lambda wildcards, input: ','.join(input.concordance_output_discordance_filtered_per_sample),
    output:
        discordance_dogs_full = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_dogs_full.png',
        discordance_dogs = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_dogs.png',
        discordance_wolves_full = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_wolves_full.png',
        discordance_wolves = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_wolves.png',
    shell:
        '''
        Rscript {params.path_script}/genotype_discordance.R \
        {params.discordance_phased} \
        {input.concordance_metadata} \
        {output.discordance_dogs_full} \
        {output.discordance_dogs} \
        {output.discordance_wolves_full} \
        {output.discordance_wolves}
        '''