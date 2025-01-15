#######################################################################
# Preperation for and running GLIMPSE concordance filtered validation #
#######################################################################

rule prepare_validation_samples_filt:
    """
    Take target bams (initial coverage) and call genotypes on the same site as the filtered reference panel, then filter sites based on Sousa da Mota 2022
    """
    input:
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam',
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz',
        ref_panel_sites_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.tsv.gz',
        ref_fasta_chr = '{path}/output/reference_genome/CanFam31_{chrom_con}.fasta',
    output:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.bcf',
        validation_sample_filt_csi = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.tsv'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools mpileup -f {input.ref_fasta_chr} \
        -I -E -a 'FORMAT/DP' \
        -T {input.ref_panel_sites_vcf} \
        -r {wildcards.chrom_con} \
        -q 30 -Q 20 \
        {input.target_bams_chr} -Ou --threads {threads} | bcftools call -Aim -C alleles \
        -T {input.ref_panel_sites_tsv} --threads {threads} -Ou | \
        bcftools filter -e "%QUAL<30" \
        -Ob -o {output.validation_sample_filt} 2> {log}
        
        bcftools index -f {output.validation_sample_filt}
        '''

rule prepare_concordance_lst_info_score_filt:
    """
    Prepare the lst files required to run GLIMPSE_concordance using the filtered dataset !!!! I might have to filter out the sites in the imputed bcf with the same sites as in the validation bcf!!!
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.bcf',
        info_filt_ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}.bcf'
    output:
        concordance_lst_info_score_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.lst'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter} {input.validation_sample_filt} {input.info_filt_ligated_bcf} > {output.concordance_lst_info_score_filt}
        '''

rule glimpse_concordance_info_score_filt:
    """
    Run GLIMPSE concordance using the filtered validation specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}.txt'
    output:
        concordance_output_info_score_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.rsquare.grp.txt.gz',
        concordance_output_discordance_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt'
    log:
        '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.tsv'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        {glimpse_concordance} \
        --input {input.concordance_lst_info_score_filt} \
        --minDP 8 \
        --output {params.prefix} \
        --minPROB 0.9 \
        --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
        --sample {input.sm_samples} \
        --af-tag AF \
        --thread {threads} 2> {log}
        '''

rule plot_rsquare_accuracy_filt:
    input:
        concordance_output_info_score_1_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_0.8_filt.rsquare.grp.txt.gz',
        concordance_output_info_score_2_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_0.9_filt.rsquare.grp.txt.gz',
        concordance_output_info_score_3_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_0.95_filt.rsquare.grp.txt.gz',
        concordance_output_info_score_4_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_0.0_filt.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x_filt.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = '{coverage_val}'
    script:
        "../scripts/rsquare_accuracy.R"


rule merge_plot_rsquare_accuracy_filt:
    input:
        plot = expand('{path}/output/GLIMPSE_concordance/plots/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x_filt.png', sample=SAMPLE, chrom_con=CHROM_CON, allow_missing=True)
    output:
        plots_per_coverage = '{path}/output/GLIMPSE_concordance/plots/rsquare_accuracy_{chrom_con}_{coverage_val}x_filt.pdf'
    shell:
        '''
        convert {input.plot} {output.plots_per_coverage}
        '''
        
rule prepare_concordance_lst_info_score_sample_filtered:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.bcf',
        info_imputed_sample_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}.bcf'
    output:
        concordance_lst_info_score_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.lst'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter} {input.validation_sample_filt} {input.info_imputed_sample_info} > {output.concordance_lst_info_score_sample_filtered}
        '''

rule GLIMPSE_concordance_sample_info_score_filtered:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}.txt'
    output:
        concordance_output_info_score_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz',
        concordance_output_discordance_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered'
    log:
        '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        {glimpse_concordance} \
        --input {input.concordance_lst_info_score_sample_filtered} \
        --minDP 8 \
        --output {params.prefix} \
        --minPROB 0.9 \
        --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
        --sample {input.sm_samples} \
        --af-tag AF \
        --thread {threads} 2> {log}
        '''

rule plot_rsquare_accuracy_sample_filtered:
    input:
        concordance_output_info_score_1 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_0.8_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_2 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_0.9_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_3 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_0.95_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_4 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_0.0_filtered.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x-only_sample_filtered.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = '{coverage_val}'
    script:
        "../scripts/rsquare_accuracy.R"

rule merge_plot_rsquare_accuracy_sample_filt:
    input:
        plot = expand('{path}/output/GLIMPSE_concordance/plots/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x-only_sample_filtered.png', sample=SAMPLE, chrom_con=CHROM_CON, allow_missing=True)
    output:
        plots_per_coverage = '{path}/output/GLIMPSE_concordance/plots/rsquare_accuracy_{chrom_con}_{coverage_val}x-only_sample_filtered.pdf'
    shell:
        '''
        convert {input.plot} {output.plots_per_coverage}
        '''



rule prepare_concordance_output_filt:
    input:
        concordance_output_discordance_filt = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.error.spl.txt.gz'
    output:
        concordance_output_discordance_filt_temp = temp('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/temp_concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.txt'),
        concordance_output_discordance_filt_prep = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.txt'
    shell:
        '''
        zcat {input.concordance_output_discordance_filt} | sed -n '3p' >> {output.concordance_output_discordance_filt_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filt_temp} > {output.concordance_output_discordance_filt_prep}
        '''

rule merge_concordance_output_filt:
    input:
        concordance_output_discordance_filt_prep = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_filt.txt', chrom=CHROM, coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        concordance_output_discordance_filt_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filt.txt'
    shell:
        '''
        cat {input.concordance_output_discordance_filt_prep} > {output.concordance_output_discordance_filt_per_sample}
        '''

rule plot_discordance_filt:
    input:
        concordance_output_discordance_filt_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filt.txt'
    params:
        chr = '{chrom_con}',
        name = '{sample}'
    output:
        concordance_output_discordance_filt_per_sample_plot = '{path}/output/GLIMPSE_concordance/plots/concordance_{sample}_{chrom_con}_filt.png'
    script:
        "../scripts/genotype_discordance.R"


rule prepare_concordance_output_filt_sample_only:
    input:
        concordance_output_discordance_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    output:
        concordance_output_discordance_sample_filtered_temp = temp('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/temp_concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.txt'),
        concordance_output_discordance_sample_filtered_prep = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.txt'
    shell:
        '''
        zcat {input.concordance_output_discordance_sample_filtered} | sed -n '3p' >> {output.concordance_output_discordance_sample_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_sample_filtered_temp} > {output.concordance_output_discordance_sample_filtered_prep}
        '''

rule merge_concordance_output_filt_sample_only:
    input:
        concordance_output_discordance_sample_filtered_prep = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-only_sample-INFO_{info_cutoff}_filtered.txt', chrom=CHROM, coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        concordance_output_discordance_sample_filtered_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-only_sample_filtered.txt'
    shell:
        '''
        cat {input.concordance_output_discordance_sample_filtered_prep} > {output.concordance_output_discordance_sample_filtered_per_sample}
        '''

rule plot_discordance_sample_only_filt:
    input:
        concordance_output_discordance_sample_filtered_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-only_sample_filtered.txt'
    params:
        chr = '{chrom_con}',
        name = '{sample}'
    output:
        concordance_output_discordance_sample_per_sample_plot = '{path}/output/GLIMPSE_concordance/plots/concordance_{sample}_{chrom_con}-only_sample_filtered.png'
    script:
        "../scripts/genotype_discordance.R"





















