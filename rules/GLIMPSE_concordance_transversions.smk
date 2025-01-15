###############################
# Running GLIMPSE concordance #
###############################

rule transversions_validation_concordance:
    """
    Only take transversions from filtered validation data
    """
    input:
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf'
    output:
        tranversion_sites_allelic = '{path}/output/GLIMPSE_concordance/validation_bams_transversions/{sample}_{chrom_con}_validation_filt_transversions.tsv.gz',
        validation_transversions_allelic = '{path}/output/GLIMPSE_concordance/validation_bams_transversions/{sample}_{chrom_con}_validation_filt_transversions.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams_transversions/{sample}_{chrom_con}_validation_filt_transversions.log'
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.validation_sample_filt_allelic} | bgzip -c > {output.tranversion_sites_allelic}

        tabix -s1 -b2 -e2 {output.tranversion_sites_allelic}

        bcftools view {input.validation_sample_filt_allelic} \
        --regions-file {output.tranversion_sites_allelic} \
        --threads {threads} \
        -Ob -o {output.validation_transversions_allelic} 2> {log}

        bcftools index -f {output.validation_transversions_allelic}
        '''

rule filter_transversions_imputed:
    input:
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf',
    output:
        tranversion_sites = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered_transversions/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_transversions.tsv.gz',
        imputed_transversions = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered_transversions/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered_transversions/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf.log'
    #conda:
    #    '../envs/environment.yaml'
    threads: 8
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.info_imputed_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.info_imputed_info} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Ob -o {output.imputed_transversions} 2> {log}

        bcftools index -f {output.imputed_transversions}
        '''

rule prepare_concordance_lst_info_score_filtered_transversions:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_transversions_allelic = '{path}/output/GLIMPSE_concordance/validation_bams_transversions/{sample}_{chrom_con}_validation_filt_transversions.bcf',
        imputed_transversions = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered_transversions/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_transversions.bcf'
    output:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.lst'
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter} {input.validation_transversions_allelic} {input.imputed_transversions} > {output.concordance_lst_info_score_filtered}
        '''

rule GLIMPSE_concordance_info_score_filtered_transversions:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}_{coverage_val}x.txt'
    output:
        concordance_output_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.rsquare.grp.txt.gz',
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions'
    log:
        '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.log'
    threads: 8
    #conda:
    #    '../envs/environment.yaml'
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

rule plot_rsquare_accuracy_filtered_transversions:
    """
    Plot accuracy 
    """
    input:
        concordance_output_info_score_1 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.8_filtered_transversions.rsquare.grp.txt.gz',
        concordance_output_info_score_2 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.9_filtered_transversions.rsquare.grp.txt.gz',
        concordance_output_info_score_3 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.95_filtered_transversions.rsquare.grp.txt.gz',
        concordance_output_info_score_4 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.0_filtered_transversions.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x_filtered_transversions.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = '{coverage_val}'
    script:
        "../scripts/rsquare_accuracy.R"


rule prepare_concordance_output_filt_transversions:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.error.spl.txt.gz'
    output:
        concordance_output_discordance_filtered_temp = temp('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/temp_concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.txt'),
        concordance_output_discordance_filtered_prep = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.txt'
    shell:
        '''
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        '''

rule merge_concordance_output_filt_transversions:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.txt', coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        concordance_output_discordance_filtered_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_filtered_transversions.txt'
    shell:
        '''
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        '''

rule plot_discordance_filt_transversions:
    """
    Plot genotype discordances
    """
    input:
        concordance_output_discordance_filtered_per_sample = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_{chrom_con}_filtered_transversions.txt', sample=SAMPLE, allow_missing=True),
        concordance_metadata = config['bam_targets']
    params:
        path_script = '{path}/scripts',
        discordance_phased=lambda wildcards, input: ','.join(input.concordance_output_discordance_filtered_per_sample),
    output:
        discordance_dogs_full = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_{chrom_con}_filtered_dogs_full_transversions.png',
        discordance_dogs = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_{chrom_con}_filtered_dogs_transversions.png',
        discordance_wolves_full = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_{chrom_con}_filtered_wolves_full_transversions.png',
        discordance_wolves = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_{chrom_con}_filtered_wolves_transversions.png',
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