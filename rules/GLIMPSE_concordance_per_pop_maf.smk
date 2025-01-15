###############################
# Running GLIMPSE concordance #
###############################


rule get_ref_samples:
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        ref_meta = config['reference_panel_metadata']
    output:    
        ref_concordance_sample_excl_filltags_filter_names_temp = temp('{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_names_temp.txt'),
        ref_concordance_sample_excl_filltags_filter_names = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_names.txt',
    params:
        dog='Dogs',
        wolf='Wolves',
        dog_pop='dogs',
        wolf_pop='wolves',
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        #get names from reference vcf
        bcftools query -l {input.ref_concordance_sample_excl_filltags_filter} > {output.ref_concordance_sample_excl_filltags_filter_names_temp}

        #match population info from ref metadata
        awk 'NR==FNR{{A[$1]=$3;next}}A[$1]{{$2=A[$1] FS $2;print}}' {input.ref_meta} FS='\t' OFS='\t' {output.ref_concordance_sample_excl_filltags_filter_names_temp} > names_temp1.txt

        #replace all dogs with "dogs" and wolves with "wolves" for the population assigment
        del={params.dog}; sed "s/[^[:blank:]]*${{del}}[^[:blank:]]*/{params.dog_pop}/g" names_temp1.txt > names_temp2.txt 
        del={params.wolf}; sed "s/[^[:blank:]]*${{del}}[^[:blank:]]*/{params.wolf_pop}/g" names_temp2.txt > {output.ref_concordance_sample_excl_filltags_filter_names}

        #remove unwanted files:
        rm names_temp1.txt names_temp2.txt 
        '''

rule bcftools_fill_tags_pop:
    input:
        ref_concordance_sample_excl_filltags_filter_names = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_names.txt',
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf'
    output:
        ref_concordance_sample_excl_filltags_filter_pop = '{path}/output/GLIMPSE_concordance2/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf'
    log:
        '{path}/output/GLIMPSE_concordance2/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf.log'
    threads: 8
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools +fill-tags {input.ref_concordance_sample_excl_filltags_filter} \
        -Ob -o {output.ref_concordance_sample_excl_filltags_filter_pop} \
        --threads {threads} \
        -- -S {input.ref_concordance_sample_excl_filltags_filter_names} \
        -t all,F_MISSING 2> {log}

        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_pop}
        '''


rule prepare_concordance_lst_info_score_filtered2:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter_pop = '{path}/output/GLIMPSE_concordance2/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf',
    output:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.lst'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter_pop} {input.validation_sample_filt_allelic} {input.info_imputed_info} > {output.concordance_lst_info_score_filtered}
        '''

rule GLIMPSE_concordance_info_score_filtered2:
    """
    Run GLIMPSE concordance specifying the target sample we want  ### PAY ATTENTION TO THE AF-TAG FIELD!
    """
    input:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}_{coverage_val}x.txt'
    output:
        concordance_output_info_score_filtered = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.rsquare.grp.txt.gz',
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}',
        af_pop = '{canid_subset}'
    log:
        '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.log'
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
        --af-tag AF_{params.af_pop} \
        --thread {threads} 2> {log}
        '''

rule plot_rsquare_accuracy_filtered2:
    """
    Plot accuracy 
    """
    input:
        concordance_output_info_score_1 = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.8_filtered_{canid_subset}.rsquare.grp.txt.gz',
        concordance_output_info_score_2 = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.9_filtered_{canid_subset}.rsquare.grp.txt.gz',
        concordance_output_info_score_3 = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.95_filtered_{canid_subset}.rsquare.grp.txt.gz',
        concordance_output_info_score_4 = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.0_filtered_{canid_subset}.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance2/plots/glimpse_concordance/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x_filtered_{canid_subset}.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = '{coverage_val}'
    script:
        "../scripts/rsquare_accuracy.R"


rule prepare_concordance_output_filt2:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.error.spl.txt.gz'
    output:
        concordance_output_discordance_filtered_temp = temp('{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/temp_concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.txt'),
        concordance_output_discordance_filtered_prep = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.txt'
    shell:
        '''
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        '''

rule merge_concordance_output_filt2:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep = expand('{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered_{canid_subset}.txt', coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        concordance_output_discordance_filtered_per_sample = '{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filtered_{canid_subset}.txt'
    shell:
        '''
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        '''

rule plot_discordance_filt2:
    """
    Plot genotype discordances
    """
    input:
        concordance_output_discordance_filtered_per_sample = expand('{path}/output/GLIMPSE_concordance2/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filtered_{canid_subset}.txt', sample=SAMPLE, allow_missing=True),
        concordance_metadata = config['bam_targets']
    params:
        path_script = '{path}/scripts',
        discordance_phased=lambda wildcards, input: ','.join(input.concordance_output_discordance_filtered_per_sample),
        #name = '{sample}',
        #info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
        #species = lambda wildcards: samples_df.loc[wildcards.sample, "Wolf_Dog_PCA"]
    output:
        discordance_dogs_full = '{path}/output/GLIMPSE_concordance2/plots/glimpse_concordance/concordance_{chrom_con}_filtered_{canid_subset}_dogs_full.png',
        discordance_dogs = '{path}/output/GLIMPSE_concordance2/plots/glimpse_concordance/concordance_{chrom_con}_filtered_{canid_subset}_dogs.png',
        discordance_wolves_full = '{path}/output/GLIMPSE_concordance2/plots/glimpse_concordance/concordance_{chrom_con}_filtered_{canid_subset}_wolves_full.png',
        discordance_wolves = '{path}/output/GLIMPSE_concordance2/plots/glimpse_concordance/concordance_{chrom_con}_filtered_{canid_subset}_wolves.png',
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