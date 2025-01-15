##############################################
# Running GLIMPSE concordance for HC imputed #
##############################################


rule compute_GLs_HC_samples_concordance:
    """
    Compute GLs of HC target bams 
    """
    input:
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam',
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz',
        ref_panel_sites_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.tsv.gz',
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta'
    output:
        GL_vcf_HC_target_bams = '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}.vcf.gz',
        GL_vcf_HC_target_bams_csi = '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}.vcf.gz.csi'
    log:
        '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/GLs_target_bams/{sample}_{chrom_con}.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools mpileup -f {input.ref_fasta_chr} -I -E -a 'FORMAT/DP' -T {input.ref_panel_sites_vcf} -r {wildcards.chrom_con} {input.target_bams_chr} -Ou | \
        bcftools call -Aim -C alleles -T {input.ref_panel_sites_tsv} -Oz -o {output.GL_vcf_HC_target_bams} --threads {threads} 2> {log}
        
        bcftools index -f {output.GL_vcf_HC_target_bams}
        '''

rule impute_HC_concordance:
    """
    Impute all samples at the same time!!!
    """
    input:
        GL_vcf_HC_target_bams = '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}.vcf.gz',
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        chunks = '{path}/output/GLIMPSE_concordance/chunks/{chrom_con}_chunks.txt'
    output:
        imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}.00.bcf',
        imputed_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}.00.bcf.csi'
    params:    
        gen_map_path = config['gen_map_path'],
        gen_map_files = config['gen_map_files'],
        prefix = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}'
    threads: 2
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}.log'
    benchmark:
        '{path}/benchmarks/GLIMPSE_imputed/{sample}_{chrom_con}.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT={params.prefix}.${{ID}}.bcf
            {glimpse_impute} \
            --input {input.GL_vcf_HC_target_bams} --reference {input.ref_concordance_sample_excl_filltags_filter} --map {params.gen_map_path}{params.gen_map_files} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} --output ${{OUT}} \
            --thread {threads}
            bcftools index -f ${{OUT}}
        done < {input.chunks} 2> {log}
        '''

rule ligate_HC_list_concordance:
    """
    Create list of imputed output files for each chunk to merge later
    """
    input:
        chunks = '{path}/output/GLIMPSE_concordance/chunks/{chrom_con}_chunks.txt',
        imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}.00.bcf'
    output:
        ligated_list = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom_con}.txt'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        '''

rule ligate_HC_concordance:
    """
    Merge all imputed chunks
    """
    input:
        ligated_list = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom_con}.txt'
    output:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.bcf',
        ligated_bcf_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        {glimpse_ligate} \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.ligated_bcf}
        '''

rule phase_HC_concordance:
    """
    Phase!!!
    """
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.bcf'
    output:
        phased_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}.bcf',
        phased_bcf_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/GLIMPSE_phased/phased.{sample}_{chrom_con}.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        {glimpse_sample} \
        --input {input.ligated_bcf} --solve --output {output.phased_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.phased_bcf}
        '''
        
rule filter_info_score_HC:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.bcf',
    output:
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}-INFO_{info_cutoff}.bcf',
        info_imputed_info_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}-INFO_{info_cutoff}.bcf.csi'
    params:
        info_val = '{info_cutoff}'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}-INFO_{info_cutoff}.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view {input.ligated_bcf} \
        --include 'INFO/INFO >= {params.info_val}' \
        --threads {threads} \
        -Ob -o {output.info_imputed_info} 2> {log}

        bcftools index -f {output.info_imputed_info}
        '''

rule get_ID_for_targets_HC:
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.bcf',
    output:
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}.txt'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools query -l {input.ligated_bcf} > {output.sm_samples}
        '''

rule prepare_HC_concordance_lst_info_score_filtered:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}-INFO_{info_cutoff}.bcf'
    output:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_{info_cutoff}_filtered.lst'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter} {input.validation_sample_filt_allelic} {input.info_imputed_info} > {output.concordance_lst_info_score_filtered}
        '''

rule GLIMPSE_concordance_HC_info_score_filtered:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_{info_cutoff}_filtered.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}.txt'
    output:
        concordance_output_info_score_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz',
        concordance_output_discordance_sample_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_{info_cutoff}_filtered'
    log:
        '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_{info_cutoff}_filtered.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
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

rule plot_rsquare_accuracy_HC_sample_filtered:
    """
    Plot accuracy 
    """
    input:
        concordance_output_info_score_1 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_0.8_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_2 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_0.9_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_3 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_0.95_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_4 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}-INFO_0.0_filtered.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_{sample}_{chrom_con}_filtered.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = 'HC_imputed'
    #conda:
        #'../envs/r4.3.1.yaml'
    script:
        "../scripts/rsquare_accuracy.R"
