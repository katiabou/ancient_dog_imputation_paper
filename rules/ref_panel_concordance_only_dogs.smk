######################################################
# Prepare reference panel for concordance imputation #
######################################################

#rule vcf_coverage:
#    """
#    Estimate depth of coverage of reference panel samples to filter out low coverage
#    """
#    input:
#        ref = config['reference_panel']
#    output:
#        ref_depth = '{path}/output/reference_panel_only_dogs/reference_panel_depth.idepth'
#    log:
#        '{path}/output/reference_panel_only_dogs/reference_panel_depth.log'
#    resources:
#        mem_mb=50*1024
#    benchmark:
#        '{path}/benchmarks/reference_panel_only_dogs/reference_panel_depth.tsv'
#    conda:
#        #'../envs/environment.yaml'
#    shell: 
#        '''
#        vcftools \
#        --gzvcf {input.ref} \
#        --depth \
#        --out {wildcards.path}/output/reference_panel_only_dogs/reference_panel_depth 2> {log}
#        '''

rule make_vcf_sample_list_only_dogs:
    """
    Make list of samples to remove from reference panel
    """
    input:
        ref_depth = '{path}/output/reference_panel/reference_panel_depth.idepth',
        reseq_names = 'sample_lists/reseq_samples.txt',
        boxer_names = 'sample_lists/boxer_samples.txt',
        non_dog_names = 'sample_lists/non_dog_samples.txt'
    output:
        sample_list_exclude = '{path}/output/reference_panel_only_dogs/remove_samples.txt'
    params:
        depth_cutoff = config["depth_cutoff"]
    shell:
        '''
        awk '{{print $1, $3}}' {input.ref_depth} | sed '1d' > tempb.txt
        awk -F" " '$2<{params.depth_cutoff}' tempb.txt | awk '{{ print $1 }}' > temp2b.txt
        cat temp2b.txt {input.reseq_names} {input.boxer_names} {input.non_dog_names} | awk '!seen[$0]++' > {output.sample_list_exclude}
        rm tempb.txt temp2b.txt
        '''

rule remove_sample_indels_multiallelic_snps_concordance_only_dogs:
    """
    Remove samples, only keep biallelic snps
    """
    input:
        ref = config['reference_panel'],
        sample_list_exclude = '{path}/output/reference_panel_only_dogs/remove_samples.txt'
    output:
        ref_sample_snp = temp('{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp.vcf.gz')
    log:
        '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -r {wildcards.chrom_con} \
        -m 2 -M 2 \
        -S ^{input.sample_list_exclude} \
        --trim-alt-alleles \
        {input.ref} -Ou | 
        bcftools filter -e "type!='snp'" -Oz -o {output.ref_sample_snp}
        '''

rule fill_tags_concordance_only_dogs:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING)
    """
    input:
        ref_sample_snp = '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp.vcf.gz'
    output:
        ref_sample_snp_filltags = temp('{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags.vcf.gz')
    log:
        '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools +fill-tags {input.ref_sample_snp} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags} \
        -- -t all,F_MISSING 
        '''

#rule get_F_missing_hist:
#    input:
#        ref_sample_snp_filltags = '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags.vcf.gz'
#    output:
#        ref_panel_f_missing_hist = temp('{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_f_missing.txt')
#    #conda:
#        '../envs/environment.yaml'
#    shell:
#        '''
#        bcftools query -f '%INFO/F_MISSING\n' {input.ref_sample_snp_filltags} > {output.ref_panel_f_missing_hist}
#        '''

#rule merge_F_missing:
#    input:
#        ref_panel_f_missing_hist = expand('{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_f_missing.txt', chrom=CHROM)
#    output:
#        ref_panel_f_missing_hist_merged = '{path}/output/reference_panel_only_dogs/ref-panel_allchrom_f_missing.txt'
#    shell:
#        '''
#        cat {input.ref_panel_f_missing_hist} > {output.ref_panel_f_missing_hist_merged}
#        '''

rule filter_sites_concordance_only_dogs:
    """
    Filter for only PASS sites and missingness (F_MISSING)
    """
    input:
        ref_sample_snp_filltags = '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags.vcf.gz'
    output:
        ref_sample_snp_filltags_filter = '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags_filter.vcf.gz',
        ref_sample_snp_filltags_filter_tbi = '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags_filter.vcf.gz.tbi'
    params:
        f_missing = config['F_MISSING']
    log:
        '{path}/output/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags_filter.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/reference_panel_only_dogs/ref-panel_{chrom_con}_sample-snp_filltags_filter.tsv'
    #conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view -i 'FILTER=="PASS" & F_MISSING<{params.f_missing}' {input.ref_sample_snp_filltags} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags_filter}

        bcftools index --tbi {output.ref_sample_snp_filltags_filter}
        '''