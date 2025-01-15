##########################################
# Prepare reference panel for imputation #
##########################################


#### NOT FILTERING CORRECLTY, LOOK AT REF_PANEL_CONCORDANCE!!!!

rule remove_sample_indels_multiallelic_snps:
    """
    Remove samples, only keep biallelic snps
    """
    input:
        ref = config['reference_panel'],
        sample_list_exclude = '{path}/output/GLIMPSE_concordance/reference_panel/remove_samples.txt'
    output:
        ref_sample_snp = temp('{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp.vcf.gz')
    log:
        '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp.tsv'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view -e "type!='snp'" --max-alleles 2 \
        -r {wildcards.chrom} \
        -S ^{input.sample_list_exclude} \
        --trim-alt-alleles \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp} {input.ref}
        '''

rule fill_tags:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING)
    """
    input:
        ref_sample_snp = '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp.vcf.gz'
    output:
        ref_sample_snp_filltags = temp('{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags.vcf.gz')
    log:
        '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags.tsv'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools +fill-tags {input.ref_sample_snp} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags} \
        -- -t all,F_MISSING 
        '''
        
rule filter_sites:
    """
    Filter for only PASS sites and missingness (F_MISSING)
    """
    input:
        ref_sample_snp_filltags = '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags.vcf.gz'
    output:
        ref_sample_snp_filltags_filter = '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        ref_sample_snp_filltags_filter_tbi = '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz.tbi'
    params:
        f_missing = config['F_MISSING']
    log:
        '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.tsv'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view -i 'FILTER=="PASS" & F_MISSING<{params.f_missing}' {input.ref_sample_snp_filltags} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags_filter}

        bcftools index --tbi {output.ref_sample_snp_filltags_filter}
        '''