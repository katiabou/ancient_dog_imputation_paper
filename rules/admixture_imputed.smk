###################################################################################
# These rules run admixture on the full imputed dataset filtered for INFO and MAF #
# merged with the modern samples                                                  #
###################################################################################



### As input, use plink format, of the merged modern and imputed dataset (can take from pca step)
### (can exclude modern samples if needed at this step)


DOCS = ['bed','bim','fam']


rule extract_modern_samples_subset: #make this temp output file once finished!!!!!
    """
    Extract either only dogs or wolves from reference panel 
    """
    input:
        ref_sample_snp_filltags_filter = '{path}/output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        modern_subset = config['modern_subset_admixture'],
    output:
        modern_subset_vcf = '{path}/output/GLIMPSE_imputation/admixture_phased/modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_subset.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/admixture_phased/modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_subset.vcf.log'
    threads: 8
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -S {input.modern_subset} \
        --trim-alt-alleles \
        --threads {threads} \
        {input.ref_sample_snp_filltags_filter} -Oz -o {output.modern_subset_vcf} 2> {log}

        bcftools index --tbi {output.modern_subset_vcf}
        '''

rule merge_phased_modern_subset:
    """
    Merge phased and modern vcfs
    """
    input:
        modern_subset_vcf = '{path}/output/GLIMPSE_imputation/admixture_phased/modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_subset.vcf.gz',
        merged_phased_vcf_maf_recalibrated_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf}_recalibrated_INFO_{info}.vcf.gz'
    output:
        phased_modern_subset_vcf = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.{chrom}_MAF_{maf}_recalibrated_INFO_{info}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.{chrom}_MAF_{maf}_recalibrated_INFO_{info}.vcf.log'
    threads: 8
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools merge \
        {input.merged_phased_vcf_maf_recalibrated_info} {input.modern_subset_vcf} \
        --threads {threads} \
        -Oz -o {output.phased_modern_subset_vcf} 2> {log}
        '''

rule prepare_merged_chr_list_admixture:
    """ 
    Prepare list to merge chromosomes
    """
    input:
        phased_modern_subset_vcf = expand('{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.{chrom}_MAF_{maf}_recalibrated_INFO_{info}.vcf.gz', chrom=CHROM, allow_missing=True)
    output:
        chr_list = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/chr_list.MAF_{maf}_recalibrated_INFO_{info}.txt'
    shell:
        '''
        ls -v {input.phased_modern_subset_vcf} >> {output.chr_list}
        ''' 

rule merge_chrom_admixture:
    """ 
    Merge all filtered phased chromosomes
    """
    input: 
        chr_list = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/chr_list.MAF_{maf}_recalibrated_INFO_{info}.txt'
    output:
        merged_vcf = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.vcf.gz.log'
    threads: 8
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools concat \
        --file-list {input.chr_list} \
        -Oz -o {output.merged_vcf} \
        --threads {threads} 2> {log}

        tabix -p vcf {output.merged_vcf}
        '''


rule make_plink_phased_recalibrated_ref_subset:
    """
    Make plink format of recalibrated filtered phased samples (say I take the recalibrated imputed vcf based on all samples, will confirm though) and subset of ref panel
    """
    input:
        merged_vcf = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}'
    log:
        '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.log'
    conda:
        '../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.merged_vcf} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''


rule admixture:
    input:
        bed = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.bed'
    output:
        q_file = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.{admix_K}.Q',
        p_file = '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.{admix_K}.P'
    log:
        '{path}/output/GLIMPSE_imputation/admixture_phased/merged_phased_modern_vcf/merged_phased_ref_panel_subset.all_chr_MAF_{maf}_recalibrated_INFO_{info}.{admix_K}.Q.log'
    threads: 8
    conda:
        '../envs/admixture1.3.yaml'
    shell:
        '''
        adxmixture -j {threads} --cv \
        {input.bed} \
        {wildcards.admix_K} 2> {log}
        '''
