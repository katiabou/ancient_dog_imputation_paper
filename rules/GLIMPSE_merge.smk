##########################################################
#   Filter for info and maf, then merge phased bam files #
##########################################################

rule prepare_merge_chr_list:
    """ 
    Prepare list to merge all filtered imputed chromosomes
    """
    input:
        phased_vcf_info_maf = expand('{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.{chrom}_INFO_{info}_MAF_{maf}.vcf.gz', chrom=CHROM, allow_missing=True)
    output:
        chr_list = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/chr_list_INFO_{info}_MAF_{maf}.txt'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/chr_list_INFO_{info}_MAF_{maf}.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        ls -v {input.phased_vcf_info_maf} >> {output.chr_list}
        '''    

rule merge_phased_bcfs_phased:
    """ 
    Merge all filtered imputed chromosomes
    """
    input: 
        chr_list = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/chr_list_INFO_{info}_MAF_{maf}.txt'
    output:
        merged_phased_vcf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.all_chr_INFO_{info}_MAF_{maf}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.all_chr_INFO_{info}_MAF_{maf}.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools concat \
        --file-list {input.chr_list} \
        -Oz -o {output.merged_phased_vcf} \
        --threads {threads}

        tabix -p vcf {output.merged_phased_vcf}
        '''

rule prepare_merge_chr_list_no_info_maf:
    """ 
    Prepare list to merge all filtered imputed chromosomes
    """
    input:
        phased_vcf_info_maf = expand('{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.{chrom}.bcf', chrom=CHROM, allow_missing=True)
    output:
        chr_list = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/chr_list.txt'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/chr_list.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        ls -v {input.phased_vcf_info_maf} >> {output.chr_list}
        '''    

rule merge_phased_bcfs_phased_no_info_maf:
    """ 
    Merge all filtered imputed chromosomes
    """
    input: 
        chr_list = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/chr_list.txt'
    output:
        merged_phased_bcf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.all_chr.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.all_chr.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools concat \
        --file-list {input.chr_list} \
        -Oz -o {output.merged_phased_bcf} \
        --threads {threads}

        tabix -p vcf {output.merged_phased_bcf}
        '''
