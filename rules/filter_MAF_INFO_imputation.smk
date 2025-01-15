##########################################################################
#  Filter imputed dataset based on MAF of reference panel and INFO score #
##########################################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']

rule merge_phased_bcfs:
    """ 
    Merge all phased samples 
    """
    input: 
        phased_bcf = expand('{path}/output/GLIMPSE_imputation/GLIMPSE_phased/{bam_imputation}_phased.{chrom}.bcf', bam_imputation = BAM, allow_missing = True)
    output:
        merged_phased_vcf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}.vcf.gz'
    threads: 10
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        bcftools merge \
        {input.phased_bcf} \
        -Oz -o {output.merged_phased_vcf} \
        --threads {threads}

        tabix -p vcf {output.merged_phased_vcf}
        '''

rule MAF_sites_ref_pan_phased:
    """
    Extract sites using a MAF filter from the reference panel 
    """
    input:
        ref_sample_snp_filltags_filter = '{path}/output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
    output:
        ref_sample_snp_filltags_filter_maf_vcf = '{path}/output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.vcf.gz',
        ref_sample_snp_filltags_filter_maf_tsv = '{path}/output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.tsv.gz'
    params:
        maf=config['maf_cutoff']
    log:
        '{path}/output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.vcf.log'
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -q {params.maf}:minor \
        {input.ref_sample_snp_filltags_filter} \
        -Oz -o {output.ref_sample_snp_filltags_filter_maf_vcf} 2> {log}
        
        bcftools index -f {output.ref_sample_snp_filltags_filter_maf_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_sample_snp_filltags_filter_maf_vcf} | \
        bgzip -c > {output.ref_sample_snp_filltags_filter_maf_tsv}

        tabix -s1 -b2 -e2 {output.ref_sample_snp_filltags_filter_maf_tsv}
        '''

rule maf_sites_phased_vcf:
    """
    Extract MAF sites from merged phased VCF
    """
    input:
        ref_sample_snp_filltags_filter_maf_tsv = '{path}/output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.tsv.gz', 
        merged_phased_vcf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}.vcf.gz'
    output:
        merged_phased_vcf_maf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}.vcf.gz.log'
    #conda:
    #    '../envs/environment.yaml'
    threads: 10
    shell:
        '''
        bcftools view {input.merged_phased_vcf} \
        --regions-file {input.ref_sample_snp_filltags_filter_maf_tsv} \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf}
        '''


rule filter_INFO_phased_vcf:
    """
    Filter for INFO based on concordance results WITHOUT re-calibrating INFO score for all samples together
    """
    input:
        merged_phased_vcf_maf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}.vcf.gz'
    output:
        merged_phased_vcf_maf_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}.vcf.gz'
    params:
        info=config['info']
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}.vcf.gz.log'
    threads: 10
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        {input.merged_phased_vcf_maf} \
        --trim-alt-alleles -Ou | \
        bcftools view \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf_info}

        tabix -p vcf {output.merged_phased_vcf_maf_info}
        '''


rule recalibrate_info_phased_vcf:
    """
    Re-calibrate INFO scores based on all samples present in the merged VCF (this is for analyses where all samples are used together, like when they're used to create PCs)
    """
    input:
        merged_phased_vcf_maf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}.vcf.gz'
    output:
        merged_phased_vcf_maf_recalibrated = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz.log'
    #conda:
    #    '../envs/environment.yaml'
    threads: 10
    shell:
        '''
        bcftools plugin impute-info \
        {input.merged_phased_vcf_maf} \
        -Ob -o {output.merged_phased_vcf_maf_recalibrated} \
        --threads {threads} 2> {log}

        bcftools index -f {output.merged_phased_vcf_maf_recalibrated}
        '''

rule filter_recalibrated_INFO_phased_vcf:
    """
    Filter for INFO after re-calibrating INFO score for all samples together
    """
    input:
        merged_phased_vcf_maf_recalibrated = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO.vcf.gz'
    output:
        merged_phased_vcf_maf_recalibrated_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz'
    params:
        info=config['info']
    log:
        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz.log'
    threads: 10
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        {input.merged_phased_vcf_maf_recalibrated} \
        --trim-alt-alleles -Ou | \
        bcftools view \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.merged_phased_vcf_maf_recalibrated_info}

        tabix -p vcf {output.merged_phased_vcf_maf_recalibrated_info}
        '''