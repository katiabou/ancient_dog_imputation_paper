###############################################################################################################
#  Filter imputed downsampled and HC datasets based on MAF of reference panel, INFO score and low cov samples #
###############################################################################################################

# define output files of make plink
DOCS = ["bed", "bim", "fam"]

### Downsampled imputed ###


rule MAF_sites_ref_pan:
    """
    Extract sites using a MAF filter from the reference panel 
    """
    input:
        ref_concordance_sample_excl_filltags_filter="{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf",
    output:
        ref_concordance_sample_excl_filltags_filter_maf_vcf="{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv="{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz",
    params:
        maf=config["maf_cutoff"],
    log:
        "{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.bcf.log",
    shell:
        """
        bcftools view \
        -q {params.maf}:minor \
        {input.ref_concordance_sample_excl_filltags_filter} \
        -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf} 2> {log}
        
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv}

        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv}
        """


# rule maf_INFO_sites_concordance_imputed:
#     """
#     Extract MAF and INFO filtered sites from imputed VCF
#     """
#     input:
#         ref_concordance_sample_excl_filltags_filter_maf_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz',
#         ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
#     output:
#         imputed_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz'
#     params:
#         info=config['info']
#     log:
#         '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz.log'
#     threads: 8
#     shell:
#         '''
#         bcftools view {input.ligated_bcf} \
#         --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv} \
#         --include 'INFO/INFO >= {params.info}' \
#         --threads {threads} \
#         -Oz -o {output.imputed_maf_info} 2> {log}

#         bcftools index -f {output.imputed_maf_info}
#         '''


rule maf_INFO_sites_concordance_phased:
    """
    Extract MAF and INFO filtered sites from phased VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv="{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz",
        phased_bcf="{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.bcf",
    output:
        phased_maf_info="{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz",
    params:
        info=config["info_cutoff"],
    log:
        "{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz.log",
    threads: 8
    shell:
        """
        bcftools view {input.phased_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.phased_maf_info} 2> {log}

        bcftools index -f {output.phased_maf_info}
        """


# rule make_lov_cov_sample_file_concordance:
#    """
#    Remove low coverage samples from imputed dataset (<0.5x)
#    Have to have a "Mean_Depth" column in metadata, along with the "Sample" column which as the bam name
#    """
#    input:
#        bams_meta ="sample_lists/bams_published_imputation_metadata_cutoff.tsv",
#    output:
#        remove_samples_imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/remove_samples_imputed_no_low_cov_{cov_cutoff}.txt'
#    params:
#        cov_imp_cutoff = config['cov_cutoff']
#    shell:
#        '''
#        f=$(head -1 {input.bams_meta} | tr '\t' '\n' | cat -n | grep "Mean_Depth" | awk '{{print $1}}')
#        j=$(head -1 {input.bams_meta} | tr '\t' '\n' | cat -n | grep "name_haplo_VCF" | awk '{{print $1}}')

#        awk -v col="$f" -F"\t" '$col<{params.cov_imp_cutoff}' {input.bams_meta} | cut -f $j > {output.remove_samples_imputed}
#        '''

# rule filter_INFO_concordance:
#    """
#    Filter for INFO based on concordance results (already filtered for MAF based on reference panel) NOT REMOVING LOW COVERAGE SINCE THERE IS ONLY ONE SAMPLE IN EACH IMPUTED VCF
#    """
#    input:
#        phased_vcf_maf = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_MAF_{maf}.vcf.gz',
# remove_samples_imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/remove_samples_imputed_no_low_cov_{cov_cutoff}.txt'
#    output:
#        phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz'
#    params:
#        info=config['info']
#    log:
#        '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz.log'
#    threads: 10
#    #conda:
#        #'../envs/environment.yaml'
#    shell:
#        '''
#        bcftools view \
#        {input.phased_vcf_maf} \
#        --trim-alt-alleles -Ou | \
#        bcftools view \
#        --include 'INFO/INFO >= {params.info}' \
#        --threads {threads} \
#        -Oz -o {output.phased_vcf_maf_info}

#        tabix -p vcf {output.phased_vcf_maf_info}
#        '''

### HC imputed ###

# rule maf_INFO_sites_concordance_HC_imputed:
#     """
#     Extract MAF sites from HC imputed VCF
#     """
#     input:
#         ref_concordance_sample_excl_filltags_filter_maf_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz',
#         ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}.bcf',
#     output:
#         imputed_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz'
#     params:
#         info=config['info']
#     log:
#         '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz.log'
#     threads: 8
#     shell:
#         '''
#         bcftools view {input.ligated_bcf} \
#         --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv} \
#         --include 'INFO/INFO >= {params.info}' \
#         --threads {threads} \
#         -Oz -o {output.imputed_maf_info} 2> {log}

#         bcftools index -f {output.imputed_maf_info}
#         '''


rule maf_INFO_sites_concordance_HC_phased:
    """
    Extract MAF sites from HC phased VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv="{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf}.phased.tsv.gz",
        phased_bcf="{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.bcf",
    output:
        phased_maf_info="{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_INFO_{info}_MAF_{maf}.vcf.gz",
    params:
        info=config["info_cutoff"],
    log:
        "{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_INFO_{info}_MAF_{maf}.vcf.gz.log",
    threads: 8
    shell:
        """
        bcftools view {input.phased_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.phased_maf_info} 2> {log}

        bcftools index -f {output.phased_maf_info}
        """


# rule filter_INFO_low_cov_HC_concordance:
#    """
#    Filter for INFO based on HC concordance results (already filtered for MAF based on reference panel) NOT REMOVING LOW COVERAGE SINCE THERE IS ONLY ONE SAMPLE IN EACH IMPUTED VCF
#    """
#    input:
#        phased_vcf_maf = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_MAF_{maf}.vcf.gz',
# remove_samples_imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/remove_samples_imputed_no_low_cov_{cov_cutoff}.txt'
#    output:
#        phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz'
#    params:
#        info=config['info']
#    log:
#        '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz.log'
#    #conda:
#        #'../envs/environment.yaml'
#    shell:
#        '''
#        bcftools view \
#        {input.phased_vcf_maf} \
#        --trim-alt-alleles -Ou | \
#        bcftools view \
#        --include 'INFO/INFO >= {params.info}' \
#        --threads {threads} \
#        -Oz -o {output.phased_vcf_maf_info}
#        tabix -p vcf {output.phased_vcf_maf_info}
#        '''
