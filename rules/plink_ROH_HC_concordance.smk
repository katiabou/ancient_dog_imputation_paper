########################################################
#  Plink ROH to test target imputed against validation #
########################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']

rule transversions_phased:
    """
    Only take transversions from phased data
    """
    input:
        phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    output:
        tranversion_sites = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions.tsv.gz',
        phased_transversions = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.phased_vcf_maf_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.phased_vcf_maf_info} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.phased_transversions} 2> {log}

        bcftools index -f {output.phased_transversions}
        '''

rule make_plink_transversions_phased:
    """
    Prepare file format for plink
    """
    input:
        phased_transversions = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_plink.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_plink'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_plink.log'
    conda:
        '../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.phased_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_transversions_phased:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_plink.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink.hom'
    params:
        prefix_in = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_plink',
        prefix_out = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink.hom.log'
    conda:
        '../envs/plink19.yaml'
    shell:
        '''
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out}
        '''

rule merge_roh_transversions_phased_prep:
    """
    Prepare ROH output files for plot
    """
    input:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink.hom'
    output:
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink-temp.hom'
    shell:
        '''
        awk -F'\t' '{{print $0 "\t" "HC_imputed"}}' {input.phased_roh} > {output.temp_roh_phased}
        sed -i '1s/HC_imputed/cov/' {output.temp_roh_phased}
        '''

rule plot_imputed_concordance_roh_transversions: 
    """
    Plot ROH transversions for imputed concordance, imputed HC, and validation
    """
    input:
        temp_roh_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink-temp.hom', chrom=CHROM, maf=config['maf'], info=config['info'], path=config['path'], allow_missing=True),
        temp_roh_concordance_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_plink-temp.hom', chrom_con=CHROM_CON, coverage_val=COVERAGE_VAL, maf=config['maf'], info=config['info'], path=config['path'], allow_missing=True),
        temp_roh_validation = expand('{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_validation_filt_transversions_hom_win_het_{hom_win_het}_plink-temp.hom', chrom_con=CHROM_CON, path=config['path'], allow_missing=True),
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome'
    output:
        plot = '{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_transversions_hom_win_het_{hom_win_het}_ROH.png'
    params:
        files_concordance_phased=lambda wildcards, input: ','.join(input.temp_roh_concordance_phased),
        files_validation=lambda wildcards, input: input.temp_roh_validation,
        files_phased=lambda wildcards, input: input.temp_roh_phased,
        name = '{sample}',
        path_script = '{path}/scripts'
    shell:
        "Rscript {params.path_script}/ROH_plotting.R {params.files_concordance_phased} {params.files_validation} {params.files_phased} {params.name} {input.ref_fasta_chr_size} {output.plot}"
        
        
####################################################        
#### Run ROHs for transversions and transitions ####
####################################################      

rule make_plink_all_sites_phased:
    """
    Prepare file format for plink
    """
    input:
        phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_plink.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_plink'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_plink.log'
    conda:
        '../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.phased_vcf_maf_info} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_all_sites_phased:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink_hom_win_het_{hom_win_het}_plink.hom'
    params:
        prefix_in = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink',
        prefix_out = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink_hom_win_het_{hom_win_het}_plink'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink_hom_win_het_{hom_win_het}_plink.hom.log'
    conda:
        '../envs/plink19.yaml'
    shell:
        '''
        plink \
        --bfile {params.prefix_in} \
        --chr-set 38 \
        --homozyg \
        --homozyg-density {config[roh][homozyg_density]} \
        --homozyg-gap {config[roh][homozyg_gap]} \
        --homozyg-kb {config[roh][homozyg_kb]} \
        --homozyg-snp {config[roh][homozyg_snp]} \
        --homozyg-window-het {wildcards.hom_win_het} \
        --homozyg-window-missing {config[roh][homozyg_window_missing]} \
        --homozyg-window-snp {config[roh][homozyg_window_snp]} \
        --homozyg-window-threshold {config[roh][homozyg_window_threshold]} \
        --out {params.prefix_out}
        '''

rule merge_roh_all_sites_phased_prep:
    """
    Prepare ROH output files for plot
    """
    input:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink_hom_win_het_{hom_win_het}_plink.hom'
    output:
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink_hom_win_het_{hom_win_het}_plink-temp.hom'
    shell:
        '''
        awk -F'\t' '{{print $0 "\t" "HC_imputed"}}' {input.phased_roh} > {output.temp_roh_phased}
        sed -i '1s/HC_imputed/cov/' {output.temp_roh_phased}
        '''

rule plot_imputed_concordance_roh_all_sites:
    """
    Plot ROH all sites for imputed concordance, imputed HC, and validation
    """
    input:
        temp_roh_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_plink_hom_win_het_{hom_win_het}_plink-temp.hom', chrom_con=CHROM_CON, maf=config['maf'], info=config['info'], path=config['path'], allow_missing=True),
        temp_roh_concordance_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_hom_win_het_{hom_win_het}_plink-temp.hom', chrom_con=CHROM_CON, coverage_val=COVERAGE_VAL, maf=config['maf'], info=config['info'], path=config['path'], allow_missing=True),
        temp_roh_validation = expand('{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_validation_filt_all_sites_hom_win_het_{hom_win_het}_plink-temp.hom', chrom_con=CHROM_CON, path=config['path'], allow_missing=True),
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome'
    output:
        plot = '{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_all_sites_hom_win_het_{hom_win_het}_ROH.png'
    params:
        files_concordance_phased=lambda wildcards, input: ','.join(input.temp_roh_concordance_phased),
        files_validation=lambda wildcards, input: input.temp_roh_validation,
        files_phased=lambda wildcards, input: input.temp_roh_phased,
        name = '{sample}',
        path_script = '{path}/scripts'
    shell:
        "Rscript {params.path_script}/ROH_plotting.R {params.files_concordance_phased} {params.files_validation} {params.files_phased} {params.name} {input.ref_fasta_chr_size} {output.plot}"
        

