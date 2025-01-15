#######################################
#  Plink ROH for full imputed dataset #
#######################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']


rule transversions_phased_subset:
    """
    Only take transversions from phased data (non-recalibrated, since I'm focusing on each sample seperately)
    Doing this for dogs and wolves seperately
    """
    input:
        merged_phased_vcf_maf_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info}.vcf.gz',
        imputed_canid_subset = config['imputed_canid_subset']
    output:
        tranversion_sites = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.tsv.gz',
        phased_transversions = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.log'
    threads: 8
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.merged_phased_vcf_maf_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.merged_phased_vcf_maf_info} \
        -S {input.imputed_canid_subset} \
        --trim-alt-alleles \
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
        phased_transversions = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.log'
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
        bim = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_roh_sum = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', 
    params:
        prefix_in = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_{canid_subset}',
        prefix_out = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
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

rule merge_chrom_ROH_phased_transversions:
    input:
        phased_roh = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
        phased_roh_sum = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
    output:
        phased_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_roh_sum_allchrom = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
    shell:
        '''
        awk 'FNR>1 || NR==1' {input.phased_roh} > {output.phased_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.phased_roh_sum} > {output.phased_roh_sum_allchrom}
        '''

        
####################################################        
#### Run ROHs for transversions and transitions ####
####################################################      

rule all_sites_phased_subset:
    """
    Take all sites from phased data (non-recalibrated, since I'm focusing on each sample seperately)
    Doing this for dogs and wolves seperately
    """
    input:
        merged_phased_vcf_maf_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_INFO_{info}.vcf.gz',
        imputed_canid_subset = config['imputed_canid_subset']
    output:
        phased_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.log'
    threads: 8
    shell:
        '''
        bcftools view {input.merged_phased_vcf_maf_info} \
        -S {input.imputed_canid_subset} \
        --trim-alt-alleles \
        --threads {threads} \
        -Oz -o {output.phased_all_sites_subset} 2> {log}

        bcftools index -f {output.phased_all_sites_subset}
        '''


rule make_plink_all_sites_phased:
    """
    Prepare file format for plink
    """
    input:
        phased_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.log'
    shell:
        '''
        plink \
        --vcf {input.phased_all_sites_subset} \
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
        bim = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_roh_ind = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        phased_roh_sum = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary'
    params:
        prefix_in = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}',
        prefix_out = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
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

rule merge_chrom_ROH_phased:
    input:
        phased_roh = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
        phased_roh_sum = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
    output:
        phased_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_roh_sum_allchrom = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
    shell:
        '''
        awk 'FNR>1 || NR==1' {input.phased_roh} > {output.phased_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.phased_roh_sum} > {output.phased_roh_sum_allchrom}
        '''


####################################################       
#### Run ROHs for transversions reference panel ####
####################################################    
        
rule transversions_ref_panel_subset:
    """
    Only take transversions from ref panel
    """
    input:
        #ref_sample_snp_filltags_filter = '{path}/output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        ref_sample_snp_filltags_filter_maf_vcf = '{path}/output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.phased.vcf.gz',
        modern_canid_subset = config['modern_canid_subset'],
    output:
        tranversion_sites = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.tsv.gz',
        ref_panel_transversions = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.log'
    threads: 8
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.ref_sample_snp_filltags_filter_maf_vcf} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.ref_sample_snp_filltags_filter_maf_vcf} \
        -S {input.modern_canid_subset} \
        --trim-alt-alleles \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.ref_panel_transversions} 2> {log}

        bcftools index -f {output.ref_panel_transversions}
        '''

rule make_plink_transversions_ref_panel:
    """
    Prepare file format for plink
    """
    input:
        ref_panel_transversions = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.log'
    shell:
        '''
        plink \
        --vcf {input.ref_panel_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_transversions_ref_panel:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_roh_sum = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        phased_roh_ind = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv'
    params:
        prefix_in = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_{canid_subset}',
        prefix_out = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
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

rule merge_chrom_ROH_transversions_ref_panel:
    input:
        modern_roh = expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
    output:
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_{canid_subset}.hom',
    shell:
        '''
        awk 'FNR>1 || NR==1' {input.modern_roh} > {output.modern_roh_allchrom}
        '''


################################################       
#### Run ROHs for all sites reference panel ####
################################################    

rule all_sites_ref_panel_subset:
    """
    Take all sites from ref panel
    """
    input:
        #ref_sample_snp_filltags_filter = '{path}/output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        ref_sample_snp_filltags_filter_maf_vcf = '{path}/output/GLIMPSE_imputation/reference_panel/{chrom}_ref_panel_filltags_filter_MAF_{maf_cutoff}.phased.vcf.gz',
        modern_canid_subset = config['modern_canid_subset'],
    output:
        ref_panel_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.log'
    threads: 8
    shell:
        '''
        bcftools view {input.ref_sample_snp_filltags_filter_maf_vcf} \
        -S {input.modern_canid_subset} \
        --trim-alt-alleles \
        --threads {threads} \
        -Oz -o {output.ref_panel_all_sites_subset} 2> {log}

        bcftools index -f {output.ref_panel_all_sites_subset}
        '''


rule make_plink_all_sites_ref_panel:
    """
    Prepare file format for plink
    """
    input:
        ref_panel_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.log'
    shell:
        '''
        plink \
        --vcf {input.ref_panel_all_sites_subset} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_all_sites_ref_panel:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_roh_sum = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        phased_roh_ind= '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv'
    params:
        prefix_in = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}',
        prefix_out = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
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

rule merge_chrom_ROH_ref_panel:
    input:
        modern_roh = expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
        modern_roh_sum = expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
    output:
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        modern_roh_allchrom_sum = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
    shell:
        '''
        awk 'FNR>1 || NR==1' {input.modern_roh} > {output.modern_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.modern_roh_sum} > {output.modern_roh_allchrom_sum}
        '''


###########################################################################################################################################
#### Estimating ROHS for merged modern and ancient (since we want the .hom.summary file which estimated based on the input file samples) 
#### This is for the heatmap for all dogs and all wolves (seperately)
#### Doing this only for all sites
###########################################################################################################################################
    
rule merge_phased_modern_all_sites:
    """
    Prepare file format for plink
    """
    input:
        phased_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz',
        #ref_panel_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_all_sites_{canid_subset}.vcf.gz'
        ref_panel_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_{chrom}_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_{canid_subset}.vcf.gz'
    output:
        phased_modern_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz',
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz.log'
    threads: 8
    shell:
        '''
        bcftools merge \
        {input.phased_all_sites_subset} {input.ref_panel_all_sites_subset} \
        --threads {threads} \
        -Oz -o {output.phased_modern_all_sites_subset} 2> {log}

        bcftools index -f {output.phased_modern_all_sites_subset}
        '''

rule make_plink_all_sites_phased_ref_panel:
    """
    Prepare file format for plink
    """
    input:
        phased_modern_all_sites_subset = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.vcf.gz',
    output:
        expand('{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.plink.log'
    shell:
        '''
        plink \
        --vcf {input.phased_modern_all_sites_subset} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_all_sites_phased_ref_panel:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}.bim'
    output:
        phased_modern_roh = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_modern_roh_ind = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        phased_modern_roh_summary = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary'
    params:
        prefix_in = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_{canid_subset}',
        prefix_out = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.log'
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

rule merge_chrom_ROH_phased_ref_panel:
    input:
        phased_modern_roh = expand('{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', chrom=CHROM, allow_missing=True),
        phased_modern_roh_sum = expand('{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.{chrom}_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary', chrom=CHROM, allow_missing=True),
    output:
        phased_modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom',
        phased_modern_roh_sum_allchrom_sum = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary'
    shell:
        '''
        awk 'FNR>1 || NR==1' {input.phased_modern_roh} > {output.phased_modern_roh_allchrom}

        awk 'FNR>1 || NR==1' {input.phased_modern_roh_sum} > {output.phased_modern_roh_sum_allchrom_sum}
        '''


##################################################################################     
#### Steps for plotting all ROH for imputed and selected modern for all sites ####
##################################################################################     

# rule estimate_chrom_size:  # this is wrong cause it has unassigned contigs 
#     """
#     Estimate chromosome sizes for plotting x axis 
#     """
#     input:
#         ref_fasta = config['ref_fasta_file']
#     output:
#         ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
#     shell:
#         '''
#         faidx {input.ref_fasta} -i chromsizes > {output.ref_fasta_size}
#         '''

rule estimate_chrom_size_per_chr:
    """
    Estimate chromosome sizes for plotting x axis 
    """
    input:
        ref_fasta_chr = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}.fasta',
    output:
        ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}_size.genome'
    shell:
        '''
        faidx {input.ref_fasta_chr} -i chromsizes > {output.ref_fasta_size}
        '''

rule estimate_allchrom_size:
    """
    Estimate all autosome size 
    """
    input:
        ref_fasta_size = expand('{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}_size.genome', chrom=CHROM, allow_missing=True)
    output:
        ref_fasta_allchr_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome'
    shell:
        '''
        cat {input.ref_fasta_size} | awk '{{Total=Total+$2}} END{{print "allchrom " Total}}' > {output.ref_fasta_allchr_size}
        '''

# rule merge_all_imputed_ROH:
#     input:
#         phased_roh_allchrom_transverions = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom', canid_subset=CANID_SUBSET, allow_missing=True),
#         phased_roh_allchrom_transverions_ind = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv', canid_subset=CANID_SUBSET, allow_missing=True),
#         phased_roh_allchrom_all_sites = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', canid_subset=CANID_SUBSET, allow_missing=True),
#         phased_roh_allchrom_all_sites_ind = expand('{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv', canid_subset=CANID_SUBSET, allow_missing=True)
#     output:
#         phased_roh_allchrom_tranversions_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}.hom',
#         phased_roh_allchrom_tranversions_merged_ind = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}.hom.indiv',
#         phased_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}.hom',
#         phased_roh_allchrom_all_sites_merged_ind = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}.hom.indiv',
#     shell:
#         '''
#         awk 'NR == 1 || FNR > 1'  {input.phased_roh_allchrom_transverions} > {output.phased_roh_allchrom_tranversions_merged}
#         awk 'NR == 1 || FNR > 1'  {input.phased_roh_allchrom_transverions_ind} > {output.phased_roh_allchrom_tranversions_merged_ind}

#         awk 'NR == 1 || FNR > 1'  {input.phased_roh_allchrom_all_sites} > {output.phased_roh_allchrom_all_sites_merged}
#         awk 'NR == 1 || FNR > 1'  {input.phased_roh_allchrom_all_sites_ind} > {output.phased_roh_allchrom_all_sites_merged_ind}
#         '''

# rule merge_all_modern_ROH:
#     input:
#         ref_panel_roh_allchrom_transversions = expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_{canid_subset}.hom', canid_subset=CANID_SUBSET, allow_missing=True),
#         ref_panel_roh_allchrom_all_sites = expand('{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom', canid_subset=CANID_SUBSET, allow_missing=True),
#     output:
#         ref_panel_roh_allchrom_transversions_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}.hom',
#         ref_panel_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}.hom',
#     shell:
#         '''
#         awk 'NR == 1 || FNR > 1'  {input.ref_panel_roh_allchrom_transversions} > {output.ref_panel_roh_allchrom_transversions_merged}
#         awk 'NR == 1 || FNR > 1'  {input.ref_panel_roh_allchrom_all_sites} > {output.ref_panel_roh_allchrom_all_sites_merged}
#         '''


rule plot_transversions_ROH_all_chr:
    """
    Using only the dogwolf files for these plots (contains all imputed samples)
    """
    input:
        bam_metadata = config['bam_imputation_meta'],
        ref_metadata = config['reference_panel_metadata'],
        phased_roh_allchrom_tranversions_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_panel_roh_allchrom_transversions_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_dogwolf.hom',
    output:
        plot_dogs = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_only.png',
        plot_wolves = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_only.png'
    params:
        sites = 'transversions'
    script:
        "../scripts/ROH_bands_all_chr.R"

rule plot_all_sites_ROH_all_chr:
    """
    Using only the dogwolf files for these plots (contains all imputed samples)
    """
    input:
        bam_metadata = config['bam_imputation_meta'],
        ref_metadata = config['reference_panel_metadata'],
        phased_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_panel_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
    output:
        plot_dogs = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_only.png',
        plot_wolves = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_only.png'
    params:
        sites = 'all sites'
    script:
        "../scripts/ROH_bands_all_chr.R"


rule plot_all_sites_ROH_count_length_dogs:
    input:
        bam_metadata = config['bam_imputation_meta'],
        ref_metadata = config['reference_panel_metadata'],
        phased_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_panel_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome'
    output:
        plot_dogs_length_count = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs.png',
        plot_dogs_length_count_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs-labelled.png',
        plot_dogs_coeff = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs.png',
        plot_dogs_coeff_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs-labelled.png',
        plot_dogs_length_count_long = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs.png',
        plot_dogs_length_count_long_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs-labelled.png',
        plot_dogs_coeff_long = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs.png',
        plot_dogs_coeff_long_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs-labelled.png',
        plot_dogs_length_count_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs.png',
        plot_dogs_length_count_short_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs-labelled.png',
        plot_dogs_coeff_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs.png',
        plot_dogs_coeff_short_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs-labelled.png',
        froh_test = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_tests.tsv',
        roh_results_all = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_all_roh_results.tsv',
        plot_dogs_coeff_map_main = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs_map.png',
        plot_dogs_coeff_long_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs.png',
        plot_dogs_coeff_boxplot = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs_boxplot.png',
    script:
        "../scripts/ROH_count_length_coeff_dogs.R"


rule plot_transversions_ROH_count_length_dogs:
    input:
        bam_metadata = config['bam_imputation_meta'],
        ref_metadata = config['reference_panel_metadata'],
        phased_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_panel_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome'
    output:
        plot_dogs_length_count = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs.png',
        plot_dogs_length_count_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_all_ROHs-labelled.png',
        plot_dogs_coeff = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs.png',
        plot_dogs_coeff_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs-labelled.png',
        plot_dogs_length_count_long = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs.png',
        plot_dogs_length_count_long_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_long_ROHs-labelled.png',
        plot_dogs_coeff_long = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs.png',
        plot_dogs_coeff_long_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_ROHs-labelled.png',
        plot_dogs_length_count_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs.png',
        plot_dogs_length_count_short_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_count_length_short_ROHs-labelled.png',
        plot_dogs_coeff_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs.png',
        plot_dogs_coeff_short_laebelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_short_ROHs-labelled.png',
        froh_test_dogs = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_tests.tsv',
        roh_results_all = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_all_roh_results.tsv',
        plot_dogs_coeff_map_main = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_all_ROHs_map.png',
        plot_dogs_coeff_long_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs.png',
        plot_dogs_coeff_boxplot = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogs_coeff_long_short_ROHs_boxplot.png',
    script:
        "../scripts/ROH_count_length_coeff_dogs.R"



rule plot_all_sites_ROH_count_length_wolves:
    input:
        bam_metadata = config['bam_imputation_meta'],
        ref_metadata = config['reference_panel_metadata'],
        phased_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_panel_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome'
    output:
        plot_wolves_length_count = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_all_ROHs.png',
        plot_wolves_length_count_long = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_long_ROHs.png',
        plot_wolves_length_count_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_count_length_short_ROHs.png',
        plot_wolves_coeff = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs.png',
        plot_wolves_coeff_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs-labelled.png',
        #plot_wolves_coeff_boxplot_modern = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_boxplot.png',
        plot_wolves_coeff_boxplot_modern_ancient = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_ancient_boxplot.png',
        plot_pleistocene_wolves_coeff = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene.png',
        plot_pleistocene_wolves_coeff_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene-labelled.png',
        roh_results_all = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_all_roh_results.tsv',
    script:
        "../scripts/ROH_count_length_coeff_wolves.R"


rule plot_transversions_ROH_count_length_wolves:
    input:
        bam_metadata = config['bam_imputation_meta'],
        ref_metadata = config['reference_panel_metadata'],
        phased_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_panel_roh_allchrom_all_sites_merged = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_hom_win_het_{hom_win_het}_dogwolf.hom',
        modern_roh_allchrom = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_transversions_MAF_{maf_cutoff}_hom_win_het_{hom_win_het}_dogwolf.hom',
        #ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchr_size.genome'
        ref_fasta_allchr_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_allchrom_size.genome'
    output:
        plot_wolves_length_count = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_count_length_all_ROHs.png',
        plot_wolves_length_count_long = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_count_length_long_ROHs.png',
        plot_wolves_length_count_short = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_count_length_short_ROHs.png',
        plot_wolves_coeff = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs.png',
        plot_wolves_coeff_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs-labelled.png',
        #plot_wolves_coeff_boxplot_modern = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_boxplot.png',
        plot_wolves_coeff_boxplot_modern_ancient = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_modern_ancient_boxplot.png',
        plot_pleistocene_wolves_coeff = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene.png',
        plot_pleistocene_wolves_coeff_labelled = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_coeff_all_long_short_ROHs_pleistocene-labelled.png',
        roh_results_all = '{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_{hom_win_het}_wolves_all_roh_results.tsv',
    script:
        "../scripts/ROH_count_length_coeff_wolves.R"


##################################################################################     
#### Steps for plotting all ROH for imputed and selected modern for all sites ####
##################################################################################     


rule estimate_ROH_windows:
    """
    Get file with 500KB windows, to estimate coverage after. These will be the same, since the genome size is the same for all samples. 
    I'll further estimate coverages using either only the ancient wolf bams or only the ancient dog bams, or both dogs and wolves (seperate window coverage masks)
    """
    input:
        ref_fasta_size = '{path}/output/GLIMPSE_imputation/reference_genome/CanFam31_{chrom}_size.genome'
    output:
        ROH_500_kb_windows = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/CanFam31_{chrom}_500_kb_windows.txt',
    script:
        "../scripts/ROH_500_kb_window.R"

rule bamlist_subset:
    """
    Create bam list for either dogs or wolves or dogs/wolves
    """
    input:
        bam_meta =  config['bam_imputation'],
        imputed_canid_subset = config['imputed_canid_subset'],
    output:
        bam_list_subset = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/bam_list_{canid_subset}_window_cov.txt'
    shell:
        '''
        awk -F '\t' 'NR==FNR{{a[$1]; next}} FNR==1 || $1 in a' {input.imputed_canid_subset} {input.bam_meta} | awk -F'\t' 'FNR>1{{print $3}}' > {output.bam_list_subset}
        '''

rule window_coverage_estimate:
    """
    Esimate window coverage for imputed bams
    """
    input:
        bam_list_subset = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/bam_list_{canid_subset}_window_cov.txt',
        ROH_500_kb_windows = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/CanFam31_{chrom}_500_kb_windows.txt',
    output:
        windows_cov_subset = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_{chrom}_windows_cov_500kb.txt'
    shell:
        '''
        while read -r bam_line
        do
        echo $bam_line > temp_{wildcards.chrom}_{wildcards.canid_subset}.txt 
        while read -r line
        do
        samtools coverage -r $line --bam-list temp_{wildcards.chrom}_{wildcards.canid_subset}.txt    | grep -v '#' >> {output.windows_cov_subset}
        done < {input.ROH_500_kb_windows}
        done < {input.bam_list_subset}

        rm temp_{wildcards.chrom}_{wildcards.canid_subset}.txt 
        '''

rule merge_chrom_window_coverage_estimate:
    """
    Merge all chromosomes 
    """
    input:
        windows_cov_subset = expand('{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_{chrom}_windows_cov_500kb.txt', chrom=CHROM, allow_missing=True),
    output:
        windows_cov_subset_allchrom = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_allchrom_windows_cov_500kb.txt'
    shell:
        '''
        cat {input.windows_cov_subset} >> {output.windows_cov_subset_allchrom}
        '''

rule plot_ROH_window_prevelance_depth:
    """
    Plotting all dogs, all wolves, and dogs/wolves
    Getting overlapping windows
    """
    input:
        phased_roh_sum_allchrom = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        phased_roh_sum_ind = '{path}/output/GLIMPSE_imputation/ROH_phased/merged_phased.chr1_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        modern_roh_allchrom_sum = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        modern_roh_ind = '{path}/output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_chr1_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        phased_modern_roh_sum_allchrom_sum = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.summary',
        phased_modern_roh_sum_ind = '{path}/output/GLIMPSE_imputation/ROH_phased_modern/merged_phased_modern.chr1_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.hom.indiv',
        windows_cov_subset_allchrom = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/{canid_subset}_allchrom_windows_cov_500kb.txt'
    output:
        windows_prev_depth_plot = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/windows_prevelance_cov_500kb_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        heatmap_imputed_modern = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_modern_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        heatmap_density_imputed_modern = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_modern_heatmap_density_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        heatmap_imputed = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        heatmap_density_imputed = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/imputed_heatmap_density_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        heatmap_modern = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        heatmap_density_modern = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_heatmap_density_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.png',
        ROH_windows_modern_imputed = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_imputed_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.txt',
        ROH_windows_imputed = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.txt',
        ROH_windows_modern = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}.txt',
        imputed_windows_bed_islands = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed',
        imputed_windows_bed_deserts = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed',
        modern_windows_bed_islands = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed',
        modern_windows_bed_deserts = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed',
        imputed_modern_windows_bed_islands = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed',
        imputed_modern_windows_bed_deserts = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed',
        top_go_terms = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts_GO_terms.txt',
        main_figure = '{path}/output/GLIMPSE_imputation/plots/ROH_islands_deserts/modern_ancient_heatmap_ROH_500kb_windows_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_main.png',
    script:
        "../scripts/ROH_deserts_islands.R"



rule get_canfam31_gene_annotations:
    """
    Get gene annotations for ROH islands and deserts in common between modern and ancient samples 
    """
    input:
        cf31_ann = config['canfam31_annotation'],
        imputed_windows_bed_islands = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed',
        imputed_windows_bed_deserts = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed',
        modern_windows_bed_islands = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed',
        modern_windows_bed_deserts = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed',
        imputed_modern_windows_bed_islands = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_islands.bed',
        imputed_modern_windows_bed_deserts = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}_deserts.bed',
    output:
        genes_temp = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/cf31_annotation_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-temp.bed',
        imputed_genes_desert = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes_deserts.bed',
        imputed_genes_island = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes_islands.bed',
        imputed_genes_unique_desert = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique_deserts.bed',
        imputed_genes_unique_island = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique_islands.bed',
        modern_genes_desert = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes_deserts.bed',
        modern_genes_island = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes_islands.bed',
        modern_genes_unique_desert = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique_deserts.bed',
        modern_genes_unique_island = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique_islands.bed',
        imputed_modern_genes_desert = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes_deserts.bed',
        imputed_modern_genes_island = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes_islands.bed',
        imputed_modern_genes_unique_desert = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique_deserts.bed',
        imputed_modern_genes_unique_island = '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique_islands.bed',
    log:
        '{path}/output/GLIMPSE_imputation/ROH_islands_deserts/modern_imputed_overlap_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_{canid_subset}-genes-unique.log'
    threads: 8
    shell:
        '''
        awk '$3=="gene"' {input.cf31_ann} > {output.genes_temp}

        bedtools intersect \
        -b {input.imputed_windows_bed_islands} \
        -a {output.genes_temp} \
        -wa -wb >> {output.imputed_genes_island}

        awk -F ';' '{{print $1}}' {output.imputed_genes_island} | uniq | awk '{{ print $10 }}' | tr -d '"' >> {output.imputed_genes_unique_island}

        bedtools intersect \
        -b {input.imputed_windows_bed_deserts} \
        -a {output.genes_temp} \
        -wa -wb >> {output.imputed_genes_desert}

        awk -F ';' '{{print $1}}' {output.imputed_genes_desert} | uniq | awk '{{ print $10 }}' | tr -d '"' >> {output.imputed_genes_unique_desert}

        bedtools intersect \
        -b {input.modern_windows_bed_islands} \
        -a {output.genes_temp} \
        -wa -wb >> {output.modern_genes_island}

        awk -F ';' '{{print $1}}' {output.modern_genes_island} | uniq | awk '{{ print $10 }}' | tr -d '"' >> {output.modern_genes_unique_island}

        bedtools intersect \
        -b {input.modern_windows_bed_deserts} \
        -a {output.genes_temp} \
        -wa -wb >> {output.modern_genes_desert}

        awk -F ';' '{{print $1}}' {output.modern_genes_desert} | uniq | awk '{{ print $10 }}' | tr -d '"' >> {output.modern_genes_unique_desert}

        bedtools intersect \
        -b {input.imputed_modern_windows_bed_islands} \
        -a {output.genes_temp} \
        -wa -wb >> {output.imputed_modern_genes_island}

        awk -F ';' '{{print $1}}' {output.imputed_modern_genes_island} | uniq | awk '{{ print $10 }}' | tr -d '"' >> {output.imputed_modern_genes_unique_island}

        bedtools intersect \
        -b {input.imputed_modern_windows_bed_deserts} \
        -a {output.genes_temp} \
        -wa -wb >> {output.imputed_modern_genes_desert}

        awk -F ';' '{{print $1}}' {output.imputed_modern_genes_desert} | uniq | awk '{{ print $10 }}' | tr -d '"' >> {output.imputed_modern_genes_unique_desert}
        '''


