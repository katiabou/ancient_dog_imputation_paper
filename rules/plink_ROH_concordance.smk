########################################################
#  Plink ROH to test target imputed against validation #
########################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']


rule transversions_phased_concordance:
    """
    Only take transversions from filtered phased concordance data
    """
    input:
        phased_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        tranversion_sites = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions.tsv.gz',
        phased_transversions = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions.vcf.gz.log'
    #conda:
    #    '../envs/environment.yaml'
    threads: 10
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.phased_maf_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.phased_maf_info} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.phased_transversions} 2> {log}

        bcftools index -f {output.phased_transversions}
        '''

rule make_plink_transversions_phased_concordance:
    """
    Prepare file format for plink
    """
    input:
        phased_transversions = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions'
    #conda:
        #'../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.phased_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 
        '''

rule roh_transversions_phased_concordance:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}.hom'
    params:
        prefix_in = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions',
        prefix_out = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}.hom.log'
    #conda:
        #'../envs/plink19.yaml'
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

rule merge_roh_transversions_phased_concordance_prep:
    """
    Prepare ROH output files for plot
    """
    input:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}.hom'
    output:
        temp1_roh_phased = temp('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp1.hom'),
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp.hom'
    shell:
        '''
        scp {input.phased_roh} {output.temp1_roh_phased}
        grep -qxF '{wildcards.sample}' {output.temp1_roh_phased} || printf "{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output.temp1_roh_phased}
        awk -F'\t' '{{print $0 "\t" "{wildcards.coverage_val}x"}}' {output.temp1_roh_phased} > {output.temp_roh_phased}
        sed -i '1s/{wildcards.coverage_val}x/cov/' {output.temp_roh_phased}
        '''

rule estimate_chr_size:
    """
    Estimate chromosome size for plotting x axis 
    """
    input:
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta'
    output:
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome'
    shell:
        '''
        faidx {input.ref_fasta_chr} -i chromsizes > {output.ref_fasta_chr_size}
        '''
   

#### Run ROHs for HC concordance transversions only ####

rule transversions_phased_HC:
    """
    Only take transversions from phased data
    """
    input:
        phased_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        tranversion_sites = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.tsv.gz',
        phased_transversions = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.vcf.gz.log'
    threads: 10
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools query -e 'REF="A" && ALT="G" || REF="G" && ALT="A" || REF="C" && ALT="T" || REF="T" && ALT="C"' \
        -f'%CHROM\t%POS\n' {input.phased_maf_info} | bgzip -c > {output.tranversion_sites}

        tabix -s1 -b2 -e2 {output.tranversion_sites}

        bcftools view {input.phased_maf_info} \
        --regions-file {output.tranversion_sites} \
        --threads {threads} \
        -Oz -o {output.phased_transversions} 2> {log}

        bcftools index -f {output.phased_transversions}
        '''

rule make_plink_transversions_phased_HC:
    """
    Prepare file format for plink
    """
    input:
        phased_transversions = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.log'
    #conda:
        #'../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.phased_transversions} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_transversions_phased_HC:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}.hom'
    params:
        prefix_in = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions',
        prefix_out = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}.hom.log'
    #conda:
        #'../envs/plink19.yaml'
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

rule merge_roh_transversions_phased_prep_HC:
    """
    Prepare ROH output files for plot
    """
    input:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}.hom'
    output:
        temp1_roh_phased = temp('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp1.hom'),
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp.hom'
    shell:
        '''
        scp {input.phased_roh} {output.temp1_roh_phased}
        grep -qxF '{wildcards.sample}' {output.temp1_roh_phased} || printf "{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output.temp1_roh_phased}
        awk -F'\t' '{{print $0 "\t" "HC_imputed"}}' {output.temp1_roh_phased} > {output.temp_roh_phased}
        sed -i '1s/HC_imputed/cov/' {output.temp_roh_phased}
        '''

####################################################        
#### Run ROHs for transversions and transitions ####
####################################################        

rule make_plink_all_sites_phased_concordance:
    """
    Prepare file format for plink
    """
    input:
        phased_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites.log'
    #conda:
        #'../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.phased_maf_info} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_all_sites_phased_concordance:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}.hom'
    params:
        prefix_in = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites',
        prefix_out = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}.hom.log'
    #conda:
        #'../envs/plink19.yaml'
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

rule merge_roh_all_sites_phased_concordance_prep:
    """
    Prepare ROH output files for plot
    """
    input:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}.hom'
    output:
        temp1_roh_phased = temp('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp1.hom'),
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp.hom'
    shell:
        '''
        scp {input.phased_roh} {output.temp1_roh_phased}
        grep -qxF '{wildcards.sample}' {output.temp1_roh_phased} || printf "{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output.temp1_roh_phased}
        awk -F'\t' '{{print $0 "\t" "{wildcards.coverage_val}x"}}' {output.temp1_roh_phased} > {output.temp_roh_phased}
        sed -i '1s/{wildcards.coverage_val}x/cov/' {output.temp_roh_phased}
        '''


#### Run ROHs for HC concordance all sites ####

rule make_plink_all_sites_phased_HC:
    """
    Prepare file format for plink
    """
    input:
        phased_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites.log'
    #conda:
        #'../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.phased_maf_info} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule roh_all_sites_phased_HC:
    """
    Run ROH estimation with plink
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites.bim'
    output:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}.hom'
    params:
        prefix_in = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites',
        prefix_out = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}'
    log:
        '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}.hom.log'
    #conda:
        #'../envs/plink19.yaml'
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

rule merge_roh_all_sites_phased_prep_HC:
    """
    Prepare ROH output files for plot
    """
    input:
        phased_roh = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}.hom'
    output:
        temp1_roh_phased = temp('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp1.hom'),
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp.hom'
    shell:
        '''
        scp {input.phased_roh} {output.temp1_roh_phased}
        grep -qxF '{wildcards.sample}' {output.temp1_roh_phased} || printf "{wildcards.sample}\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" >> {output.temp1_roh_phased}
        awk -F'\t' '{{print $0 "\t" "HC_imputed"}}' {output.temp1_roh_phased} > {output.temp_roh_phased}
        sed -i '1s/HC_imputed/cov/' {output.temp_roh_phased}
        '''

################################################################
##### Plot concordance downsampled, HC and validation ROHs #####
################################################################

rule plot_imputed_concordance_roh_transversions: 
    """
    Plot ROH transversions for imputed concordance, imputed HC, and validation
    """
    input:
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp.hom',
        temp_roh_concordance_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp.hom', coverage_val=COVERAGE_VAL, allow_missing=True),
        temp_roh_validation = '{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink-temp.hom',
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_ROH.png'
    params:
        files_concordance_phased=lambda wildcards, input: ','.join(input.temp_roh_concordance_phased),
        files_validation=lambda wildcards, input: input.temp_roh_validation,
        files_phased=lambda wildcards, input: input.temp_roh_phased,
        name = '{sample}',
        path_script = '{path}/scripts',
        chrom = '{chrom_con}',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
        site_type = 'Transversions'
    shell:
        "Rscript {params.path_script}/ROH_plotting.R {params.files_concordance_phased} {params.files_validation} {params.files_phased} {params.name} {params.chrom} {params.cov_sample} {params.info_sample} {params.site_type} {output.plot} {input.ref_fasta_chr_size}"

    

rule plot_imputed_concordance_roh_all_sites:
    """
    Plot ROH all sites for imputed concordance, imputed HC, and validation
    """
    input:
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp.hom',
        temp_roh_concordance_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp.hom', coverage_val=COVERAGE_VAL, allow_missing=True),
        temp_roh_validation = '{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink-temp.hom',
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome',
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_ROH.png'
    params:
        files_concordance_phased=lambda wildcards, input: ','.join(input.temp_roh_concordance_phased),
        files_validation=lambda wildcards, input: input.temp_roh_validation,
        files_phased=lambda wildcards, input: input.temp_roh_phased,
        name = '{sample}',
        path_script = '{path}/scripts',
        chrom = '{chrom_con}',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
        site_type = 'Transversions+transitions'
    shell:
        "Rscript {params.path_script}/ROH_plotting.R {params.files_concordance_phased} {params.files_validation} {params.files_phased} {params.name} {params.chrom} {params.cov_sample} {params.info_sample} {params.site_type} {output.plot} {input.ref_fasta_chr_size}"
        

################################################################
#####      Estimate F1 and MCC stats on TP, FP, TN, FN     #####
################################################################

rule roh_accuracy_downsampled_imputed_transversions:
    """
    Estimate accuracy metrics for downsampled imputed against the validation HC for transversions
    """
    input:
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp.hom',
        temp_roh_concordance_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}-temp.hom', coverage_val=COVERAGE_VAL, allow_missing=True),
        temp_roh_validation = '{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_validation_qual_dp_ab_filt_transversions_hom_win_het_{hom_win_het}_plink-temp.hom',
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome'
    output:
        accuracy_seg = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_accuracy_segment.tsv',
        accuracy_len = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_accuracy_length.tsv',
        accuracy_plot = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_accuracy.png',
        accuracy_plot2 = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_accuracy_FDR_sensitivity_specificity.png',
        accuracy_plot3 = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_transversions_hom_win_het_{hom_win_het}_accuracy_sensitivity_specificity.png'
    params:
        files_concordance_phased=lambda wildcards, input: ','.join(input.temp_roh_concordance_phased),
        files_validation=lambda wildcards, input: input.temp_roh_validation,
        files_phased=lambda wildcards, input: input.temp_roh_phased,
        name = '{sample}',
        path_script = '{path}/scripts',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
        site_type = 'Transversions'
    shell:
        '''
        Rscript {params.path_script}/ROH_concordance_accuracy_copy.R \
        {params.files_concordance_phased} \
        {params.files_validation} \
        {params.files_phased} \
        {params.name} \
        {params.cov_sample} \
        {params.info_sample} \
        {params.site_type} \
        {output.accuracy_seg} \
        {output.accuracy_len} \
        {output.accuracy_plot} \
        {output.accuracy_plot2} \
        {output.accuracy_plot3} \
        {input.ref_fasta_chr_size}
        '''


rule roh_accuracy_downsampled_imputed_all_sites:
    """
    Estimate accuracy metrics for downsampled imputed against the validation HC for all sites
    """
    input:
        temp_roh_phased = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp.hom',
        temp_roh_concordance_phased = expand('{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}-temp.hom', coverage_val=COVERAGE_VAL, allow_missing=True),
        temp_roh_validation = '{path}/output/GLIMPSE_concordance/ROH_validation/{sample}_{chrom_con}_validation_filt_qual_dp_ab_all_sites_hom_win_het_{hom_win_het}_plink-temp.hom',
        ref_fasta_chr_size = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}_size.genome'
    output:
        accuracy_seg = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy_segment.tsv',
        accuracy_len = '{path}/output/GLIMPSE_concordance/ROH_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy_length.tsv',
        accuracy_plot = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy.png',
        accuracy_plot2 = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy_FDR_sensitivity_specificity.png',
        accuracy_plot3 = '{path}/output/GLIMPSE_concordance/plots/ROH_concordance/{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy_sensitivity_specificity.png'
    params:
        files_concordance_phased=lambda wildcards, input: ','.join(input.temp_roh_concordance_phased),
        files_validation=lambda wildcards, input: input.temp_roh_validation,
        files_phased=lambda wildcards, input: input.temp_roh_phased,
        name = '{sample}',
        path_script = '{path}/scripts',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
        site_type = 'Transversions+transitions'
    shell:
        '''
       Rscript {params.path_script}/ROH_concordance_accuracy_copy.R \
        {params.files_concordance_phased} \
        {params.files_validation} \
        {params.files_phased} \
        {params.name} \
        {params.cov_sample} \
        {params.info_sample} \
        {params.site_type} \
        {output.accuracy_seg} \
        {output.accuracy_len} \
        {output.accuracy_plot} \
        {output.accuracy_plot2} \
        {output.accuracy_plot3} \
        {input.ref_fasta_chr_size}
        '''
