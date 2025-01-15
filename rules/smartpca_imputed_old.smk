###########################################################################################################
# These rules merge the imputed vcf with the selection of modern dogs and wolves from the reference panel #
# Then smartpca is carried out                                                                            #
###########################################################################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']

####################################### HAVE TO CORRECT FOR HOW MAF IS FILTERED!

rule extract_imputed_sites:
    """
    Extract sites from imputed VCF, filtered for MAF and INFO
    """
    input:
        phased_vcf_info_maf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.{chrom}_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        imputed_sites_temp = temp('{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/imputed_sites_{chrom}_INFO_{info}_MAF_{maf}-temp.tsv'),
        imputed_sites = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/imputed_sites_{chrom}_INFO_{info}_MAF_{maf}.tsv'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/imputed_sites_{chrom}_INFO_{info}_MAF_{maf}.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools query \
        -f '%CHROM %POS\n' \
        {input.phased_vcf_info_maf} > {output.imputed_sites_temp} 2> {log}

        cat {output.imputed_sites_temp} | sed 's/ \+/\t/g' > {output.imputed_sites}
        '''

rule extract_modern_samples_imputed_sites:
    """
    Extract imputed sites and selected samples from FILTERED reference panel 
    """
    input:
        imputed_sites = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/imputed_sites_{chrom}_INFO_{info}_MAF_{maf}.tsv',
        modern_canid_subset = config['modern_canid_subset'],
        ref_sample_snp_filltags_filter = '{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz'
    output:
        modern_subset_vcf = temp('{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_INFO_{info}_MAF_{maf}_{canid_subset}.vcf.gz')
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/ref-panel_{chrom}_imputed_sites_filt_{canid_subset}_INFO_{info}_MAF_{maf}.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -S {input.modern_canid_subset} \
        --regions-file {input.imputed_sites} \
        --threads {threads} \
        {input.ref_sample_snp_filltags_filter} -Oz -o {output.modern_subset_vcf} 2> {log}

        bcftools index --tbi {output.modern_subset_vcf}
        '''

rule extract_imputed_dogs_wolves:
    """
    Extract either only dogs or dogs and wolves from imputed vcf
    """
    input:
        phased_vcf_info_maf = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.{chrom}_INFO_{info}_MAF_{maf}.vcf.gz',
        imputed_canid_subset = config['imputed_canid_subset']
    output:
        imputed_subset_vcf = temp('{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_phased.{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}.vcf.gz')
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_phased.{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -S {input.imputed_canid_subset} \
        --threads {threads} \
        {input.phased_vcf_info_maf} -Oz -o {output.imputed_subset_vcf} 2> {log}

        bcftools index --tbi {output.imputed_subset_vcf}
        '''

rule merge_imputed_modern_vcfs:
    """
    Merge imputed and modern vcf
    """
    input:
        modern_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_INFO_{info}_MAF_{maf}_{canid_subset}.vcf.gz',
        imputed_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_phased.{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}.vcf.gz'
    output:
        imputed_modern_subset_vcf = temp('{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}.vcf.gz')
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools merge \
        {input.imputed_subset_vcf} {input.modern_subset_vcf} \
        --threads {threads} \
        -Oz -o {output.imputed_modern_subset_vcf} 2> {log}
        '''

rule imputed_modern_vcfs_filltags:
    """
    Fill info fields in order to filter based on MAF
    """
    input:
        imputed_modern_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}.vcf.gz'
    output:
        imputed_modern_merged_vcf_filltags = temp('{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}-filltags.vcf.gz')
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}-filltags.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools +fill-tags {input.imputed_modern_subset_vcf} \
        -Oz -o {output.imputed_modern_merged_vcf_filltags} \
        --threads {threads} \
        -- -t all,F_MISSING 2> {log}
        '''

rule imputed_modern_vcfs_filltags_minmaf:
    """
    Filter based on MAF (!!! HAVE NOT FILTERED FOR F_MISSING )
    """
    input:
        imputed_modern_merged_vcf_filltags = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}-filltags.vcf.gz'
    output:
        imputed_modern_merged_vcf_filltags_minmaf = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}-filltags-minmaf.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}-filltags-minmaf.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -q {wildcards.maf}:minor \
        {input.imputed_modern_merged_vcf_filltags} \
        --threads {threads} \
        -Oz -o {output.imputed_modern_merged_vcf_filltags_minmaf} 2> {log}

        bcftools index --tbi {output.imputed_modern_merged_vcf_filltags_minmaf}
        '''

rule prepare_merged_chr_list:
    """ 
    Prepare list to merge chromosomes
    """
    input:
        imputed_modern_merged_vcf_filltags_minmaf = expand('{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_{chrom}_INFO_{info}_MAF_{maf}_{canid_subset}-filltags-minmaf.vcf.gz', chrom=CHROM, allow_missing=True)
    output:
        chr_list = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/chr_list_INFO_{info}_MAF_{maf}_{canid_subset}.txt'
    shell:
        '''
        ls -v {input.imputed_modern_merged_vcf_filltags_minmaf} >> {output.chr_list}
        '''    

rule merge_chrom:
    """ 
    Merge all filtered imputed chromosomes
    """
    input: 
        chr_list = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/chr_list_INFO_{info}_MAF_{maf}_{canid_subset}.txt'
    output:
        merged_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-filltags-minmaf.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-filltags-minmaf.log'
    threads: 10
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

rule make_plink:
    """
    Create plink files
    """
    input:
        merged_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_imputed_modern_vcf/merged_ref-panel_imputed_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-filltags-minmaf.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.log'
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

rule fix_chr_column:
    """
    Fill info fields in order to filter based on MAF
    """
    input:
        bim = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.bim'
    output:
        temp_mod_bim = temp('{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-mod-temp.bim'),
        mod_bim = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-mod.bim'
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.bim} > {output.temp_mod_bim}
        awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' {output.temp_mod_bim} > {output.mod_bim}
        '''

rule fix_fam_file:
    """
    Modify the last column of the .fam file from 0 to 2 (otherwise if it's 0, it does not recognize individuals properly)
    If the sample names are the path to the file, I have to add the fix_fam_file rule from before
    """
    input:
        fam = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.fam'
    output:
        fam_mod = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-mod.fam'
    shell:
        '''
        awk '{{$6=2 ; print ; }}' {input.fam} > {output.fam_mod}
        '''

rule prepare_convertf_parfile:
    """
    Prepare convertf file to convert to eigenstrat format
    """
    input:
        bed = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.bed',
        bim_mod = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-mod.bim',
        fam_mod = '{path}/output/GLIMPSE_imputation/PCA/plink/plink_all_chr_INFO_{info}_MAF_{maf}_{canid_subset}-mod.fam'
    output:
        convertf_file = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}_convertf_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}'
    shell:
        '''
        echo "
        genotypename: {input.bed}
        snpname:      {input.bim_mod}
        indivname:    {input.fam_mod}
        outputformat:    EIGENSTRAT                                     
        genooutfilename:   {params.prefix}.eigenstratgeno
        snpoutfilename:    {params.prefix}.snp
        indoutfilename:    {params.prefix}.ind
        familynames:       NO
        " >> {output.convertf_file}
        '''

rule convertf:
    """
    Convert to eigenstrat format
    """
    input:
        convertf_file = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}_convertf_parfile'
    output:
        geno = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.eigenstratgeno',
        ind = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.ind',
        snp = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.snp'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.log'
    conda:
        '../envs/eigensoft8.yaml'
    shell:
        '''
        convertf \
        -p {input.convertf_file} 2> {log}
        '''

rule prepare_smartpca_parfile:
    """
    Prepare input file for smartpca
    """
    input:
        geno = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.eigenstratgeno',
        ind = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.ind',
        snp = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}.snp',
    output:
        smartpca_parfile = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}_smartpca_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}'
    shell:
        '''
        echo "
        genotypename: {input.geno}
        snpname:      {input.snp}
        indivname:    {input.ind}
        evecoutname:  {params.prefix}_eigenvec_output
        evaloutname:  {params.prefix}_eigenval_output
        numchrom: 38
        numoutlieriter: 0
        " >> {output.smartpca_parfile}
        '''

rule smartpca:
    """
    Run smartpca
    """
    input:
        smartpca_parfile = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}_smartpca_parfile'
    output:
        smartpca_log = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/all_chr_INFO_{info}_MAF_{maf}_{canid_subset}_smartpca.log'
    conda:
        '../envs/eigensoft8.yaml'
    shell:
        '''
        smartpca \
        -p {input.smartpca_parfile} > {output.smartpca_log}
        '''