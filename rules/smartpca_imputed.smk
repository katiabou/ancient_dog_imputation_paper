###########################################################################################################
# These rules merge the phased vcf with the selection of modern dogs and wolves from the reference panel #
# Then smartpca is carried out                                                                            #
###########################################################################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']


#rule extract_phased_sites:
#    """
#    Extract sites from phased VCF, filtered for MAF and INFO 
#    """
#    input:
#        #phased_vcf_maf_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased/merged_phased.{chrom}_INFO_{info}_MAF_{maf_cutoff}_no_low_cov_{cov_cutoff}.vcf.gz',
#        merged_phased_vcf_maf_recalibrated_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz'
#    output:
#        imputed_sites_temp = temp('{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_INFO_{info}_MAF_{maf_cutoff}_no_low_cov_{cov_cutoff}-temp.tsv'),
#        imputed_sites = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_INFO_{info}_MAF_{maf_cutoff}_no_low_cov_{cov_cutoff}.tsv'
#    log:
#        '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_INFO_{info}_MAF_{maf_cutoff}_no_low_cov_{cov_cutoff}.log'
#    #conda:
#        #'../envs/environment.yaml'
#    shell:
#        '''
#        bcftools query \
#        -f '%CHROM %POS\n' \
#        {input.merged_phased_vcf_maf_recalibrated_info} > {output.imputed_sites_temp} 2> {log}

#        cat {output.imputed_sites_temp} | sed 's/ \+/\t/g' > {output.imputed_sites}
#        '''

rule extract_modern_samples_canid_subset: #make this temp output file once finished!!!!!
    """
    Extract either only dogs or wolves from reference panel 
    """
    input:
        ref_sample_snp_filltags_filter = '{path}/output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.vcf.gz',
        modern_canid_subset = config['modern_canid_subset'],
    output:
        modern_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_{canid_subset}.vcf.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -S {input.modern_canid_subset} \
        --trim-alt-alleles \
        --threads {threads} \
        {input.ref_sample_snp_filltags_filter} -Oz -o {output.modern_subset_vcf} 2> {log}

        bcftools index --tbi {output.modern_subset_vcf}
        '''

rule extract_imputed_dogs_wolves: #make this temp output file once finished!!!!!
    """
    Extract either only dogs or wolves from imputed vcf
    """
    input:
        merged_phased_vcf_maf_recalibrated_info = '{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz',
        imputed_canid_subset = config['imputed_canid_subset']
    output:
        phased_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_vcf/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_phased_vcf/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -S {input.imputed_canid_subset} \
        --trim-alt-alleles \
        --threads {threads} \
        {input.merged_phased_vcf_maf_recalibrated_info} -Oz -o {output.phased_subset_vcf} 2> {log}

        bcftools index --tbi {output.phased_subset_vcf}
        '''

rule merge_phased_modern_vcfs:
    """
    Merge phased and modern vcfs
    """
    input:
        modern_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/modern_vcf/ref-panel_{chrom}_sample-snp_filltags_filter_{canid_subset}.vcf.gz',
        phased_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_vcf/merged_phased.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz'
    output:
        phased_modern_subset_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_phased_ref_panel.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_phased_ref_panel.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools merge \
        {input.phased_subset_vcf} {input.modern_subset_vcf} \
        --threads {threads} \
        -Oz -o {output.phased_modern_subset_vcf} 2> {log}
        '''

rule prepare_merged_chr_list:
    """ 
    Prepare list to merge chromosomes
    """
    input:
        phased_modern_subset_vcf = expand('{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_phased_ref_panel.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz', chrom=CHROM, allow_missing=True)
    output:
        chr_list = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/chr_list.MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.txt'
    shell:
        '''
        ls -v {input.phased_modern_subset_vcf} >> {output.chr_list}
        '''    

rule merge_chrom:
    """ 
    Merge all filtered phased chromosomes
    """
    input: 
        chr_list = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/chr_list.MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.txt'
    output:
        merged_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
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
        merged_vcf = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.log'
    #conda:
        #'../envs/plink19.yaml'
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
        bim = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.bim'
    output:
        temp_mod_bim = temp('{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod-temp.bim'),
        mod_bim = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.bim'
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
        fam = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.fam'
    output:
        fam_mod = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.fam'
    shell:
        '''
        awk '{{$6=2 ; print ; }}' {input.fam} > {output.fam_mod}
        '''

rule prepare_convertf_parfile:
    """
    Prepare convertf file to convert to eigenstrat format
    """
    input:
        bed = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.bed',
        bim_mod = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.bim',
        fam_mod = '{path}/output/GLIMPSE_imputation/PCA/plink/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.fam'
    output:
        convertf_file = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_convertf_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}'
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
        convertf_file = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_convertf_parfile'
    output:
        geno = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.eigenstratgeno',
        ind = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.ind',
        snp = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.snp'
    log:
        '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.log'
    #conda:
        #'../envs/eigensoft8.yaml'
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
        geno = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.eigenstratgeno',
        ind = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.ind',
        snp = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.snp',
    output:
        smartpca_parfile = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}'
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
        smartpca_parfile = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_parfile'
    output:
        smartpca_log = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca.log'
    #conda:
        #'../envs/eigensoft8.yaml'
    shell:
        '''
        smartpca \
        -p {input.smartpca_parfile} > {output.smartpca_log}
        '''


rule plot_smartpca:
    """
    plot smartpca
    """
    input:
        ind_dog = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_dogs.ind',
        smartpca_eigenvalue_dog = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_dogs_eigenval_output', 
        smartpca_eigenvector_dog = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_dogs_eigenvec_output', 
        ind_wolf = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_wolves.ind',
        smartpca_eigenvalue_wolf = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_wolves_eigenval_output', 
        smartpca_eigenvector_wolf = '{path}/output/GLIMPSE_imputation/PCA/eigensoft/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_wolves_eigenvec_output', 
        ref_metadata = config['reference_panel_metadata'],
        bam_metadata = config['bam_imputation_meta']
    output:
        smartpca_plot_dog = '{path}/output/GLIMPSE_imputation/plots/PCA/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_dogs.png',
        smartpca_plot_wolf = '{path}/output/GLIMPSE_imputation/plots/PCA/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_wolves.png'
    script:
        "../scripts/smartpca_phased_dataset.R"
