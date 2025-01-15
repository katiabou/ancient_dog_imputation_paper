########################################################################################################################
# These rules merge the imputed downsampled vcfs with the selection of modern dogs and wolves from the reference panel #
# Then smartpca is carried out                                                                                         #
########################################################################################################################

#define output files of make plink
DOCS = ['bed', 'bim', 'fam']


rule extract_imputed_samples_concordance:
    """
    Extract target samples from phased VCFs and put into new VCFs
    """
    input:
        phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz',
    output:
        sample_phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -s {wildcards.sample} \
        {input.phased_vcf_maf_info} \
        -Oz -o {output.sample_phased_vcf_maf_info} 2> {log}
        '''

rule rename_sample_concordance:
    """
    Rename sample in VCF header, to include coverage info
    """
    input:
        sample_phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz',
    output:
        new_name_file = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_names.txt',
        sample_phased_new_name = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo '{wildcards.sample} {wildcards.sample}_{wildcards.coverage_val}x_imputed' > {output.new_name_file}

        bcftools reheader \
        -s {output.new_name_file} \
        {input.sample_phased_vcf_maf_info} \
        -o {output.sample_phased_new_name}

        tabix -p vcf {output.sample_phased_new_name}
        '''

rule extract_HC_imputed_samples_concordance:
    """
    Extract target HC samples from phased VCFs and put into new VCFs
    """
    input:
        phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    output:
        sample_phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -s {wildcards.sample} \
        {input.phased_vcf_maf_info} \
        -Oz -o {output.sample_phased_vcf_maf_info} 2> {log}
        '''

rule rename_HC_imputed_sample:
    """
    Rename sample in VCF header of HC_imputed, to include info
    """
    input:
        sample_phased_vcf_maf_info = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_only.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    output:
        new_name_file_HC = temp('{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.txt'),
        sample_phased_new_name_HC = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo '{wildcards.sample} {wildcards.sample}_HC_imputed' > {output.new_name_file_HC}

        bcftools reheader \
        -s {output.new_name_file_HC} \
        {input.sample_phased_vcf_maf_info} \
        -o {output.sample_phased_new_name_HC} 2> {log}

        tabix -p vcf {output.sample_phased_new_name_HC}
        '''

# RENAME DOWNSAMPLED GENOTYPED VALIDATION
rule rename_downsampled_validation_sample:
    """
    Rename sample in header of downsampled validation samples
    """
    input:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_{coverage_val}x_validation_filt_qual.bcf',
    output:
        new_name_file = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_{coverage_val}x_validation_filt_qual_names.txt',
        sample_validation_new_name = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_{coverage_val}x_validation_filt_qual_renamed.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_{coverage_val}x_validation_filt_qual_renamed.bcf.log'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        echo '{wildcards.sample} {wildcards.sample}_{wildcards.coverage_val}x' > {output.new_name_file}

        bcftools reheader \
        -s {output.new_name_file} \
        {input.validation_sample_filt} \
        -o {output.sample_validation_new_name}

        bcftools index -f {output.sample_validation_new_name}
        '''


rule prepare_merge_imputed_genotyped_sample: # BE CAREFUL WITH VALIDATION SAMPLE, WHICH FILTERING HAS BEEN APPLIED ON IT
    """
    Prepare file with VCF files to merge
    """
    input:
        sample_phased_new_name = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz', coverage_val=COVERAGE_VAL, allow_missing=True),
        sample_phased_new_name_HC = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_renamed_sample.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz',
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        sample_validation_new_name = expand('{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_{coverage_val}x_validation_filt_qual_renamed.bcf', coverage_val=COVERAGE_VAL, allow_missing=True)
    output:
        vcf_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_list_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.txt'
    shell:
        '''
        ls -v {input.sample_phased_new_name} {input.sample_phased_new_name_HC} {input.validation_sample_filt_allelic} {input.sample_validation_new_name} >> {output.vcf_list}
        '''    


rule merge_renamed_imputed_genotyped_sample:
    """
    Merge all downsampled and HC version of a target sample into the same VCF
    """
    input:
        vcf_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/phased_sample_list_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.txt'
    output:
        merged_phased_new_name = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/validation_phased_renamed_sample.merged.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/validation_phased_renamed_sample.merged.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz.log'
    conda:
        '../envs/environment.yaml'
    threads: 10
    shell:
        '''
        bcftools merge \
        --file-list {input.vcf_list} \
        -Oz -o {output.merged_phased_new_name} \
        --threads {threads} 2> {log}

        tabix -p vcf {output.merged_phased_new_name}
        '''

#################### PCA steps from here on ########################

rule extract_modern_samples_canid_subset_concordance:
    """
    Extract selected samples from FILTERED reference panel (filtered regarding initial samples, sites and latest MAF filtering)
    """
    input:
        modern_canid_subset = config['modern_canid_subset'],
        #ref_concordance_sample_excl_filltags_filter_maf_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_MAF_{maf}.vcf.gz',
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf'
    output:
        #ref_concordance_sample_excl_filltags_filter_vcf_pca = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_MAF_{maf}_{canid_subset}.vcf.gz'
        ref_concordance_sample_excl_filltags_filter_vcf_pca = '{path}/output/GLIMPSE_concordance/reference_panel_rerun/{chrom_con}_ref_panel_filltags_filter_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/reference_panel_rerun/{chrom_con}_ref_panel_filltags_filter_{canid_subset}.vcf.gz.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        -S {input.modern_canid_subset} \
        --threads {threads} \
        {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_vcf_pca} 2> {log}

        bcftools index --tbi {output.ref_concordance_sample_excl_filltags_filter_vcf_pca}
        '''

rule merge_imputed_modern_vcfs_concordance:
    """
    Merge imputed and modern vcf 
    """
    input:
        ref_concordance_sample_excl_filltags_filter_vcf_pca = '{path}/output/GLIMPSE_concordance/reference_panel_rerun/{chrom_con}_ref_panel_filltags_filter_{canid_subset}.vcf.gz',
        merged_phased_new_name = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/sample_VCFs/validation_phased_renamed_sample.merged.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}.vcf.gz'
    output:
        merged_modern_imputed_concordance = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.vcf.gz'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.vcf.gz.log'
    threads: 10
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        bcftools merge \
        {input.ref_concordance_sample_excl_filltags_filter_vcf_pca} {input.merged_phased_new_name} \
        --threads {threads} \
        -Oz -o {output.merged_modern_imputed_concordance} 2> {log} 
        
        bcftools index --tbi {output.merged_modern_imputed_concordance} 
        '''

rule make_plink_con:
    """
    Create plink files
    """
    input:
        merged_modern_imputed_concordance = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.vcf.gz'
    output:
        expand('{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.log'
    conda:
        '../envs/plink19.yaml'
    shell:
        '''
        plink \
        --vcf {input.merged_modern_imputed_concordance} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        '''

rule fix_chr_column_con:
    """
    Fill info fields in order to filter based on MAF
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.bim'
    output:
        temp_mod_bim = temp('{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}-mod-temp.bim'),
        mod_bim = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}-mod.bim'
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.bim} > {output.temp_mod_bim}
        awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' {output.temp_mod_bim} > {output.mod_bim}
        '''

rule fix_fam_file_con:
    """
    Modify the last column of the .fam file from 0 to 2 (otherwise if it's 0, it does not recognize individuals properly)
    If the sample names are the path to the file, I have to add the fix_fam_file rule from before
    """
    input:
        fam = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.fam'
    output:
        fam_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}-mod.fam'
    shell:
        '''
        awk '{{$6=2 ; print ; }}' {input.fam} > {output.fam_mod}
        '''

rule prepare_convertf_parfile_con:
    """
    Prepare convertf file to convert to eigenstrat format
    """
    input:
        bed = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.bed',
        bim_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}-mod.bim',
        fam_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}-mod.fam'
    output:
        convertf_file = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_convertf_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}'
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

rule convertf_con:
    """
    Convert to eigenstrat format
    """
    input:
        convertf_file = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_convertf_parfile'
    output:
        geno = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.eigenstratgeno',
        ind = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.ind',
        snp = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.snp'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.log'
    conda:
        '../envs/eigensoft8.yaml'
    shell:
        '''
        convertf \
        -p {input.convertf_file} 2> {log}
        '''

rule get_number_of_modern: #Best to fix this rule better (in order to set the high and low coverage for the projection option)
    input:
       # ref_concordance_sample_excl_filltags_filter_maf_vcf_pca = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter_MAF_{maf}_{canid_subset}.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_vcf_pca = '{path}/output/GLIMPSE_concordance/reference_panel_rerun/{chrom_con}_ref_panel_filltags_filter_{canid_subset}.vcf.gz',
        ind = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.ind',
    output:
        ind_temp = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_temp.ind',
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_mod.ind',
        poplist = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_poplist.txt'
    conda:
        '../envs/environment.yaml'
    shell:
        '''
        mod_samples=$(bcftools query -l {input.ref_concordance_sample_excl_filltags_filter_vcf_pca} | wc -l)
        awk '{{for (i=1;i<=NF;i++) if ($i=="Case" && n++<'$mod_samples') $i="high"; print}}' {input.ind} > {output.ind_temp}
        awk '{{for (i=1;i<=NF;i++) if ($i=="Case") $i="low"; print}}' {output.ind_temp} > {output.ind_mod}
        echo "high" > {output.poplist}
        '''

rule prepare_smartpca_parfile_con:
    """
    Prepare input file for smartpca
    """
    input:
        geno = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.eigenstratgeno',
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_mod.ind',
        snp = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.snp',
        poplist = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_poplist.txt'
    output:
        smartpca_parfile = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_smartpca_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}'
    shell:
        '''
        echo "
        genotypename: {input.geno}
        snpname:      {input.snp}
        indivname:    {input.ind_mod}
        evecoutname:  {params.prefix}_eigenvec_output
        evaloutname:  {params.prefix}_eigenval_output
        numchrom: 38
        lsqproject: YES
        poplistname: {input.poplist}
        numoutlieriter: 0
        " >> {output.smartpca_parfile}
        '''

rule smartpca_con:
    """
    Run smartpca
    """
    input:
        smartpca_parfile = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_smartpca_parfile'
    output:
        smartpca_log = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_smartpca.log',
        smartpca_eigenvalue = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_eigenval_output',
        smartpca_eigenvector = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_eigenvec_output'
    conda:
        '../envs/eigensoft8.yaml'
    shell:
        '''
        smartpca \
        -p {input.smartpca_parfile} > {output.smartpca_log}
        '''

rule plot_pca_concordance:
    input:
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_mod.ind',
        smartpca_eigenvalue = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_eigenval_output',
        smartpca_eigenvector = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/modern_imputed_VCFs/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}_eigenvec_output',
        ref_metadata = config['reference_panel_metadata'],
        concordance_metadata = config['bam_targets']
    output:
        smartpca_plot = '{path}/output/GLIMPSE_concordance/PCA_concordance_re_run/plots/merged_ref-panel_validation_imputed_concordance_{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_no_low_cov_{cov_cutoff}_{canid_subset}.png'
    params:
        sample = '{sample}',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"]
    script:
        "../scripts/smartpca_ph.R"