##################################################################################
# These rules run smarpca on the whole imputed dataset and the pseudohaploid dataset #
##################################################################################

global CHROM

# define output files of make plink
DOCS = ["bed", "bim", "fam"]


rule sites_phased:
    """
    Extract sites from the filtered recalibrated phased data
    """
    input:
        merged_phased_vcf_maf_recalibrated_info="{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz",
    output:
        merged_phased_vcf_maf_recalibrated_info_sites="{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.tsv.gz",
    log:
        "{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.tsv.gz.log",
    shell:
        """
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {input.merged_phased_vcf_maf_recalibrated_info} 2> {log} | \
        bgzip -c > {output.merged_phased_vcf_maf_recalibrated_info_sites}

        tabix -s1 -b2 -e2 {output.merged_phased_vcf_maf_recalibrated_info_sites}
        """


rule extract_sites_phased_ref_panel:
    """
    Extract filtered phased sites from ref panel
    """
    input:
        merged_phased_vcf_maf_recalibrated_info_sites="{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.tsv.gz",
        ref_panel_phased="{path}/output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
    output:
        ref_panel_imputed_sites="{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter-imputed_sites_MAF_{maf_cutoff}_recalibrated_INFO_{info}.phased.vcf.gz",
    log:
        "{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter-imputed_sites_MAF_{maf_cutoff}_recalibrated_INFO_{info}.phased.vcf.gz.log",
    threads: 4
    shell:
        """
        bcftools view {input.ref_panel_phased} \
        --regions-file {input.merged_phased_vcf_maf_recalibrated_info_sites} \
        --threads {threads} \
        -Oz -o {output.ref_panel_imputed_sites} 2> {log}

        bcftools index -f {output.ref_panel_imputed_sites}
        """


rule merge_maf_ref_panel_imputed:
    """
    Merge phased filt ref panel with filtered imputed
    """
    input:
        merged_phased_vcf_maf_recalibrated_info="{path}/output/GLIMPSE_imputation/GLIMPSE_phased_merged/merged_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz",
        ref_panel_imputed_sites="{path}/output/GLIMPSE_imputation/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter-imputed_sites_MAF_{maf_cutoff}_recalibrated_INFO_{info}.phased.vcf.gz",
    output:
        modern_imputed="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz",
    shell:
        """
        bcftools merge \
        {input.ref_panel_imputed_sites} \
        {input.merged_phased_vcf_maf_recalibrated_info} \
        -Oz -o {output.modern_imputed} \
        --threads {threads}

        bcftools index -f {output.modern_imputed}
        """


rule prep_canid_list:
    """
    Prepare canid group list of samples (both modern and ancient)
    """
    input:
        modern_canid_subset="sample_lists/ref_panel_filt_{canid_subset}.tsv",
        imputed_canid_subset="sample_lists/names_imputed_{canid_subset}.tsv",
    output:
        modern_imputed_canid_subset="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/modern_imputed_{canid_subset}.txt",
    shell:
        """
        cat {input.modern_canid_subset} {input.imputed_canid_subset} > {output.modern_imputed_canid_subset}
        """


rule canid_subset_merged_vcf:
    """
    Make canid specific groups
    """
    input:
        modern_imputed="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}.vcf.gz",
        modern_imputed_canid_subset="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/modern_imputed_{canid_subset}.txt",
    output:
        modern_imputed_subset="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz",
    threads: 4
    shell:
        """
        bcftools view \
        -S {input.modern_imputed_canid_subset} \
        --threads {threads} \
        {input.modern_imputed} -Oz -o {output.modern_imputed_subset}

        bcftools index -f {output.modern_imputed_subset}
        """


# rule canid_subset_merged_vcf_recalculate:
#     """
#     Re-calculate fields
#     """
#     input:
#         modern_imputed_subset = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz',
#     output:
#         modern_imputed_subset_filltags = temp('{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-filltags.vcf.gz'),
#     threads: 10
#     shell:
#         '''
#         bcftools +fill-tags {input.modern_imputed_subset} \
#         --threads {threads} \
#         -Oz -o {output.modern_imputed_subset_filltags} \
#         -- -t all,F_MISSING
#         '''

# rule canid_subset_merged_vcf_filter_no_missing:
#     """
#     Filter for no missing sites (F_MISSING<0.01)
#     """
#     input:
#         modern_imputed_subset_filltags = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-filltags.vcf.gz',
#     output:
#         modern_imputed_subset_no_missing = '{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_no-missing.vcf.gz',
#     threads: 10
#     shell:
#         '''
#         bcftools view -i 'F_MISSING<0.01' {input.modern_imputed_subset_filltags} \
#         --threads {threads} \
#         -Oz -o {output.modern_imputed_subset_no_missing}

#         tabix -p vcf {output.modern_imputed_subset_no_missing}
#         '''


rule canid_subset_prepare_merged_chr_list:
    """ 
    Prepare list to merge chromosomes
    """
    input:
        modern_imputed_subset=expand(
            "{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.{chrom}_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/chr_list.MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.txt",
    shell:
        """
        ls -v {input.modern_imputed_subset} >> {output.chr_list}
        """


rule merge_chrom:
    """ 
    Merge all chromosomes
    """
    input:
        chr_list="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/chr_list.MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.txt",
    output:
        merged_vcf="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz",
    threads: 4
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Oz -o {output.merged_vcf} \
        --threads {threads}

        bcftools index -f {output.merged_vcf}
        """


rule make_plink:
    """
    Create plink files
    """
    input:
        merged_vcf="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.vcf.gz",
    output:
        expand(
            "{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}",
    log:
        "{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.log",
    shell:
        """
        plink \
        --vcf {input.merged_vcf} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38 2> {log}
        """


rule fix_chr_column:
    """
    Fix chromosome and snps columns
    """
    input:
        bim="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.bim",
    output:
        temp_mod_bim=temp(
            "{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod-temp.bim"
        ),
        mod_bim="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.bim",
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.bim} > {output.temp_mod_bim}
        awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' {output.temp_mod_bim} > {output.mod_bim}
        """


rule fix_fam_file:
    """
    Modify the last column of the .fam file from 0 to 2 (otherwise if it's 0, it does not recognize individuals properly)
    If the sample names are the path to the file, I have to add the fix_fam_file rule from before
    """
    input:
        fam="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.fam",
    output:
        fam_mod="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.fam",
    shell:
        """
        awk '{{$6=2 ; print ; }}' {input.fam} > {output.fam_mod}
        """


rule prepare_convertf_parfile:
    """
    Prepare convertf file to convert to eigenstrat format
    """
    input:
        bed="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.bed",
        bim_mod="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.bim",
        fam_mod="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.fam",
    output:
        convertf_file="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_convertf_parfile",
    params:
        prefix="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}",
    shell:
        """
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
        """


rule convertf:
    """
    Convert to eigenstrat format
    """
    input:
        convertf_file="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_convertf_parfile",
    output:
        geno="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.eigenstratgeno",
        ind="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.ind",
        snp="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.snp",
    log:
        "{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.log",
    shell:
        """
        convertf \
        -p {input.convertf_file} 2> {log}
        """


rule prepare_smartpca_parfile:
    """
    Prepare input file for smartpca
    """
    input:
        geno="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.eigenstratgeno",
        ind="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.ind",
        snp="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}.snp",
    output:
        smartpca_parfile="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_parfile",
    params:
        prefix="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}",
    shell:
        """
        echo "
        genotypename: {input.geno}
        snpname:      {input.snp}
        indivname:    {input.ind}
        evecoutname:  {params.prefix}_eigenvec_output
        evaloutname:  {params.prefix}_eigenval_output
        numchrom: 38
        numoutlieriter: 0
        " >> {output.smartpca_parfile}
        """


rule smartpca:
    """
    Run smartpca
    """
    input:
        smartpca_parfile="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_parfile",
    output:
        smartpca_log="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca.log",
        smartpca_eigenvec="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_eigenvec_output",
        smartpca_eigenval="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_eigenval_output",
    shell:
        """
        smartpca \
        -p {input.smartpca_parfile} > {output.smartpca_log}
        """


rule plot_smartpca:
    """
    Plot smartpca 
    """
    input:
        smartpca_eigenval="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_eigenval_output",
        smartpca_eigenvec="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_eigenvec_output",
        fam_mod="{path}/output/GLIMPSE_imputation/PCA/merged_phased_modern_vcf/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}-mod.fam",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        bam_metadata="sample_lists/Dog_Wolf_aDNA_WG-Master.tsv",
    params:
        canid_group="{canid_subset}",
    output:
        smartpca_P1P2="{path}/output/GLIMPSE_imputation/plots/PCA/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_P1P2.png",
        smartpca_P1P2_labelled="{path}/output/GLIMPSE_imputation/plots/PCA/merged_modern_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_{canid_subset}_smartpca_P1P2_labelled.png",
    script:
        "../scripts/smartpca_phased_dataset.R"
