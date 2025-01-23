#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

########################################################################################################################
# These rules merge the imputed downsampled vcfs and pseudohaploid downsampled with the selection of modern dogs and   #
# wolves from the reference panel, imputed are not filtered for MAF and INFO                                           #
# Then smartpca is carried out                                                                                         #
########################################################################################################################

global COVERAGE_VAL, CHROM, samples_df

# define output files of haplotoplink
DOCS_hc = ["tped", "tfam"]
DOCS = ["bed", "bim", "fam"]


###########################################################################################################################################
#
# Prepare the imputed HC and downsampled for plink format (imputed is re-calibrated for INFO score per sample and >=0.8 with MAF > 0.01)
#
###########################################################################################################################################


rule rename_sample_concordance_PH_no_filt:
    """
    Rename sample in VCF header, to include coverage info
    """
    input:
        #phased_maf_info = 'output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz'
        phased_bcf="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x.bcf",
    output:
        new_name_file=temp(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x_names.txt"
        ),
        sample_imputed_new_name="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x.bcf",
    log:
        "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x.bcf",
    shell:
        """
        echo '{wildcards.sample} {wildcards.sample}_{wildcards.coverage_val}x_imputed' > {output.new_name_file}

        bcftools reheader \
        -s {output.new_name_file} \
        {input.phased_bcf} \
        -o {output.sample_imputed_new_name} 2> {log}

        bcftools index -f {output.sample_imputed_new_name}
        """


rule rename_HC_imputed_sample_PH_no_filt:
    """
    Rename sample in VCF header of HC_imputed, to include info
    """
    input:
        #phased_maf_info = 'output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}.vcf.gz'
        phased_bcf="output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}.bcf",
    output:
        new_name_file_HC=temp(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}.txt"
        ),
        sample_imputed_new_name_HC="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}.bcf",
    log:
        "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}.bcf.log",
    shell:
        """
        echo '{wildcards.sample} {wildcards.sample}_HC_imputed' > {output.new_name_file_HC}

        bcftools reheader \
        -s {output.new_name_file_HC} \
        {input.phased_bcf} \
        -o {output.sample_imputed_new_name_HC} 2> {log}

        bcftools index -f {output.sample_imputed_new_name_HC}
        """


rule prepare_merge_imputed_genotyped_sample_PH_no_filt:
    """
    Prepare file with VCF files to merge
    """
    input:
        sample_imputed_new_name=expand(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x.bcf",
            coverage_val=COVERAGE_VAL,
            allow_missing=True,
        ),
        sample_imputed_new_name_HC="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}.bcf",
    output:
        vcf_list="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample_list.{sample}_{chrom}.txt",
    shell:
        """
        ls -v {input.sample_imputed_new_name} {input.sample_imputed_new_name_HC} >> {output.vcf_list}
        """


rule merge_renamed_imputed_genotyped_sample_PH_no_filt:
    """
    Merge all downsampled and HC version of a target sample into the same VCF
    """
    input:
        vcf_list="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample_list.{sample}_{chrom}.txt",
    output:
        merged_imputed_new_name="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.bcf",
    log:
        "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.bcf.log",
    threads: 10
    shell:
        """
        bcftools merge \
        --file-list {input.vcf_list} \
        -Ob -o {output.merged_imputed_new_name} \
        --threads {threads} 2> {log}

        bcftools index -f {output.merged_imputed_new_name}
        """


rule imputed_to_plink_no_filt:
    """
    Take imputed merged VCF and create plink format
    """
    input:
        merged_imputed_new_name="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.bcf",
    output:
        expand(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}",
    shell:
        """
        plink \
        --bcf {input.merged_imputed_new_name} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38
        """


rule imputed_prepare_correct_format_no_filt:
    """
    Fix bim and fam file columns
    """
    input:
        imputed_bim="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.bim",
        imputed_bed="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.bed",
        imputed_fam="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}.fam",
    output:
        imputed_bim_corr="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr.bim",
        imputed_fam_corr="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr.fam",
        imputed_bed_corr="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr.bed",
    shell:
        """
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.imputed_bim} | awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' > {output.imputed_bim_corr}

        awk '{{$6=2 ; print ; }}' {input.imputed_fam} > {output.imputed_fam_corr}

        scp {input.imputed_bed} {output.imputed_bed_corr}
        """


rule create_sample_list_imputed_no_filt:
    """
    Create list of imputed files to merge with reference and PH validation files from above
    """
    input:
        imputed_bed_corr="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr.bed",
    output:
        HC_prefix_list_imputed="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr.txt",
    params:
        prefix="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr",
    shell:
        """
        echo "{params.prefix}" >> {output.HC_prefix_list_imputed}
        """


rule merge_downsampled_reference_imputed_no_filt:
    """
    Merge reference and PH validation files with imputed files
    """
    input:
        #ref_bed_corr = 'output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom_con}_validation_no_missnp_all_{canid_subset}.bed',
        merged="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_validation_no_missnp_all_{canid_subset}.bed",
        HC_prefix_list_imputed="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_corr.txt",
    output:
        merged=expand(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_corr_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix_ref="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_validation_no_missnp_all_{canid_subset}",
        prefix_output="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_corr_{canid_subset}",
    shell:
        """
        plink \
        --bfile {params.prefix_ref} \
        --make-bed \
        --merge-list {input.HC_prefix_list_imputed} \
        --allow-no-sex \
        --chr-set 38 \
        --out {params.prefix_output} 
        """


rule merge_allchrom_downsampled_reference_imputed_list_no_filt:
    """ 
    Prepare list to merge chromosomes 
    """
    input:
        bed=expand(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_corr_{canid_subset}.bed",
            chrom=CHROM,
            allow_missing=True,
        ),
    params:
        prefix_in=expand(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_corr_{canid_subset}",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_list.txt",
    shell:
        """
        echo -e {params.prefix_in} | tr " " "\n" | sed '1d' >> {output.chr_list}
        """


rule merge_allchrom_downsampled_reference_imputed_no_filt:
    """
    Merge all chromosomes
    """
    input:
        chr_list="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_list.txt",
    output:
        merged=expand(
            "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.{doc}",
            doc=DOCS,
            allow_missing=True,
        ),
    params:
        prefix_chr1="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_chr1_corr_{canid_subset}",
        prefix_output="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}",
    shell:
        """
        plink \
        --bfile {params.prefix_chr1} \
        --merge-list {input.chr_list} \
        --make-bed \
        --allow-no-sex \
        --chr-set 38 \
        --out {params.prefix_output} 
        """


#### smartpca part


rule prepare_convertf_parfile_con_PH_no_filt:
    """
    Prepare convertf file to convert to eigenstrat format
    """
    input:
        bed="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.bed",
        bim_mod="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.bim",
        fam_mod="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.fam",
    output:
        convertf_file="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_convertf_parfile",
    params:
        prefix="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}",
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


rule convertf_con_PH_no_filt:
    """
    Convert to eigenstrat format
    """
    input:
        convertf_file="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_convertf_parfile",
    output:
        geno="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.eigenstratgeno",
        ind="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.ind",
        snp="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.snp",
    log:
        "output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.log",
    shell:
        """
        convertf \
        -p {input.convertf_file} 2> {log}
        """


rule set_high_low_coverage_no_filt:
    """
    Set high and low coverage categories for projection
    """
    input:
        ind="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.ind",
    output:
        ind_mod="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_mod.ind",
        poplist="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_poplist.txt",
    shell:
        """
        awk 'BEGIN {{IFS = OFS = "\t"}} {{
        if ($1 ~ "{wildcards.sample}")
            $3 = "low";
        else 
            $3 = "high";
        print
        }}' {input.ind} > {output.ind_mod}

        echo "high" > {output.poplist}
        """


rule prepare_smartpca_parfile_con_PH_no_filt:
    """
    Prepare input file for smartpca
    """
    input:
        geno="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.eigenstratgeno",
        ind_mod="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_mod.ind",
        snp="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.snp",
        poplist="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_poplist.txt",
    output:
        smartpca_parfile="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_smartpca_parfile",
    params:
        prefix="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}",
    shell:
        """
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
        """


rule smartpca_con_PH_no_filt:
    """
    Run smartpca
    """
    input:
        smartpca_parfile="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_smartpca_parfile",
    output:
        smartpca_log="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_smartpca_parfile.log",
        smartpca_eigenvalue="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_eigenval_output",
        smartpca_eigenvector="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_eigenvec_output",
    shell:
        """
        smartpca \
        -p {input.smartpca_parfile} > {output.smartpca_log}
        """


rule plot_pca_concordance_ph_no_filt:
    """
    Plot downsampled/HC imputed and non imputed samples PCA
    """
    input:
        ind_mod="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_mod.ind",
        smartpca_eigenvalue="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_eigenval_output",
        smartpca_eigenvector="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_eigenvec_output",
        ref_metadata="sample_lists/Dog_Wolf_aDNA_WG-Modern.tsv",
        concordance_metadata="sample_lists/concordance_bams_published.tsv",
    output:
        smartpca_plot="output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}.png",
    params:
        sample="{sample}",
        cov_sample=lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample=lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
    script:
        "../scripts/smartpca_ph.R"


#########################################################
#####      Estimate Sum of weighted PC distance     #####
#########################################################


rule pca_accuracy_downsampled_imputed_no_filt:
    """
    Estimate pca distances for downsampled imputed against the validation HC 
    """
    input:
        ind_mod="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_mod.ind",
        smartpca_eigenvalue="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_eigenval_output",
        smartpca_eigenvector="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_eigenvec_output",
    output:
        accuracy_plot_HC="output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_accuracy_HC.png",
        accuracy_ratio="output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ratio_pseudohaploid_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_corr_{canid_subset}_accuracy_HC.tsv",
    params:
        name="{sample}",
        cov_sample=lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample=lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
    script:
        "../scripts/smartpca_ph_distances.R"
