########################################################################################################################
# These rules merge the imputed downsampled vcfs and pseudohaploid downsampled with the selection of modern dogs and   #
# wolves from the reference panel                                                                                      #                             
# Then smartpca is carried out                                                                                         #
########################################################################################################################

#define output files of haplotoplink
DOCS_hc = ['tped', 'tfam']
DOCS = ['bed','bim','fam']


rule merge_reference_HC_genotyped:
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom}_ref_panel_filltags_filter.phased.bcf',
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_validation_filt_qual_dp_ab.bcf',
    output:
        ref_panel_HC_gen = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.bcf.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.bcf.tsv'
    shell:
        '''
        bcftools merge \
        {input.ref_concordance_sample_excl_filltags_filter} {input.validation_sample_filt_allelic} \
        --threads {threads} \
        -Ob -o {output.ref_panel_HC_gen} 2> {log}

        bcftools index -f {output.ref_panel_HC_gen}
        '''


rule reference_HC_gen_to_plink:
    """
    Take merged reference and HC genotyped and create plink format
    """
    input:
        ref_panel_HC_gen = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.bcf'
    output:
        expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab'
    shell:
        '''
        plink \
        --bcf {input.ref_panel_HC_gen} \
        --make-bed \
        --allow-extra-chr \
        --out {params.prefix} \
        --double-id \
        --chr-set 38
        '''

rule reference_prepare_correct_format:
    """
    Fix bim and fam file columns
    """
    input:
        ref_bim = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.bim',
        ref_bed = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.bed',
        ref_fam = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab.fam'
    output:
        ref_bim_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr.bim',
        ref_fam_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr.fam',
        ref_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr.bed'
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.ref_bim} | awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' > {output.ref_bim_corr}

        awk '{{$6=2 ; print ; }}' {input.ref_fam} > {output.ref_fam_corr}

        scp {input.ref_bed} {output.ref_bed_corr}
        '''


rule haplo_to_plink:
    """
    Convert the pseudohaploid sites from angsd to plink (tped and tfam)
    """
    input:
        validation_sample_filt_pseudohaploid_sites = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_{coverage_val}x_validation.haplo.gz',
    output:
        expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation.{doc_hc}', doc_hc=DOCS_hc, allow_missing=True),
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation',
    shell:
        '''
        haploToPlink {input.validation_sample_filt_pseudohaploid_sites} {params.prefix}
        '''

rule tped_to_bed:
    """
    Turn into bed, bim, fam for reference panel sites
    """
    input:
        tped = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation.tped',
        tfam = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation.tfam',
        ref_bim_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr.bim',
        #ref_bim_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom_con}_ref_panel_filltags_filter_corr.bim'
    output:
        expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink'
    shell:
        '''
        plink \
        --make-bed \
        --extract {input.ref_bim_corr} \
        --tped {input.tped} \
        --tfam {input.tfam} \
        --missing-genotype N \
        --output-missing-genotype 0 \
        --out {params.prefix} \
        --chr-set 38 
        '''

rule ph_called_prepare_correct_format:
    """
    Fix bim and fam columns
    """
    input:
        bim = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink.bim',
        fam = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink.fam',
        bed = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink.bed'
    output:
        bim_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.bim',
        fam_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.fam',
        bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.bed'
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.bim} > {output.bim_corr}

        sed 's/ind0/{wildcards.sample}_{wildcards.coverage_val}x/g' {input.fam} | awk '{{$6=2 ; print ; }}' > {output.fam_corr}

        scp {input.bed} {output.bed_corr}
        '''

rule create_sample_list_downsampled:  
    """
    Create list containing downsampled prefixes for following merge step
    """
    input:
        bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.bed'
    output:
        HC_prefix_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/list_{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.txt'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr'
    shell:
        '''
        echo {params.prefix} > {output.HC_prefix_list}
        '''

rule merge_downsampled_reference_HC_gen_missnp:
    """
    First attempt to merge downsampled ph sites with reference panel and HC genotyped, to get the missnp file (in order to avoid poly-allelic sites. The || true option avoids errors occuring from bash strict mode). 
    The touch option makes sure that there is a missnp file produced even if all the sites are ok within that sample
    """
    input:
        ref_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr.bed',
        HC_prefix_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/list_{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.txt'
    output:
        missnps = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_{coverage_val}x_validation.missnp'
    params:
        prefix_ref = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr',
        prefix_output = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_{coverage_val}x_validation'
    shell:
        '''
        plink \
        --bfile {params.prefix_ref} \
        --merge-list {input.HC_prefix_list} \
        --make-just-bim \
        --chr-set 38 \
        --out {params.prefix_output} || true

        touch {params.prefix_output}.missnp
        '''

rule plink_exclude_missnp:
    """
    Exclude the missnp sites 
    """
    input:
        missnp = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_{coverage_val}x_validation.missnp',
        bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.bed'
    output:
        expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{sample}_{chrom}_{coverage_val}x_validation_no_missnp.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix_input = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr',
        prefix_output = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{sample}_{chrom}_{coverage_val}x_validation_no_missnp'
    shell:
        '''
        plink \
        --make-bed \
        --exclude {input.missnp} \
        --bfile {params.prefix_input} \
        --out {params.prefix_output} \
        --chr-set 38
        '''

rule check_percentage_of_dropped_sites:
    input:
        missnp = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_{coverage_val}x_validation.missnp',
        bim_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ph_called_plink/{sample}_{chrom}_{coverage_val}x_validation_ref_sites_plink_corr.bim'
    output:
        perc = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/percentage_missnp/merged_ph_called_modern-{sample}_{chrom}_{coverage_val}x_validation.missnp'
    shell:
        '''
        test=$(wc -l < {input.missnp})
        test2=$(wc -l < {input.bim_corr})
        frac=$(echo "$test $test2" | awk '{{print $1/$2}}')
        echo $frac > {output.perc}
        '''

rule create_sample_list_no_missnp:
    """
    Create list containing downsampled prefixes without missnps for following merge step
    """
    input:
        bed_downsampled = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{sample}_{chrom}_{coverage_val}x_validation_no_missnp.bed', coverage_val=COVERAGE_VAL, allow_missing=True),
    output:
        HC_prefix_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/list_{sample}_{chrom}_validation_no_missnp.txt'
    params:
        prefix = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{sample}_{chrom}_{coverage_val}x_validation_no_missnp', coverage_val=COVERAGE_VAL, allow_missing=True),
    shell:
        '''
        echo "{params.prefix}" | tr " " "\n" >> {output.HC_prefix_list}
        '''


#### rules for splitting reference into canid subset ###

rule prepare_canid_subset_file:
    input:
        modern_canid_subset = config['modern_canid_subset'],
    output:
        modern_canid_subset_plink = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/ref_panel_filt_{canid_subset}_plink.txt'
    shell:
        '''
        awk -F'\t' '{{ print $1 "  " $1 }}' {input.modern_canid_subset} > {output.modern_canid_subset_plink}
        '''

rule prepare_canid_subset_fileb:
    input:
        modern_canid_subset_plink = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/ref_panel_filt_{canid_subset}_plink.txt'
    output:
        sample_temp = temp('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/ref_panel_filt_{sample}_{canid_subset}_plink-temp.txt'),
        modern_canid_subset_plink_sample = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/ref_panel_filt_{sample}_{canid_subset}_plink.txt'
    params:
        sample_name =  '{sample}'
    shell:
        '''
        echo '{params.sample_name}\t{params.sample_name}' > {output.sample_temp}
        cat {input.modern_canid_subset_plink} {output.sample_temp} >> {output.modern_canid_subset_plink_sample}
        '''

rule select_canid_subset_from_plink:
    """
    Extract selected samples from FILTERED/FILLTAGS reference panel 
    """
    input:
        modern_canid_subset_plink_sample = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/ref_panel_filt_{sample}_{canid_subset}_plink.txt',
        ref_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr.bed',
        #ref_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom_con}_ref_panel_filltags_filter_corr.bed'
    output:
        ref_bed_corr_subset = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom}_ref_panel_filltags_filter_corr_{sample}_{canid_subset}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix_input = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/{chrom}_ref_panel_{sample}_validation_filt_qual_dp_ab_corr',
        prefix_output = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom}_ref_panel_filltags_filter_corr_{sample}_{canid_subset}'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom}_ref_panel_filltags_filter_corr_{sample}_{canid_subset}.log'
    shell:
        '''
        plink \
        -bfile {params.prefix_input} \
        --make-bed \
        --keep {input.modern_canid_subset_plink_sample} \
        --allow-no-sex \
        --chr-set 38 \
        --out {params.prefix_output} 
        '''

rule merge_downsampled_reference:
    """
    Merge reference panel subsets with the downsampled PH validation files
    """
    input:
        ref_bed_corr_subset = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom}_ref_panel_filltags_filter_corr_{sample}_{canid_subset}.bed',
        HC_prefix_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/list_{sample}_{chrom}_validation_no_missnp.txt'
        #HC_prefix_list_all = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/list_{sample}_{chrom_con}_validation_no_missnp_all.txt'
    output:
        merged = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_validation_no_missnp_all_{canid_subset}.{doc}', doc=DOCS, allow_missing=True) 
    params:
        prefix_ref = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/modern_plink/{chrom}_ref_panel_filltags_filter_corr_{sample}_{canid_subset}',
        prefix_output = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_validation_no_missnp_all_{canid_subset}'
    shell:
        '''
        plink \
        --bfile {params.prefix_ref} \
        --make-bed \
        --merge-list {input.HC_prefix_list} \
        --allow-no-sex \
        --chr-set 38 \
        --out {params.prefix_output} 
        '''



###########################################################################################################################################
#
# Prepare the imputed HC and downsampled for plink format (imputed is re-calibrated for INFO score per sample and >=0.8 with MAF > 0.01)
#
###########################################################################################################################################

rule rename_sample_concordance_PH:
    """
    Rename sample in VCF header, to include coverage info
    """
    input:
       # imputed_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-only_sample_INFO_{info}_MAF_{maf}.vcf.gz'
        phased_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        new_name_file = temp('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}_names.txt'),
        sample_imputed_new_name = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.bcf'
    shell:
        '''
        echo '{wildcards.sample} {wildcards.sample}_{wildcards.coverage_val}x_imputed' > {output.new_name_file}

        bcftools reheader \
        -s {output.new_name_file} \
        {input.phased_maf_info} \
        -o {output.sample_imputed_new_name}

        bcftools index -f {output.sample_imputed_new_name}
        '''


rule rename_HC_imputed_sample_PH:
    """
    Rename sample in VCF header of HC_imputed, to include info
    """
    input:
        phased_maf_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom}_INFO_{info}_MAF_{maf}.vcf.gz'
    output:
        new_name_file_HC = temp('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_INFO_{info}_MAF_{maf}.txt'),
        sample_imputed_new_name_HC = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bcf.log'
    shell:
        '''
        echo '{wildcards.sample} {wildcards.sample}_HC_imputed' > {output.new_name_file_HC}

        bcftools reheader \
        -s {output.new_name_file_HC} \
        {input.phased_maf_info} \
        -o {output.sample_imputed_new_name_HC} 2> {log}

        bcftools index -f {output.sample_imputed_new_name_HC}
        '''


rule prepare_merge_imputed_genotyped_sample_PH: 
    """
    Prepare file with VCF files to merge
    """
    input:
        sample_imputed_new_name = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_{coverage_val}x_INFO_{info}_MAF_{maf}.bcf', coverage_val=COVERAGE_VAL, allow_missing=True),
        sample_imputed_new_name_HC = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bcf'
    output:
        vcf_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample_list.{sample}_{chrom}_INFO_{info}_MAF_{maf}.txt'
    shell:
        '''
        ls -v {input.sample_imputed_new_name} {input.sample_imputed_new_name_HC} >> {output.vcf_list}
        '''    

rule merge_renamed_imputed_genotyped_sample_PH:
    """
    Merge all downsampled and HC version of a target sample into the same VCF
    """
    input:
        vcf_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_sample_list.{sample}_{chrom}_INFO_{info}_MAF_{maf}.txt'
    output:
        merged_imputed_new_name = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bcf'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bcf.log'
    threads: 10
    shell:
        '''
        bcftools merge \
        --file-list {input.vcf_list} \
        -Ob -o {output.merged_imputed_new_name} \
        --threads {threads} 2> {log}

        bcftools index -f {output.merged_imputed_new_name}
        '''

rule imputed_to_plink:
    """
    Take imputed merged VCF and create plink format
    """
    input:
        merged_imputed_new_name = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bcf'
    output:
        expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.{doc}', doc=DOCS, allow_missing=True)
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}'
    shell:
        '''
        plink \
        --bcf {input.merged_imputed_new_name} \
        --make-bed \
        --out {params.prefix} \
        --double-id \
        --chr-set 38
        '''

rule imputed_prepare_correct_format:
    """
    Fix bim and fam file columns
    """
    input:
        imputed_bim = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bim',
        imputed_bed = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.bed',
        imputed_fam = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}.fam'
    output:
        imputed_bim_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr.bim',
        imputed_fam_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr.fam',
        imputed_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr.bed'
    shell:
        '''
        awk 'BEGIN{{OFS="\t"}}$1="chr"$1' {input.imputed_bim} | awk 'BEGIN{{OFS="\t"}}$2=$1"_"$4' > {output.imputed_bim_corr}

        awk '{{$6=2 ; print ; }}' {input.imputed_fam} > {output.imputed_fam_corr}

        scp {input.imputed_bed} {output.imputed_bed_corr}
        '''

rule create_sample_list_imputed:
    """
    Create list of imputed files to merge with reference and PH validation files from above
    """
    input:
        imputed_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr.bed'
    output:
        HC_prefix_list_imputed = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr.txt',
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr'
    shell:
        '''
        echo "{params.prefix}" >> {output.HC_prefix_list_imputed}
        '''

rule merge_downsampled_reference_imputed:
    """
    Merge reference and PH validation files with imputed files
    """
    input:
        #ref_bed_corr = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom_con}_validation_no_missnp_all_{canid_subset}.bed',
        merged = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_validation_no_missnp_all_{canid_subset}.bed', 
        HC_prefix_list_imputed = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/sample_VCFs/phased_renamed_merged.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr.txt',
    output:
        merged = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr_{canid_subset}.{doc}', doc=DOCS, allow_missing=True) 
    params:
        prefix_ref = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ancient_modern/merged_ph_called_modern-{sample}_{chrom}_validation_no_missnp_all_{canid_subset}',
        prefix_output = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr_{canid_subset}'
    shell:
        '''
        plink \
        --bfile {params.prefix_ref} \
        --make-bed \
        --merge-list {input.HC_prefix_list_imputed} \
        --allow-no-sex \
        --chr-set 38 \
        --out {params.prefix_output} 
        '''

rule merge_allchrom_downsampled_reference_imputed_list:
    """ 
    Prepare list to merge chromosomes 
    """
    input:
        bed = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr_{canid_subset}.bed', chrom=CHROM, allow_missing=True),
       # bim_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr_{canid_subset}.bim',
       # fam_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr_{canid_subset}.fam'
    params:
        prefix_in = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_{chrom}_INFO_{info}_MAF_{maf}_corr_{canid_subset}', chrom=CHROM, allow_missing=True)
    output:
        chr_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_list.txt'
    shell:
        '''
        echo -e {params.prefix_in} | tr " " "\n" | sed '1d' >> {output.chr_list}
        '''  

rule merge_allchrom_downsampled_reference_imputed:
    """
    Merge all chromosomes
    """
    input:
        chr_list = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_list.txt'
    output:
        merged = expand('{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.{doc}', doc=DOCS, allow_missing=True) 
    params:
        prefix_chr1 = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_chr1_INFO_{info}_MAF_{maf}_corr_{canid_subset}',
        prefix_output = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}'
    shell:
        '''
        plink \
        --bfile {params.prefix_chr1} \
        --merge-list {input.chr_list} \
        --make-bed \
        --allow-no-sex \
        --chr-set 38 \
        --out {params.prefix_output} 
        '''

#### smartpca part

rule prepare_convertf_parfile_con_PH:
    """
    Prepare convertf file to convert to eigenstrat format
    """
    input:
        bed = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.bed',
        bim_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.bim',
        fam_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.fam'
    output:
        convertf_file = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_convertf_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}'
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

rule convertf_con_PH:
    """
    Convert to eigenstrat format
    """
    input:
        convertf_file = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_convertf_parfile'
    output:
        geno = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.eigenstratgeno',
        ind = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.ind',
        snp = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.snp'
    log:
        '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.log'
    shell:
        '''
        convertf \
        -p {input.convertf_file} 2> {log}
        '''

rule set_high_low_coverage: 
    """
    Set high and low coverage categories for projection
    """
    input:
        ind = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.ind',
    output:
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_mod.ind',
        poplist = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_poplist.txt'
    shell:
        '''
        awk 'BEGIN {{IFS = OFS = "\t"}} {{
        if ($1 ~ "{wildcards.sample}")
            $3 = "low";
        else 
            $3 = "high";
        print
        }}' {input.ind} > {output.ind_mod}

        echo "high" > {output.poplist}
        '''

rule prepare_smartpca_parfile_con_PH:
    """
    Prepare input file for smartpca
    """
    input:
        geno = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.eigenstratgeno',
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_mod.ind',
        snp = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.snp',
        poplist = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_poplist.txt'
    output:
        smartpca_parfile = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_smartpca_parfile'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}'
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

rule smartpca_con_PH:
    """
    Run smartpca
    """
    input:
        smartpca_parfile = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_smartpca_parfile'
    output:
        smartpca_log = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_smartpca_parfile.log',
        smartpca_eigenvalue = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_eigenval_output',
        smartpca_eigenvector = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_eigenvec_output'
    shell:
        '''
        smartpca \
        -p {input.smartpca_parfile} > {output.smartpca_log}
        '''

rule plot_pca_concordance_ph:
    """
    Plot downsampled/HC imputed and non imputed samples PCA
    """
    input:
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_mod.ind',
        smartpca_eigenvalue = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_eigenval_output',
        smartpca_eigenvector = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_eigenvec_output',
        ref_metadata = config['reference_panel_metadata'],
        concordance_metadata = config['bam_targets']
    output:
        smartpca_plot = '{path}/output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}.png'
    params:
        sample = '{sample}',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"]
    script:
        "../scripts/smartpca_ph.R"


#########################################################
#####      Estimate Sum of weighted PC distance     #####
#########################################################

rule pca_accuracy_downsampled_imputed:
    """
    Estimate pca distances for downsampled imputed against the validation HC 
    """
    input:
        ind_mod = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_mod.ind',
        smartpca_eigenvalue = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_eigenval_output',
        smartpca_eigenvector = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/smartpca/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_eigenvec_output',
    output:
        accuracy_plot_HC = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ratio_pseudohaploid_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_accuracy_HC.png',
        accuracy_ratio = '{path}/output/GLIMPSE_concordance/PCA_concordance_PH_HC_genotyped/ratio_pseudohaploid_imputed/merged_ph_called_modern_imputed.{sample}_allchrom_INFO_{info}_MAF_{maf}_corr_{canid_subset}_accuracy_HC.tsv'
    params:
        name = '{sample}',
        cov_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Coverage"],
        info_sample = lambda wildcards: samples_df.loc[wildcards.sample, "Info"],
    script:
        "../scripts/smartpca_ph_distances.R"

