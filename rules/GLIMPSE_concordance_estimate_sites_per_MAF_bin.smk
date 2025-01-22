######################################################
#  Estimate sites per MAF bin for each target sample #
######################################################

global CHROM, INFO_CUTOFF


rule MAF_bins_sites_ref_pan:
    """
    Extract sites per MAF bin for the reference panel
    """
    input:
        #ref_concordance_sample_excl_filltags_filter = 'output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.phased.bcf'
        ref_concordance_sample_excl_filltags_filter_allchrom="output/GLIMPSE_concordance/reference_panel/allchrom_ref_panel_filltags_filter.phased.bcf",
    output:
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_0_001.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_001_0_002.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_002_0_005.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_005_0_01.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_01_0_05.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_05_0_1.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_1_0_2.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_2_0_5.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz",
    shell:
        """
        bcftools view -i 'MAF > 0 & MAF <=0.001' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001}
        bcftools view -i 'MAF > 0.001 & MAF <=0.002' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002}
        bcftools view -i 'MAF > 0.002 & MAF <=0.005' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005}
        bcftools view -i 'MAF > 0.005 & MAF <=0.01' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01}
        bcftools view -i 'MAF > 0.01 & MAF <=0.05' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05}
        bcftools view -i 'MAF > 0.05 & MAF <=0.1' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1}
        bcftools view -i 'MAF > 0.1 & MAF <=0.2' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2}
        bcftools view -i 'MAF > 0.2 & MAF <=0.5' {input.ref_concordance_sample_excl_filltags_filter_allchrom} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5}

        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5}

        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5}
        """


rule prepare_merged_chr_list_ligated:
    """ 
    Prepare list to merge chromosomes for imputed unfiltered
    """
    input:
        ligated_bcf=expand(
            "output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_ligated/chr_list.{sample}_{coverage_val}x.txt",
    shell:
        """
        ls -v {input.ligated_bcf} >> {output.chr_list}
        """


rule merge_chr_concordance_ligated:
    """
    Merge all chrom imputed unfiltered
    """
    input:
        chr_list="output/GLIMPSE_concordance/GLIMPSE_ligated/chr_list.{sample}_{coverage_val}x.txt",
    output:
        imputed_allchrom="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.bcf",
        imputed_allchrom_csi="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.bcf.csi",
    log:
        "output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.log",
    threads: 8
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Ob -o {output.imputed_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.imputed_allchrom}
        """


rule MAF_bins_sites_ref_pan_concordance_imputed:
    """
    Extract MAF bin sites and INFO filtered sites from imputed VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5="output/GLIMPSE_concordance/reference_panel_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz",
        #ligated_bcf = 'output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf',
        imputed_allchrom="output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.bcf",
    output:
        imputed_maf_info_0_0_001="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz",
        imputed_maf_info_0_001_0_002="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz",
        imputed_maf_info_0_002_0_005="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz",
        imputed_maf_info_0_005_0_01="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz",
        imputed_maf_info_0_01_0_05="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz",
        imputed_maf_info_0_05_0_1="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz",
        imputed_maf_info_0_1_0_2="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz",
        imputed_maf_info_0_2_0_5="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz",
    params:
        info="{info_cutoff}",
    threads: 8
    shell:
        """
        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_0_001} 

        bcftools index -f {output.imputed_maf_info_0_0_001} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_001_0_002} 

        bcftools index -f {output.imputed_maf_info_0_001_0_002} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_002_0_005} 

        bcftools index -f {output.imputed_maf_info_0_002_0_005} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_005_0_01} 

        bcftools index -f {output.imputed_maf_info_0_005_0_01} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_01_0_05} 

        bcftools index -f {output.imputed_maf_info_0_01_0_05} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_05_0_1} 

        bcftools index -f {output.imputed_maf_info_0_05_0_1} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_1_0_2} 

        bcftools index -f {output.imputed_maf_info_0_1_0_2} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_2_0_5} 

        bcftools index -f {output.imputed_maf_info_0_2_0_5} 
        """


rule count_sites:
    """
    Count number of sites per bin
    """
    input:
        imputed_maf_info_0_0_001="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz",
        imputed_maf_info_0_001_0_002="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz",
        imputed_maf_info_0_002_0_005="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz",
        imputed_maf_info_0_005_0_01="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz",
        imputed_maf_info_0_01_0_05="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz",
        imputed_maf_info_0_05_0_1="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz",
        imputed_maf_info_0_1_0_2="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz",
        imputed_maf_info_0_2_0_5="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz",
    output:
        sample_file="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}.txt",
    threads: 8
    shell:
        """
        SITES=$(bcftools view -H {input.imputed_maf_info_0_0_001} | wc -l)
        echo {input.imputed_maf_info_0_0_001} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES > {output.sample_file}
        
        SITES=$(bcftools view -H {input.imputed_maf_info_0_001_0_002} | wc -l)
        echo {input.imputed_maf_info_0_001_0_002} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
       
        SITES=$(bcftools view -H {input.imputed_maf_info_0_002_0_005} | wc -l)
        echo {input.imputed_maf_info_0_002_0_005} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_005_0_01} | wc -l)
        echo {input.imputed_maf_info_0_005_0_01} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_01_0_05} | wc -l)
        echo {input.imputed_maf_info_0_01_0_05} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_05_0_1} | wc -l)
        echo {input.imputed_maf_info_0_05_0_1} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_1_0_2} | wc -l)
        echo {input.imputed_maf_info_0_1_0_2} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_2_0_5} | wc -l)
        echo {input.imputed_maf_info_0_2_0_5} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
        """


rule merge_counted_sites:
    """
    Count number of sites per bin
    """
    input:
        sample_file=expand(
            "output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}.txt",
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        sample_file_merge="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_all.txt",
    shell:
        """
        cat {input.sample_file} >> {output.sample_file_merge}
        """


###################################
#### DOGS ONLY REFERENCE PANEL ####
###################################


rule MAF_bins_sites_ref_pan_dogs_only:
    """
    Extract sites per MAF bin for the reference panel
    """
    input:
        ref_concordance_sample_excl_filltags_filter="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/allchrom_ref_panel_filltags_filter.phased.bcf",
    output:
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_0_001.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_001_0_002.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_002_0_005.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_005_0_01.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_01_0_05.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_05_0_1.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_1_0_2.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_2_0_5.vcf.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz",
    shell:
        """
        bcftools view -i 'MAF > 0 & MAF <=0.001' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001}
        bcftools view -i 'MAF > 0.001 & MAF <=0.002' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002}
        bcftools view -i 'MAF > 0.002 & MAF <=0.005' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005}
        bcftools view -i 'MAF > 0.005 & MAF <=0.01' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01}
        bcftools view -i 'MAF > 0.01 & MAF <=0.05' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05}
        bcftools view -i 'MAF > 0.05 & MAF <=0.1' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1}
        bcftools view -i 'MAF > 0.1 & MAF <=0.2' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2}
        bcftools view -i 'MAF > 0.2 & MAF <=0.5' {input.ref_concordance_sample_excl_filltags_filter} -Oz -o {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5}

        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2}
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5} | \
        bgzip -c > {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5}

        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2}
        tabix -s1 -b2 -e2 {output.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5}
        """


rule prepare_merged_chr_list_ligated_only_dogs:
    """ 
    Prepare list to merge chromosomes for imputed unfiltered
    """
    input:
        ligated_bcf=expand(
            "output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf",
            chrom=CHROM,
            allow_missing=True,
        ),
    output:
        chr_list="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/chr_list.{sample}_{coverage_val}x.txt",
    shell:
        """
        ls -v {input.ligated_bcf} >> {output.chr_list}
        """


rule merge_chr_concordance_ligated_only_dogs:
    """
    Merge all chrom imputed unfiltered
    """
    input:
        chr_list="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/chr_list.{sample}_{coverage_val}x.txt",
    output:
        imputed_allchrom="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.bcf",
        imputed_allchrom_csi="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.bcf.csi",
    log:
        "output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.log",
    threads: 8
    shell:
        """
        bcftools concat \
        --file-list {input.chr_list} \
        -Ob -o {output.imputed_allchrom} \
        --threads {threads} 2> {log}

        bcftools index -f {output.imputed_allchrom}
        """


rule MAF_bins_sites_ref_pan_dogs_only_concordance_imputed:
    """
    Extract MAF bin sites and INFO filtered sites from imputed VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz",
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5="output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/allchrom_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz",
        #ligated_bcf = 'output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_{chrom}_{coverage_val}x.bcf',
        imputed_allchrom="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_allchrom_{coverage_val}x.bcf",
    output:
        imputed_maf_info_0_0_001="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz",
        imputed_maf_info_0_001_0_002="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz",
        imputed_maf_info_0_002_0_005="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz",
        imputed_maf_info_0_005_0_01="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz",
        imputed_maf_info_0_01_0_05="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz",
        imputed_maf_info_0_05_0_1="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz",
        imputed_maf_info_0_1_0_2="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz",
        imputed_maf_info_0_2_0_5="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz",
    params:
        info="{info_cutoff}",
    threads: 8
    shell:
        """
        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_0_001} 

        bcftools index -f {output.imputed_maf_info_0_0_001} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_001_0_002} 

        bcftools index -f {output.imputed_maf_info_0_001_0_002} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_002_0_005} 

        bcftools index -f {output.imputed_maf_info_0_002_0_005} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_005_0_01} 

        bcftools index -f {output.imputed_maf_info_0_005_0_01} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_01_0_05} 

        bcftools index -f {output.imputed_maf_info_0_01_0_05} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_05_0_1} 

        bcftools index -f {output.imputed_maf_info_0_05_0_1} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_1_0_2} 

        bcftools index -f {output.imputed_maf_info_0_1_0_2} 

        bcftools view {input.imputed_allchrom} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_2_0_5} 

        bcftools index -f {output.imputed_maf_info_0_2_0_5} 
        """


rule count_sites_dogs_only:
    """
    Count number of sites per bin
    """
    input:
        imputed_maf_info_0_0_001="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz",
        imputed_maf_info_0_001_0_002="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz",
        imputed_maf_info_0_002_0_005="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz",
        imputed_maf_info_0_005_0_01="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz",
        imputed_maf_info_0_01_0_05="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz",
        imputed_maf_info_0_05_0_1="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz",
        imputed_maf_info_0_1_0_2="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz",
        imputed_maf_info_0_2_0_5="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz",
    output:
        sample_file="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}.txt",
    threads: 8
    shell:
        """
        SITES=$(bcftools view -H {input.imputed_maf_info_0_0_001} | wc -l)
        echo {input.imputed_maf_info_0_0_001} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES > {output.sample_file}
        
        SITES=$(bcftools view -H {input.imputed_maf_info_0_001_0_002} | wc -l)
        echo {input.imputed_maf_info_0_001_0_002} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
       
        SITES=$(bcftools view -H {input.imputed_maf_info_0_002_0_005} | wc -l)
        echo {input.imputed_maf_info_0_002_0_005} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_005_0_01} | wc -l)
        echo {input.imputed_maf_info_0_005_0_01} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_01_0_05} | wc -l)
        echo {input.imputed_maf_info_0_01_0_05} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_05_0_1} | wc -l)
        echo {input.imputed_maf_info_0_05_0_1} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_1_0_2} | wc -l)
        echo {input.imputed_maf_info_0_1_0_2} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_2_0_5} | wc -l)
        echo {input.imputed_maf_info_0_2_0_5} {wildcards.sample} 'all_autosomes' {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
        """


rule merge_counted_sites_dogs_only:
    """
    Count number of sites per bin
    """
    input:
        sample_file=expand(
            "output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_{info_cutoff}.txt",
            info_cutoff=INFO_CUTOFF,
            allow_missing=True,
        ),
    output:
        sample_file_merge="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_allchrom_{coverage_val}x_INFO_all.txt",
    shell:
        """
        cat {input.sample_file} >> {output.sample_file_merge}
        """


rule plot_sites_ref_panels:
    """
    Count number of sites per bin
    """
    input:
        sample_file_merge_1="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.NGDG_allchrom_0.5x_INFO_all.txt",
        sample_file_merge_2="output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.CGG32_allchrom_1x_INFO_all.txt",
        sample_file_merge_3="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.NGDG_allchrom_0.5x_INFO_all.txt",
        sample_file_merge_4="output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.CGG32_allchrom_1x_INFO_all.txt",
    output:
        plot_sites_ref_panels="output/GLIMPSE_concordance/plots/glimpse_concordance_MAF_bins_reference_panel/merged_ligated.NGDG_CGG32_allchrom_INFO_all.png",
    script:
        "../scripts/reference_panel_MAF_sites_comparison.R"
