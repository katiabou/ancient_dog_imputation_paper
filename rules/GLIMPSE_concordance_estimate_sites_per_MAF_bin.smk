######################################################
#  Estimate sites per MAF bin for each target sample #
######################################################

rule MAF_bins_sites_ref_pan:
    """
    Extract sites per MAF bin for the reference panel
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf'
    output:
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_0_001.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_001_0_002.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_002_0_005.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_005_0_01.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_01_0_05.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_05_0_1.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_1_0_2.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_2_0_5.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz',
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
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
        '''


rule MAF_bins_sites_ref_pan_concordance_imputed:
    """
    Extract MAF bin sites and INFO filtered sites from imputed VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5 = '{path}/output/GLIMPSE_concordance/reference_panel_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz',
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
    output:
        imputed_maf_info_0_0_001 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz',
        imputed_maf_info_0_001_0_002 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz',
        imputed_maf_info_0_002_0_005 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz',
        imputed_maf_info_0_005_0_01 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz',
        imputed_maf_info_0_01_0_05 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz',
        imputed_maf_info_0_05_0_1 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz',
        imputed_maf_info_0_1_0_2 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz',
        imputed_maf_info_0_2_0_5 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz'
    params:
        info='{info_cutoff}'
    #conda:
    #    '../envs/environment.yaml'
    threads: 8
    shell:
        '''
        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_0_001} 

        bcftools index -f {output.imputed_maf_info_0_0_001} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_001_0_002} 

        bcftools index -f {output.imputed_maf_info_0_001_0_002} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_002_0_005} 

        bcftools index -f {output.imputed_maf_info_0_002_0_005} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_005_0_01} 

        bcftools index -f {output.imputed_maf_info_0_005_0_01} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_01_0_05} 

        bcftools index -f {output.imputed_maf_info_0_01_0_05} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_05_0_1} 

        bcftools index -f {output.imputed_maf_info_0_05_0_1} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_1_0_2} 

        bcftools index -f {output.imputed_maf_info_0_1_0_2} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_2_0_5} 

        bcftools index -f {output.imputed_maf_info_0_2_0_5} 
        '''
    

rule count_sites:
    """
    Count number of sites per bin
    """
    input:
        imputed_maf_info_0_0_001 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz',
        imputed_maf_info_0_001_0_002 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz',
        imputed_maf_info_0_002_0_005 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz',
        imputed_maf_info_0_005_0_01 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz',
        imputed_maf_info_0_01_0_05 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz',
        imputed_maf_info_0_05_0_1 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz',
        imputed_maf_info_0_1_0_2 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz',
        imputed_maf_info_0_2_0_5 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz'
    output:
        sample_file = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}.txt'
    #conda:
    #    '../envs/environment.yaml'
    threads: 8
    shell:
        '''
        SITES=$(bcftools view -H {input.imputed_maf_info_0_0_001} | wc -l)
        echo {input.imputed_maf_info_0_0_001} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES > {output.sample_file}
        
        SITES=$(bcftools view -H {input.imputed_maf_info_0_001_0_002} | wc -l)
        echo {input.imputed_maf_info_0_001_0_002} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
       
        SITES=$(bcftools view -H {input.imputed_maf_info_0_002_0_005} | wc -l)
        echo {input.imputed_maf_info_0_002_0_005} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_005_0_01} | wc -l)
        echo {input.imputed_maf_info_0_005_0_01} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_01_0_05} | wc -l)
        echo {input.imputed_maf_info_0_01_0_05} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_05_0_1} | wc -l)
        echo {input.imputed_maf_info_0_05_0_1} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_1_0_2} | wc -l)
        echo {input.imputed_maf_info_0_1_0_2} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_2_0_5} | wc -l)
        echo {input.imputed_maf_info_0_2_0_5} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
        '''

rule merge_counted_sites:
    """
    Count number of sites per bin
    """
    input:
        sample_file = expand('{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}.txt', info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        sample_file_merge = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_all.txt'
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
        cat {input.sample_file} >> {output.sample_file_merge}
        '''

###################################
#### DOGS ONLY REFERENCE PANEL ####
###################################


rule MAF_bins_sites_ref_pan_dogs_only:
    """
    Extract sites per MAF bin for the reference panel
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs/{chrom_con}_ref_panel_filltags_filter.bcf'
    output:
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_0_001 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_0_001.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_001_0_002 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_001_0_002.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_002_0_005 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_002_0_005.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_005_0_01 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_005_0_01.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_01_0_05 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_01_0_05.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_05_0_1 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_05_0_1.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_1_0_2 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_1_0_2.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_vcf_0_2_0_5 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_2_0_5.vcf.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz',
    #conda:
    #    '../envs/environment.yaml'
    shell:
        '''
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
        '''


rule MAF_bins_sites_ref_pan_dogs_only_concordance_imputed:
    """
    Extract MAF bin sites and INFO filtered sites from imputed VCF
    """
    input:
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_0_001.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_001_0_002.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_002_0_005.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_005_0_01.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_01_0_05.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_05_0_1.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_1_0_2.tsv.gz',
        ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5 = '{path}/output/GLIMPSE_concordance_only_dogs/reference_panel_only_dogs_MAF_bins/{chrom_con}_ref_panel_filltags_filter_MAF_0_2_0_5.tsv.gz',
        ligated_bcf = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
    output:
        imputed_maf_info_0_0_001 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz',
        imputed_maf_info_0_001_0_002 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz',
        imputed_maf_info_0_002_0_005 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz',
        imputed_maf_info_0_005_0_01 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz',
        imputed_maf_info_0_01_0_05 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz',
        imputed_maf_info_0_05_0_1 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz',
        imputed_maf_info_0_1_0_2 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz',
        imputed_maf_info_0_2_0_5 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz'
    params:
        info='{info_cutoff}'
    #conda:
    #    '../envs/environment.yaml'
    threads: 8
    shell:
        '''
        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_0_001} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_0_001} 

        bcftools index -f {output.imputed_maf_info_0_0_001} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_001_0_002} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_001_0_002} 

        bcftools index -f {output.imputed_maf_info_0_001_0_002} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_002_0_005} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_002_0_005} 

        bcftools index -f {output.imputed_maf_info_0_002_0_005} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_005_0_01} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_005_0_01} 

        bcftools index -f {output.imputed_maf_info_0_005_0_01} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_01_0_05} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_01_0_05} 

        bcftools index -f {output.imputed_maf_info_0_01_0_05} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_05_0_1} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_05_0_1} 

        bcftools index -f {output.imputed_maf_info_0_05_0_1} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_1_0_2} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_1_0_2} 

        bcftools index -f {output.imputed_maf_info_0_1_0_2} 

        bcftools view {input.ligated_bcf} \
        --regions-file {input.ref_concordance_sample_excl_filltags_filter_maf_tsv_0_2_0_5} \
        --include 'INFO/INFO >= {params.info}' \
        --threads {threads} \
        -Oz -o {output.imputed_maf_info_0_2_0_5} 

        bcftools index -f {output.imputed_maf_info_0_2_0_5} 
        '''


rule count_sites_dogs_only:
    """
    Count number of sites per bin
    """
    input:
        imputed_maf_info_0_0_001 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_0_001.vcf.gz',
        imputed_maf_info_0_001_0_002 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_001_0_002.vcf.gz',
        imputed_maf_info_0_002_0_005 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_002_0_005.vcf.gz',
        imputed_maf_info_0_005_0_01 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_005_0_01.vcf.gz',
        imputed_maf_info_0_01_0_05 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_01_0_05.vcf.gz',
        imputed_maf_info_0_05_0_1 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_05_0_1.vcf.gz',
        imputed_maf_info_0_1_0_2 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_1_0_2.vcf.gz',
        imputed_maf_info_0_2_0_5 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}_MAF_0_2_0_5.vcf.gz'
    output:
        sample_file = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}.txt'
    #conda:
    #    '../envs/environment.yaml'
    threads: 8
    shell:
        '''
        SITES=$(bcftools view -H {input.imputed_maf_info_0_0_001} | wc -l)
        echo {input.imputed_maf_info_0_0_001} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES > {output.sample_file}
        
        SITES=$(bcftools view -H {input.imputed_maf_info_0_001_0_002} | wc -l)
        echo {input.imputed_maf_info_0_001_0_002} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
       
        SITES=$(bcftools view -H {input.imputed_maf_info_0_002_0_005} | wc -l)
        echo {input.imputed_maf_info_0_002_0_005} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_005_0_01} | wc -l)
        echo {input.imputed_maf_info_0_005_0_01} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_01_0_05} | wc -l)
        echo {input.imputed_maf_info_0_01_0_05} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_05_0_1} | wc -l)
        echo {input.imputed_maf_info_0_05_0_1} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_1_0_2} | wc -l)
        echo {input.imputed_maf_info_0_1_0_2} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}

        SITES=$(bcftools view -H {input.imputed_maf_info_0_2_0_5} | wc -l)
        echo {input.imputed_maf_info_0_2_0_5} {wildcards.sample} {wildcards.chrom_con} {wildcards.coverage_val} {wildcards.info_cutoff} $SITES >> {output.sample_file}
        '''

rule merge_counted_sites_dogs_only:
    """
    Count number of sites per bin
    """
    input:
        sample_file = expand('{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_{info_cutoff}.txt', info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        sample_file_merge = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.{sample}_{chrom_con}_{coverage_val}x_INFO_all.txt'
    shell:
        '''
        cat {input.sample_file} >> {output.sample_file_merge}
        '''


rule plot_sites_ref_panels:
    """
    Count number of sites per bin
    """
    input:
        sample_file_merge_1 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.NGDG_{chrom_con}_0.5x_INFO_all.txt',
        sample_file_merge_2 = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.CGG32_{chrom_con}_1x_INFO_all.txt',
        sample_file_merge_3 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.NGDG_{chrom_con}_0.5x_INFO_all.txt',
        sample_file_merge_4 = '{path}/output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.CGG32_{chrom_con}_1x_INFO_all.txt'
    output:
        plot_sites_ref_panels = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_MAF_bins_reference_panel/merged_ligated.NGDG_CGG32_{chrom_con}_INFO_all.png'
    #conda:
    #    '../envs/r4.3.1.yaml'
    script:
        "../scripts/reference_panel_MAF_sites_comparison.R"
   