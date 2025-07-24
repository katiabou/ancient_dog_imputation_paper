#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

###################################################################################################
#  Exctract common CNV windows from dog10k data and remove ROH desert windows with high frequency #
#  Then run GO analysis excluding these from the ROH deserts, also test effect of DLA region      #
###################################################################################################

rule get_cnv_files:
    """
    Download CNV files from dog10k
    """
    output:
        windows_depth='data/CNV_windows/dog10k_no_outliers_up_to_chrx.tsv',
        windows='data/CNV_windows/window_coords.txt',
        dog10k_names='data/CNV_windows/dog10k_good_names.txt',
        dog10k_meta='data/CNV_windows/dog10k_metadata.txt',
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/CNV_windows/dog10k_no_outliers_up_to_chrx.tsv > {output.windows_depth}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/CNV_windows/window_coords.txt > {output.windows}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/CNV_windows/dog10k_good_names.txt > {output.dog10k_names}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/CNV_windows/dog10k_metadata.txt > {output.dog10k_meta}
        """

rule per_chr_windows:
    """
    Extract read depth window file from dog10k for each chromosome (easier to process downstream)
    """
    input:
        windows_depth='data/CNV_windows/dog10k_no_outliers_up_to_chrx.tsv',
        windows='data/CNV_windows/window_coords.txt',
    params:
        chrom='{chrom}'
    output:
        windows_depth_chr='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/dog10k_no_outliers_{chrom}.tsv',
    threads: 24
    script:
        "../scripts/create_chr_subset.R"

rule prepare_depth_file:
    """
    Prepare read depth window file from dog10k for dogs and wolves 
    Good are windows which have coverage =2
    Bad are windows which have coverage <2 or >2 
    """
    input:
        windows_depth_chr='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/dog10k_no_outliers_{chrom}.tsv',
        windows='data/CNV_windows/window_coords.txt',
        dog10k_names='data/CNV_windows/dog10k_good_names.txt',
        dog10k_meta='data/CNV_windows/dog10k_metadata.txt',
    output:
        good_wolves='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/good_windows_wolves_{chrom}.bed',
        bad_wolves='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bad_windows_wolves_{chrom}.bed',
        good_dogs='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/good_windows_dogs_{chrom}.bed',
        bad_dogs='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bad_windows_dogs_{chrom}.bed',
    params:
        chrom='{chrom}'
    script:
        "../scripts/prep_cnv_window_depth.R"

rule download_liftover_chain:
    """
    Download the canFam4 to canfam3.1 liftover chain
    """
    output:
        chain="data/liftover_chain/canFam4ToCanFam3.over.chain.gz",
    shell:
        "wget --quiet -O - https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/liftOver/canFam4ToCanFam3.over.chain.gz > {output}"


rule liftover_good_windows:
    """
    Liftover good regions from canfam4 to canfam3.1
    """
    input:
        good_wolves='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/good_windows_wolves_{chrom}.bed',
        good_dogs='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/good_windows_dogs_{chrom}.bed',
        chain="data/liftover_chain/canFam4ToCanFam3.over.chain.gz",
    output:
        lifted_wolf_temp = temp('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_good_wolf_lifted_{chrom}-temp.bed'),
        unlifted_wolf= 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_good_wolf_unlifted_{chrom}.bed',
        lifted_wolf = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_good_wolf_lifted_{chrom}.bed',
        lifted_dog_temp = temp('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_good_dog_lifted_{chrom}-temp.bed'),
        unlifted_dog= 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_good_dog_unlifted_{chrom}.bed',
        lifted_dog = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_good_dog_lifted_{chrom}.bed',
    shell:
        '''
        liftOver \
        {input.good_wolves} \
        {input.chain} \
        {output.lifted_wolf_temp} \
        {output.unlifted_wolf}

        sed 's/chr//g' {output.lifted_wolf_temp} > {output.lifted_wolf}

        liftOver \
        {input.good_dogs} \
        {input.chain} \
        {output.lifted_dog_temp} \
        {output.unlifted_dog}

        sed 's/chr//g' {output.lifted_dog_temp} > {output.lifted_dog}
        '''

rule liftover_bad_windows:
    """
    Liftover bad regions from canfam4 to canfam3.1
    """
    input:
        bad_wolves='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bad_windows_wolves_{chrom}.bed',
        bad_dogs='output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bad_windows_dogs_{chrom}.bed',
        chain="data/liftover_chain/canFam4ToCanFam3.over.chain.gz",
    output:
        lifted_wolf_temp = temp('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_{chrom}-temp.bed'),
        unlifted_wolf= 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_unlifted_{chrom}.bed',
        lifted_wolf = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_{chrom}.bed',
        lifted_dog_temp = temp('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_{chrom}-temp.bed'),
        unlifted_dog= 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_unlifted_{chrom}.bed',
        lifted_dog = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_{chrom}.bed',
    shell:
        '''
        liftOver \
        {input.bad_wolves} \
        {input.chain} \
        {output.lifted_wolf_temp} \
        {output.unlifted_wolf}

        sed 's/chr//g' {output.lifted_wolf_temp} > {output.lifted_wolf}

        liftOver \
        {input.bad_dogs} \
        {input.chain} \
        {output.lifted_dog_temp} \
        {output.unlifted_dog}

        sed 's/chr//g' {output.lifted_dog_temp} > {output.lifted_dog}
        '''

rule prep_chrom_bed:
    """
    Prepare format for only autosomes
    """
    input:
        lifted_dog = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_{chrom}.bed',
        lifted_wolf = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_{chrom}.bed',
    output:
        lifted_dog_mod = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_{chrom}-mod.bed',
        lifted_wolf_mod = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_{chrom}-mod.bed',
    shell:
        '''
        awk 'BEGIN {{OFS="\t"}} {{print "chr"$1,$2,$3}}' {input.lifted_dog} | awk '$1 == "{wildcards.chrom}"' > {output.lifted_dog_mod}
        awk 'BEGIN {{OFS="\t"}} {{print "chr"$1,$2,$3}}' {input.lifted_wolf} | awk '$1 == "{wildcards.chrom}"' > {output.lifted_wolf_mod}
        '''

rule merge_chrom_bed:
    """
    Merge chromosomes
    """
    input:
        lifted_dog = expand('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_{chrom}-mod.bed', chrom=CHROM, allow_missing=True),
        lifted_wolf = expand('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_{chrom}-mod.bed', chrom=CHROM, allow_missing=True),
    output:
        allchrom_dogs = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_allchrom-mod.bed',
        allchrom_wolves = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_allchrom-mod.bed',
    shell:
        '''
        cat {input.lifted_dog} > {output.allchrom_dogs}
        cat {input.lifted_wolf} > {output.allchrom_wolves}
        '''

rule merge_chrom_genome_windows:
    """
    Merge chromosomes for 500kb windows used in ROH analysis
    """
    input:
        ROH_500_kb_windows = expand('output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/CanFam31_{chrom}_500_kb_windows.txt', chrom=CHROM, allow_missing=True),
    output:
        allchrom_genome_windows_temp = temp('output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/CanFam31_allchrom_500_kb_windows-temp.bed'),
        allchrom_genome_windows = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/CanFam31_allchrom_500_kb_windows.bed',
    shell:
        '''
        cat {input.ROH_500_kb_windows} > {output.allchrom_genome_windows_temp}
        awk -F '[:-]' '{{print $1 "\t" $2 "\t" $3}}' {output.allchrom_genome_windows_temp} > {output.allchrom_genome_windows}
        '''

rule count_occurences:
    """
    Estimate frequency of each "bad" window in dogs and wolves
    """
    input:
        allchrom_dogs = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_dog_lifted_allchrom-mod.bed',
        allchrom_wolves = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/liftover/canfam31_bad_wolf_lifted_allchrom-mod.bed',
        allchrom_genome_windows = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/CanFam31_allchrom_500_kb_windows.bed',
        dog10k_names='data/CNV_windows/dog10k_good_names.txt',
        dog10k_meta='data/CNV_windows/dog10k_metadata.txt',
    output:
        allchrom_dogs_count = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/canfam31_bad_dog_lifted_allchrom-count.bed',
        allchrom_wolves_count = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/canfam31_bad_wolf_lifted_allchrom-count.bed',
        allchrom_genome_windows_mod = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/CanFam31_allchrom_500_kb_windows-mod.bed',
    script:
        "../scripts/count_bad_window.R"


rule bedtools_intersect_dogs:
    """
    Get intersection of the "bad" regions with the whole genome 500kb windows
    """
    input:
        allchrom_dogs_count = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/canfam31_bad_dog_lifted_allchrom-count.bed',
        allchrom_genome_windows_mod = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/CanFam31_allchrom_500_kb_windows-mod.bed',
    output:
        all_chrom_genome_windows_intersect_dog = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/genome_bad_windows_dog_allchrom.bed',
    shell:
        '''
        bedtools intersect \
        -a {input.allchrom_genome_windows_mod} \
        -b {input.allchrom_dogs_count} \
        -wo > {output.all_chrom_genome_windows_intersect_dog}
        '''

rule bedtools_intersect_wolves:
    """
    Get intersection of the "bad" regions with the whole genome 500Kb windows
    """
    input:
        allchrom_wolves_count = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/canfam31_bad_wolf_lifted_allchrom-count.bed',
        allchrom_genome_windows_mod = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/CanFam31_allchrom_500_kb_windows-mod.bed',
    output:
        all_chrom_genome_windows_intersect_wolf = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/genome_bad_windows_wolf_allchrom.bed',
    shell:
        '''
        bedtools intersect \
        -a {input.allchrom_genome_windows_mod} \
        -b {input.allchrom_wolves_count} \
        -wo > {output.all_chrom_genome_windows_intersect_wolf}
        '''

rule go_analysis_dogs:
    """
    Look into ROH deserts and GO analysis, removing windows with high frequency of CNVs
    Also test results removing DLA region
    """
    input:
        imputed_modern_windows_bed_deserts = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts.bed',
        all_chrom_genome_windows_intersect_dog = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/genome_bad_windows_dog_allchrom.bed',
        phased_roh_sum_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs.hom.summary',
        phased_roh_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.chr1_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs.hom.indiv',
        modern_roh_allchrom_sum = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogs.hom.summary',
        modern_roh_ind = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_chr1_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_dogs.hom.indiv',
        windows_cov_subset_allchrom = 'output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/dogs_allchrom_windows_cov_500kb.txt',
        DLA_regions = 'data/dla_regions.tsv'
    output:
        cnv_freq_density = 'output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/cnv_window_freq_density_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs.png',
        ROH_prev_cnv_freq = 'output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/ROH_prevelance_cnv_freq_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs.png',
        top_go_terms = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_GO_terms.txt',
        candidate_genes = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_candidate_genes.txt',
        top_go_terms_noDLA = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_GO_terms_noDLA.txt',
        candidate_genes_nDLA = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_dogs_deserts_candidate_genes_noDLA.txt',
    script:
        "../scripts/GO_analysis_dogs.R"


rule go_analysis_wolves:
    """
    Look into ROH deserts and GO analysis, removing windows with high frequency of CNVs
    Also test results removing DLA region
    """
    input:
        imputed_modern_windows_bed_deserts = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts.bed',
        all_chrom_genome_windows_intersect_wolf = 'output/GLIMPSE_imputation/ROH_islands_deserts/per_chrom_window_depth/bedtools_intersect/genome_bad_windows_wolf_allchrom.bed',
        phased_roh_sum_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.allchrom_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves.hom.summary',
        phased_roh_allchrom = 'output/GLIMPSE_imputation/ROH_phased/merged_phased_annotated.chr1_MAF_{maf_cutoff}_recalibrated_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves.hom.indiv',
        modern_roh_allchrom_sum = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_allchrom_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_wolves.hom.summary',
        modern_roh_ind = 'output/GLIMPSE_imputation/ROH_ref_panel/ref-panel_chr1_sample-snp_filltags_filter_MAF_{maf_cutoff}_all_sites_hom_win_het_{hom_win_het}_wolves.hom.indiv',
        windows_cov_subset_allchrom = 'output/GLIMPSE_imputation/ROH_islands_deserts/window_depth/dogs_allchrom_windows_cov_500kb.txt',
        DLA_regions = 'data/dla_regions.tsv'
    output:
        cnv_freq_density = 'output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/cnv_window_freq_density_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves.png',
        ROH_prev_cnv_freq = 'output/GLIMPSE_imputation/ROH_islands_deserts/CNV_windows/ROH_prevelance_cnv_freq_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves.png',
        top_go_terms = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_GO_terms.txt',
        candidate_genes = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_candidate_genes.txt',
        top_go_terms_noDLA = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_GO_terms_noDLA.txt',
        candidate_genes_nDLA = 'output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_{hom_win_het}_wolves_deserts_candidate_genes_noDLA.txt',
    script:
        "../scripts/GO_analysis_wolves.R"

