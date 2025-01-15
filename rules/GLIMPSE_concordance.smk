##############################################
# Running imputation for GLIMPSE concordance #
##############################################

rule extract_chrom_ref_fast_concordance:
    """
    Extract chromosome from fasta reference file
    """
    input:
        ref_fasta = config['ref_fasta_file']
    output:
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta',
        ref_fasta_chr_fai = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta.fai'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        samtools faidx {input.ref_fasta} {wildcards.chrom_con} > {output.ref_fasta_chr}
        samtools faidx {output.ref_fasta_chr}
        '''

rule extract_chr_target_bams:
    """
    Extract chromosomes from target bam files 
    """
    input:
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta',
        target_bams = lambda wildcards: samples_df.loc[wildcards.sample, "bam_path"]
    output:
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam',
        target_bams_chr_bai = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam.bai'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        samtools view -T {input.ref_fasta_chr} \
        -bo {output.target_bams_chr} \
        {input.target_bams} \
        {wildcards.chrom_con}

        samtools index {output.target_bams_chr}
        '''

rule estimate_coverage_fraction:
    """
    Estimate coverage fraction to be used for downsampling the target bams to lower coverages (specified in Snakefile)
    """
    input:
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta',
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam'
    output:
        seed_frac = temp('{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}_{coverage_val}.txt')
    params:
        coverage_val = '{coverage_val}'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        LENGTH=$(cut -f1,2 {input.ref_fasta_chr}.fai | awk -F ' ' '{{print $2}}')
        COV=$(samtools depth -a {input.target_bams_chr}  |  awk -v var="$LENGTH" '{{sum+=$3}} END {{ print sum/var}}')
        NUM_READS=$(samtools view -c {input.target_bams_chr})
        coverage_val={params.coverage_val}
        reads_down=$(echo "$coverage_val $NUM_READS $COV" | awk '{{print ($1 * $2)/$3}}')
        frac=$( samtools idxstats {input.target_bams_chr} | cut -f3 | awk -v var="$reads_down" 'BEGIN {{total=0}} {{total += $1}} END {{frac=var/total; if (frac > 1) {{print 1}} else {{print frac}}}}' )
        seed_frac=$(echo "$frac 1" | awk '{{print $1 + $2}}')
        echo $seed_frac > {output.seed_frac}
        '''

rule downsample_target_bam:
    """
    Downsample target bams to lower coverages
    """
    input:
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta',
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam',
        seed_frac = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}_{coverage_val}.txt'
    output:
        downsampled_bam = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}_{coverage_val}x.bam',
        downsampled_bam_bai = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}_{coverage_val}x.bam.bai'
    log:
        '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}_{coverage_val}.log'
    benchmark:
        '{path}/benchmarks/target_bams/{sample}_{chrom_con}_{coverage_val}.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        s=$(cat {input.seed_frac}) 

        samtools view -T {input.ref_fasta_chr} \
        -s $s \
        -bo {output.downsampled_bam} \
        {input.target_bams_chr} 2> {log}

        samtools index {output.downsampled_bam}
        '''

rule prepare_list_validation_reference_removal:
    """
    Create a list of target samples which are also in the reference VCF but under another name (have to provide file with two columns for each sample)
    Only needed if samples overlap in target bams and VCF panel
    """
    input:
        ref_val_samples = config['ref_validation_list']
    output:
        ref_val_sample_file = '{path}/output/GLIMPSE_concordance/reference_panel/ref_val_sample.txt'
    shell:
        '''
        cat {input.ref_val_samples} | awk -F '\t' '{{print $2}}' >> {output.ref_val_sample_file}
        '''

rule prepare_ref_panel:
    """
    Removes overlapping target and reference panel samples from the reference panel
    """
    input:
        ref_sample_snp_filltags_filter = '{path}/output/reference_panel/ref-panel_{chrom_con}_sample-snp_filltags_filter.vcf.gz',
        ref_val_sample_file = '{path}/output/GLIMPSE_concordance/reference_panel/ref_val_sample.txt' 
    output:
        ref_concordance_sample_excl = temp('{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel.vcf.gz')
    log:
        '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/reference_panel/{chrom_con}_ref_panel.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view \
        {input.ref_sample_snp_filltags_filter} \
        -m 2 -M 2 \
        -S ^{input.ref_val_sample_file} \
        --trim-alt-alleles \
        --threads {threads} \
        -Ou | bcftools filter -e "type!='snp'" -Oz -o {output.ref_concordance_sample_excl} 2> {log}
  
        bcftools index -f {output.ref_concordance_sample_excl}
        '''

rule fill_tags_ref_sample_excl:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING)
    """
    input:
        ref_concordance_sample_excl = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel.vcf.gz'
    output:
        ref_concordance_sample_excl_filltags = temp('{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags.vcf.gz')
    log:
        '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/reference_panel/{chrom_con}_ref_panel_filltags.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools +fill-tags {input.ref_concordance_sample_excl} \
        -Oz -o {output.ref_concordance_sample_excl_filltags} \
        --threads {threads} \
        -- -t all,F_MISSING 2> {log}
        '''

rule filter_sites_ref_sample_excl:
    """
    Filter for missingness (F_MISSING) again, since we removed individuals
    """
    input:
        ref_concordance_sample_excl_filltags = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags.vcf.gz'
    output:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        ref_concordance_sample_excl_filltags_filter_csi = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf.csi'
    params:
        f_missing = config['F_MISSING']
    log:
        '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/reference_panel/{chrom_con}_ref_panel_filltags_filter.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view -i 'F_MISSING<{params.f_missing}' \
        {input.ref_concordance_sample_excl_filltags} \
        --threads {threads} \
        -Ob -o {output.ref_concordance_sample_excl_filltags_filter} 2> {log}
    
        bcftools index -f {output.ref_concordance_sample_excl_filltags_filter}
        '''

rule extract_var_pos_concordance:
    """
    Extract sites from reference panel to use for GLs estimation from bams afterwards
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf'
    output:
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz',
        ref_panel_sites_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.tsv.gz'
    log:
        '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/reference_panel/{chrom_con}_ref_panel_sites.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view -G -m 2 -M 2 -v snps \
        {input.ref_concordance_sample_excl_filltags_filter} \
        --threads {threads} \
        -Oz -o {output.ref_panel_sites_vcf} 2> {log}

        bcftools index -f {output.ref_panel_sites_vcf}

        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' {output.ref_panel_sites_vcf} | \
        bgzip -c > {output.ref_panel_sites_tsv}
        
        tabix -s1 -b2 -e2 {output.ref_panel_sites_tsv}
        '''

rule compute_GLs_downsampled_samples_concordance:
    """
    Compute GLs of target bams 
    """
    input:
        downsampled_bam = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}_{coverage_val}x.bam',
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz',
        ref_panel_sites_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.tsv.gz',
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta'
    output:
        GL_vcf_target_bams = '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}_{coverage_val}x.vcf.gz',
        GL_vcf_target_bams_csi = '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}_{coverage_val}x.vcf.gz.csi'
    log:
        '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}_{coverage_val}x.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/GLs_target_bams/{sample}_{chrom_con}_{coverage_val}x.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools mpileup -f {input.ref_fasta_chr} -I -E -a 'FORMAT/DP' -T {input.ref_panel_sites_vcf} -r {wildcards.chrom_con} {input.downsampled_bam} -Ou | \
        bcftools call -Aim -C alleles -T {input.ref_panel_sites_tsv} -Oz -o {output.GL_vcf_target_bams} --threads {threads} 2> {log}
        
        bcftools index -f {output.GL_vcf_target_bams}
        '''

rule chunk_spliting_concordance:
    """
    Split chromosome into chunks for imputation
    """
    input:
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz'
    output:
        chunks = '{path}/output/GLIMPSE_concordance/chunks/{chrom_con}_chunks.txt'
    params:
        window_size = config['window_size'],
        buffer_size = config['buffer_size']
    log:
        '{path}/output/GLIMPSE_concordance/chunks/{chrom_con}_chunks.log'
    benchmark:
        '{path}/benchmarks/chunks/{chrom_con}_chunks.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        {glimpse_chunk} \
        --input {input.ref_panel_sites_vcf} \
        --region {wildcards.chrom_con} \
        --window-size {params.window_size} --buffer-size {params.buffer_size} \
        --thread 5 \
        --output {output.chunks} 2> {log}
        '''

rule impute_concordance:
    """
    Impute each sample seperately
    """
    input:
        GL_vcf_target_bams = '{path}/output/GLIMPSE_concordance/GLs_target_bams/{sample}_{chrom_con}_{coverage_val}x.vcf.gz',
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        chunks = '{path}/output/GLIMPSE_concordance/chunks/{chrom_con}_chunks.txt'
    output:
        imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x.00.bcf',
        imputed_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x.00.bcf.csi'
    params:    
        gen_map_path = config['gen_map_path'],
        gen_map_files = config['gen_map_files'],
        prefix = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x'
    threads: 2
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x.log'
    benchmark:
        '{path}/benchmarks/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            IRG=$(echo $LINE | cut -d" " -f3)
            ORG=$(echo $LINE | cut -d" " -f4)
            OUT={params.prefix}.${{ID}}.bcf
            {glimpse_impute} \
            --input {input.GL_vcf_target_bams} --reference {input.ref_concordance_sample_excl_filltags_filter} --map {params.gen_map_path}{params.gen_map_files} \
            --input-region ${{IRG}} \
            --output-region ${{ORG}} --output ${{OUT}} \
            --thread {threads}
            bcftools index -f ${{OUT}}
        done < {input.chunks} 2> {log}
        '''

rule ligate_list_concordance:
    """
    Create list of imputed output files for each chunk to merge later
    """
    input:
        chunks = '{path}/output/GLIMPSE_concordance/chunks/{chrom_con}_chunks.txt',
        imputed = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x.00.bcf'
    output:
        ligated_list = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom_con}_{coverage_val}x.txt'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/GLIMPSE_imputed/{sample}_{chrom_con}_{coverage_val}x'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        while IFS="" read -r LINE || [ -n "$LINE" ];
        do
            printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
            ls {params.prefix}.${{ID}}.bcf >> {output.ligated_list}
        done < {input.chunks}
        '''

rule ligate_concordance:
    """
    Merge all imputed chunks
    """
    input:
        ligated_list = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/ligated_list_{sample}_{chrom_con}_{coverage_val}x.txt'
    output:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
        ligated_bcf_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        {glimpse_ligate} \
        --input {input.ligated_list} \
        --output {output.ligated_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.ligated_bcf}
        '''

rule phase_concordance:
    """
    Phase!!!
    """
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf'
    output:
        phased_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x.bcf',
        phased_bcf_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x.log'
    threads: 8
    benchmark:
        '{path}/benchmarks/GLIMPSE_phased/phased.{sample}_{chrom_con}_{coverage_val}x.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        {glimpse_sample} \
        --input {input.ligated_bcf} --solve --output {output.phased_bcf} \
        --thread {threads} 2> {log}

        bcftools index -f {output.phased_bcf}
        '''
        


#####################################################
#                                                   #
#   Getting validation genotypes (the ground truth) #
#                                                   #
#####################################################


rule estimate_coverage_target_sample:
    input:
        target_bams = lambda wildcards: samples_df.loc[wildcards.sample, "bam_path"]
    output:
        cov_depth_cutoff = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_cov_depth_cutoff.txt'
    log:
        '{path}/output/GLIMPSE_concordance/target_bams/{sample}_genome_coverage.txt.log'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        #estimate the length of the regions within the bam file 
        LENGTH_GENOME=$(samtools view -H {input.target_bams} \
        | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{{sum+=$1}} END {{print sum}}')

        #estimate sample coverage
        COV_GENOME=$(samtools depth -a {input.target_bams}  \
        | awk -v var="$LENGTH_GENOME" '{{sum+=$3}} END {{ print sum/var}}')

        #avoid "invalid arithmetic operator" error when trying to add variables in Bash script
        C=$(echo $COV_GENOME/3 | bc)
        D=$(echo $COV_GENOME*2 | bc)

        #find max of cov/3 and 8 and add twice avg coverage to the same file
        echo $(($C>8 ? $C : 8)) > {output.cov_depth_cutoff}
        echo $D >> {output.cov_depth_cutoff}
        '''

rule prepare_validation_samples_filt: 
    """
    Take target bams (initial coverage) and call genotypes on the same site as the filtered reference panel, then filter sites based on Sousa da Mota 2022
    """
    input:
        target_bams_chr = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_{chrom_con}.bam',
        ref_panel_sites_vcf = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.vcf.gz',
        ref_panel_sites_tsv = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_sites.tsv.gz',
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta',
        cov_depth_cutoff = '{path}/output/GLIMPSE_concordance/target_bams/{sample}_cov_depth_cutoff.txt'
    output:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.bcf',
        validation_sample_filt_csi = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.log'
    threads: 10
    benchmark:
        '{path}/benchmarks/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.tsv'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        min=$(sed -n '1p' {input.cov_depth_cutoff})
        max=$(sed -n '2p' {input.cov_depth_cutoff})
        
        bcftools mpileup -f {input.ref_fasta_chr} \
        -I -E \
        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
        -T {input.ref_panel_sites_vcf} \
        -r {wildcards.chrom_con} \
        -q 30 -Q 30 -C 50 \
        --threads {threads} \
        {input.target_bams_chr} -Ou | \
        bcftools call -Aim -C alleles \
        -T {input.ref_panel_sites_tsv} --threads {threads} -Ou | \
        bcftools filter -i "%QUAL>=30 && FORMAT/DP<$max && FORMAT/DP>$min" \
        -Ob -o {output.validation_sample_filt} 2> {log}
        
        bcftools index -f {output.validation_sample_filt}
        '''

rule plot_ad_dp_prep:
    input:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt.bcf'
    output:
        validation_sample_filt_het = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.bcf',
        validation_sample_filt_het_sites = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.txt'
    shell:
        '''
        bcftools view --include 'GT="het"' {input.validation_sample_filt} -Ob -o {output.validation_sample_filt_het}
        bcftools query -f '%CHROM %POS %INFO/DP %INFO/AD{{1}}\n' {output.validation_sample_filt_het} > {output.validation_sample_filt_het_sites}
        '''

rule plot_ad_dp:
    input:
        validation_sample_filt_het_sites = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.txt'
    output:
        validation_sample_filt_het_plot = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_het.png'
    script:
        "../scripts/allelic_depth_read_depth_plot.R"
        

rule prepare_validation_samples_filt_allelic: 
    """
    Filter sites for allelic imbalance
    """
    input:
        validation_sample_filt = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp.bcf',
        ref_fasta_chr = '{path}/output/GLIMPSE_concordance/reference_genome/CanFam31_{chrom_con}.fasta' 
    output:
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        validation_sample_filt_csi_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf.csi'
    log:
        '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.log'
    threads: 10
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view --exclude 'GT="het" && ((INFO/AD[1] / INFO/DP < 0.15) || (INFO/AD[1] / INFO/DP > 0.85))' \
        {input.validation_sample_filt} \
        --threads {threads} \
        -Ob -o {output.validation_sample_filt_allelic} 2> {log}
        
        bcftools index -f {output.validation_sample_filt_allelic}
        '''


#############################
#                           #
#  Run GLIMPSE concordance  #
#                           #
#############################

rule filter_info_score:
    """
    Filter sites based on different INFO score cutoffs
    """
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
    output:
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf',
        info_imputed_info_csi = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf.csi'
    params:
        info_val = '{info_cutoff}'
    log:
        '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools view {input.ligated_bcf} \
        --include 'INFO/INFO >= {params.info_val}' \
        --threads {threads} \
        -Ob -o {output.info_imputed_info} 2> {log}

        bcftools index -f {output.info_imputed_info}
        '''

rule get_ID_for_targets:
    input:
        ligated_bcf = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated/merged_ligated.{sample}_{chrom_con}_{coverage_val}x.bcf',
    output:
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}_{coverage_val}x.txt'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        bcftools query -l {input.ligated_bcf} > {output.sm_samples}
        '''

rule prepare_concordance_lst_info_score_filtered:
    """
    Prepare the lst files required to run GLIMPSE_concordance
    """
    input:
        ref_concordance_sample_excl_filltags_filter = '{path}/output/GLIMPSE_concordance/reference_panel/{chrom_con}_ref_panel_filltags_filter.bcf',
        validation_sample_filt_allelic = '{path}/output/GLIMPSE_concordance/validation_bams/{sample}_{chrom_con}_validation_filt_qual_dp_ab.bcf',
        info_imputed_info = '{path}/output/GLIMPSE_concordance/GLIMPSE_ligated_INFO_filtered/merged_ligated.{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}.bcf',
    output:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.lst'
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        echo {wildcards.chrom_con} {input.ref_concordance_sample_excl_filltags_filter} {input.validation_sample_filt_allelic} {input.info_imputed_info} > {output.concordance_lst_info_score_filtered}
        '''

rule GLIMPSE_concordance_info_score_filtered:
    """
    Run GLIMPSE concordance specifying the target sample we want
    """
    input:
        concordance_lst_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.lst',
        sm_samples = '{path}/output/GLIMPSE_concordance/validation_bams/sm_{sample}_{chrom_con}_{coverage_val}x.txt'
    output:
        concordance_output_info_score_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.rsquare.grp.txt.gz',
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    params:
        prefix = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered'
    log:
        '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.log'
    threads: 8
    #conda:
        #'../envs/environment.yaml'
    shell:
        '''
        {glimpse_concordance} \
        --input {input.concordance_lst_info_score_filtered} \
        --minDP 8 \
        --output {params.prefix} \
        --minPROB 0.9 \
        --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
        --sample {input.sm_samples} \
        --af-tag AF \
        --thread {threads} 2> {log}
        '''

rule plot_rsquare_accuracy_filtered:
    """
    Plot accuracy 
    """
    input:
        concordance_output_info_score_1 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.8_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_2 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.9_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_3 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.95_filtered.rsquare.grp.txt.gz',
        concordance_output_info_score_4 = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_0.0_filtered.rsquare.grp.txt.gz'
    output:
        plot = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_{sample}_{chrom_con}_{coverage_val}x_filtered.png'
    params:
        chr = '{chrom_con}',
        name = '{sample}',
        cov = '{coverage_val}'
    #conda:
        #'../envs/r4.3.1.yaml'
    script:
        "../scripts/rsquare_accuracy.R"


rule prepare_concordance_output_filt:
    """
    Prepare files for genotype discordance plot
    """
    input:
        concordance_output_discordance_filtered = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.error.spl.txt.gz'
    output:
        concordance_output_discordance_filtered_temp = temp('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/temp_concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.txt'),
        concordance_output_discordance_filtered_prep = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.txt'
    shell:
        '''
        zcat {input.concordance_output_discordance_filtered} | sed -n '3p' >> {output.concordance_output_discordance_filtered_temp}
        awk '{{print "{wildcards.coverage_val}  {wildcards.info_cutoff}   "$0}}' {output.concordance_output_discordance_filtered_temp} > {output.concordance_output_discordance_filtered_prep}
        '''

rule merge_concordance_output_filt:
    """
    Merge all coverages and INFO per sample
    """
    input:
        concordance_output_discordance_filtered_prep = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_{coverage_val}x-INFO_{info_cutoff}_filtered.txt', coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF, allow_missing=True)
    output:
        concordance_output_discordance_filtered_per_sample = '{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filtered.txt'
    shell:
        '''
        cat {input.concordance_output_discordance_filtered_prep} > {output.concordance_output_discordance_filtered_per_sample}
        '''

rule plot_discordance_filt:
    """
    Plot genotype discordances
    """
    input:
        concordance_output_discordance_filtered_per_sample = expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered/concordance_{sample}_{chrom_con}_filtered.txt', sample=SAMPLE, allow_missing=True),
        concordance_metadata = config['bam_targets']
    params:
        path_script = '{path}/scripts',
        discordance_phased=lambda wildcards, input: ','.join(input.concordance_output_discordance_filtered_per_sample),
    output:
        discordance_dogs_full = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_dogs_full.png',
        discordance_dogs = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_dogs.png',
        discordance_wolves_full = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_wolves_full.png',
        discordance_wolves = '{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/concordance_{chrom_con}_filtered_wolves.png',
    #conda:
        #'../envs/r4.3.1.yaml'
    shell:
        '''
        Rscript {params.path_script}/genotype_discordance.R \
        {params.discordance_phased} \
        {input.concordance_metadata} \
        {output.discordance_dogs_full} \
        {output.discordance_dogs} \
        {output.discordance_wolves_full} \
        {output.discordance_wolves}
        '''