##########################################################
# Call pseudohaploid sites on downsampled and HC for PCA #
##########################################################


rule create_bamlist_per_sample_downsampled:
    """
    Create list of bam files for downsampled target bams before PH calling
    """
    input:
        downsampled_bam="output/GLIMPSE_concordance/target_bams/{sample}_{chrom}_{coverage_val}x.bam",
    output:
        downsampled_bam_list="output/GLIMPSE_concordance/target_bams/list_{sample}_{chrom}_{coverage_val}x.txt",
    shell:
        """
        ls  {input.downsampled_bam} > {output.downsampled_bam_list}
        """


rule pseudohaploid_calling_downsampled:
    """
    PH calling for downsampled target samples on ALL SITES and not the reference panel sites, since it's crazy slow if you specify regions
    """
    input:
        downsampled_bam_list="output/GLIMPSE_concordance/target_bams/list_{sample}_{chrom}_{coverage_val}x.txt",
    output:
        validation_sample_filt_pseudohaploid_sites="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_{coverage_val}x_validation.haplo.gz",
    params:
        prefix="output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_{coverage_val}x_validation",
    threads: 2
    log:
        "output/GLIMPSE_concordance/validation_bams/{sample}_{chrom}_{coverage_val}x_validation.haplo.gz.log",
    shell:
        """
        angsd -b {input.downsampled_bam_list} \
        -doHaploCall 1 \
        -doCounts 1 \
        -minMapQ 30 -minQ 30 \
        -trim 5 \
        -noTrans 1 \
        -P {threads} \
        -checkBamHeaders 0 \
        -out {params.prefix} 2> {log}
        """
