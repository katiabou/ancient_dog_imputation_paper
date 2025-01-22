######################################################
# Prepare reference panel for concordance imputation #
######################################################


rule remove_sample_indels_multiallelic_snps_concordance_only_dogs:
    """
    Remove samples, trim-alt alleles created, normalise indels, only keep biallelic snps
    """
    input:
        ref_panel_phased="output/reference_panel/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        non_dog_names="sample_lists/non_dog_samples.tsv",
    output:
        ref_sample_snp=temp(
            "output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp.phased.vcf.gz"
        ),
    log:
        "output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp.phased.tsv"
    shell:
        """
        (
        bcftools view \
        -r {wildcards.chrom} \
        -S ^{input.non_dog_names} \
        --trim-alt-alleles \
        {input.ref_panel_phased} -Ou | \
        bcftools norm -a -Ou | \
        bcftools view -m 2 -M 2 -v snps  \
        -Oz -o {output.ref_sample_snp}
        ) 2> {log}
        """


rule fill_tags_concordance_only_dogs:
    """
    Fill tags to re-estimate fields after sample removal (have to specify F_MISSING, which is the fraction of missing genotypes)
    """
    input:
        ref_sample_snp="output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp.phased.vcf.gz",
    output:
        ref_sample_snp_filltags=temp(
            "output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags.phased.vcf.gz"
        ),
    log:
        "output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags.phased.tsv"
    shell:
        """
        (
        bcftools +fill-tags {input.ref_sample_snp} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags} \
        -- -t all,F_MISSING
        ) 2> {log} 
        """


rule filter_sites_concordance_only_dogs:
    """
    Filter for only PASS sites and missingness (F_MISSING)
    """
    input:
        ref_sample_snp_filltags="output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags.phased.vcf.gz",
    output:
        ref_sample_snp_filltags_filter="output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz",
        ref_sample_snp_filltags_filter_tbi="output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags_filter.phased.vcf.gz.tbi",
    params:
        f_missing=config["F_MISSING"],
    log:
        "output/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags_filter.phased.log",
    threads: 4
    benchmark:
        "benchmarks/reference_panel_only_dogs/ref-panel_{chrom}_sample-snp_filltags_filter.phased.tsv"
    shell:
        """
        bcftools view -i 'FILTER=="PASS" & F_MISSING<{params.f_missing}' {input.ref_sample_snp_filltags} \
        --threads {threads} \
        -Oz -o {output.ref_sample_snp_filltags_filter} 2> {log}

        bcftools index --tbi {output.ref_sample_snp_filltags_filter}
        """
