#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

############################################################
#  Download and merge all the publicly available BAM files #
############################################################

import pandas as pd


def get_accession_ftp(wildcards):
    """
    Get the FTP URL for an ENA accession
    """
    ena_metadata = pd.read_table(
        f"sample_lists/ena/filereport_read_run_{wildcards.prj_code}.tsv"
    ).set_index("run_accession", drop=False)
    bam_ftp = ena_metadata.loc[wildcards.accession, "bam_ftp"]

    return bam_ftp


rule download_ena_accession:
    """
    Download a published BAM file from the ENA 
    """
    params:
        url=get_accession_ftp,
    resources:
        ftp=1,
    output:
        bam=temp("samples/bam/ena/{prj_code}_{accession}.bam"),
    shell:
        """
        wget --quiet -O {output.bam} -o /dev/null ftp://{params.url}
        """


def list_sample_accesions(wildcards):
    """
    Fetch a list of all the ENA accession codes for a given sample
    """
    samples = pd.read_table(
        "sample_lists/bams_published_imputation_metadata_cutoff.tsv"
    ).set_index("Sample", drop=False)

    # get the ENA project code for the given sample
    prj_code = samples.loc[wildcards.sample, "ena_prj"]

    # load the accession metadata for that PRJ code
    ena_metadata = pd.read_table(
        f"sample_lists/ena/filereport_read_run_{prj_code}.tsv"
    ).set_index("sample_id", drop=False)
    accessions = ena_metadata.loc[wildcards.sample, "run_accession"]

    return expand(
        "samples/bam/ena/{prj_code}_{accession}.bam",
        prj_code=prj_code,
        accession=accessions,
    )


rule merge_ena_accessions:
    """
    Merge all individual BAM files for a published sample. 
    """
    input:
        bams=list_sample_accesions,
    output:
        bam="samples/bam/{sample}_merged.bam",
        bai="samples/bam/{sample}_merged.bai",
    log:
        bam="samples/bam/{sample}_merged.log",
    threads: 4
    shell:
        """
        bcftools merge {input.bams} --threads {threads} -Ob -o {output.bam} 2> {log}
        bcftools index -f {output.bam}
        """
