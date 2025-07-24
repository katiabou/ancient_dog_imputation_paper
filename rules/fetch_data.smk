#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Katia Bougiouri"
__copyright__ = "Copyright 2025, University of Copenhagen"
__email__ = "katia.bougiouri@gmail.com"
__license__ = "MIT"

###########################################
#  Download reference panel and BAM files #
###########################################

rule get_reference_panel:
    """
    Download VCF to prepare reference panel
    """
    output:
        vcf = "data/reference_panel_vcf/1697g_WildSled.SNP.INDEL.chrAll.newID.vcf.gz",
        tbi = "data/reference_panel_vcf/1697g_WildSled.SNP.INDEL.chrAll.newID.vcf.gz.tbi",
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/reference_panel_vcf/1697g_WildSled.SNP.INDEL.chrAll.newID.vcf.gz > {output.vcf}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/reference_panel_vcf/1697g_WildSled.SNP.INDEL.chrAll.newID.vcf.gz.tbi > {output.tbi}
        """

rule get_ref_fasta:
    """
    Download canfam3.1 reference fasta file
    """
    output:
        fa = "data/reference_fasta/canFam3_withY.fa",
        fai = "data/reference_fasta/canFam3_withY.fa.fai",
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/reference_fasta/canFam3_withY.fa > {output.fa}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/reference_fasta/canFam3_withY.fa.fai > {output.fai}
        """

rule get_bams_concordance:
    """
    Download BAM files concordance
    """
    output:
        bam="data/samples/bams_concordance/{bam_con}",
        bai="data/samples/bams_concordance/{bam_con}.bai",
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/samples/bams_concordance/{wildcards.bam_con} > {output.bam}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/samples/bams_concordance/{wildcards.bam_con}.bai > {output.bai}
        """

rule get_bams:
    """
    Download BAM files for full dataset
    """
    output:
        bam="data/samples/bams_full_dataset/{bam}",
        bai="data/samples/bams_full_dataset/{bam}.bai",
    shell:
        """
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/samples/bams_full_dataset/{wildcards.bam} > {output.bam}
        wget --quiet -O - https://sid.erda.dk/share_redirect/Hvjfs9cMxP/samples/bams_full_dataset/{wildcards.bam}.bai > {output.bai}
        """