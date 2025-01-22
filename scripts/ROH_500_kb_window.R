#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(windowscanr)

options(scipen = 999)

# take chr size file:
ref_fasta_chr_size <- read.delim(snakemake@input[[1]], header = FALSE)

colnames(ref_fasta_chr_size) <- c("CHR", "BP")
ref_fasta_chr_size$CHR <- as.numeric(gsub("chr", "", ref_fasta_chr_size$CHR))

hom_sum <- ref_fasta_chr_size %>%
    arrange(CHR, BP) %>%
    head(38) %>%
    mutate(index = 1:nrow(.))

# dummy column to make the winscan run:
hom_sum$UNAFF <- 0

running_roh <- winScan(
    x = hom_sum,
    groups = "CHR",
    position = "BP",
    values = "UNAFF",
    win_size = 500000,
    win_step = 500000,
    funs = c("mean"),
    cores = 8
)

# take for df:
running_roh$region <- paste("chr", running_roh$CHR, ":", running_roh$win_start, "-", running_roh$win_end, sep = "")

# export sites to estimate coverage:
write.table(running_roh$region, snakemake@output[[1]], col.names = FALSE, row.names = FALSE, quote = FALSE)
