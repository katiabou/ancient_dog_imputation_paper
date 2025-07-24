#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(dplyr)
library(tidyr)
library(data.table)
library(readr)

# target chromosome:
chr <- snakemake@params[[1]]

# read depth window file
cnv <- fread(snakemake@input[[1]])

print("Imported cnv file")

# window coordinates
window_coords <- read.table(snakemake@input[[2]])
colnames(window_coords) <- c("region")

# make columns rows and rows columns
cnv_t <- as.data.frame(t(cnv))

print("Transposed dataframe")

# merge window names to depth file
merged_cnv <- cbind(window_coords, cnv_t)
merged_cnv <- as.data.table(merged_cnv)

# extract chromosome
filtered_dt <- merged_cnv[grepl(paste0(chr,":"), region)]

# round read depth to nearest integer (do not include first column with region name)
cnv_round <- as.data.table(filtered_dt[,-1])
cnv_round[] <- lapply(cnv_round, round)

print("Rounded read depths to nearest integer")

# then if x=2 make it binary, give it a 1
# if x>2 or x<2, give it a 0
cnv_round_binary <- ifelse(cnv_round==2,1,0)
cnv_round_binary <- as.data.frame(cnv_round_binary)

print("Turned into binary format")

# export specific dataframe
write_delim(cnv_round_binary, snakemake@output[[1]], col_names = FALSE)

