#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(readr)
library(dplyr)

# import DP estimates
DP_list <- read.table(snakemake@input[[1]], quote = "\"", comment.char = "")

# change column names
colnames(DP_list) <- c("chrom", "pos", "DP")

# remove missing sites (I think)
DP_list_n <- DP_list %>% filter(DP != ".")

DP_list_n$DP <- as.numeric(DP_list_n$DP)

# group by chrom and pos
test <- DP_list_n %>%
    group_by(chrom, pos) %>%
    summarise_at(vars(DP), list(name = mean))

# group by chrom
test_1 <- test %>%
    group_by(chrom) %>%
    summarise_at(vars(name), list(name = mean))

colnames(test_1) <- c("chrom", "DoC")

# find the maximum of doc/3 and 8 to use as a cutoff
doc_min <- max((test_1[1, 2] / 3), 8)

# set the maximum doc (twice the average)
doc_max <- test_1[1, 2] * 2

# combine min and max
final <- cbind(doc_min, doc_max)

# add column names
colnames(final) <- c("min_doc", "max_doc")

write.table(final, snakemake@output[[1]], sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
