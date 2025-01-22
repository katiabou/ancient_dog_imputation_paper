#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)

# import data
initial <- read.table(snakemake@input[[1]], quote = "\"", comment.char = "")
colnames(initial) <- c("chr", "pos", "dp", "ad_alt")

df2 <- subset(initial, dp != ".")

df2$dp <- as.numeric(df2$dp)
df2$ad_alt <- as.numeric(df2$ad_alt)

df2$ratio <- df2$ad_alt / df2$dp

options(scipen = 999)

png(snakemake@output[[1]])
ggplot(df2, aes(x = ratio)) +
    geom_histogram() +
    scale_x_continuous(n.breaks = 10) +
    geom_vline(xintercept = 0.15, colour = "red") +
    geom_vline(xintercept = 0.85, colour = "blue") +
    xlab("AD alt / DP")

dev.off()
