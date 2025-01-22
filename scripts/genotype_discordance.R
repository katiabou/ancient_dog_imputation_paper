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
library("MetBrewer")
library(stringr)
library(data.table)

# import data
args <- commandArgs(trailingOnly = TRUE)

discordance_phased <- args[1]
meta <- read.csv(args[2], sep = "")

a <- str_split(discordance_phased, pattern = ",")

all <- c()

for (i in 1:length(a[[1]])) {
    tmp <- read.table(a[[1]][i], quote = "\"", comment.char = "")[, c(1:2, 5, 12:15)]
    all <- rbind(all, tmp)
}

colnames(all) <- c("cov", "info", "sample", "RR", "RA", "AA", "NRD")

mdata <- reshape2::melt(all, id = c("sample", "cov", "info"))
mdata$cov <- as.character(mdata$cov)

# add into to all df from meta df
mdata$species <- meta$Wolf_Dog_PCA[match(mdata$sample, meta$Original_ID)]
mdata$info_sample <- meta$Info[match(mdata$sample, meta$Original_ID)]
mdata$info_sample_name <- paste(mdata$sample, mdata$info_sample, sep = " - ")
mdata$info_sample_name <- gsub("_", " ", mdata$info_sample_name)

info_names <- c(
    "0" = "No INFO cutoff",
    "0.8" = "0.8 INFO cutoff",
    "0.9" = "0.9 INFO cutoff",
    "0.95" = "0.95 INFO cutoff",
    "RR" = "RR",
    "RA" = "RA",
    "AA" = "AA",
    "NRD" = "NRD"
)


# subset dogs
mdata_dogs <- mdata %>% filter(species == "Dogs")

# order alphabetically
mdata_dogs <- mdata_dogs[order(mdata_dogs$sample), ]

# specify colour palettes
cols <- met.brewer(name = "Cross", n = 7, type = "discrete")

# plot with all dog samples and all INFO/cov
png(args[3], width = 10, height = 7, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
m <- ggplot(data = mdata_dogs, aes(cov, value, col = info_sample_name))
m + geom_point(aes(col = info_sample_name), size = 3, alpha = 0.8) +
    labs(y = "Genotyping error %", x = "Coverage") +
    scale_colour_manual(values = cols, name = "Dog samples") +
    scale_y_continuous(
        breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    ) +
    scale_x_discrete(
        labels = c("0.05x", "0.1x", "0.2x", "0.5x", "1x", "2x")
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, size = 10, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 10, vjust = 0.5)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 12)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white")) +
    geom_hline(
        yintercept = 10, linetype = "dashed",
        color = "grey4", linewidth = 0.5
    ) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    guides(colour = guide_legend(byrow = TRUE)) +
    facet_grid(info ~ variable, labeller = as_labeller(info_names))
dev.off()


# plot with all dog samples and 0.8 INFO, above 0.5x only
mdata_dogs_sub <- mdata_dogs %>% filter(info == 0.8 & cov >= 0.5)

png(args[4], width = 12, height = 3, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
m <- ggplot(data = mdata_dogs_sub, aes(cov, value, col = info_sample_name))
m + geom_point(aes(col = info_sample_name), size = 3, alpha = 0.8) +
    labs(y = "Genotyping error %", x = "Coverage") +
    scale_colour_manual(values = cols, name = "Dog samples") +
    scale_y_continuous(
        limits = c(0, 30),
        breaks = c(0, 5, 10, 15, 20, 25, 30)
    ) +
    scale_x_discrete(
        labels = c("0.5x", "1x", "2x")
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, size = 10, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 10, vjust = 0.5)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white")) +
    geom_hline(
        yintercept = 10, linetype = "dashed",
        color = "grey4", size = 0.5
    ) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    guides(colour = guide_legend(byrow = TRUE)) +
    facet_grid(info ~ variable, labeller = as_labeller(info_names))
dev.off()


# plot with all wolf samples and all INFO/cov

# subset dogs
mdata_wolves <- mdata %>% filter(species == "Pleistocene_wolf")

# order alphabetically
mdata_wolves <- mdata_wolves[order(mdata_wolves$sample), ]

cols <- met.brewer(name = "Hokusai3", n = 3, type = "discrete")

# plot with all dog samples and all INFO/cov
png(args[5], width = 10, height = 7, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
m <- ggplot(data = mdata_wolves, aes(cov, value, col = info_sample_name))
m + geom_point(aes(col = info_sample_name), size = 3, alpha = 0.8) +
    labs(y = "Genotyping error %", x = "Coverage") +
    scale_colour_manual(values = cols, name = "Wolf samples") +
    scale_y_continuous(
        breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    ) +
    scale_x_discrete(
        labels = c("0.05x", "0.1x", "0.2", "0.5x", "1x", "2x")
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, size = 10, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 10, vjust = 0.5)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 12)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white")) +
    geom_hline(
        yintercept = 10, linetype = "dashed",
        color = "grey4", linewidth = 0.5
    ) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    guides(colour = guide_legend(byrow = TRUE)) +
    facet_grid(info ~ variable, labeller = as_labeller(info_names))
dev.off()


# plot with all wolf samples and 0.8 INFO, above 0.5x only
mdata__wolves_sub <- mdata_wolves %>% filter(info == 0.8 & cov >= 0.5)

png(args[6], width = 12, height = 3, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
m <- ggplot(data = mdata__wolves_sub, aes(cov, value, col = info_sample_name))
m + geom_point(aes(col = info_sample_name), size = 3, alpha = 0.8) +
    labs(y = "Genotyping error %", x = "Coverage") +
    scale_colour_manual(values = cols, name = "Wolf samples") +
    scale_y_continuous(
        limits = c(0, 30),
        breaks = c(0, 5, 10, 15, 20, 25, 30)
    ) +
    scale_x_discrete(
        labels = c("0.5x", "1x", "2x")
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, size = 10, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 10, vjust = 0.5)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 10)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white")) +
    geom_hline(
        yintercept = 10, linetype = "dashed",
        color = "grey4", size = 0.5
    ) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    guides(colour = guide_legend(byrow = TRUE)) +
    facet_grid(info ~ variable, labeller = as_labeller(info_names))
dev.off()
