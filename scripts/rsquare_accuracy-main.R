#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(stringr)
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# import list of files
files <- args[1]

a <- str_split(files, pattern = ",")

# merge all flare dfs
all <- c()

for (i in 1:length(a[[1]])) {
    tmp <- read.table(a[[1]][i], quote = "\"", comment.char = "")
    all <- rbind(all, tmp)
}

# colnames(all) <- c('V1','V2','V3','V4','V5','V6','V7','V8')

# s1 <- read.table("~/Downloads/concordance_NGDG_allchrom_0.5x-INFO_0.8_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s2 <- read.table("~/Downloads/concordance_NGDG_allchrom_0.5x-INFO_0.9_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s3 <- read.table("~/Downloads/concordance_NGDG_allchrom_1x-INFO_0.8_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s4 <- read.table("~/Downloads/concordance_NGDG_allchrom_1x-INFO_0.9_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s5 <- read.table("~/Downloads/concordance_NGDG_allchrom_2x-INFO_0.8_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s6 <- read.table("~/Downloads/concordance_NGDG_allchrom_2x-INFO_0.9_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s7 <- read.table("~/Downloads/concordance_PortauChoix_allchrom_0.5x-INFO_0.8_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s8 <- read.table("~/Downloads/concordance_PortauChoix_allchrom_0.5x-INFO_0.9_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s9 <- read.table("~/Downloads/concordance_PortauChoix_allchrom_1x-INFO_0.8_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s10 <- read.table("~/Downloads/concordance_PortauChoix_allchrom_1x-INFO_0.9_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s11 <- read.table("~/Downloads/concordance_PortauChoix_allchrom_2x-INFO_0.8_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
# s12 <- read.table("~/Downloads/concordance_PortauChoix_allchrom_2x-INFO_0.9_filtered.rsquare-mod.grp.txt.gz", quote="\"", comment.char="")
#
# all <- rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)


cols <- c("#CC6677", "#009E73", "#DDCC77", "#0072B2")


# get break points for x axis, the average of the MAF values in the df
breaks <- all %>%
    group_by(V1) %>%
    summarize(Avg = log(mean(V3)))

sample_info <- c(
    "NGDG" = "Neolithic European dog\nNGDG",
    "PortauChoix" = "North American pre-contact dog\nPortauChoix",
    "TRF.05.05" = "Iron age Siberian dog\nTRF.05.05",
    "CGG33" = "Pleistocene wolf\nCGG33"
)

cov_info <- c(
    "0.5" = "0.5x",
    "1" = "1x",
    "2" = "2x"
)

# define order of facet
all$V8_b <- factor(all$V8, levels = c("NGDG", "PortauChoix", "TRF.05.05", "CGG33"))

png(args[2], width = 25, height = 16, units = "in", res = 250, pointsize = 4)
ggplot(all, aes(log(V3), V4, colour = as.character(V7))) +
    geom_rect(
        xmin = breaks$Avg[5], xmax = max(breaks$Avg), ymin = -0.5, ymax = 1.5,
        fill = "gray80", colour = "white", alpha = 0.01
    ) +
    geom_line(linetype = "solid", size = 1.5) +
    geom_point(size = 3) +
    geom_vline(xintercept = breaks$Avg[5], linetype = "dashed") +
    scale_colour_manual(values = cols, name = "INFO cutoff", labels = c("No cutoff", "0.8", "0.9", "0.95")) +
    # scale_colour_manual(values = cols, name = "INFO cutoff", labels = c("0.8", "0.9"))+
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    scale_x_continuous(
        breaks = breaks$Avg,
        labels = c("0-0.001", "0.001-0.002", "0.002-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.5")
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 30, size = 15, vjust = 0.5),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white"),
        strip.text.y = element_text(angle = 0, size = 18),
        strip.text.x = element_text(size = 18),
        panel.grid.minor = element_blank()
    ) +
    labs(y = bquote("r"^2), x = "Minor allele frequency") +
    facet_grid(V6 ~ V8_b, labeller = labeller(V8_b = sample_info, V6 = cov_info), scales = "free_x")
dev.off()
