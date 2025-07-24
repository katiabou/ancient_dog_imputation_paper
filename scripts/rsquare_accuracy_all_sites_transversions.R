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
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

# import list of files all sites
files <- args[1]

a <- str_split(files, pattern = ",")

# merge all flare dfs
all <- c()

for (i in 1:length(a[[1]])) {
    tmp <- read.table(a[[1]][i], quote = "\"", comment.char = "")
    all <- rbind(all, tmp)
}

colnames(all) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")

# import list of files all sites
files_trans <- args[2]

b <- str_split(files_trans, pattern = ",")

# merge all flare dfs
trans <- c()

for (i in 1:length(b[[1]])) {
    tmp <- read.table(b[[1]][i], quote = "\"", comment.char = "")
    trans <- rbind(trans, tmp)
}

colnames(trans) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")


cols <- c("#CC6677", "#009E73", "#DDCC77", "#0072B2")


# get break points for x axis, the average of the MAF values in the df
breaks <- all %>%
    group_by(V1) %>%
    summarize(Avg = log(mean(V3)))

# import metadata file
meta <- read.csv(args[3], sep = "")

# add all sites or trans info to dataframe
all$type <- "Transversions & Transitions"
trans$type <- "Transversions"

# merge two datasets
df_both <- rbind(all, trans)

# add into to all df from meta df
df_both$species <- meta$Wolf_Dog_PCA[match(df_both$V8, meta$Original_ID)]
df_both$info_sample <- meta$Info[match(df_both$V8, meta$Original_ID)]
df_both$info_sample_name <- paste(df_both$V8, df_both$info_sample, sep = " - ")
df_both$info_sample_name <- gsub("_", " ", df_both$info_sample_name)


cov_info <- c(
    "0.05" = "0.05x",
    "0.1" = "0.1x",
    "0.2" = "0.2x",
    "0.5" = "0.5x",
    "1" = "1x",
    "2" = "2x"
)


# define order of facet
df_both$type_b <- factor(df_both$type, levels = c("Transversions & Transitions", "Transversions"))
df_both$V8_b <- factor(df_both$V8, levels = c(
    "NGDG", "SOTN01_merged", "PortauChoix", "TRF.02.53", "TRF.05.05",
    "FAMICHN00012", "FAMINGR00004", "CGG32", "CGG33", "WolfHead"
))

p1 <- df_both %>%
    filter(species == "Dogs" & V6 == 0.5) %>%
    ggplot(aes(log(V3), V4, colour = as.character(V7))) +
    geom_line(linetype = "solid", size = 1.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = cols, name = "INFO cutoff", labels = c("No cutoff", "0.8", "0.9", "0.95")) +
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
        strip.text.y = element_text(angle = -90, size = 18),
        strip.text.x = element_text(size = 18),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    labs(y = bquote("r"^2), x = "Minor allele frequency") +
    facet_grid(V8_b ~ type_b)


p2 <- df_both %>%
    filter(species == "Pleistocene_wolf" & V6 == 1) %>%
    ggplot(aes(log(V3), V4, colour = as.character(V7))) +
    geom_line(linetype = "solid", size = 1.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = cols, name = "INFO cutoff", labels = c("No cutoff", "0.8", "0.9", "0.95")) +
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
        strip.text.y = element_text(angle = -90, size = 18),
        strip.text.x = element_text(size = 18),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    labs(y = bquote("r"^2), x = "Minor allele frequency") +
    facet_grid(V8 ~ type_b)


png(args[4], width = 21, height = 24, units = "in", res = 200, pointsize = 4)
plot1 <- ggarrange(p1,
    ggarrange(p2, nrow = 2.5, labels = c("b"), font.label = list(size = 20)),
    ncol = 2,
    labels = "a", font.label = list(size = 20)
)
plot1
dev.off()
