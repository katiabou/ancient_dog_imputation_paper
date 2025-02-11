#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# load libraries
library(data.table)
library(tidyverse)
library(windowscanr)
library(cowplot)
library(grid)
library(ggplotify)
library(patchwork)
library(viridis)
library(scales)
library(ggridges)
library(lme4)
library(broom)
library(ggpubr)

options(scipen = 999)

# imputed dogs
# hom_sum <-fread("~/Desktop/merged_phased.allchrom_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary")
# ind_file <- fread("~/Desktop/merged_phased.chr1_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.indiv")
hom_sum <- fread(snakemake@input[[1]])
ind_file <- fread(snakemake@input[[2]])

# modern dogs
# hom_sum_modern <-fread("~/Desktop/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.summary")
# ind_file_modern <- fread("~/Desktop/ref-panel_chr1_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogs.hom.indiv")
hom_sum_modern <- fread(snakemake@input[[3]])
ind_file_modern <- read.csv(snakemake@input[[4]], sep = "")

# imputed and modern dogs
# hom_sum_alldogs <-fread("~/Desktop/merged_phased_modern.allchrom_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.summary")
# ind_file_alldogs <- fread("~/Desktop/merged_phased_modern.chr1_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogs.hom.indiv")
hom_sum_alldogs <- fread(snakemake@input[[5]])
ind_file_alldogs <- fread(snakemake@input[[6]])

# per window coverage estimate file
# all_cov_window <- read.delim("~/Desktop/dogs_allchrom_windows_cov_500kb.txt", header=FALSE)
all_cov_window <- read.delim(snakemake@input[[7]], header = FALSE)



# Estimating the window cutoffs based on imputed dogs or imputed wolves (from the bams)

# create index column
hom_sum <- hom_sum %>%
    mutate(index = 1:nrow(.))

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

# remove windows without snps
running_roh %>%
    mutate(UNAFF_mean = UNAFF_mean / nrow(ind_file)) %>%
    filter(UNAFF_n > 0) -> running_roh_p

# take for df:
running_roh_p$region <- paste(running_roh_p$CHR, ":", running_roh_p$win_start, "-", running_roh_p$win_end, sep = "")


# prepare window coverage file
all_cov_window$V1 <- str_replace_all(all_cov_window$V1, "chr", "")
all_cov_window$region <- paste(all_cov_window$V1, ":", all_cov_window$V2, "-", all_cov_window$V3, sep = "")

# count number of bam samples used for the window coverage calculation (same number as the imputed dogs or wolves)
num_ind <- nrow(ind_file)

# Group by region and estimate across all samples the average mean depth for each window:
all_cov_window %>%
    group_by(region) %>%
    summarize(Mean = sum(V7) / num_ind) -> test


# merge the window coverage with the ROH data (either imputed dogs or imputed wolves)
joinn <- left_join(running_roh_p, test, by = c("region"))

# get mean and std
meann <- mean(joinn$Mean)
stdd <- 2 * sd(joinn$Mean)
low <- meann - stdd
high <- meann + stdd


png(snakemake@output[[1]], width = 9, height = 6, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
ggplot(joinn, aes(x = Mean, UNAFF_mean)) +
    geom_point(data = . %>% filter(Mean > low & Mean < high), alpha = 0.6, colour = "cyan3") +
    geom_point(data = . %>% filter(Mean < low | Mean > high), alpha = 0.6, colour = "grey40") +
    geom_smooth(data = . %>% filter(Mean > low & Mean < high), method = "lm", colour = "cyan4") +
    geom_vline(
        xintercept = low, linetype = "dashed",
        color = "red3", size = 0.8
    ) +
    geom_vline(
        xintercept = high, linetype = "dashed",
        color = "red3", size = 0.8
    ) +
    theme_bw() +
    ylim(0, 1) +
    labs(x = "Mean depth per window (500 KB)", y = "Mean ROH prevelance per window (500 KB)")
dev.off()



######################################
######### HEAT MAPS ##################
######################################

####################################### All (imputed and modern)  ######################################

hom_sum_all <- hom_sum_alldogs %>%
    mutate(
        MB = BP / 1000000,
        KB = BP / 1000,
        index = 1:nrow(.)
    )

# count ROH in running windows of 500 Kb
# UNAFF	is the number of non-cases with a ROH including this SNP
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh_all <- winScan(
    x = hom_sum_all,
    groups = "CHR",
    position = "BP",
    values = "UNAFF",
    win_size = 500000,
    win_step = 500000,
    funs = c("mean"),
    cores = 8
)

# prepare region colomn:
running_roh_all$region <- paste(running_roh_all$CHR, ":", running_roh_all$win_start, "-", running_roh_all$win_end, sep = "")


# remove windows without snps
running_roh_all %>%
    mutate(UNAFF_mean = UNAFF_mean / nrow(ind_file_alldogs)) %>%
    filter(UNAFF_n > 0) -> running_roh_p_all


# prepare colour scale:
fill_cols <- viridis(20, option = "A")
qn <- scales::rescale(quantile(running_roh_p_all$UNAFF_mean,
    probs = seq(0, 1, length.out = length(fill_cols))
))


# windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean < low | Mean > high)


# remove windows with extreme depth estimates (may have CNVs):
running_roh_p_all_nogrey <- running_roh_p_all %>%
    filter(!region %in% final_window_remove$region)

p1 <- ggplot(running_roh_p_all_nogrey, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) +
    # geom_tile(color = "grey", size = 0) +
    geom_tile() +
    geom_tile(data = final_window_remove, aes(x = win_start, y = 0.5), fill = "grey") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(0, 125000000, by = 10000000),
        labels = as.character(seq(0, 125000, 10000) / 1000)
    ) +
    ylab("Chromosome") +
    scale_fill_gradientn("% of present-day and ancient samples with ROH",
        colors = rev(fill_cols),
        # values = qn,
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50),
        na.value = "grey"
    ) +
    facet_grid(CHR ~ ., switch = "both") +
    xlab("Position in Mb") +
    theme_minimal(
        base_family = "Helvetica",
        base_size = 13
    ) +
    theme(
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 5), size = 15),
        axis.title.y = element_text(margin = margin(r = 5), size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 14),
        strip.text.y.left = element_text(size = 14, angle = 0),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.75, 0.13),
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        # axis.line.x = element_line(size = 0.3)
    ) +
    ggtitle("Ancient and present-day") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
    guides(fill = guide_colourbar(
        title.position = "bottom",
        barwidth = 12, barheight = 0.8
    ))
p1

ggsave(snakemake@output[[2]], p1, width = 7, height = 6, bg = "transparent")


# density plot (add manually to previous plot)
x <- running_roh_p_all_nogrey$UNAFF_mean
y <- density(x, n = 2^12)

p2 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) +
    geom_line() +
    geom_segment(aes(xend = x, yend = 0, color = x)) +
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_gradientn(
        colors = rev(fill_cols),
        # values = qn,
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50)
    ) +
    theme_minimal(
        base_family = "Helvetica",
        base_size = 13
    ) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    ylab("Density")

p2

ggsave(snakemake@output[[3]], p2, width = 5, height = 2.5)

# get 5% of rows:
# perc <- round(nrow(running_roh_p_all_nogrey)*0.05)

# get top and bottom 5% of sites
# extremes_all_high <- running_roh_p_all_nogrey %>%
#  arrange(UNAFF_mean) %>%
#  tail(n=perc)

# extremes_all_low<- running_roh_p_all_nogrey %>%
#  arrange(UNAFF_mean) %>%
#  head(n=perc)

# extremes_all <- rbind(extremes_all_high, extremes_all_low)





####################################### Imputed dogs  ######################################

### only ancient dog heatmap

hom_sum <- hom_sum %>%
    mutate(
        MB = BP / 1000000,
        KB = BP / 1000,
        index = 1:nrow(.)
    )

# count ROH in running windows of 500 Kb
# UNAFF	is the number of non-cases with a ROH including this SNP
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh_imputed <- winScan(
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
running_roh_imputed$region <- paste(running_roh_imputed$CHR, ":", running_roh_imputed$win_start, "-", running_roh_imputed$win_end, sep = "")


# remove windows without snps
running_roh_imputed %>%
    mutate(UNAFF_mean = UNAFF_mean / nrow(ind_file)) %>%
    filter(UNAFF_n > 0) -> running_roh_p_imputed

# prepare colour scale:
fill_cols <- viridis(20, option = "A")
qn <- scales::rescale(quantile(running_roh_p_imputed$UNAFF_mean,
    probs = seq(0, 1, length.out = length(fill_cols))
))

# windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean < low | Mean > high)

# remove windows with extreme depth estimates (may have CNVs):
running_roh_p_imputed_nogrey <- running_roh_p_imputed %>%
    filter(!region %in% final_window_remove$region)

p3 <- ggplot(running_roh_p_imputed_nogrey, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) +
    # geom_tile(color = "grey", size = 0) +
    geom_tile() +
    geom_tile(data = final_window_remove, aes(x = win_start, y = 0.5), fill = "grey") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(0, 125000000, by = 10000000),
        labels = as.character(seq(0, 125000, 10000) / 1000)
    ) +
    ylab("Chromosome") +
    scale_fill_gradientn("% of ancient samples with ROH",
        colors = rev(fill_cols),
        # values = qn,
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50),
        na.value = "grey"
    ) +
    facet_grid(CHR ~ ., switch = "both") +
    xlab("Position in Mb") +
    theme_minimal(
        base_family = "Helvetica",
        base_size = 13
    ) +
    theme(
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 5), size = 15),
        axis.title.y = element_text(margin = margin(r = 5), size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 14),
        strip.text.y.left = element_text(size = 14, angle = 0),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.75, 0.13),
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        # axis.line.x = element_line(size = 0.3)
    ) +
    ggtitle("Ancient") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
    guides(fill = guide_colourbar(
        title.position = "bottom",
        barwidth = 12, barheight = 0.8
    ))
p3

ggsave(snakemake@output[[4]], p3, width = 7, height = 6, bg = "transparent")


# density plot (add manually to previous plot)
x <- running_roh_p_imputed_nogrey$UNAFF_mean
y <- density(x, n = 2^12)

p4 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) +
    geom_line() +
    geom_segment(aes(xend = x, yend = 0, color = x)) +
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50)
    ) +
    scale_y_continuous(expand = c(0, 0),
        breaks = c(2,4,6,8,10,12)
    ) +
    scale_color_gradientn(
        colors = rev(fill_cols),
        # values = qn,
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50)
    ) +
    theme_minimal(
        base_family = "Helvetica",
        base_size = 13
    ) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    ylab("Density")

p4

ggsave(snakemake@output[[5]], p4, width = 5, height = 2.5)
# ggsave('~/Desktop/den_imp.png', p4, width = 5, height = 2.5)

# get 5% of rows:
# perc <- round(nrow(running_roh_p_imputed_nogrey)*0.05)

# get top and bottom 5% of sites
# extremes_imputed_high <- running_roh_p_imputed_nogrey %>%
#  arrange(UNAFF_mean) %>%
#  tail(n=perc)

# extremes_imputed_low<- running_roh_p_imputed_nogrey %>%
#  arrange(UNAFF_mean) %>%
#  head(n=perc)

# extremes_imputed <- rbind(extremes_imputed_high, extremes_imputed_low)


####################################### Modern dogs  ######################################

### only ancient dog heatmap
hom_sum_modern <- hom_sum_modern %>%
    mutate(
        MB = BP / 1000000,
        KB = BP / 1000,
        index = 1:nrow(.)
    )

# count ROH in running windows of 500 Kb
# UNAFF	is the number of non-cases with a ROH including this SNP
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh_modern <- winScan(
    x = hom_sum_modern,
    groups = "CHR",
    position = "BP",
    values = "UNAFF",
    win_size = 500000,
    win_step = 500000,
    funs = c("mean"),
    cores = 8
)

# take for df:
running_roh_modern$region <- paste(running_roh_modern$CHR, ":", running_roh_modern$win_start, "-", running_roh_modern$win_end, sep = "")


# remove windows without snps
running_roh_modern %>%
    mutate(UNAFF_mean = UNAFF_mean / nrow(ind_file_modern)) %>%
    filter(UNAFF_n > 0) -> running_roh_p_modern

# prepare colour scale:
fill_cols <- viridis(20, option = "A")
qn <- scales::rescale(quantile(running_roh_p_modern$UNAFF_mean,
    probs = seq(0, 1, length.out = length(fill_cols))
))

# windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean < low | Mean > high)


# remove windows with extreme depth estimates (may have CNVs):
running_roh_p_modern_nogrey <- running_roh_p_modern %>%
    filter(!region %in% final_window_remove$region)

p5 <- ggplot(running_roh_p_modern_nogrey, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) +
    # geom_tile(color = "grey", size = 0) +
    geom_tile() +
    geom_tile(data = final_window_remove, aes(x = win_start, y = 0.5), fill = "grey") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(0, 125000000, by = 10000000),
        labels = as.character(seq(0, 125000, 10000) / 1000)
    ) +
    ylab("Chromosome") +
    scale_fill_gradientn("% of present-day samples with ROH",
        colors = rev(fill_cols),
        # values = qn,
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50),
        na.value = "grey"
    ) +
    facet_grid(CHR ~ ., switch = "both") +
    xlab("Position in Mb") +
    # theme_simple(base_size = 13, grid_lines = FALSE, base_family = "Helvetica") +
    theme_minimal(
        base_family = "Helvetica",
        base_size = 13
    ) +
    theme(
        panel.spacing.y = unit(0.1, "lines"),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(t = 5), size = 15),
        axis.title.y = element_text(margin = margin(r = 5), size = 15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 14),
        strip.text.y.left = element_text(size = 14, angle = 0),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.75, 0.13),
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        # axis.line.x = element_line(size = 0.3)
    ) +
    ggtitle("Present-day") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)) +
    guides(fill = guide_colourbar(
        title.position = "bottom",
        barwidth = 12, barheight = 0.8
    ))
p5

ggsave(snakemake@output[[6]], p5, width = 7, height = 6, bg = "transparent")


# density plot (add manually to previous plot)
x <- running_roh_p_modern_nogrey$UNAFF_mean
y <- density(x, n = 2^12)

p6 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) +
    geom_line() +
    geom_segment(aes(xend = x, yend = 0, color = x)) +
    scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_color_gradientn(
        colors = rev(fill_cols),
        # values = qn,
        limits = c(0, 0.55),
        breaks = c(0.1, 0.3, 0.5),
        labels = c(10, 30, 50)
    ) +
    # theme_simple(axis_lines = TRUE, grid_lines = FALSE, base_family = "Helvetica") +
    theme_minimal(
        base_family = "Helvetica",
        base_size = 13
    ) +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()
    ) +
    ylab("Density")

p6

ggsave(snakemake@output[[7]], p6, width = 5, height = 2.5)
# ggsave('~/Desktop/den_mod.png', p6, width = 5, height = 2.5)


# get 5% of rows:
# perc <- round(nrow(running_roh_p_modern_nogrey)*0.05)

# get top and bottom 5% of sites
# extremes_modern_high <- running_roh_p_modern_nogrey %>%
#             arrange(UNAFF_mean) %>%
#             tail(n=perc)

# extremes_modern_low<- running_roh_p_modern_nogrey %>%
#           arrange(UNAFF_mean) %>%
#           head(n=perc)

# extremes_modern <- rbind(extremes_modern_high, extremes_modern_low)



# Export window files:
write.table(running_roh_p_all_nogrey, snakemake@output[[8]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(running_roh_p_imputed_nogrey, snakemake@output[[9]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(running_roh_p_modern_nogrey, snakemake@output[[10]], sep = "\t", row.names = FALSE, quote = FALSE)



####################################### Gene ontology  ######################################

########### Ancient imputed dogs:

#### ISLANDS
imputed_high <- running_roh_p_imputed_nogrey %>%
    filter(UNAFF_mean >= 0.05) %>%
    select(region)

# prepare bed format for bedtools intersect:
imputed_high_bed <- imputed_high %>%
    separate_wider_delim(region, ":", names = c("chr", "region")) %>%
    separate_wider_delim(region, "-", names = c("start", "end"))

imputed_high_bed$chr <- as.numeric(imputed_high_bed$chr)
imputed_high_bed$start <- as.numeric(imputed_high_bed$start)
imputed_high_bed$end <- as.numeric(imputed_high_bed$end)

# sort
imputed_high_bed_sorted <- imputed_high_bed %>%
    arrange(chr, start)

#### DESERTS
imputed_low <- running_roh_p_imputed_nogrey %>%
    filter(UNAFF_mean < 0.05) %>%
    select(region)

# prepare bed format for bedtools intersect:
imputed_low_bed <- imputed_low %>%
    separate_wider_delim(region, ":", names = c("chr", "region")) %>%
    separate_wider_delim(region, "-", names = c("start", "end"))

imputed_low_bed$chr <- as.numeric(imputed_low_bed$chr)
imputed_low_bed$start <- as.numeric(imputed_low_bed$start)
imputed_low_bed$end <- as.numeric(imputed_low_bed$end)

# sort
imputed_low_bed_sorted <- imputed_low_bed %>%
    arrange(chr, start)


########### Modern  dogs:

#### ISLANDS
modern_high <- running_roh_p_modern_nogrey %>%
    filter(UNAFF_mean >= 0.05) %>%
    select(region)

# prepare bed format for bedtools intersect:
modern_high_bed <- modern_high %>%
    separate_wider_delim(region, ":", names = c("chr", "region")) %>%
    separate_wider_delim(region, "-", names = c("start", "end"))

modern_high_bed$chr <- as.numeric(modern_high_bed$chr)
modern_high_bed$start <- as.numeric(modern_high_bed$start)
modern_high_bed$end <- as.numeric(modern_high_bed$end)

# sort
modern_high_bed_sorted <- modern_high_bed %>%
    arrange(chr, start)


#### DESERTS
modern_low <- running_roh_p_modern_nogrey %>%
    filter(UNAFF_mean < 0.05) %>%
    select(region)

# prepare bed format for bedtools intersect:
modern_low_bed <- modern_low %>%
    separate_wider_delim(region, ":", names = c("chr", "region")) %>%
    separate_wider_delim(region, "-", names = c("start", "end"))

modern_low_bed$chr <- as.numeric(modern_low_bed$chr)
modern_low_bed$start <- as.numeric(modern_low_bed$start)
modern_low_bed$end <- as.numeric(modern_low_bed$end)

# sort
modern_low_bed_sorted <- modern_low_bed %>%
    arrange(chr, start)

# export bed files of windows for gene annotation:
write.table(imputed_high_bed_sorted, snakemake@output[[11]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(imputed_low_bed_sorted, snakemake@output[[12]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(modern_high_bed_sorted, snakemake@output[[13]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(modern_low_bed_sorted, snakemake@output[[14]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





##############################################################
# FIND COMMON OVERLAPS BETWEEN ANCIENT AND MODERN for ISLANDS:

modern_imputed_overlap_islands <- as.data.frame(intersect(modern_high$region, imputed_high$region))
colnames(modern_imputed_overlap_islands) <- "Common_regions"

# prepare bed format for bedtools intersect:
modern_imputed_overlap_islands_bed <- modern_imputed_overlap_islands %>%
    separate_wider_delim(Common_regions, ":", names = c("chr", "region")) %>%
    separate_wider_delim(region, "-", names = c("start", "end"))

modern_imputed_overlap_islands_bed$chr <- as.numeric(modern_imputed_overlap_islands_bed$chr)
modern_imputed_overlap_islands_bed$start <- as.numeric(modern_imputed_overlap_islands_bed$start)
modern_imputed_overlap_islands_bed$end <- as.numeric(modern_imputed_overlap_islands_bed$end)

# sort
modern_imputed_overlap_islands_bed_sorted <- modern_imputed_overlap_islands_bed %>%
    arrange(chr, start)

# export windows for gene annotation:
write.table(modern_imputed_overlap_islands_bed_sorted, snakemake@output[[15]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)




##############################################################
# FIND COMMON OVERLAPS BETWEEN ANCIENT AND MODERN for DESERTS:

modern_imputed_overlap_deserts <- as.data.frame(intersect(modern_low$region, imputed_low$region))
colnames(modern_imputed_overlap_deserts) <- "Common_regions"

# prepare bed format for bedtools intersect:
modern_imputed_overlap_deserts_bed <- modern_imputed_overlap_deserts %>%
    separate_wider_delim(Common_regions, ":", names = c("chr", "region")) %>%
    separate_wider_delim(region, "-", names = c("start", "end"))

modern_imputed_overlap_deserts_bed$chr <- as.numeric(modern_imputed_overlap_deserts_bed$chr)
modern_imputed_overlap_deserts_bed$start <- as.numeric(modern_imputed_overlap_deserts_bed$start)
modern_imputed_overlap_deserts_bed$end <- as.numeric(modern_imputed_overlap_deserts_bed$end)

# sort
modern_imputed_overlap_deserts_bed_sorted <- modern_imputed_overlap_deserts_bed %>%
    arrange(chr, start)

# export windows for gene annotation:
write.table(modern_imputed_overlap_deserts_bed_sorted, snakemake@output[[16]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



##############################################################
# RUN GOfunR for gene ontology analysis:

library(GOfuncR)

# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("TxDb.Cfamiliaris.UCSC.canFam3.refGene")
# BiocManager::install("org.Cf.eg.db")

# windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean < low | Mean > high)

# grey regions to not use as background
final_window_keep <- test %>%
    filter(!region %in% final_window_remove$region) %>%
    dplyr::select(region)


# desert regions
imputed_modern_regions <- modern_imputed_overlap_deserts_bed_sorted
imputed_modern_regions$range <- paste(imputed_modern_regions$start, imputed_modern_regions$end, sep = "-")
imputed_modern_regions$region <- paste(imputed_modern_regions$chr, imputed_modern_regions$range, sep = ":")
imputed_modern_regions_sub <- imputed_modern_regions %>% dplyr::select(region)


# prepare the background and target regions
is_candidate_back <- data.frame(final_window_keep$region, is_candidate = c(rep(0, nrow(final_window_keep))))
is_candidate_target <- data.frame(imputed_modern_regions_sub$region, is_candidate = c(rep(1, nrow(imputed_modern_regions_sub))))
colnames(is_candidate_back) <- c("regions", "is_candidate")
colnames(is_candidate_target) <- c("regions", "is_candidate")
is_candidate_all <- rbind(is_candidate_back, is_candidate_target)

# Running a hypergeometric test with correction for gene length using regions
# I have to set background regions when I use the regions option:
res_hyper1 <- go_enrich(is_candidate_all,
    test = "hyper", n_randsets = 1000,
    orgDb = "org.Cf.eg.db",
    txDb = "TxDb.Cfamiliaris.UCSC.canFam3.refGene", gene_len = TRUE,
    regions = TRUE
)


## first element of go_enrich result has the stats
stats1 <- res_hyper1[[1]]
## top-GO categories
# head(stats1)

## top GO-categories per domain
# by(stats1, stats1$ontology, head, n=3)

## all valid input genes
# head(res_hyper1[[2]])

## annotation package used (default='Homo.sapiens') and GO-graph version
# res_hyper1[[3]]

## minimum p-values from randomsets
# head(res_hyper1[[4]])

## hypergeometric test
# top_gos_hyper1 = res_hyper1[[1]][1:10, 'node_id']
# GO-categories with a high proportion of candidate genes
# top_gos_hyper1

# plot_anno_scores(res_hyper1, top_gos_hyper1)


# output table with GO terms:
write.table(stats1, snakemake@output[[17]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)





####################################### merge_plots  ######################################

# supp figure heatmap:
png(snakemake@output[[18]], width = 17, height = 10, units = "in", res = 200, pointsize = 4)
# png('~/Desktop/heat.png', width=17, height=10, units='in', res=200, pointsize=4)
ggarrange(p3, p5,
    labels = c("a", "b"),
    ncol = 2, nrow = 1, font.label = list(size = 20)
)
dev.off()
