#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

#############################################################
# ROH counts and length on imputed and modern dogs and wolves#
#############################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(readr)
library("MetBrewer")
library(scales)
library("cowplot")
library(ggpubr)


# import imputed ROH data and info
info <- read.delim(snakemake@input[[1]])
roh <- read.csv(snakemake@input[[3]], sep = "")

colnames(info)[colnames(info) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] <- "Dog_PCA"
colnames(info)[colnames(info) == "Wolf.Dog_PCA"] <- "Wolf_Dog_PCA"
colnames(info)[1] <- "Sample"
colnames(roh)[1] <- "Sample"

# fix Tumat (Tumat2) and WolfHead (Wolf_head_IN18-016)
roh$Sample[roh$Sample == "Tumat"] <- "Tumat2"
roh$Sample[roh$Sample == "WolfHead"] <- "Wolf_head_IN18-016"

# merge two based on specific columns
final <- left_join(roh, info %>% dplyr::select(Sample, Wolf_Dog_PCA, Dog_PCA, Meta.Population, Species, Age_Mean_BP), "Sample")

# import modern ROH data and info
info_ref <- read.delim(snakemake@input[[2]])
roh_ref <- read.csv(snakemake@input[[4]], sep = "")
colnames(info_ref)[colnames(info_ref) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] <- "Dog_PCA"
colnames(info_ref)[colnames(info_ref) == "Wolf.Dog_PCA"] <- "Wolf_Dog_PCA"

# make new column without number in IDs
colnames(info_ref)[1] <- "Sample"
colnames(roh_ref)[1] <- "Sample"

# fix sample names which lost part of name in previous step:
roh_ref$Sample[roh_ref$Sample == "Wolf107"] <- "RWJR007"
roh_ref$Sample[roh_ref$Sample == "Wolf108"] <- "RWJR016"
roh_ref$Sample[roh_ref$Sample == "Wolf109"] <- "RWJR012"
roh_ref$Sample[roh_ref$Sample == "Wolf110"] <- "RWJR003"
roh_ref$Sample[roh_ref$Sample == "Wolf92"] <- "RKW7639"
roh_ref$Sample[roh_ref$Sample == "Wolf95"] <- "RKW7619"

# fix Wolf08 metapopulation:
info_ref$Meta.Population[info_ref$Sample == "Wolf08"] <- "Western_Eurasian_Wolves"

# remove duplicate samples from metadata (which were reseq)
# info_ref_no_duplicate = info_ref[!duplicated(info_ref$Sample),]

# merge two based on specific columns
final_ref <- left_join(roh_ref, info_ref %>% dplyr::select(Sample, Wolf_Dog_PCA, Dog_PCA, Meta.Population, Species, Age_Mean_BP), "Sample")
final_ref$Age_Mean_BP <- 0


# import genome sizes per chromosome (for Froh estimation)
sizes_autosomes <- read.table(snakemake@input[[5]], quote = "\"", comment.char = "")
total_genome_size <- sum(sizes_autosomes$V2)

# add type column
final$type <- "imputed"
final_ref$type <- "modern"

# extract only common columns from imputed and modern
common_cols <- intersect(colnames(final), colnames(final_ref))
all <- rbind(
    subset(final, select = common_cols),
    subset(final_ref, select = common_cols)
)

##### ######### ####
##### ALL ROHS #####
##### ######### ####

# sum total ROH count per sample:
a <- all %>%
    group_by(Sample) %>%
    mutate(n = n()) %>%
    distinct(Sample, .keep_all = TRUE)

# sum total ROH length per sample:
b <- all %>%
    group_by(Sample) %>%
    summarise(ROH_tol = sum(KB))

ab <- merge(a, b[c(1, 2)], by = "Sample")

# estimate Froh for all ROH sizes (have to multiply by 1000 since the ROHs are given in KBs):
ab$froh <- (ab$ROH_tol * 1000) / total_genome_size

# Get unique samples and add dummy values in columns for next step
t.first <- all[match(unique(all$Sample), all$Sample), ]
t.first$n <- "0"
t.first$ROH_tol <- "0"
t.first$froh <- "0"

# add row with missing sample to dataframe
for (i in 1:nrow(t.first)) {
    if (!(t.first$Sample[i] %in% ab$Sample)) {
        ab[nrow(ab) + 1, ] <- t.first[i, ]
    }
}

##########  PLOTTING  ##########

#############
# Wolves only
#############

my_breaks <- c(5000, 10000, 20000, 50000, 100000)
my_labels <- c("5,000", "10,000", "20,000", "50,000", "100,000")

# filter out wolves and unwanted outgroup populations:
wolves <- ab %>% filter(Dog_PCA == "Wolves" & Wolf_Dog_PCA != "Outgroup" & Wolf_Dog_PCA != "Wolves")

wolves$Age_Mean_BP <- as.numeric(wolves$Age_Mean_BP)
wolves$ROH_tol <- as.numeric(wolves$ROH_tol)
wolves$n <- as.numeric(wolves$n)
wolves$froh <- as.numeric(wolves$froh)

# re-scale x axis and age
wolves <- wolves %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )


wolves$Meta.Population[wolves$Sample == "Wolf08"] <- "Western_Eurasian_Wolves"


# make modern and ancient layers for plotting:
df_layer_1 <- wolves[wolves$type == "modern", ]
df_layer_2 <- wolves[wolves$type != "modern", ]

group_names <- c(
    "Pleistocene_Wolves" = "Pleistocene",
    "Eastern_Eurasian_Wolves" = "East Eurasia",
    "Western_Eurasian_Wolves" = "West Eurasia",
    "North_American_Wolves" = "North America"
)


# ROH count against size
png(snakemake@output[[1]], width = 11, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 4, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 6, shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "F") +
    labs(x = "Total ROH length (Mb)", y = "Total # ROH") +
    labs(fill = "Sample age (kya)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("Pleistocene_Wolves", "Eastern_Eurasian_Wolves", "Western_Eurasian_Wolves", "North_American_Wolves")), labeller = as_labeller(group_names)) +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()





########################
# only pleistocene wolves

# filter out samples:
wolves_pleist <- ab %>% filter(Wolf_Dog_PCA == "Pleistocene_Wolves")

wolves_pleist$Age_Mean_BP <- as.numeric(wolves_pleist$Age_Mean_BP)
wolves_pleist$ROH_tol <- as.numeric(wolves_pleist$ROH_tol)
wolves_pleist$n <- as.numeric(wolves_pleist$n)
wolves_pleist$froh <- as.numeric(wolves_pleist$froh)

# re-scale x axis and age
wolves_pleist <- wolves_pleist %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )




##### ######### #### #####
##### Long ROHS >=1.6Mb #####
##### ######### #### #####

# sum total long ROH count per sample:

a_long <- all %>%
    filter(KB >= 1600) %>%
    group_by(Sample) %>%
    mutate(n = n()) %>%
    distinct(Sample, .keep_all = TRUE)


b_long <- all %>%
    filter(KB >= 1600) %>%
    group_by(Sample) %>%
    summarise(ROH_tol = sum(KB))

ab_long <- merge(a_long, b_long[c(1, 2)], by = "Sample")

# estimate Froh for all ROH sizes:
ab_long$froh <- (ab_long$ROH_tol * 1000) / total_genome_size

# Get unique samples and add dummy values in columns for next step
t.first <- all[match(unique(all$Sample), all$Sample), ]
t.first$n <- "0"
t.first$ROH_tol <- "0"
t.first$froh <- "0"

# add row with missing sample to dataframe
for (i in 1:nrow(t.first)) {
    if (!(t.first$Sample[i] %in% ab_long$Sample)) {
        ab_long[nrow(ab_long) + 1, ] <- t.first[i, ]
    }
}

##########  PLOTTING  ##########

#############
# Wolves only
#############

my_breaks <- c(5000, 10000, 20000, 50000, 100000)
my_labels <- c("5,000", "10,000", "20,000", "50,000", "100,000")

# filter out wolves and unwanted outgroup populations:
wolves_long <- ab_long %>% filter(Dog_PCA == "Wolves" & Wolf_Dog_PCA != "Outgroup" & Wolf_Dog_PCA != "Wolves")

wolves_long$Age_Mean_BP <- as.numeric(wolves_long$Age_Mean_BP)
wolves_long$ROH_tol <- as.numeric(wolves_long$ROH_tol)
wolves_long$n <- as.numeric(wolves_long$n)
wolves_long$froh <- as.numeric(wolves_long$froh)

# re-scale x axis and age
wolves_long <- wolves_long %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )


wolves_long$Meta.Population[wolves_long$Sample == "Wolf08"] <- "Western_Eurasian_Wolves"


# make modern and ancient layers for plotting:
df_layer_1 <- wolves_long[wolves_long$type == "modern", ]
df_layer_2 <- wolves_long[wolves_long$type != "modern", ]

group_names <- c(
    "Pleistocene_Wolves" = "Pleistocene",
    "Eastern_Eurasian_Wolves" = "East Eurasia",
    "Western_Eurasian_Wolves" = "West Eurasia",
    "North_American_Wolves" = "North America"
)


png(snakemake@output[[2]], width = 11, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(wolves_long, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 4, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 6, shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "F") +
    labs(x = "Total ROH length (Mb) (ROH >= 1.6Mb)", y = "Total # ROH (ROH >= 1.6Mb)") +
    labs(fill = "Sample age (kya)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("Pleistocene_Wolves", "Eastern_Eurasian_Wolves", "Western_Eurasian_Wolves", "North_American_Wolves")), labeller = as_labeller(group_names)) +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()





########################
# only pleistocene wolves

# filter out samples:
wolves_pleist_long <- ab_long %>% filter(Wolf_Dog_PCA == "Pleistocene_Wolves")

wolves_pleist_long$Age_Mean_BP <- as.numeric(wolves_pleist_long$Age_Mean_BP)
wolves_pleist_long$ROH_tol <- as.numeric(wolves_pleist_long$ROH_tol)
wolves_pleist_long$n <- as.numeric(wolves_pleist_long$n)
wolves_pleist_long$froh <- as.numeric(wolves_pleist_long$froh)


# re-scale x axis and age
wolves_pleist_long <- wolves_pleist_long %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )




##### ######### #### #####
##### Short ROHS <1.6Mb #####
##### ######### #### #####

# sum total short ROH count per sample:

a_short <- all %>%
    filter(KB < 1600) %>%
    group_by(Sample) %>%
    mutate(n = n()) %>%
    distinct(Sample, .keep_all = TRUE)


b_short <- all %>%
    filter(KB < 1600) %>%
    group_by(Sample) %>%
    summarise(ROH_tol = sum(KB))

ab_short <- merge(a_short, b_short[c(1, 2)], by = "Sample")

# estimate Froh for all ROH sizes:
ab_short$froh <- (ab_short$ROH_tol * 1000) / total_genome_size

# Get unique samples and add dummy values in columns for next step
t.first <- all[match(unique(all$Sample), all$Sample), ]
t.first$n <- "0"
t.first$ROH_tol <- "0"
t.first$froh <- "0"

# add row with missing sample to dataframe
for (i in 1:nrow(t.first)) {
    if (!(t.first$Sample[i] %in% ab_short$Sample)) {
        ab_short[nrow(ab_short) + 1, ] <- t.first[i, ]
    }
}

##########  PLOTTING  ##########

#############
# Wolves only
#############

my_breaks <- c(5000, 10000, 20000, 50000, 100000)
my_labels <- c("5,000", "10,000", "20,000", "50,000", "100,000")

# filter out wolves and unwanted outgroup populations:
wolves_short <- ab_short %>% filter(Dog_PCA == "Wolves" & Wolf_Dog_PCA != "Outgroup" & Wolf_Dog_PCA != "Wolves")

wolves_short$Age_Mean_BP <- as.numeric(wolves_short$Age_Mean_BP)
wolves_short$ROH_tol <- as.numeric(wolves_short$ROH_tol)
wolves_short$n <- as.numeric(wolves_short$n)
wolves_short$froh <- as.numeric(wolves_short$froh)

# re-scale x axis and age
wolves_short <- wolves_short %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )


wolves_short$Meta.Population[wolves_short$Sample == "Wolf08"] <- "Western_Eurasian_Wolves"


# make modern and ancient layers for plotting:
df_layer_1 <- wolves_short[wolves_short$type == "modern", ]
df_layer_2 <- wolves_short[wolves_short$type != "modern", ]

group_names <- c(
    "Pleistocene_Wolves" = "Pleistocene",
    "Eastern_Eurasian_Wolves" = "East Eurasia",
    "Western_Eurasian_Wolves" = "West Eurasia",
    "North_American_Wolves" = "North America"
)


png(snakemake@output[[3]], width = 11, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(wolves_short, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 4, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 6, shape = 21, alpha = 0.7) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "F") +
    labs(x = "Total ROH length (Mb) (ROH < 1.6Mb)", y = "Total # ROH (ROH < 1.6Mb)") +
    labs(fill = "Sample age (kya)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("Pleistocene_Wolves", "Eastern_Eurasian_Wolves", "Western_Eurasian_Wolves", "North_American_Wolves")), labeller = as_labeller(group_names)) +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()





########################
# only pleistocene wolves

# filter out samples:
wolves_pleist_short <- ab_short %>% filter(Wolf_Dog_PCA == "Pleistocene_Wolves")

wolves_pleist_short$Age_Mean_BP <- as.numeric(wolves_pleist_short$Age_Mean_BP)
wolves_pleist_short$ROH_tol <- as.numeric(wolves_pleist_short$ROH_tol)
wolves_pleist_short$n <- as.numeric(wolves_pleist_short$n)
wolves_pleist_short$froh <- as.numeric(wolves_pleist_short$froh)

# re-scale x axis and age
wolves_pleist_short <- wolves_pleist_short %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )





###### ###### ###### ######
###### COMBINED PLOTS #####
###### ###### ###### ######

wolves$category <- "all_ROHs"
wolves_short$category <- "short_ROHs"
wolves_long$category <- "long_ROHs"

wolves_all_ROHs <- rbind(wolves, wolves_short, wolves_long)

# make modern and ancient layers for plotting:
df_layer_1 <- wolves_all_ROHs[wolves_all_ROHs$type == "modern", ]
df_layer_2 <- wolves_all_ROHs[wolves_all_ROHs$type != "modern", ]

group_names <- c(
    "Pleistocene_Wolves" = "Pleistocene",
    "Eastern_Eurasian_Wolves" = "East Eurasia",
    "Western_Eurasian_Wolves" = "West Eurasia",
    "North_American_Wolves" = "North America"
)

category_names <- c(
    "all_ROHs" = "All ROH",
    "long_ROHs" = "ROH >= 1.6Mb",
    "short_ROHs" = "ROH < 1.6Mb"
)

cols <- c("royalblue4", "darkturquoise", "yellow3")

png(snakemake@output[[4]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
wolves_coeff <- ggplot(wolves_all_ROHs, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 2, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 3, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols, labels = group_names, name = "Population") +
    scale_x_continuous(breaks = seq(round(min(wolves_short$Age_Mean_KBP + 1)), 0, 10)) +
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH]))) +
    theme_bw() +
    facet_grid(factor(category) ~ ., labeller = as_labeller(category_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
wolves_coeff
dev.off()


png(snakemake@output[[5]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
wolves_coeff <- ggplot(wolves_all_ROHs, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 2, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 3, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols, labels = group_names, name = "Population") +
    scale_x_continuous(breaks = seq(round(min(wolves_short$Age_Mean_KBP + 1)), 0, 10)) +
    geom_label_repel(
        data = wolves_all_ROHs %>% filter(type == "imputed"),
        aes(x = Age_Mean_KBP, y = froh, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH]))) +
    theme_bw() +
    facet_grid(factor(category) ~ ., labeller = as_labeller(category_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
wolves_coeff
dev.off()

################
#### Box plots
################

cols <- c("royalblue4", "green4", "yellow3")

modern_wolves_box <- wolves_all_ROHs %>%
    filter(type == "modern") %>%
    ggplot(aes(x = Meta.Population, y = froh)) +
    geom_boxplot(aes(fill = Meta.Population), outlier.shape = NA) +
    geom_jitter(colour = "gray28", size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = cols) +
    scale_x_discrete(labels = c(
        "Eastern_Eurasian_Wolves" = "East Eurasia", "North_American_Wolves" = "North America",
        "Western_Eurasian_Wolves" = "West Eurasia"
    )) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(y = expression(paste(italic("F")[ROH]))) +
    facet_grid(. ~ category, labeller = as_labeller(category_names))

cols <- c("royalblue4", "darkturquoise", "yellow3")

ancient_wolves_sub <- wolves_all_ROHs %>% filter(type == "imputed")
ancient_wolves_sub$Meta.Population[ancient_wolves_sub$Meta.Population == "Eastern_Eurasian_Wolves"] <- "Holocene"
ancient_wolves_sub$Meta.Population[ancient_wolves_sub$Meta.Population == "Western_Eurasian_Wolves"] <- "Holocene"

ancient_wolves_box <- ancient_wolves_sub %>%
    ggplot(aes(x = Meta.Population, y = froh)) +
    geom_boxplot(aes(fill = Meta.Population), outlier.shape = NA) +
    geom_jitter(colour = "gray28", size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = cols) +
    scale_x_discrete(labels = c("Holocene" = "Holocene", "Pleistocene_Wolves" = "Pleistocene")) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(y = expression(paste(italic("F")[ROH]))) +
    facet_grid(. ~ category, labeller = as_labeller(category_names))

####################################### merge_plots  ######################################

# supp figure boxplot modern ancient:
png(snakemake@output[[6]], width = 14, height = 10, units = "in", res = 200, pointsize = 4)
ggarrange(ancient_wolves_box, modern_wolves_box,
    labels = c("a", "b"),
    ncol = 1, nrow = 2, font.label = list(size = 15)
)
dev.off()


#### #### #### #### #### ####
#### PLEISTOCENE WOLVES ONLY
#### #### #### #### #### ####

wolves_pleist$category <- "all_ROHs"
wolves_pleist_short$category <- "short_ROHs"
wolves_pleist_long$category <- "long_ROHs"


wolves_pleist_all_ROHs <- rbind(wolves_pleist, wolves_pleist_short, wolves_pleist_long)


category_names <- c(
    "all_ROHs" = "All ROH",
    "long_ROHs" = "ROH >= 1.6Mb",
    "short_ROHs" = "ROH < 1.6Mb"
)
cols <- "darkturquoise"

png(snakemake@output[[7]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
wolves_pleist_coeff <- ggplot(wolves_pleist_all_ROHs, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(aes(fill = Meta.Population), size = 3, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(wolves_pleist_all_ROHs$Age_Mean_KBP + 1)), 0, 10)) +
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH]))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(category) ~ ., labeller = as_labeller(category_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
wolves_pleist_coeff
dev.off()

png(snakemake@output[[8]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
wolves_pleist_coeff <- ggplot(wolves_pleist_all_ROHs, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(aes(fill = Meta.Population), size = 3, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(wolves_pleist_all_ROHs$Age_Mean_KBP + 1)), 0, 10)) +
    geom_label_repel(
        data = wolves_pleist_all_ROHs %>% filter(type != "modern"),
        aes(x = Age_Mean_KBP, y = froh, label = Sample), size = 3.5, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH]))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(category) ~ ., labeller = as_labeller(category_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 11),
        axis.text.y = element_text(size = 11),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
wolves_pleist_coeff
dev.off()


#### #### #### #### #### ####
#### ALL ROH DATA WOLVES
#### #### #### #### #### ####

# export FROH values etc:
ab$category <- "all_rohs"
ab_long$category <- "long_rohs"
ab_short$category <- "short_rohs"
all_roh_results <- rbind(ab, ab_long, ab_short)

write.table(all_roh_results, file = snakemake@output[[9]], quote = FALSE, sep = "\t", row.names = FALSE)
