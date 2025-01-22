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


# import imputed ROH data and info
# info <- read.delim('~/Downloads/roh_check/Dog_Wolf_aDNA_WG-Master.tsv')
# roh <- read.csv('~/Downloads/roh_check/merged_phased.allchrom_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogwolf.hom', sep="")
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
# info_ref <- read.delim('~/Downloads/roh_check/Dog_Wolf_aDNA_WG-Modern.tsv')
# roh_ref <- read.csv('~/Downloads/roh_check/ref-panel_allchrom_sample-snp_filltags_filter_MAF_0.01_all_sites_hom_win_het_1_dogwolf.hom', sep="")
info_ref <- read.delim(snakemake@input[[2]])
roh_ref <- read.csv(snakemake@input[[4]], sep = "")
colnames(info_ref)[colnames(info_ref) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] <- "Dog_PCA"
colnames(info_ref)[colnames(info_ref) == "Wolf.Dog_PCA"] <- "Wolf_Dog_PCA"

# make new column without number in IDs
# roh_ref$Sample <- sub("_[^_]+$", "", roh_ref$FID)
# colnames(info_ref)[1] <- "Sample"
colnames(info_ref)[3] <- "Sample"
colnames(roh_ref)[1] <- "Sample"

# fix sample names which lost part of name in previous step:
# roh_ref$Sample[roh_ref$Sample == 'Bern'] <- 'Bern_AlpineDachsbracke'
# roh_ref$Sample[roh_ref$Sample == 'CatahoulaLeopardDog01_Reseq'] <- 'CatahoulaLeopardDog01'
# roh_ref$Sample[roh_ref$Sample == 'MIX'] <- 'MIX_Dachshund01'
# roh_ref$Sample[roh_ref$Sample == 'MIX_AmericanCocker'] <- 'MIX_AmericanCocker_Beagle01'
# roh_ref$Sample[roh_ref$Sample == 'MIX_KerryBlueTerrier'] <- 'MIX_KerryBlueTerrier_Beagle01'
# roh_ref$Sample[roh_ref$Sample == 'MIX_MiniatureSchnauzer'] <- 'MIX_MiniatureSchnauzer_Beagle01'
# roh_ref$Sample[roh_ref$Sample == 'VillDog'] <- 'VillDog_Australia01'
# roh_ref$Sample[roh_ref$Sample == 'Wolf_WO001'] <- 'Wolf_WO001_895'
# roh_ref$Sample[roh_ref$Sample == 'Wolf_WO002'] <- 'Wolf_WO002_732'
# roh_ref$Sample[roh_ref$Sample == 'Wolf_WO003'] <- 'Wolf_WO003_636'

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
# sizes_autosomes <- read.delim('~/Downloads/CanFam31_allchr_size.genome', header=FALSE)
# sizes_autosomes <- read.delim(snakemake@input[[5]], header=FALSE)
sizes_autosomes <- read.table(snakemake@input[[5]], quote = "\"", comment.char = "")
# sizes_autosomes <- read.table('~/Downloads/roh_check/CanFam31_allchrom_size.genome', quote="\"", comment.char="")
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

# export table with Froh and ROH data
write.table(ab, file = snakemake@output[[]], quote = FALSE, sep = "\t", row.names = FALSE)



##########  PLOTTING  ##########

#############
# Dogs only
#############

my_breaks <- c(0, 100, 200, 500, 1000, 2000, 5000, 10000)
my_labels <- c("0", "100", "200", "500", "1,000", "2,000", "5,000", "10,000")

# filter out wolves and unwanted dog populations:
no_wolves <- ab %>% filter(Dog_PCA != "Wolves" & Dog_PCA != "Americas_Dogs" & Dog_PCA != "East_Asian_Dogs")

# make port au choix arctic:
no_wolves$Dog_PCA[no_wolves$Dog_PCA == "preContact_Dogs"] <- "Arctic_Dogs"
no_wolves$Meta.Population[no_wolves$Meta.Population == "preContact_Dogs"] <- "Arctic_Dogs"
no_wolves$Meta.Population[no_wolves$Meta.Population == "American_European_Dogs"] <- "European_Dogs"

no_wolves$Age_Mean_BP <- as.numeric(no_wolves$Age_Mean_BP)
no_wolves$ROH_tol <- as.numeric(no_wolves$ROH_tol)
no_wolves$n <- as.numeric(no_wolves$n)
no_wolves$froh <- as.numeric(no_wolves$froh)

# re-scale x axis and age
no_wolves <- no_wolves %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )

# make modern and ancient layers for plotting:
df_layer_1 <- no_wolves[no_wolves$type == "modern", ]
df_layer_2 <- no_wolves[no_wolves$type != "modern", ]

group_names <- c(
    "African_NearEast_India_Dogs" = "Near East",
    "Arctic_Dogs" = "Arctic",
    "European_Dogs" = "Europe"
)

# chosen_dogs <- c('AL2022_Turkey','AL2571_Iran', 'ASHQ06_P8903', 'C89_Ajvide','C90_Ajvide','SOTN01_merged',
#                 'TRF.04.09','CGG6','PortauChoix','TRF.05.14','TRF.02.25','OL4061_Veretye', 'KT0056')

# ROH count against size
png(snakemake@output[[1]], width = 7, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "G") +
    geom_smooth(method = "lm", se = FALSE, color = "gray28", size = 0.5, alpha = 0.8) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=ROH_tol_2, y=n, label=Sample),size=3.5, box.padding = 3, max.overlaps = Inf)+
    labs(x = "Total ROH length (Mb)", y = "Total # ROH") +
    labs(fill = "Sample age (ybp)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")), labeller = as_labeller(group_names), ncol = 1, strip.position = "right") +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

# ROH count against size labelled
png(snakemake@output[[2]], width = 7, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "G") +
    geom_smooth(method = "lm", se = FALSE, color = "gray28", size = 0.5, alpha = 0.8) +
    geom_label_repel(
        data = no_wolves %>% filter(type == "imputed"),
        aes(x = ROH_tol_2, y = n, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Total ROH length (Mb)", y = "Total # ROH") +
    labs(fill = "Sample age (ybp)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")), labeller = as_labeller(group_names), ncol = 1, strip.position = "right") +
    theme_bw() +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

# Froh against time per population
cols <- c("#ce4441", "#62929a", "#ffbb44")
png(snakemake@output[[3]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(no_wolves$Age_Mean_KBP)), 0, 1)) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3, box.padding = 3, max.overlaps = Inf)+
    # geom_label_repel(data = no_wolves %>% filter(type=='imputed'),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH]))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")) ~ ., labeller = as_labeller(group_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

# Froh against time per population labelled
cols <- c("#ce4441", "#62929a", "#ffbb44")
png(snakemake@output[[4]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(no_wolves$Age_Mean_KBP)), 0, 1)) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3, box.padding = 3, max.overlaps = Inf)+
    geom_label_repel(
        data = no_wolves %>% filter(type == "imputed"),
        aes(x = Age_Mean_KBP, y = froh, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH]))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")) ~ ., labeller = as_labeller(group_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

#### Mann–Whitney U test to check if Froh is different between the 3 dog pops

# perform the Mann Whitney U test for all three possible combos between pops
merge_froh <- function(pop1, pop2, pop3) {
    a <- wilcox.test(pop1$froh, pop2$froh)[c(1, 3)]
    aname <- paste(deparse(substitute(pop1)), deparse(substitute(pop2)), sep = ",")
    b <- wilcox.test(pop1$froh, pop3$froh)[c(1, 3)]
    bname <- paste(deparse(substitute(pop1)), deparse(substitute(pop3)), sep = ",")
    c <- wilcox.test(pop2$froh, pop3$froh)[c(1, 3)]
    cname <- paste(deparse(substitute(pop2)), deparse(substitute(pop3)), sep = ",")
    all <- rbind(a, b, c)
    all <- as.data.frame(all)
    all$group <- NA
    all$group[1] <- aname
    all$group[2] <- bname
    all$group[3] <- cname
    return(all)
}

# only on ancient samples
near_east_ancient <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & type == "imputed")
europe_ancient <- no_wolves %>% filter(Meta.Population == "European_Dogs" & type == "imputed")
arctic_ancient <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & type == "imputed")

df_froh_1 <- merge_froh(near_east_ancient, europe_ancient, arctic_ancient)

# on modern and ancient samples
near_east_all <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs")
europe_all <- no_wolves %>% filter(Meta.Population == "European_Dogs")
arctic_all <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs")

df_froh_2 <- merge_froh(near_east_all, europe_all, arctic_all)


# only modern samples
near_east_modern <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & type == "modern")
europe_modern <- no_wolves %>% filter(Meta.Population == "European_Dogs" & type == "modern")
arctic_modern <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & type == "modern")

df_froh_3 <- merge_froh(near_east_modern, europe_modern, arctic_modern)


# test on samples pre and post 1500 years ago to see if there is a difference (between pops and within pops)
near_east_early <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & Age_Mean_BP < 1500)
europe_early <- no_wolves %>% filter(Meta.Population == "European_Dogs" & Age_Mean_BP < 1500)
arctic_early <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & Age_Mean_BP < 1500)

df_froh_4 <- merge_froh(near_east_early, europe_early, arctic_early)

near_east_late <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & Age_Mean_BP >= 1500)
europe_late <- no_wolves %>% filter(Meta.Population == "European_Dogs" & Age_Mean_BP >= 1500)
arctic_late <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & Age_Mean_BP >= 1500)

df_froh_5 <- merge_froh(near_east_late, europe_late, arctic_late)


# perform the Mann Whitney U test for all three possible combos within pops
merge_froh_within <- function(pop1, pop2) {
    a <- wilcox.test(pop1$froh, pop2$froh)[c(1, 3)]
    aname <- paste(deparse(substitute(pop1)), deparse(substitute(pop2)), sep = ",")
    a <- as.data.frame(a)
    a$group <- NA
    a$group[1] <- aname
    return(a)
}

df_froh_6 <- merge_froh_within(near_east_early, near_east_late)
df_froh_7 <- merge_froh_within(europe_early, europe_late)
df_froh_8 <- merge_froh_within(arctic_early, arctic_late)

df_froh_9 <- merge_froh_within(near_east_ancient, near_east_modern)
df_froh_10 <- merge_froh_within(europe_ancient, europe_modern)
df_froh_11 <- merge_froh_within(arctic_ancient, arctic_modern)

# summarize results into same file:
all_froh <- rbind(
    df_froh_1, df_froh_2, df_froh_3, df_froh_4, df_froh_5,
    df_froh_6, df_froh_7, df_froh_8, df_froh_9, df_froh_10, df_froh_11
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
# Dogs only
#############

my_breaks <- c(0, 100, 200, 500, 1000, 2000, 5000, 10000)
my_labels <- c("0", "100", "200", "500", "1,000", "2,000", "5,000", "10,000")

# filter out wolves and unwanted dog populations:
no_wolves <- ab_long %>% filter(Dog_PCA != "Wolves" & Dog_PCA != "Americas_Dogs" & Dog_PCA != "East_Asian_Dogs")

# make port au choix arctic:
no_wolves$Dog_PCA[no_wolves$Dog_PCA == "preContact_Dogs"] <- "Arctic_Dogs"
no_wolves$Meta.Population[no_wolves$Meta.Population == "preContact_Dogs"] <- "Arctic_Dogs"
no_wolves$Meta.Population[no_wolves$Meta.Population == "American_European_Dogs"] <- "European_Dogs"

no_wolves$Age_Mean_BP <- as.numeric(no_wolves$Age_Mean_BP)
no_wolves$ROH_tol <- as.numeric(no_wolves$ROH_tol)
no_wolves$n <- as.numeric(no_wolves$n)
no_wolves$froh <- as.numeric(no_wolves$froh)

# re-scale x axis and age
no_wolves <- no_wolves %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )

# make modern and ancient layers for plotting:
df_layer_1 <- no_wolves[no_wolves$type == "modern", ]
df_layer_2 <- no_wolves[no_wolves$type != "modern", ]

group_names <- c(
    "African_NearEast_India_Dogs" = "Near East",
    "Arctic_Dogs" = "Arctic",
    "European_Dogs" = "Europe"
)

# ROH count against size
png(snakemake@output[[5]], width = 7, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "G") +
    geom_smooth(method = "lm", se = FALSE, color = "gray28", size = 0.5, alpha = 0.8) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=ROH_tol_2, y=n, label=Sample),size=3.5, box.padding = 3, max.overlaps = Inf)+
    labs(x = "Total ROH length (Mb) (ROH >= 1.6Mb)", y = "Total # ROH (ROH >= 1.6Mb)") +
    labs(fill = "Sample age (ybp)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")), labeller = as_labeller(group_names), ncol = 1, strip.position = "right") +
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


# ROH count against size labelled
png(snakemake@output[[6]], width = 7, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "G") +
    geom_smooth(method = "lm", se = FALSE, color = "gray28", size = 0.5, alpha = 0.8) +
    geom_label_repel(
        data = no_wolves %>% filter(type == "imputed"),
        aes(x = ROH_tol_2, y = n, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Total ROH length (Mb) (ROH >= 1.6Mb)", y = "Total # ROH (ROH >= 1.6Mb)") +
    labs(fill = "Sample age (ybp)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")), labeller = as_labeller(group_names), ncol = 1, strip.position = "right") +
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

# Froh against time per population
cols <- c("#ce4441", "#62929a", "#ffbb44")
png(snakemake@output[[7]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(no_wolves$Age_Mean_KBP)), 0, 1)) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3, box.padding = 3, max.overlaps = Inf)+
    # geom_label_repel(data = no_wolves %>% filter(type=='imputed'),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH], " (ROH >= 1.6Mb)", sep = ""))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")) ~ ., labeller = as_labeller(group_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

# Froh against time per population labelled
cols <- c("#ce4441", "#62929a", "#ffbb44")
png(snakemake@output[[8]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(no_wolves$Age_Mean_KBP)), 0, 1)) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3, box.padding = 3, max.overlaps = Inf)+
    geom_label_repel(
        data = no_wolves %>% filter(type == "imputed"),
        aes(x = Age_Mean_KBP, y = froh, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH], " (ROH >= 1.6Mb)", sep = ""))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")) ~ ., labeller = as_labeller(group_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

#### Mann–Whitney U test to check if Froh is different between the 3 dog pops

# perform the Mann Whitney U test for all three possible combos between pops
merge_froh <- function(pop1, pop2, pop3) {
    a <- wilcox.test(pop1$froh, pop2$froh)[c(1, 3)]
    aname <- paste(deparse(substitute(pop1)), deparse(substitute(pop2)), sep = ",")
    b <- wilcox.test(pop1$froh, pop3$froh)[c(1, 3)]
    bname <- paste(deparse(substitute(pop1)), deparse(substitute(pop3)), sep = ",")
    c <- wilcox.test(pop2$froh, pop3$froh)[c(1, 3)]
    cname <- paste(deparse(substitute(pop2)), deparse(substitute(pop3)), sep = ",")
    all <- rbind(a, b, c)
    all <- as.data.frame(all)
    all$group <- NA
    all$group[1] <- aname
    all$group[2] <- bname
    all$group[3] <- cname
    return(all)
}

# only on ancient samples
near_east_ancient <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & type == "imputed")
europe_ancient <- no_wolves %>% filter(Meta.Population == "European_Dogs" & type == "imputed")
arctic_ancient <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & type == "imputed")

df_froh_1 <- merge_froh(near_east_ancient, europe_ancient, arctic_ancient)

# on modern and ancient samples
near_east_all <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs")
europe_all <- no_wolves %>% filter(Meta.Population == "European_Dogs")
arctic_all <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs")

df_froh_2 <- merge_froh(near_east_all, europe_all, arctic_all)


# only modern samples
near_east_modern <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & type == "modern")
europe_modern <- no_wolves %>% filter(Meta.Population == "European_Dogs" & type == "modern")
arctic_modern <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & type == "modern")

df_froh_3 <- merge_froh(near_east_modern, europe_modern, arctic_modern)


# test on samples pre and post 1500 years ago to see if there is a difference (between pops and within pops)
near_east_early <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & Age_Mean_BP < 1500)
europe_early <- no_wolves %>% filter(Meta.Population == "European_Dogs" & Age_Mean_BP < 1500)
arctic_early <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & Age_Mean_BP < 1500)

df_froh_4 <- merge_froh(near_east_early, europe_early, arctic_early)

near_east_late <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & Age_Mean_BP >= 1500)
europe_late <- no_wolves %>% filter(Meta.Population == "European_Dogs" & Age_Mean_BP >= 1500)
arctic_late <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & Age_Mean_BP >= 1500)

df_froh_5 <- merge_froh(near_east_late, europe_late, arctic_late)


# perform the Mann Whitney U test for all three possible combos within pops
merge_froh_within <- function(pop1, pop2) {
    a <- wilcox.test(pop1$froh, pop2$froh)[c(1, 3)]
    aname <- paste(deparse(substitute(pop1)), deparse(substitute(pop2)), sep = ",")
    a <- as.data.frame(a)
    a$group <- NA
    a$group[1] <- aname
    return(a)
}

df_froh_6 <- merge_froh_within(near_east_early, near_east_late)
df_froh_7 <- merge_froh_within(europe_early, europe_late)
df_froh_8 <- merge_froh_within(arctic_early, arctic_late)

df_froh_9 <- merge_froh_within(near_east_ancient, near_east_modern)
df_froh_10 <- merge_froh_within(europe_ancient, europe_modern)
df_froh_11 <- merge_froh_within(arctic_ancient, arctic_modern)

# summarize results into same file:
all_froh_long <- rbind(
    df_froh_1, df_froh_2, df_froh_3, df_froh_4, df_froh_5,
    df_froh_6, df_froh_7, df_froh_8, df_froh_9, df_froh_10, df_froh_11
)



##### ######### #### #####
##### Short ROHS <1.6Mb #####
##### ######### #### #####

# sum total long ROH count per sample:

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
# Dogs only
#############

my_breaks <- c(0, 100, 200, 500, 1000, 2000, 5000, 10000)
my_labels <- c("0", "100", "200", "500", "1,000", "2,000", "5,000", "10,000")

# filter out wolves and unwanted dog populations:
no_wolves <- ab_short %>% filter(Dog_PCA != "Wolves" & Dog_PCA != "Americas_Dogs" & Dog_PCA != "East_Asian_Dogs")

# make port au choix arctic:
no_wolves$Dog_PCA[no_wolves$Dog_PCA == "preContact_Dogs"] <- "Arctic_Dogs"
no_wolves$Meta.Population[no_wolves$Meta.Population == "preContact_Dogs"] <- "Arctic_Dogs"
no_wolves$Meta.Population[no_wolves$Meta.Population == "American_European_Dogs"] <- "European_Dogs"

no_wolves$Age_Mean_BP <- as.numeric(no_wolves$Age_Mean_BP)
no_wolves$ROH_tol <- as.numeric(no_wolves$ROH_tol)
no_wolves$n <- as.numeric(no_wolves$n)
no_wolves$froh <- as.numeric(no_wolves$froh)

# re-scale x axis and age
no_wolves <- no_wolves %>%
    mutate(
        ROH_tol_2 = ROH_tol / 1e+3,
        Age_Mean_KBP = (Age_Mean_BP / 1000) * (-1),
        Age_Mean_KBP_2 = (Age_Mean_BP / 1000)
    )

# make modern and ancient layers for plotting:
df_layer_1 <- no_wolves[no_wolves$type == "modern", ]
df_layer_2 <- no_wolves[no_wolves$type != "modern", ]

group_names <- c(
    "African_NearEast_India_Dogs" = "Near East",
    "Arctic_Dogs" = "Arctic",
    "European_Dogs" = "Europe"
)

# ROH count against size
png(snakemake@output[[9]], width = 7, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "G") +
    geom_smooth(method = "lm", se = FALSE, color = "gray28", size = 0.5, alpha = 0.8) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=ROH_tol_2, y=n, label=Sample),size=3.5, box.padding = 3, max.overlaps = Inf)+
    labs(x = "Total ROH length (Mb) (ROH < 1.6Mb)", y = "Total # ROH (ROH < 1.6Mb)") +
    labs(fill = "Sample age (ybp)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")), labeller = as_labeller(group_names), ncol = 1, strip.position = "right") +
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

# ROH count against size labelled
png(snakemake@output[[10]], width = 7, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = ROH_tol_2, y = n)) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Age_Mean_BP), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_viridis_c(trans = "log", breaks = my_breaks, labels = my_labels, option = "G") +
    geom_smooth(method = "lm", se = FALSE, color = "gray28", size = 0.5, alpha = 0.8) +
    geom_label_repel(
        data = no_wolves %>% filter(type == "imputed"),
        aes(x = ROH_tol_2, y = n, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    labs(x = "Total ROH length (Mb) (ROH < 1.6Mb)", y = "Total # ROH (ROH < 1.6Mb)") +
    labs(fill = "Sample age (ybp)") +
    facet_wrap(. ~ factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")), labeller = as_labeller(group_names), ncol = 1, strip.position = "right") +
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

# Froh against time per population
cols <- c("#ce4441", "#62929a", "#ffbb44")
png(snakemake@output[[11]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(no_wolves$Age_Mean_KBP)), 0, 1)) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3, box.padding = 3, max.overlaps = Inf)+
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH], " (ROH < 1.6Mb)", sep = ""))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")) ~ ., labeller = as_labeller(group_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()

# Froh against time per population labelled
cols <- c("#ce4441", "#62929a", "#ffbb44")
png(snakemake@output[[12]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
options(scipen = 999)
ggplot(no_wolves, aes(x = Age_Mean_KBP, y = froh)) +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span = 1, aes(fill = Meta.Population), colour = "black", size = 0.5) +
    geom_point(data = df_layer_1, size = 3, alpha = 0.6, colour = "grey") +
    geom_point(data = df_layer_2, aes(fill = Meta.Population), size = 4, shape = 21, alpha = 0.8) +
    scale_fill_manual(values = cols) +
    scale_x_continuous(breaks = seq(round(min(no_wolves$Age_Mean_KBP)), 0, 1)) +
    geom_label_repel(
        data = no_wolves %>% filter(type == "imputed"),
        aes(x = Age_Mean_KBP, y = froh, label = Sample), size = 2, box.padding = 1, max.overlaps = Inf
    ) +
    # geom_label_repel(data = no_wolves %>% filter(Sample %in% chosen_dogs),
    #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3, box.padding = 3, max.overlaps = Inf)+
    labs(x = "Time (kya)", y = expression(paste(italic("F")[ROH], " (ROH < 1.6Mb)", sep = ""))) +
    theme_bw() +
    theme(legend.position = "none") +
    facet_grid(factor(Meta.Population, levels = c("African_NearEast_India_Dogs", "European_Dogs", "Arctic_Dogs")) ~ ., labeller = as_labeller(group_names)) +
    theme(
        strip.background = element_rect(fill = "gray28"),
        strip.text = element_text(colour = "white", size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )
dev.off()


#### Mann–Whitney U test to check if Froh is different between the 3 dog pops

# perform the Mann Whitney U test for all three possible combos between pops
merge_froh <- function(pop1, pop2, pop3) {
    a <- wilcox.test(pop1$froh, pop2$froh)[c(1, 3)]
    aname <- paste(deparse(substitute(pop1)), deparse(substitute(pop2)), sep = ",")
    b <- wilcox.test(pop1$froh, pop3$froh)[c(1, 3)]
    bname <- paste(deparse(substitute(pop1)), deparse(substitute(pop3)), sep = ",")
    c <- wilcox.test(pop2$froh, pop3$froh)[c(1, 3)]
    cname <- paste(deparse(substitute(pop2)), deparse(substitute(pop3)), sep = ",")
    all <- rbind(a, b, c)
    all <- as.data.frame(all)
    all$group <- NA
    all$group[1] <- aname
    all$group[2] <- bname
    all$group[3] <- cname
    return(all)
}

# only on ancient samples
near_east_ancient <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & type == "imputed")
europe_ancient <- no_wolves %>% filter(Meta.Population == "European_Dogs" & type == "imputed")
arctic_ancient <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & type == "imputed")

df_froh_1 <- merge_froh(near_east_ancient, europe_ancient, arctic_ancient)

# on modern and ancient samples
near_east_all <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs")
europe_all <- no_wolves %>% filter(Meta.Population == "European_Dogs")
arctic_all <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs")

df_froh_2 <- merge_froh(near_east_all, europe_all, arctic_all)


# only modern samples
near_east_modern <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & type == "modern")
europe_modern <- no_wolves %>% filter(Meta.Population == "European_Dogs" & type == "modern")
arctic_modern <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & type == "modern")

df_froh_3 <- merge_froh(near_east_modern, europe_modern, arctic_modern)


# test on samples pre and post 1500 years ago to see if there is a difference (between pops and within pops)
near_east_early <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & Age_Mean_BP < 1500)
europe_early <- no_wolves %>% filter(Meta.Population == "European_Dogs" & Age_Mean_BP < 1500)
arctic_early <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & Age_Mean_BP < 1500)

df_froh_4 <- merge_froh(near_east_early, europe_early, arctic_early)

near_east_late <- no_wolves %>% filter(Meta.Population == "African_NearEast_India_Dogs" & Age_Mean_BP >= 1500)
europe_late <- no_wolves %>% filter(Meta.Population == "European_Dogs" & Age_Mean_BP >= 1500)
arctic_late <- no_wolves %>% filter(Meta.Population == "Arctic_Dogs" & Age_Mean_BP >= 1500)

df_froh_5 <- merge_froh(near_east_late, europe_late, arctic_late)


# perform the Mann Whitney U test for all three possible combos within pops
merge_froh_within <- function(pop1, pop2) {
    a <- wilcox.test(pop1$froh, pop2$froh)[c(1, 3)]
    aname <- paste(deparse(substitute(pop1)), deparse(substitute(pop2)), sep = ",")
    a <- as.data.frame(a)
    a$group <- NA
    a$group[1] <- aname
    return(a)
}

df_froh_6 <- merge_froh_within(near_east_early, near_east_late)
df_froh_7 <- merge_froh_within(europe_early, europe_late)
df_froh_8 <- merge_froh_within(arctic_early, arctic_late)

df_froh_9 <- merge_froh_within(near_east_ancient, near_east_modern)
df_froh_10 <- merge_froh_within(europe_ancient, europe_modern)
df_froh_11 <- merge_froh_within(arctic_ancient, arctic_modern)

# summarize results into same file:
all_froh_short <- rbind(
    df_froh_1, df_froh_2, df_froh_3, df_froh_4, df_froh_5,
    df_froh_6, df_froh_7, df_froh_8, df_froh_9, df_froh_10, df_froh_11
)



# put all into one file
all_froh$type <- "all_ROHs"
all_froh_long$type <- "long_ROHs"
all_froh_short$type <- "short_ROHs"

froh_tests <- rbind(all_froh, all_froh_long, all_froh_short)
froh_tests <- apply(froh_tests, 2, as.character)

write.table(froh_tests, file = snakemake@output[[13]], quote = FALSE, sep = "\t", row.names = FALSE)



###############################
#### Map with Froh for dogs ###
###############################
library(readr)
library(maps)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(scales)
library(viridis)

# Load dog dataset:
options(scipen = 999)

# all samples, no cutoff
samples <- read.delim(snakemake@input[[1]])
# samples <- read.delim('~/Downloads/roh_check/Dog_Wolf_aDNA_WG-Master.tsv')


sites <- samples %>% select("name_haplo_VCF", "Lat", "Long", "Meta.Population", "Age_Mean_BP", "Other_ID")
colnames(sites)[1] <- "Sample"

# fix Tumat2 and Wolfhead
ab$Sample[ab$Sample == "Wolf_head_IN18-016"] <- "WolfHead"
ab$Sample[ab$Sample == "Tumat2"] <- "Tumat"

final_map <- left_join(ab, sites, "Sample")
final_map_imputed_dogs <- final_map %>% filter(Species == "familiaris" & type == "imputed")

# give approximate values for Lithuanian/Latvian FIX THIS LATER!!!
final_map_imputed_dogs$Lat[is.na(final_map_imputed_dogs$Lat)] <- 54.901821
final_map_imputed_dogs$Long[is.na(final_map_imputed_dogs$Long)] <- 23.924860

# Start plotting:
worldmap <- map_data("world")
set.seed(24)

# add noise to the data to avoid overlapping points and add size variance for date_bp
sites_2 <- final_map_imputed_dogs
sites_2$Lat <- sites_2$Lat + rnorm(length(sites_2$Lat), 0, 0.7)
sites_2$Long <- sites_2$Long + rnorm(length(sites_2$Long), 0, 0.7)

sites_2$Long[sites_2$Sample == "TRF.02.49"] <- 189

# add column with name and age
sites_2$name_age <- paste(sites_2$Sample, " (", sites_2$Age_Mean_BP.x, " BP)", sep = "")

png(snakemake@output[[14]], width = 8, height = 8, units = "in", res = 200, pointsize = 4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab = 2)
ggplot() +
    geom_polygon(
        data = worldmap,
        aes(x = long, y = lat, group = group),
        fill = "gray85", color = "gray85"
    ) +
    coord_fixed(ratio = 1, xlim = c(-5, 190), ylim = c(20, 80)) +
    geom_point(
        data = sites_2,
        aes(
            x = as.numeric(Long),
            y = as.numeric(Lat), fill = froh
        ), alpha = .7, size = 3.5, shape = 21
    ) +
    scale_fill_viridis_c(option = "A", trans = "log") +
    geom_label_repel(
        data = sites_2 %>% filter(froh > 0.1),
        aes(x = as.numeric(Long), y = as.numeric(Lat), label = name_age), size = 3, box.padding = 1, max.overlaps = Inf
    ) +
    theme(legend.position = "right", legend.direction = "vertical") +
    theme(title = element_text(size = 12)) +
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        legend.position = c(0.88, 0.23),
        legend.direction = "vertical"
    ) +
    labs(fill = bquote("" ~ F[ROH] * ""))
dev.off()


# Box plots for dogs, ancient vs present-day:




### Export tables with ROH data per sample:

# merge all rohs, long and short:
ab$category <- "all_rohs"
ab_long$category <- "long_rohs"
ab_short$category <- "short_rohs"
all_roh_results <- rbind(ab, ab_long, ab_short)

write.table(all_roh_results, file = snakemake@output[[15]], quote = FALSE, sep = "\t", row.names = FALSE)
