#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(ggpubr)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

list_of_concordance_phased <- args[1]
list_of_validation <- args[2]
list_of_phased <- args[3]
rohan <- read.delim(args[4])
chr_length <- read.delim(args[5], header = FALSE)
name <- args[6]
chrom <- args[7]
cov_hc <- args[8]

# name <- 'CGG32'
# chrom <- 'chr1'
# cov_hc <- 14.5

name_title <- args[9]
name_title_2 <- gsub("_", " ", name_title)

# name_title <- 'Pleistocene_wolf'
# name_title_2 <- gsub('_',' ',name_title)

site_type <- args[10]
# site_type <- 'Transversions+transitions'

a <- str_split(list_of_concordance_phased, pattern = ",")
b <- str_split(list_of_validation, pattern = ",")
f <- str_split(list_of_phased, pattern = ",")

a[[1]][7] <- NA
a[[1]][7] <- b[[1]]
a[[1]][8] <- NA
a[[1]][8] <- f[[1]]

d <- c()

for (i in 1:length(a[[1]])) {
    tmp <- read.csv(a[[1]][i], sep = "")
    d <- rbind(d, tmp)
}

# import imputed and HC ROH data
# b1 <- read.csv("~/Downloads/phased.CGG32_chr1_0.05x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b2 <- read.csv("~/Downloads/phased.CGG32_chr1_0.1x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b3 <- read.csv("~/Downloads/phased.CGG32_chr1_0.2x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b4 <- read.csv("~/Downloads/phased.CGG32_chr1_0.5x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b5 <- read.csv("~/Downloads/phased.CGG32_chr1_1x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b6 <- read.csv("~/Downloads/phased.CGG32_chr1_2x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b7 <- read.csv("~/Downloads/phased.CGG32_chr1_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b8 <- read.csv("~/Downloads/CGG32_chr1_validation_filt_qual_dp_ab_all_sites_hom_win_het_1_plink-temp.hom", sep="")
#
# d <- rbind(b1, b2, b3, b4, b5, b6, b7, b8)

# import ROHan results
# rohan <- read.delim("~/Downloads/ROHan_4e-5_KatiaPipeline.txt")

# chr length
# chr_length <- read.delim("~/Downloads/CanFam31_chr1_size.genome", header=FALSE)
size_chr <- chr_length[1, 2]

# put same cov for imputed and genotyped HC:
d$name[d$cov == "HC_imputed"] <- "HC"
d$name[d$cov == "HC_genotyped"] <- "HC"
d$name <- ifelse(is.na(d$name), d$cov, d$name)

# replace name
d$cov <- gsub("x", "x imputed", d$cov)
d$cov <- gsub("HC_genotyped", paste("HC (", cov_hc, "x)", sep = ""), d$cov)
d$cov <- gsub("HC_imputed", paste("HC imputed (", cov_hc, "x)", sep = ""), d$cov)

name_HC_imputed <- paste("HC imputed (", cov_hc, "x)", sep = "")
name_HC <- paste("HC (", cov_hc, "x)", sep = "")


##### ADD ROHAN DATA
cov_list <- c("0.05x", "0.1x", "0.2x", "0.5x", "1x", "2x", "ori")

# play with min, mid and max to see how it changes
rohan$FID[rohan$FID == "SOTN01"] <- "SOTN01_merged"
rohan$IID[rohan$IID == "SOTN01"] <- "SOTN01_merged"

b9 <- rohan %>%
    filter(hEst == "mid") %>%
    select(FID, IID, CHR, POS1, POS2, KB, NSNP, cov) %>%
    filter(FID == name, CHR == 1)

# add row for missing ROHan coverages (where no ROHs were found, but still need to plot them)
for (i in cov_list) {
    b9[nrow(b9) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, i)
}

# format coverage and name columns appropriately
b9$name <- b9$cov
b9$cov <- paste(b9$cov, "ROHan", sep = " ")
b9$cov[b9$cov == "ori ROHan"] <- paste("HC ROHan ", "(", cov_hc, "x)", sep = "")
b9$name[b9$name == "ori"] <- paste("HC")


# merge ROHan with imputed and HC
d2 <- rbind.fill(d, b9)

HC_category <- paste("HC ROHan ", "(", cov_hc, "x)", sep = "")

# reorder names
d2$cov <- factor(d2$cov, levels = c("0.05x ROHan", "0.05x imputed", "0.1x ROHan", "0.1x imputed", "0.2x ROHan", "0.2x imputed", "0.5x ROHan", "0.5x imputed", "1x ROHan", "1x imputed", "2x ROHan", "2x imputed", HC_category, name_HC_imputed, name_HC))

d2$POS1 <- as.numeric(d2$POS1)
d2$POS2 <- as.numeric(d2$POS2)
d2$KB <- as.numeric(d2$KB)

# rescale x axis
new_sub_final <- d2 %>%
    mutate(
        POS1 = POS1 / 1e+6,
        POS2 = POS2 / 1e+6,
        MB = KB / 1000
    )

size_chr_mb <- size_chr / 1e+6

# name_title_final <- paste(name_title_2, ' - ', name,sep='')
name_title_final <- paste(name_title_2, " - ", name, " (", site_type, ")", sep = "")

# define label based on sample name
new_sub_final$label <- ifelse(new_sub_final$FID == "NGDG", "a",
    ifelse(new_sub_final$FID == "SOTN01_merged", "b",
        ifelse(new_sub_final$FID == "PortauChoix", "c",
            ifelse(new_sub_final$FID == "TRF.02.53", "d",
                ifelse(new_sub_final$FID == "TRF.05.05", "e",
                    ifelse(new_sub_final$FID == "FAMICHN00012", "f",
                        ifelse(new_sub_final$FID == "FAMINGR00004", "g",
                            ifelse(new_sub_final$FID == "CGG32", "h",
                                ifelse(new_sub_final$FID == "CGG33", "i",
                                    ifelse(new_sub_final$FID == "WolfHead", "j", "NA")
                                )
                            )
                        )
                    )
                )
            )
        )
    )
)

unique_label <- unique(new_sub_final$label)

# plot
# png(args[9], width=16, height=4, units='in', res=250, pointsize=4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
options(scipen = 10000)
a1 <- ggplot(data = new_sub_final) +
    geom_hline(aes(yintercept = cov), color = "#d8dee9", size = 0.4) +
    geom_segment(aes(y = cov, yend = cov, x = POS1, xend = POS2, colour = name), linewidth = 10) +
    scale_colour_viridis(discrete = TRUE, option = "D") +
    scale_x_continuous(breaks = seq(0, size_chr_mb, 10)) +
    # theme_void()+
    theme_bw() +
    xlab(paste("Genomic position ", chrom, " (Mb)", sep = "")) +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.45, size = 18, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "gray60")
    ) +
    labs(colour = "Coverage") +
    ggtitle(name_title_final)
a1
# dev.off


#############################################
#                                           #
#             ACCURACY GRAPHS.              #
#                                           #
#############################################

list_of_concordance_phased <- args[11]
list_of_validation <- args[12]
list_of_phased <- args[13]
ref_fasta_chr_size <- read.table(args[14], quote = "\"", comment.char = "")
# ref_fasta_chr_size <- read.table("~/Downloads/CanFam31_allchrom_size.genome", quote="\"", comment.char="")

aa <- str_split(list_of_concordance_phased, pattern = ",")
bb <- str_split(list_of_validation, pattern = ",")
ff <- str_split(list_of_phased, pattern = ",")


a <- read.csv(aa[[1]][1], sep = "")
b <- read.csv(aa[[1]][2], sep = "")
c <- read.csv(aa[[1]][3], sep = "")
d <- read.csv(aa[[1]][4], sep = "")
e <- read.csv(aa[[1]][5], sep = "")
f <- read.csv(aa[[1]][6], sep = "")

g <- read.csv(ff[[1]][1], sep = "")

validation <- read.csv(bb[[1]][1], sep = "")


# a <- read.csv("~/Downloads/phased.CGG32_allchrom_0.05x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# b <- read.csv("~/Downloads/phased.CGG32_allchrom_0.1x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# c <- read.csv("~/Downloads/phased.CGG32_allchrom_0.2x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# d <- read.csv("~/Downloads/phased.CGG32_allchrom_0.5x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# e <- read.csv("~/Downloads/phased.CGG32_allchrom_1x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# f <- read.csv("~/Downloads/phased.CGG32_allchrom_2x_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# g <- read.csv("~/Downloads/phased.CGG32_allchrom_INFO_0.8_MAF_0.01_all_sites_hom_win_het_1-temp.hom", sep="")
# validation <- read.csv("~/Downloads/CGG32_allchrom_validation_filt_qual_dp_ab_all_sites_hom_win_het_1_plink-temp.hom", sep="")

# import ROHan results
# rohan <- read.delim("~/Downloads/ROHan_4e-5_KatiaPipeline.txt")

# split ROHAn per sample and per coverage for all chromosomes
rohan2 <- rohan %>%
    filter(hEst == "mid" & FID == name) %>%
    select(FID, IID, CHR, POS1, POS2, KB, NSNP, cov)

a_roh <- rohan2 %>% filter(cov == "0.05x")
b_roh <- rohan2 %>% filter(cov == "0.1x")
c_roh <- rohan2 %>% filter(cov == "0.2x")
d_roh <- rohan2 %>% filter(cov == "0.5x")
e_roh <- rohan2 %>% filter(cov == "1x")
f_roh <- rohan2 %>% filter(cov == "2x")
g_roh <- rohan2 %>% filter(cov == "ori")

# enter empty row for cases where no ROH was found
a_roh[nrow(a_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "0.05x")
b_roh[nrow(b_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "0.1x")
c_roh[nrow(c_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "0.2x")
d_roh[nrow(d_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "0.5x")
e_roh[nrow(e_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "1x")
f_roh[nrow(f_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "2x")
g_roh[nrow(g_roh) + 1, ] <- c(name, name, "1", 0, 0, 0, 0, "ori")


# remove extra line used for plotting (in case no ROH was found) and use GRanges to create approprite format
get_gr <- function(roh) {
    test <- roh %>%
        filter(IID != 0 & IID == name) %>%
        mutate(range = paste(POS1, "-", POS2, sep = ""))
    GRanges(
        seqnames = test$CHR,
        ranges = test$range,
        sample = test$FID,
        cov = test$cov
    )
}

get_val <- get_gr(validation)
get_a <- get_gr(a)
get_b <- get_gr(b)
get_c <- get_gr(c)
get_d <- get_gr(d)
get_e <- get_gr(e)
get_f <- get_gr(f)
get_g <- get_gr(g)

get_a_roh <- get_gr(a_roh)
get_b_roh <- get_gr(b_roh)
get_c_roh <- get_gr(c_roh)
get_d_roh <- get_gr(d_roh)
get_e_roh <- get_gr(e_roh)
get_f_roh <- get_gr(f_roh)
get_g_roh <- get_gr(g_roh)


# find overlaps between the two samples
hits_a <- findOverlaps(get_val, get_a)
hits_b <- findOverlaps(get_val, get_b)
hits_c <- findOverlaps(get_val, get_c)
hits_d <- findOverlaps(get_val, get_d)
hits_e <- findOverlaps(get_val, get_e)
hits_f <- findOverlaps(get_val, get_f)
hits_g <- findOverlaps(get_val, get_g)

hits_a_roh <- findOverlaps(get_val, get_a_roh)
hits_b_roh <- findOverlaps(get_val, get_b_roh)
hits_c_roh <- findOverlaps(get_val, get_c_roh)
hits_d_roh <- findOverlaps(get_val, get_d_roh)
hits_e_roh <- findOverlaps(get_val, get_e_roh)
hits_f_roh <- findOverlaps(get_val, get_f_roh)
hits_g_roh <- findOverlaps(get_val, get_g_roh)



# get overlaps
get_overlap <- function(query, subject, hits) {
    overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
    overlaps$proportion_subject <- width(overlaps) / width(subject[subjectHits(hits)]) # this is the proportion of the subject that is contained within the query
    overlaps$proportion_query <- width(overlaps) / width(query[queryHits(hits)]) # this is the proportion of the query that is overlapping with the subject

    return(overlaps)
}

# ab <- get_overlap(gr, gr1, hits)
a_overlap <- get_overlap(get_val, get_a, hits_a)
b_overlap <- get_overlap(get_val, get_b, hits_b)
c_overlap <- get_overlap(get_val, get_c, hits_c)
d_overlap <- get_overlap(get_val, get_d, hits_d)
e_overlap <- get_overlap(get_val, get_e, hits_e)
f_overlap <- get_overlap(get_val, get_f, hits_f)
g_overlap <- get_overlap(get_val, get_g, hits_g)


a_roh_overlap <- get_overlap(get_val, get_a_roh, hits_a_roh)
b_roh_overlap <- get_overlap(get_val, get_b_roh, hits_b_roh)
c_roh_overlap <- get_overlap(get_val, get_c_roh, hits_c_roh)
d_roh_overlap <- get_overlap(get_val, get_d_roh, hits_d_roh)
e_roh_overlap <- get_overlap(get_val, get_e_roh, hits_e_roh)
f_roh_overlap <- get_overlap(get_val, get_f_roh, hits_f_roh)
g_roh_overlap <- get_overlap(get_val, get_g_roh, hits_g_roh)


############ SEGMENT BASED ESTIMATION ###########

## Now get stats for segment based overlaps:

get_seg_stats <- function(hits, get_phased, get_val, df) {
    # ROHs identified in the phased_target that are also present in the val (so are actually true)
    if (length(unique(queryHits(hits))) != 0) {
        TP <- length(unique(queryHits(hits)))
    } else {
        TP <- 1e-16
    }

    # false positives will be a ROH found in the phased_target that was not present in the val
    if (length(get_phased[!(get_phased %in% get_phased[subjectHits(hits)])]) != 0) {
        FP <- length(get_phased[!(get_phased %in% get_phased[subjectHits(hits)])])
    } else {
        FP <- 1e-16
    }
    # false negatives will be a ROH not found in the phased_target that was present in the val
    if (length(get_val[!(get_val %in% get_val[queryHits(hits)])]) != 0) {
        FN <- length(get_val[!(get_val %in% get_val[queryHits(hits)])])
    } else {
        FN <- 1e-16
    }

    ## Now get F1 stat for segment based overlaps:
    recall_seg <- TP / (TP + FN)
    precision_seg <- TP / (TP + FP)
    F1_seg <- 2 * (precision_seg * recall_seg / (precision_seg + recall_seg))

    # False discovery rate:
    FDR <- FP / (TP + FP)

    # summarize:
    all_stats_segment <- data.frame(
        sample = name,
        # cov=unique(get_phased$cov),
        cov = unique(df$cov),
        TP = TP,
        FP = FP,
        FN = FN,
        sensitivity = recall_seg,
        precision = precision_seg,
        F1 = F1_seg,
        FDR = FDR
    )
}

get_seg_stats_a <- get_seg_stats(hits_a, get_a, get_val, a)
get_seg_stats_b <- get_seg_stats(hits_b, get_b, get_val, b)
get_seg_stats_c <- get_seg_stats(hits_c, get_c, get_val, c)
get_seg_stats_d <- get_seg_stats(hits_d, get_d, get_val, d)
get_seg_stats_e <- get_seg_stats(hits_e, get_e, get_val, e)
get_seg_stats_f <- get_seg_stats(hits_f, get_f, get_val, f)
get_seg_stats_g <- get_seg_stats(hits_g, get_g, get_val, g)

get_seg_stats_a_roh <- get_seg_stats(hits_a_roh, get_a_roh, get_val, a_roh)
get_seg_stats_b_roh <- get_seg_stats(hits_b_roh, get_b_roh, get_val, b_roh)
get_seg_stats_c_roh <- get_seg_stats(hits_c_roh, get_c_roh, get_val, c_roh)
get_seg_stats_d_roh <- get_seg_stats(hits_d_roh, get_d_roh, get_val, d_roh)
get_seg_stats_e_roh <- get_seg_stats(hits_e_roh, get_e_roh, get_val, e_roh)
get_seg_stats_f_roh <- get_seg_stats(hits_f_roh, get_f_roh, get_val, f_roh)
get_seg_stats_g_roh <- get_seg_stats(hits_g_roh, get_g_roh, get_val, g_roh)

all_seg_stats <- rbind(
    get_seg_stats_a,
    get_seg_stats_b,
    get_seg_stats_c,
    get_seg_stats_d,
    get_seg_stats_e,
    get_seg_stats_f,
    get_seg_stats_g
)

all_seg_stats_roh <- rbind(
    get_seg_stats_a_roh,
    get_seg_stats_b_roh,
    get_seg_stats_c_roh,
    get_seg_stats_d_roh,
    get_seg_stats_e_roh,
    get_seg_stats_f_roh,
    get_seg_stats_g_roh
)


all_seg_stats$Data <- "Imputed"
all_seg_stats_roh$Data <- "Non-imputed"

# merge Imputed and ROHAN
all_seg_stats_all <- rbind(all_seg_stats, all_seg_stats_roh)

############ TOTAL LENGTH BASED ESTIMATION ###########
genome_size <- ref_fasta_chr_size[1, 2]

library(mltools)

get_length_stats <- function(overlap, gen_size, get_phased, hits, get_val, df) {
    # TP should be the total ROH of the subject overlapping the validation
    TP_overlaps_sum <- sum(width(overlap))
    if (TP <- TP_overlaps_sum / gen_size != 0) {
        TP <- TP_overlaps_sum / gen_size
    } else {
        TP <- 1e-16
    }

    # FP should be the total ROH of the subject minus the TP:
    FP_overlaps_sum <- sum(width(get_phased)) - sum(width(overlap))
    if (FP <- FP_overlaps_sum / gen_size != 0) {
        FP <- FP_overlaps_sum / gen_size
    } else {
        FP <- 1e-16
    }

    # FN should be the total ROH of the query minus the TP:
    FN_overlaps_sum <- sum(width(get_val)) - sum(width(overlap))
    if (FN <- FN_overlaps_sum / gen_size != 0) {
        FN <- FN_overlaps_sum / gen_size
    } else {
        FN <- 1e-16
    }

    # TN should be the total genome size minus the TP, FP, FN:
    TN_sum <- gen_size - TP_overlaps_sum - FP_overlaps_sum - FN_overlaps_sum
    TN <- TN_sum / gen_size

    recall_length <- TP / (TP + FN)
    precision_length <- TP / (TP + FP)
    F1_length <- 2 * (precision_length * recall_length / (precision_length + recall_length))

    mcc <- mcc(TP = TP, FP = FP, TN = TN, FN = FN)

    mcc_norm <- (mcc + 1) / 2 # IF 0 until <0.5 it's worse than random, 0.5 random and >0.5 is better estimations

    # accuracy (worst value: 0; best value: 1):
    accuracy <- (TP + TN) / (TP + TN + FP + FN)

    # Specificity
    specificity <- TN / (TN + FP)

    # False discovery rate:
    FDR <- FP / (TP + FP)

    all_stats_length <- data.frame(
        sample = name,
        # cov=unique(get_phased$cov),
        cov = unique(df$cov),
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        sensitivity = recall_length,
        precision = precision_length,
        F1 = F1_length,
        FDR = FDR,
        mcc = mcc,
        mcc_norm = mcc_norm,
        accuracy = accuracy,
        specificity = specificity
    )
}

get_length_stats_a <- get_length_stats(a_overlap, genome_size, get_a, hits_a, get_val, a)
get_length_stats_b <- get_length_stats(b_overlap, genome_size, get_b, hits_b, get_val, b)
get_length_stats_c <- get_length_stats(c_overlap, genome_size, get_c, hits_c, get_val, c)
get_length_stats_d <- get_length_stats(d_overlap, genome_size, get_d, hits_d, get_val, d)
get_length_stats_e <- get_length_stats(e_overlap, genome_size, get_e, hits_e, get_val, e)
get_length_stats_f <- get_length_stats(f_overlap, genome_size, get_f, hits_f, get_val, f)
get_length_stats_g <- get_length_stats(g_overlap, genome_size, get_g, hits_g, get_val, g)

get_length_stats_a_roh <- get_length_stats(a_roh_overlap, genome_size, get_a_roh, hits_a_roh, get_val, a_roh)
get_length_stats_b_roh <- get_length_stats(b_roh_overlap, genome_size, get_b_roh, hits_b_roh, get_val, b_roh)
get_length_stats_c_roh <- get_length_stats(c_roh_overlap, genome_size, get_c_roh, hits_c_roh, get_val, c_roh)
get_length_stats_d_roh <- get_length_stats(d_roh_overlap, genome_size, get_d_roh, hits_d_roh, get_val, d_roh)
get_length_stats_e_roh <- get_length_stats(e_roh_overlap, genome_size, get_e_roh, hits_e_roh, get_val, e_roh)
get_length_stats_f_roh <- get_length_stats(f_roh_overlap, genome_size, get_f_roh, hits_f_roh, get_val, f_roh)
get_length_stats_g_roh <- get_length_stats(g_roh_overlap, genome_size, get_g_roh, hits_g_roh, get_val, g_roh)

all_length_stats <- rbind(
    get_length_stats_a,
    get_length_stats_b,
    get_length_stats_c,
    get_length_stats_d,
    get_length_stats_e,
    get_length_stats_f,
    get_length_stats_g
)

all_length_stats_roh <- rbind(
    get_length_stats_a_roh,
    get_length_stats_b_roh,
    get_length_stats_c_roh,
    get_length_stats_d_roh,
    get_length_stats_e_roh,
    get_length_stats_f_roh,
    get_length_stats_g_roh
)

all_length_stats$Data <- "Imputed"
all_length_stats_roh$Data <- "Non-imputed"

# merge Imputed and ROHAN
all_length_stats_all <- rbind(all_length_stats, all_length_stats_roh)

# export both files:
write.table(all_seg_stats_all, file = args[15], quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(all_length_stats_all, file = args[16], quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


# merge both df together (segment and length) based on common columns (mcc is only for length)
all_seg_stats_all$type <- "Segments"
all_seg_stats_all$mcc_norm <- NA
all_length_stats_all$type <- "Length"

common_cols <- intersect(colnames(all_seg_stats_all), colnames(all_length_stats_all))
all_stats <- rbind(
    subset(all_seg_stats_all, select = common_cols),
    subset(all_length_stats_all, select = common_cols)
)

# melt df
all_stats_final <- all_stats %>%
    select(sample, cov, F1, mcc_norm, type, Data) %>%
    gather(metric, values, -sample, -cov, -type, -Data)

all_stats_final <- na.omit(all_stats_final)



## plotting

name_HC_imputed <- paste("HC (", cov_hc, "x)", sep = "")
all_stats_final$cov <- gsub("HC_imputed", name_HC_imputed, all_stats_final$cov)
all_stats_final$metric <- gsub("mcc_norm", "nMCC", all_stats_final$metric)

# change ori to HC
all_stats_final$cov <- gsub("ori", name_HC_imputed, all_stats_final$cov)

# name_title_2 <- gsub('_',' ',name_title)
# name_title_final <- paste(name_title_2, ' - ', name, " (", site_type, ")", sep='')


# F1 and nMCC (segment and length together)
# png(args[10], width=9, height=6, units='in', res=200, pointsize=4)
# par(
#   mar      = c(5, 5, 2, 2),
#   xaxs     = "i",
#   yaxs     = "i",
#   cex.axis = 2,
#   cex.lab  = 2)
# b1 <- ggplot(all_stats_final, aes(x=cov, y=values, colour=Data, shape=type))+
#   geom_line(aes(group = interaction(type, metric, Data), linetype=metric), linewidth = 1.2)+
#   geom_point(size=4)+
#   scale_color_manual(values = c("steelblue","orange"))+
#   ylim(0,1)+
#   labs(x = 'Coverage', colour='Data', linetype= 'Metric', shape='Count type') +
#   #ggtitle(name_title_final) +
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 30, size = 16, vjust = 0.5),
#         axis.text.y=element_text(size=16),
#         axis.title.x=element_text(size=18),
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16),
#         axis.title.y = element_blank(),
#         plot.title = element_text(hjust = 0.5, size=18),
#         panel.grid.minor = element_blank(),
#         legend.key.width= unit(1.5, 'cm'))
# b1
# dev.off()



### Looking into other accuracy metrics for length:

# melt all_stats_final df
all_stats_final <- all_length_stats_all %>%
    select(sample, cov, type, FDR, specificity, sensitivity, Data) %>%
    gather(metric, values, -sample, -cov, -type, -Data)

all_stats_final <- na.omit(all_stats_final)
all_stats_final$cov <- gsub("HC_imputed", name_HC_imputed, all_stats_final$cov)
all_stats_final$cov <- gsub("ori", name_HC_imputed, all_stats_final$cov)

## plotting FDR, specificity and sensitivity
# png(args[11], width=9, height=6, units='in', res=200, pointsize=4)
# par(
#   mar      = c(5, 5, 2, 2),
#   xaxs     = "i",
#   yaxs     = "i",
#   cex.axis = 2,
#   cex.lab  = 2)
# c1 <- ggplot(all_stats_final, aes(x=cov, y=values, colour=metric, linetype=Data))+
#   geom_line(aes(group = interaction(type, metric, Data)), linewidth = 1.2)+
#   geom_point(size=3)+
#   scale_color_manual(values = c("gold2","darkolivegreen4", "lightpink3"))+
#   ylim(0,1)+
#   labs(x = 'Coverage', colour='Metric') +
#   #ggtitle(name_title_final) +
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 30, size = 16, vjust = 0.5),
#         axis.text.y=element_text(size=16),
#         axis.title.x=element_text(size=18),
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16),
#         axis.title.y = element_blank(),
#         plot.title = element_text(hjust = 0.5, size=18),
#         panel.grid.minor = element_blank(),
#         legend.key.width= unit(1.5, 'cm'))
# c1
# dev.off()


## plotting specificity against sensitivity
all_length_stats_all$cov <- gsub("HC_imputed", name_HC_imputed, all_length_stats_all$cov)
all_length_stats_all$cov <- gsub("ori", name_HC_imputed, all_length_stats_all$cov)

# png(args[12], width=9, height=6, units='in', res=200, pointsize=4)
# par(
#   mar      = c(5, 5, 2, 2),
#   xaxs     = "i",
#   yaxs     = "i",
#   cex.axis = 2,
#   cex.lab  = 2)
# d1 <- ggplot(all_length_stats_all, aes(x=specificity, y=sensitivity, colour=cov, shape=Data))+
#   geom_point(size=4)+
#   scale_colour_viridis(discrete = TRUE, option='D') +
#   ylim(0,1)+
#   xlim(0,1)+
#   labs(x = 'Specificity', y='Sensitivity', colour='Coverage') +
#   #ggtitle(name_title_final) +
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 30, size = 16, vjust = 0.5),
#         axis.text.y=element_text(size=16),
#         axis.title.x=element_text(size=18),
#         axis.title.y=element_text(size=18),
#         legend.text = element_text(size=16),
#         legend.title = element_text(size=16),
#         plot.title = element_text(hjust = 0.5, size=18),
#         panel.grid.minor = element_blank(),
#         legend.key.width= unit(1.5, 'cm'))
# d1
# dev.off()



#### merge all plots per sample together:

# all1 <- ggarrange(a1,
#                   labels = c("a"),
#                   ncol = 1, nrow = 1,
#                   font.label=list(size=20))
# #all1


# all2 <- ggarrange(b1,c1,d1,
#                   labels=c("b","c","d"),
#                   ncol = 3, nrow = 1,
#                   font.label=list(size=20),
#                   vjust=-0.1)
# #all2


# all <- ggarrange(all1, all2,
#                  ncol = 1, nrow = 2,
#                  font.label=list(size=20))
# #all



# png(args[[17]], width=31, height=15, units='in', res=200, pointsize=4)
# #png('~/Downloads/testf.png', width=31, height=15, units='in', res=200, pointsize=4)
# alll <- ggarrange(all1, NULL, all2,
#                   ncol = 1, nrow = 3,
#                   heights = c(0.45, 0.05, 0.5))

# alll
# dev.off()




#######################################
# comparing imputed vc non-imputed ROHs
#######################################


### plot nMCC for imputed and non-imputed ROHs
m1 <- ggplot(all_length_stats_all, aes(x = cov, y = mcc_norm, linetype = Data)) +
    geom_line(aes(group = interaction(Data)), linewidth = 1.2, colour = "lightseagreen") +
    geom_point(size = 3, colour = "lightseagreen") +
    # scale_color_manual(values = c("lightseagreen"))+
    ylim(0, 1) +
    labs(x = "Coverage", y = "nMCC", colour = "Data") +
    # ggtitle(name_title_final) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 30, size = 16, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        # axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.5, "cm")
    ) +
    guides(linetype = guide_legend(override.aes = list(color = "black")))


### plot sensitivity, specificity and FDR for imputed and non-imputed ROHs

met1 <- all_stats_final %>%
    filter(metric == "specificity") %>%
    ggplot(aes(x = cov, y = values, linetype = Data)) +
    geom_line(aes(group = interaction(type, Data)), linewidth = 1.2, colour = "lightpink3") +
    geom_point(size = 3, colour = "lightpink3") +
    # scale_color_manual(values = c("lightpink3"))+
    ylim(0, 1) +
    labs(x = "Coverage", y = "Specificity", colour = "Metric") +
    # ggtitle(name_title_final) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 30, size = 16, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        # axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.5, "cm")
    ) +
    guides(linetype = guide_legend(override.aes = list(color = "black")))

met2 <- all_stats_final %>%
    filter(metric == "sensitivity") %>%
    ggplot(aes(x = cov, y = values, linetype = Data)) +
    geom_line(aes(group = interaction(type, Data)), linewidth = 1.2, colour = "darkolivegreen4") +
    geom_point(size = 3, colour = "darkolivegreen4") +
    # scale_color_manual(values = c("darkolivegreen4"))+
    ylim(0, 1) +
    labs(x = "Coverage", y = "Sensitivity", colour = "Metric") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 30, size = 16, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        # axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.5, "cm")
    ) +
    guides(linetype = guide_legend(override.aes = list(color = "black")))


met3 <- all_stats_final %>%
    filter(metric == "FDR") %>%
    ggplot(aes(x = cov, y = values, linetype = Data)) +
    geom_line(aes(group = interaction(type, Data)), linewidth = 1.2, colour = "gold2") +
    geom_point(size = 3, colour = "gold2") +
    # scale_color_manual(values = c("gold2"))+
    ylim(0, 1) +
    labs(x = "Coverage", y = "False Discovery Rate", colour = "Metric") +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 30, size = 16, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        # axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank(),
        legend.key.width = unit(1.5, "cm")
    ) +
    guides(linetype = guide_legend(override.aes = list(color = "black")))


# combine all together:
met_all <- ggarrange(a1,
    ggarrange(m1, met1, met2, met3, ncol = 4, nrow = 1, common.legend = TRUE, legend = "bottom", labels = c("b", "c", "d", "e"), font.label = list(size = 20)),
    labels = c("a"),
    # labels = unique_label,
    ncol = 1, nrow = 2,
    font.label = list(size = 20),
    heights = c(0.5, 0.45)
)

png(args[[17]], width = 21, height = 12, units = "in", res = 200, pointsize = 4)
# png('~/Downloads/test.png', width=22, height=12, units='in', res=200, pointsize=4)
met_all
dev.off()
