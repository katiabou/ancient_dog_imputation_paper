#! /usr/bin/env Rscript

# import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)


# import data
phased <- read.csv(snakemake@input[[1]], sep = "")
validation <- read.csv(snakemake@input[[2]], sep = "")
ref_fasta_chr_size <- read.delim(snakemake@input[[3]], header = FALSE)
sample_name <- as.character(snakemake@params[["name"]])
cov_imputed <- as.character(snakemake@params[["cov"]])

# remove extra line used for plotting from validation
val <- validation %>%
    filter(IID != 0) %>%
    mutate(range = paste(POS1, "-", POS2, sep = ""))

# subset for one sample phased
phased_target <- phased %>%
    filter(IID != 0 & IID == sample_name) %>%
    mutate(range = paste(POS1, "-", POS2, sep = ""))

# use GRanges to create approprite format
gr <- GRanges(
    seqnames = val$CHR,
    ranges = val$range,
    sample = val$FID,
    cov = val$cov
)
# gr

gr1 <- GRanges(
    seqnames = phased_target$CHR,
    ranges = phased_target$range,
    sample = phased_target$FID,
    cov = phased_target$cov
)
# gr1

# length(gr1)
# range(gr)
# gaps(gr)
# width(gr)

# find overlaps between the two samples
hits <- findOverlaps(gr, gr1)

# width(overlaps)
# pintersect(gr[queryHits(hits)], gr1[subjectHits(hits)])
# queryHits(hits)
# subjectHits(hits)

# shows the gr ranges that were overlapped with gr1
# gr[gr %over% gr1]
# gr[queryHits(hits)]

# shows the gr1 ranges that were overlapped with gr
# gr1[subjectHits(hits)]

# get overlaps
get_overlap <- function(query, subject, hits) {
    overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
    overlaps$proportion_subject <- width(overlaps) / width(subject[subjectHits(hits)]) # this is the proportion of the subject that is contained within the query
    overlaps$proportion_query <- width(overlaps) / width(query[queryHits(hits)]) # this is the proportion of the query that is overlapping with the subject

    return(overlaps)
}

ab <- get_overlap(gr, gr1, hits)


############ SEGMENT BASED ESTIMATION ###########

# so in this example, I have 35 ROHs identified in the phased_target that are also present in the val (so are actually true)
TP <- length(unique(queryHits(hits)))

# false positives will be a ROH found in the phased_target that was not present in the val
FP <- length(gr1[!(gr1 %in% gr1[subjectHits(hits)])])

# false negatives will be a ROH not found in the phased_target that was present in the val
FN <- length(gr[!(gr %in% gr[queryHits(hits)])])


## Now get F1 stat for segment based overlaps:
recall_seg <- TP / (TP + FN)
precision_seg <- TP / (TP + FP)
F1_seg <- 2 * (precision_seg * recall_seg / (precision_seg + recall_seg))

# False discovery rate:
FDR <- FP / (TP + FP)
# FDR <- 1 - precision_seg

# summarize
df_seg <- data.frame(
    sample = sample_name,
    cov = cov_imputed,
    TP = TP,
    FP = FP,
    FN = FN,
    sensitivity = recall_seg,
    precision = precision_seg,
    F1 = F1_seg,
    FDR = FDR
)


############ TOTAL LENGTH BASED ESTIMATION ###########
genome_size <- ref_fasta_chr_size[1, 2]

# TP: sum up the overlaps from the ab dataframe and divide by the whole chr size
TP_overlaps_sum <- sum(width(ab))
TP <- TP_overlaps_sum / genome_size

# FP: ROHs found in the downsampled imputed that do not overlap a ROH in the HC
FP_overlaps_sum <- sum(width(gr1[!(gr1 %in% gr1[subjectHits(hits)])]))
FP <- FP_overlaps_sum / genome_size

# false negatives will be a ROH not found in the phased_target that was present in the val
FN_overlaps_sum <- sum(width(gr[!(gr %in% gr[queryHits(hits)])]))
FN <- FN_overlaps_sum / genome_size

# True negative is the total chr length minus the total ROH length from the validation
TN_sum <- genome_size - sum(width(gr))
TN <- TN_sum / genome_size

## Now get F1 stat for length based overlaps (worst value: 0; best value: 1):
recall_length <- TP / (TP + FN)
precision_length <- TP / (TP + FP)
F1_length <- 2 * (precision_length * recall_length / (precision_length + recall_length))

library(mltools)
mcc <- mcc(TP = TP, FP = FP, TN = TN, FN = FN)

# MCC <- ((TP*TN) - (FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

# unnormilized is from -1 to 1 (worst value: –1; best value: +1)
# normalize if needed: is from 0 to 1
mcc_norm <- (mcc + 1) / 2 # IF 0 until <0.5 it's worse than random, 0.5 random and >0.5 is better estimations

# other metrics to keep in mind:

# accuracy (worst value: 0; best value: 1):
accuracy <- (TP + TN) / (TP + TN + FP + FN)

# Specificity
specificity <- TN / (TN + FP)

# False discovery rate:
FDR <- FP / (TP + FP)
FDR <- 1 - precision_length

# summarize
df_length <- data.frame(
    sample = sample_name,
    cov = cov_imputed,
    TP = TP,
    FP = FP,
    FN = FN,
    sensitivity = recall_length,
    precision = precision_length,
    F1 = F1_length,
    FDR = FDR,
    mcc = mcc,
    mcc_norm = mcc_norm,
    accuracy = accuracy,
    specificity = specificity
)


# export both files:
write.table(df_seg, file = snakemake@output[[1]], quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
write.table(df_length, file = snakemake@output[[2]], quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
