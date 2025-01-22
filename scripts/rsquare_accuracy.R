#! /usr/bin/env Rscript

# import libraries
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)

# import data
# ALL_INFO <- read.table("~/Downloads/concordance_NGDG_chr1_0.5x-only_sample-INFO_0.0_filtered.rsquare.grp.txt.gz", quote="\"", comment.char="")
# INFO_0.8 <- read.table("~/Downloads/concordance_NGDG_chr1_0.5x-only_sample-INFO_0.8_filtered.rsquare.grp.txt.gz", quote="\"", comment.char="")
# INFO_0.9 <- read.table("~/Downloads/concordance_NGDG_chr1_0.5x-only_sample-INFO_0.9_filtered.rsquare.grp.txt.gz", quote="\"", comment.char="")
# INFO_0.95 <- read.table("~/Downloads/concordance_NGDG_chr1_0.5x-only_sample-INFO_0.95_filtered.rsquare.grp.txt.gz", quote="\"", comment.char="")

# import data
INFO_0.8 <- read.table(snakemake@input[[1]])
INFO_0.9 <- read.table(snakemake@input[[2]])
INFO_0.95 <- read.table(snakemake@input[[3]])
ALL_INFO <- read.table(snakemake@input[[4]])

chr <- as.character(snakemake@params[["chr"]])
# chr <- 'chr1'
name <- as.character(snakemake@params[["name"]])
# name <- 'NGDG'
cov <- as.character(snakemake@params[["cov"]])
# cov <- '0.5'

x <- "x"
covx <- paste(cov, x, sep = "")

target <- paste("Target sample:", name, sep = " ")
down_cov <- paste("Coverage:", covx, sep = " ")

name_chr_cov <- paste(target, down_cov, sep = " | ")

names(INFO_0.8)[c(3, 5)] <- c("x", "y")
names(INFO_0.9)[c(3, 5)] <- c("x", "y")
names(INFO_0.95)[c(3, 5)] <- c("x", "y")
names(ALL_INFO)[c(3, 5)] <- c("x", "y")

newData <- melt(list(
    INFO_0.8 = INFO_0.8[c(3, 5)], INFO_0.9 = INFO_0.9[c(3, 5)],
    INFO_0.95 = INFO_0.95[c(3, 5)], ALL_INFO = ALL_INFO[c(3, 5)]
), id.vars = "x")

cols <- c("#CC6677", "#009E73", "#DDCC77", "#0072B2")

test <- newData[order(newData$x), ]
GroupLabels <- 0:(nrow(test) - 1) %/% 4
test$Group <- GroupLabels
Avgs <- test %>%
    group_by(Group) %>%
    summarize(Avg = mean(x))
breaks <- log(Avgs[, 2])

png(snakemake@output[[1]], width = 8, height = 6, units = "in", res = 250, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
ggplot(newData, aes(log(x), value, colour = L1)) +
    geom_line(linetype = "solid", size = 1.5) +
    geom_point(size = 3) +
    scale_colour_manual(values = cols, name = "INFO cutoff", labels = c("No cutoff", "0.8", "0.9", "0.95")) +
    # ylim(0,1)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
    labs(y = bquote("r"^2), x = "Minor allele frequency") +
    theme_bw() +
    scale_x_continuous(
        breaks = breaks$Avg,
        labels = c("0-0.001", "0.001-0.002", "0.002-0.005", "0.005-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.2", "0.2-0.5")
    ) +
    theme(
        axis.text.x = element_text(angle = 30, size = 16, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor = element_blank()
    ) +
    ggtitle(name_chr_cov)
dev.off()
