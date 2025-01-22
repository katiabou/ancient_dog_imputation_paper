#! /usr/bin/env Rscript

# import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(data.table)
library(viridis)

args <- commandArgs(trailingOnly = TRUE)

list_of_concordance_phased <- args[1]
list_of_validation <- args[2]
list_of_phased <- args[3]
name <- args[4]
chrom <- args[5]
cov_hc <- args[6]

name_title <- args[7]
name_title_2 <- gsub("_", " ", name_title)

site_type <- args[8]

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

# chr length
chr_length <- read.delim(args[10], header = FALSE)
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

# reorder names
d$cov <- factor(d$cov, levels = c("0.05x imputed", "0.1x imputed", "0.2x imputed", "0.5x imputed", "1x imputed", "2x imputed", name_HC_imputed, name_HC))

# rescale x axis
new_sub_final <- d %>%
    mutate(
        POS1 = POS1 / 1e+6,
        POS2 = POS2 / 1e+6,
        MB = KB / 1000
    )

size_chr_mb <- size_chr / 1e+6

# name_title_final <- paste(name_title_2, ' - ', name,sep='')
name_title_final <- paste(name_title_2, " - ", name, " (", site_type, ")", sep = "")

# plot
png(args[9], width = 16, height = 4, units = "in", res = 250, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
options(scipen = 10000)
ggplot(data = new_sub_final) +
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
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "gray60")
    ) +
    labs(colour = "Coverage") +
    ggtitle(name_title_final)
dev.off()
