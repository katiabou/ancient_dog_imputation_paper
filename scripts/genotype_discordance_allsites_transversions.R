#! /usr/bin/env Rscript

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
discordance_phased_trans <- args[2]
meta <- read.csv(args[3], sep = "")

a <- str_split(discordance_phased, pattern = ",")
b <- str_split(discordance_phased_trans, pattern = ",")

all <- c()
all_trans <- c()

for (i in 1:length(a[[1]])) {
    tmp <- read.table(a[[1]][i], quote = "\"", comment.char = "")[, c(1:2, 5, 12:15)]
    all <- rbind(all, tmp)
}

for (i in 1:length(b[[1]])) {
    tmp <- read.table(b[[1]][i], quote = "\"", comment.char = "")[, c(1:2, 5, 12:15)]
    all_trans <- rbind(all_trans, tmp)
}

colnames(all) <- c("cov", "info", "sample", "RR", "RA", "AA", "NRD")
colnames(all_trans) <- c("cov", "info", "sample", "RR", "RA", "AA", "NRD")

mdata1 <- reshape2::melt(all, id = c("sample", "cov", "info"))
mdata2 <- reshape2::melt(all_trans, id = c("sample", "cov", "info"))
colnames(mdata2)[5] <- "value_trans"

df_merge <- merge(mdata1, mdata2, by = c("cov", "info", "sample", "variable"))
mdata <- df_merge
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
    "0.95" = "0.95 INFO cutoff"
)


# subset dogs
mdata_dogs <- mdata %>% filter(species == "Dogs")

# order alphabetically
mdata_dogs <- mdata_dogs[order(mdata_dogs$sample), ]

# specify colour palettes
# cols = met.brewer(name="Juarez", n=6, type="discrete")
cols <- met.brewer(name = "Archambault", n = 6, type = "discrete")
# cols = met.brewer(name="Derain", n=6, type="discrete")

png(args[4], width = 12, height = 19, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
m <- ggplot(data = mdata_dogs, aes(value, value_trans, col = cov, shape = variable))
m + geom_point(aes(col = cov), size = 3, alpha = 0.8) +
    scale_colour_manual(values = cols, name = "Coverage", labels = c("0.05x", "0.1x", "0.2x", "0.5x", "1x", "2x")) +
    scale_shape_manual(values = c(15, 16, 17, 8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, size = 10, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 10, vjust = 0.5)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 12)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white")) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    guides(colour = guide_legend(byrow = TRUE)) +
    labs(x = "Genotyping error all sites %", y = "Genotyping error transversions %", shape = "Site") +
    facet_grid(info_sample_name ~ info, labeller = labeller(info = info_names))
dev.off()




# plot with all wolf samples and all INFO/cov

# subset dogs
mdata_wolves <- mdata %>% filter(species == "Pleistocene_wolf")

# order alphabetically
mdata_wolves <- mdata_wolves[order(mdata_wolves$sample), ]

cols <- met.brewer(name = "Archambault", n = 6, type = "discrete")

png(args[5], width = 12, height = 9, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
m <- ggplot(data = mdata_wolves, aes(value, value_trans, col = cov, shape = variable))
m + geom_point(aes(col = cov), size = 3, alpha = 0.8) +
    scale_colour_manual(values = cols, name = "Coverage", labels = c("0.05x", "0.1x", "0.2x", "0.5x", "1x", "2x")) +
    scale_shape_manual(values = c(15, 16, 17, 8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, size = 10, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 10, vjust = 0.5)) +
    theme(axis.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 10)) +
    theme(legend.title = element_text(size = 12)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white")) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    guides(colour = guide_legend(byrow = TRUE)) +
    labs(x = "Genotyping error all sites %", y = "Genotyping error transversions %", shape = "Site") +
    facet_grid(info_sample_name ~ info, labeller = labeller(info = info_names))
dev.off()
