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

# import list of files
files <- args[1]

a <- str_split(files, pattern = ",")

# merge all flare dfs
all <- c()

for (i in 1:length(a[[1]])) {
    tmp <- read.table(a[[1]][i], quote = "\"", comment.char = "")
    all <- rbind(all, tmp)
}

colnames(all) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")

cols <- c("#CC6677", "#009E73", "#DDCC77", "#0072B2")


# get break points for x axis, the average of the MAF values in the df
breaks <- all %>%
    group_by(V1) %>%
    summarize(Avg = log(mean(V3)))

cov_info <- c(
    "0.05" = "0.05x",
    "0.1" = "0.1x",
    "0.2" = "0.2x",
    "0.5" = "0.5x",
    "1" = "1x",
    "2" = "2x"
)

#### import metadata
name <- as.character(args[2])
name_title <- as.character(args[3])
name_title_2 <- gsub("_", " ", name_title)
name_title_final <- paste(name_title_2, " - ", name, sep = "")

# define label based on sample name
all$label <- ifelse(all$V8 == "NGDG", "a",
    ifelse(all$V8 == "SOTN01_merged", "b",
        ifelse(all$V8 == "PortauChoix", "c",
            ifelse(all$V8 == "TRF.02.53", "d",
                ifelse(all$V8 == "TRF.05.05", "e",
                    ifelse(all$V8 == "FAMICHN00012", "f",
                        ifelse(all$V8 == "FAMINGR00004", "g",
                            ifelse(all$V8 == "CGG32", "h",
                                ifelse(all$V8 == "CGG33", "i",
                                    ifelse(all$V8 == "WolfHead", "j", "NA")
                                )
                            )
                        )
                    )
                )
            )
        )
    )
)

p1 <- ggplot(all, aes(log(V3), V4, colour = as.character(V7))) +
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
        strip.text.y = element_text(angle = 0, size = 18),
        strip.text.x = element_text(size = 18),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20)
    ) +
    labs(y = bquote("r"^2), x = "Minor allele frequency") +
    facet_wrap(vars(V6), labeller = as_labeller(cov_info)) +
    ggtitle(name_title_final)


png(args[4], width = 17, height = 11, units = "in", res = 200, pointsize = 4)
plot1 <- ggarrange(p1,
    ncol = 1, nrow = 1,
    common.legend = TRUE,
    legend = "right"
)

plot1
dev.off()
