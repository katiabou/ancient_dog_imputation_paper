#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(devtools))
quiet(library(MetBrewer))
quiet(library(reshape2))
quiet(library(dplyr))
quiet(library(ggplot2))
quiet(library(qdapRegex))
quiet(library(stringr))
quiet(library(argparser))

# get the command line arguments
p <- arg_parser("Count number of sites per bin")
p <- add_argument(p, "--sample-file-merge-1", help = "", default = "output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.NGDG_allchrom_0.5x_INFO_all.txt")
p <- add_argument(p, "--sample-file-merge-2", help = "", default = "output/GLIMPSE_concordance/GLIMPSE_ligated_MAF_bins/merged_ligated.CGG32_allchrom_1x_INFO_all.txt")
p <- add_argument(p, "--sample-file-merge-3", help = "", default = "output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.NGDG_allchrom_0.5x_INFO_all.txt")
p <- add_argument(p, "--sample-file-merge-4", help = "", default = "output/GLIMPSE_concordance_only_dogs/GLIMPSE_ligated_MAF_bins/merged_ligated.CGG32_allchrom_1x_INFO_all.txt")
p <- add_argument(p, "--output", help = "", default = "output/GLIMPSE_concordance/plots/glimpse_concordance_MAF_bins_reference_panel/merged_ligated.NGDG_CGG32_allchrom_INFO_all.png")

argv <- parse_args(p)

input1 <- read.table(argv$sample_file_merge_1, quote = "\"", comment.char = "")
input2 <- read.table(argv$sample_file_merge_2, quote = "\"", comment.char = "")
input3 <- read.table(argv$sample_file_merge_3, quote = "\"", comment.char = "")
input4 <- read.table(argv$sample_file_merge_4, quote = "\"", comment.char = "")

input1$ref_panel <- "all_canids"
input2$ref_panel <- "all_canids"
input3$ref_panel <- "dogs_only"
input4$ref_panel <- "dogs_only"

input <- rbind(input1, input2, input3, input4)

# fixing the MAF_bin column
test <- sapply(strsplit(input$V1, "/\\s*"), tail, 1)
t1 <- input$V1
t2 <- ex_between(test, "MAF_", ".vcf.gz")
t3 <- do.call(rbind.data.frame, t2)
colnames(t3) <- "MAF_bin"
t3$MAF_bin <- str_replace_all(t3$MAF_bin, "_", ".")
t4 <- data.frame(sub("^([^.]*.[^.]*).", "\\1-", t3$MAF_bin))
colnames(t4) <- "MAF_bin"
t4$MAF_bin[t4$MAF_bin == "0.0-001"] <- "0-0.001"
input$V8 <- t4$MAF_bin
input$V9 <- paste(input$V2, input$V4, sep = "_")
colnames(input) <- c("file_name", "Sample", "Chr", "Coverage", "INFO", "Sites", "Ref_panel", "MAF_bins", "Sample_cov")

input$INFO <- paste("INFO", input$INFO, sep = " ")
input$INFO[input$INFO == "INFO 0"] <- "No INFO cutoff"

# transform
input_b <- input %>% select("Sample_cov", "INFO", "Sites", "Ref_panel", "MAF_bins")
melt_data2 <- melt(input_b, id = c("MAF_bins", "Sample_cov", "Ref_panel", "INFO"))

cols <- met.brewer(name = "Demuth", n = 8, type = "discrete") # cb friendly nice

# New facet label names for reference variable
ref.labs <- c("All canid reference panel", "Dog reference panel")
names(ref.labs) <- c("all_canids", "dogs_only")

sample.labs <- c("Pleistocene wolf (CGG32 1x)", "Neolithic European dog (NGDG 0.5x)")
names(sample.labs) <- c("CGG32_1", "NGDG_0.5")

# change bar order
melt_data2$INFO <- factor(melt_data2$INFO,
    levels = c("No INFO cutoff", "INFO 0.8", "INFO 0.9", "INFO 0.95")
)

# re-scale number of sites
melt_data2 <- melt_data2 %>%
    mutate(value_2 = value / 1000)

# barplot
ggplot(data = melt_data2, aes(x = INFO, y = value_2, fill = MAF_bins, label = round(value_2, digits = 0))) +
    geom_bar(stat = "identity") +
    geom_text(size = 3.5, position = position_stack(vjust = 0.5)) +
    geom_text(
        aes(label = round(after_stat(y), digits = 0), group = INFO),
        stat = "summary", fun = sum, vjust = -0.5
    ) +
    facet_grid(Sample_cov ~ Ref_panel, labeller = labeller(Ref_panel = ref.labs, Sample_cov = sample.labs)) +
    scale_fill_manual(values = cols, name = "MAF bins") +
    ylab(bquote("Number of sites "(10^3))) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14, vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 14, vjust = 0.5)) +
    theme(axis.title = element_text(size = 14)) +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = 16)) +
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_text(size = 14)) +
    theme(strip.background = element_rect(fill = "gray28")) +
    theme(strip.text = element_text(colour = "white", size = 14)) +
    theme(legend.spacing.y = unit(0.3, "cm")) +
    theme(panel.grid.major.x = element_blank()) +
    theme(legend.key.size = unit(0.8, "cm"))

ggsave(argv$output, width = 14, height = 10)
