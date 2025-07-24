#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(dplyr)
library(tidyr)
library(data.table)
library(readr)

# import liftovered cnv_windows
dog <- fread(snakemake@input[[1]])
wolf <- fread(snakemake@input[[2]])
wg_windows <- fread(snakemake@input[[3]])

# sample list
dog10k_good_names <- read.table(snakemake@input[[4]], quote="\"", comment.char="")
colnames(dog10k_good_names) <- c("Sample.Name")

# dog10k metadata
meta <- read.delim(snakemake@input[[5]])

# merge metadata with sample list used for cnv in dog10k 
meta_merged <- merge(dog10k_good_names, meta, by="Sample.Name")[,c(1:3)]

# wolves
meta_merged_wolves <- meta_merged %>% filter(Category=='Wolf') %>% dplyr::select('Sample.Name')

# dogs
meta_merged_dogs <- meta_merged %>% filter(Category!='Wolf') %>% dplyr::select('Sample.Name')


# group per window and count frequency for dogs:
dog$V4 <- paste(dog$V1, dog$V2, dog$V3, sep="_")

all_chr_dog_counts <- dog %>%
    group_by(V1,V2,V3,V4)%>%
    summarise(V5=(n())) %>%
    mutate(V5=V5/nrow(meta_merged_dogs))


# group per window and count frequency for wolves:
wolf$V4 <- paste(wolf$V1, wolf$V2, wolf$V3, sep="_")

all_chr_wolf_counts <- wolf %>%
    group_by(V1,V2,V3,V4)%>%
    summarise(V5=(n())) %>%
    mutate(V5=V5/nrow(meta_merged_wolves))

# prepare whole genome windows as well
wg_windows$V4 <- paste(wg_windows$V1, wg_windows$V2, wg_windows$V3, sep="_")

# export 
write.table(all_chr_dog_counts, snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_chr_wolf_counts, snakemake@output[[2]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(wg_windows, snakemake@output[[3]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



