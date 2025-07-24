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

chr <- snakemake@params[[1]]

# import each chr (each has already been rounded turned to binary)
cnv_round_binary_t <- fread(snakemake@input[[1]])
print("Imported binary cnv file")

# window coordinates
window_coords <- read.table(snakemake@input[[2]])
colnames(window_coords) <- c("region")
window_coords <- as.data.table(window_coords)

# extract chromosome
window_coords_chr <- window_coords[grepl(paste0(chr,":"), region)]

# sample list
dog10k_good_names <- read.table(snakemake@input[[3]], quote="\"", comment.char="")
colnames(dog10k_good_names) <- c("Sample.Name")

# dog10k metadata
meta <- read.delim(snakemake@input[[4]])

# ROH desert 
print("Imported all input files")

# merge metadata with sample list used for cnv in dog10k 
meta_merged <- merge(dog10k_good_names, meta, by="Sample.Name")[,c(1:3)]

# wolves
meta_merged_wolves <- meta_merged %>% filter(Category=='Wolf') %>% dplyr::select('Sample.Name')

# dogs
meta_merged_dogs <- meta_merged %>% filter(Category!='Wolf') %>% dplyr::select('Sample.Name')

# add sample names to cnv file:
colnames(cnv_round_binary_t) <- dog10k_good_names$Sample.Name
cnv_round_binary_t <- as.data.frame(cnv_round_binary_t)

# extract columns which are wolves or dogs
cnv_wolves <- cnv_round_binary_t[, colnames(cnv_round_binary_t) %in% meta_merged_wolves$Sample.Name]
cnv_dogs <- cnv_round_binary_t[, colnames(cnv_round_binary_t) %in% meta_merged_dogs$Sample.Name]

print("Seperated dogs and wolves")

# merge window names with cnv df
merged_wolves <- cbind(window_coords_chr, cnv_wolves)
merged_dogs <- cbind(window_coords_chr, cnv_dogs)

print("Merged window names to read depth file")


# for each occurence of either 1 or 0, take the coordinate of that window and paste it into a new file (so from row number 1)
kmers_filter <- function(df, type_input) {
  # Convert to data.table
  dt <- as.data.table(df)
  
  # Reshape to long format and filter for specific value (0 or 1)
  dt_melt <- dt[ , data.table::melt.data.table(.SD, id.vars = "region", variable.name = "Sample", value.name = "Value", verbose = TRUE)][Value == type_input, .(region)]

  # Separate 'region' into 'chr', 'start', and 'end'
  dt_melt[, c("chr", "range") := tstrsplit(region, ":", fixed = TRUE)]
  dt_melt[, c("start", "end") := tstrsplit(range, "-", fixed = TRUE)]
  dt_melt[, range := NULL]  # Remove unnecessary column
  
  # Return final result
  return(dt_melt[, .(chr, start, end)])
}


# Example usage
good_wolves <- kmers_filter(merged_wolves,1)
print("took good kmers for wolves")

bad_wolves <- kmers_filter(merged_wolves,0)
print("took bad kmers for wolves")

good_dogs <- kmers_filter(merged_dogs,1)
print("took good kmers for dogs")

bad_dogs <- kmers_filter(merged_dogs,0)
print("took bad kmers for dogs")


# export bad and good regions 
write_delim(good_wolves, snakemake@output[[1]],col_names = FALSE)
write_delim(bad_wolves, snakemake@output[[2]],col_names = FALSE)
write_delim(good_dogs, snakemake@output[[3]],col_names = FALSE)
write_delim(bad_dogs, snakemake@output[[4]],col_names = FALSE)



