#! /usr/bin/env Rscript

# import libraries
library(readr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)

samples <- read.table(snakemake@input[[1]], quote = "\"", comment.char = "")
eigenval_output <- read.table(snakemake@input[[2]], quote = "\"", comment.char = "")
eigenvec_output <- read_table(snakemake@input[[3]], col_names = FALSE)

name <- as.character(snakemake@params[["name"]])
info_sample <- as.character(snakemake@params[["info_sample"]])
name_title_2 <- gsub("_", " ", info_sample)
cov_sample <- as.character(snakemake@params[["cov_sample"]])

# subset the dataset to only include the target sample:
target_sub1 <- eigenvec_output %>% filter(grepl(name, X1), X1 != paste(name, "_HC_imputed", sep = ""))

colnames(target_sub1) <- c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop")

# Distance for each PC

# make list with % of each PC:
mylist <- c()
for (s in eigenval_output$V1) {
    print(s / sum(eigenval_output$V1))
    a <- (s / sum(eigenval_output$V1))
    mylist <- c(mylist, a)
}


get_PC <- function(PC) {
    round(mylist[PC] * 100, digits = 2)
}

# get first 10 PCs:
PCs <- NA
for (i in 1:10) {
    PCs[i] <- get_PC(i)
}


# get first 10 PCs for high coverage genotyped:
HC_PC <- target_sub1 %>%
    filter(target_sub1$Sample == name) %>%
    select(2:11)
# HC_imputed_PC <- target_sub %>% filter(target_sub$Sample==paste(name, "_HC_imputed", sep="")) %>% select(2:11)

target_sub <- target_sub1 %>% filter(target_sub1$Sample != name)

# HC genotyped
get_PC_dist <- function(PC, df, row_num) {
    gg <- paste("PC", PC, sep = "")
    abs(HC_PC[PC] - df[row_num, gg]) * PCs[PC]
}

HC_PC_dist <- data.frame()
for (j in 1:nrow(target_sub)) {
    for (i in 1:10) {
        HC_PC_dist[j, i] <- get_PC_dist(i, target_sub, j)
    }
}

# HC imputed
# get_PC_imputed_dist <- function(PC,df,row_num){
#  gg <- paste('PC', PC, sep="")
#  abs(HC_imputed_PC[PC]-df[row_num,gg])*PCs[PC]
# }

# HC_PC_imputed_dist <- data.frame()
# for (j in 1:nrow(target_sub)){
#  for (i in 1:10){
#    HC_PC_imputed_dist[j,i] <- get_PC_imputed_dist(i, target_sub, j)
#  }
# }


# sum the weighted distances for all 10 PCs:
HC_PC_dist <- HC_PC_dist %>%
    mutate(Total = select(., PC1:PC10) %>% rowSums(na.rm = TRUE))

# HC_PC_imputed_dist <- HC_PC_imputed_dist %>%
#  mutate(Total = select(., PC1:PC10) %>% rowSums(na.rm = TRUE))

# make plot where x is the DoC and y is the sum of % across 10 PCs. Each line will be a sample and the p-value will be per sample (?)

# fix coverage column
target_sub$coverage <- NA

cov1 <- 0.05
target_sub$coverage <- ifelse(grepl(paste("_", cov1, "x", sep = ""), target_sub$Sample), paste(cov1, "x", sep = ""), target_sub$coverage)
cov2 <- 0.1
target_sub$coverage <- ifelse(grepl(paste("_", cov2, "x", sep = ""), target_sub$Sample), paste(cov2, "x", sep = ""), target_sub$coverage)
cov3 <- 0.2
target_sub$coverage <- ifelse(grepl(paste("_", cov3, "x", sep = ""), target_sub$Sample), paste(cov3, "x", sep = ""), target_sub$coverage)
cov4 <- 0.5
target_sub$coverage <- ifelse(grepl(paste("_", cov4, "x", sep = ""), target_sub$Sample), paste(cov4, "x", sep = ""), target_sub$coverage)
cov5 <- 1
target_sub$coverage <- ifelse(grepl(paste("_", cov5, "x", sep = ""), target_sub$Sample), paste(cov5, "x", sep = ""), target_sub$coverage)
cov6 <- 2
target_sub$coverage <- ifelse(grepl(paste("_", cov6, "x", sep = ""), target_sub$Sample), paste(cov6, "x", sep = ""), target_sub$coverage)

# add high coverage value:
# target_sub$coverage <- as.character(ifelse(is.na(target_sub$coverage), 'HC', target_sub$coverage))
target_sub$type <- as.character(ifelse(grepl("imputed", target_sub$Sample), "Imputed", "Pseudohaploid"))


HC_PC_dist$cov <- target_sub$coverage
# HC_PC_imputed_dist$cov <- target_sub$coverage

HC_PC_dist$type <- target_sub$type
# HC_PC_imputed_dist$type <- target_sub$type


# plot
# title

name_title_final <- paste(name_title_2, " - ", name, sep = "")
name_HC <- paste("HC (", cov_sample, "x)", sep = "")
# HC_PC_imputed_dist$cov <- gsub("HC",name_HC, HC_PC_imputed_dist$cov)
HC_PC_dist$cov <- gsub("HC", name_HC, HC_PC_dist$cov)



# distance between imputed dowsampled and PH HC
png(snakemake@output[[1]], width = 9, height = 6, units = "in", res = 200, pointsize = 4)
par(
    mar      = c(5, 5, 2, 2),
    xaxs     = "i",
    yaxs     = "i",
    cex.axis = 2,
    cex.lab  = 2
)
ggplot(HC_PC_dist, aes(x = cov, y = Total, colour = type)) +
    geom_line(aes(group = type), linewidth = 1.2) +
    geom_point() +
    scale_color_manual(values = c("lightblue3", "olivedrab4")) +
    labs(x = "Coverage", y = "Sum of weighted PC distances", colour = "Data type") +
    ylim(0, 0.5) +
    theme_bw() +
    theme(
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
ggtitle(name_title_final)
dev.off()


# get Proportion of green against blue to see how much the pseudohaploid is doing better
HC_PC_dist_prop <- HC_PC_dist %>%
    select(Total, cov, type) %>%
    group_by(cov) %>%
    summarise(ratio_pseudohaploid_imputed = Total[type == "Pseudohaploid"] / Total[type == "Imputed"]) %>%
    mutate(sample = name)

# export data table
# write.table(HC_PC_dist_prop, file='~/Downloads/test.tsv', quote=FALSE, sep='\t', row.names = FALSE)
write.table(HC_PC_dist_prop, file = snakemake@output[[2]], quote = FALSE, sep = "\t", row.names = FALSE)
