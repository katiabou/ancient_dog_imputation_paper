#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

# import libraries
library(dplyr)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GOfuncR)
library(ggpmisc)

if (!requireNamespace("windowscanr", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("tavareshugo/windowscanr")
}
library(windowscanr)

options(scipen=999)

# ROH deserts 
deserts <- read.delim(snakemake@input[[1]], header=FALSE)
deserts$region <- paste(deserts$V1,":",deserts$V2,"-",deserts$V3, sep="")

# kmer depth of coverage across whole genome windows
genome_kmer_bad_dog_allchrom <- read.delim(snakemake@input[[2]], header=FALSE)

# imputed dogs
hom_sum <-fread(snakemake@input[[3]])
ind_file <- fread(snakemake@input[[4]])

# modern dogs
hom_sum_modern <-fread(snakemake@input[[5]])
ind_file_modern <- read.csv(snakemake@input[[6]], sep="")

# per window coverage estimate file
all_cov_window <- read.delim(snakemake@input[[7]], header=FALSE)

# DLA genes
DLA_regions <- read.delim(snakemake@input[[8]], header=FALSE)



################################ Estimating the window cutoffs based on imputed dogs or imputed wolves (from the bams)  ################################
# create index column
hom_sum <- hom_sum %>%
  mutate(index = 1:nrow(.))

running_roh <- winScan(x = hom_sum,
                       groups = "CHR",
                       position = "BP",
                       values = "UNAFF",
                       win_size = 500000,
                       win_step = 500000,
                       funs = c("mean"),
                       cores = 8)

# remove windows without snps and calculate mean unaff
running_roh %>% 
  mutate(UNAFF_mean = UNAFF_mean/nrow(ind_file)) %>% 
  filter(UNAFF_n > 0) -> running_roh_p

# take for df:
running_roh_p$region <- paste(running_roh_p$CHR,":",running_roh_p$win_start,"-",running_roh_p$win_end, sep="")

# prepare window coverage file
all_cov_window$V1 <- str_replace_all(all_cov_window$V1, 'chr', '')
all_cov_window$region <- paste(all_cov_window$V1,":",all_cov_window$V2,"-",all_cov_window$V3, sep="")

# count number of bam samples used for the window coverage calculation (same number as the imputed dogs or wolves)
num_ind <- nrow(ind_file)

# Group by region and estimate across all samples the average mean depth for each window:
all_cov_window %>%
  group_by(region) %>%
  summarize(Mean = sum(V7)/num_ind) -> test

#merge the window coverage with the ROH data (either imputed dogs or imputed wolves)
joinn <- left_join(running_roh_p, test, by=c('region'))

# get mean and std
meann_wind <- mean(joinn$Mean)
stdd_wind <- 2*sd(joinn$Mean)
low_wind <- meann_wind-stdd_wind
high_wind <- meann_wind+stdd_wind


####################################### Modern dogs  ######################################

# Need the mean ROH prevelance data from the modern samples, since we have more samples there 

### only modern dog heatmap
hom_sum_modern <- hom_sum_modern %>%
  mutate(MB = BP / 1000000,
         KB = BP / 1000,
         index = 1:nrow(.))

# count ROH in running windows of 500 Kb
# UNAFF	is the number of non-cases with a ROH including this SNP
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh_modern <- winScan(x = hom_sum_modern,
                              groups = "CHR",
                              position = "BP",
                              values = "UNAFF",
                              win_size = 500000,
                              win_step = 500000,
                              funs = c("mean"),
                              cores = 8)

# take for df:
running_roh_modern$region <- paste(running_roh_modern$CHR,":",running_roh_modern$win_start,"-",running_roh_modern$win_end, sep="")

# remove windows without snps
running_roh_modern %>% 
  mutate(UNAFF_mean = UNAFF_mean/nrow(ind_file_modern)) %>% 
  filter(UNAFF_n > 0) -> running_roh_p_modern

#windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean<low_wind | Mean>high_wind)

#remove windows with extreme depth estimates (may have CNVs):
running_roh_p_modern_nogrey <- running_roh_p_modern %>%
  filter(!region %in% final_window_remove$region)


####################################### CNV depth with weighted mean frequency  ######################################

# this is what we want to weight the mean frequency of each small snv window with the overlap with our 500kb windows:
#sum(CNV_freq *(CNV_overlap / total_overlap_of_all_CNVs))

# group by unique chr, start and end and get the mean weighted frequency value 
kmer_depth_mean <- genome_kmer_bad_dog_allchrom %>%
  group_by(V1,V2,V3) %>%
  dplyr::summarize(weighted_mean = sum(V9 *(V10/sum(V10))))

# add column names
colnames(kmer_depth_mean) <- c('CHR', 'win_start', 'win_end','weighted_mean_frequency')

# replace chr in CHR column:
kmer_depth_mean$CHR <- gsub("chr", "", kmer_depth_mean$CHR)

# add region column
kmer_depth_mean$region <- paste(kmer_depth_mean$CHR,":",kmer_depth_mean$win_start,"-",kmer_depth_mean$win_end, sep="")

# remove extreme coverage windows
kmer_depth_mean_no_grey <- kmer_depth_mean %>%
  filter(!region %in% final_window_remove$region) 

# get mean and std
meann <- mean(kmer_depth_mean_no_grey$weighted_mean_frequency)
stdd <- 2*sd(kmer_depth_mean_no_grey$weighted_mean_frequency)
high <- meann+stdd

#plot distribution of mean cnv window frequency for whole genome, but without extreme coverage windows
p1 <- kmer_depth_mean_no_grey %>%
  ggplot(aes(x=weighted_mean_frequency)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  theme_bw()+
  xlab('CNV window weighted mean frequency') +
  xlim(0,max(kmer_depth_mean_no_grey$weighted_mean_frequency))

p2 <- kmer_depth_mean_no_grey %>%
  ggplot(aes(x=weighted_mean_frequency)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  theme_bw()+
  xlab('CNV window weighted mean frequency') +
  xlim(0,0.25)

png(snakemake@output[[1]], width=12, height=6, units='in', res=200, pointsize=4)

ggarrange(p1, p2,
          labels = c("a", "b"),
          ncol = 2, nrow = 1, font.label=list(size=20))
dev.off()

# filter for 500kb windows which have a cnv frequency below the selected cutoff
kmer_depth_mean_no_grey_filt <- kmer_depth_mean_no_grey %>% filter(weighted_mean_frequency<0.13) #cutoff for modern dogs



####################################### CNV frequencry against ROH prevelance  ######################################

# merge whole genome kmer with ROH prevelance from modern dogs
merged_wg_modern_mean1 <- merge(kmer_depth_mean_no_grey, running_roh_p_modern_nogrey, by=c('region', 'CHR', 'win_start', 'win_end'))
merged_wg_modern_mean <- merge(kmer_depth_mean_no_grey_filt, running_roh_p_modern_nogrey, by=c('region', 'CHR', 'win_start', 'win_end'))

# running the regression only on the deserts with no filtering cutoff
p1 <- ggplot(merged_wg_modern_mean1 %>% filter(region %in% deserts$region), aes(x=weighted_mean_frequency, y=UNAFF_mean*100))+
  theme_bw()+
  labs(x="CNV window weighted mean frequency", y="Mean ROH prevelance in ROH deserts (%)")+
  geom_smooth (data =  .%>% filter(region %in% deserts$region), method='lm', formula = "y ~ x") +
  geom_point()+
  stat_fit_glance(data =  .%>% filter(region %in% deserts$region),
    aes(label = after_stat(
      sprintf('italic(p)~"="~%s~","~italic(r)^2~"="~%.3f',
              ifelse(..p.value.. < 0.05,
                     sprintf("%.2e", ..p.value..),
                     sprintf("%.2f", ..p.value..)),
              r.squared)
    )),
    parse = TRUE,
    size = 4, color = "red", label.x = "right", label.y = "top"
  )

# running the regression only on the deserts with filtering cutoff
p2 <- ggplot(merged_wg_modern_mean %>% filter(region %in% deserts$region), aes(x=weighted_mean_frequency, y=UNAFF_mean*100))+
  theme_bw()+
  labs(x="CNV window weighted mean frequency", y="Mean ROH prevelance in ROH deserts (%)")+
  geom_smooth (data =  .%>% filter(region %in% deserts$region), method='lm', formula = "y ~ x") +
  geom_point()+
  stat_fit_glance(data =  .%>% filter(region %in% deserts$region),
    aes(label = after_stat(
      sprintf('italic(p)~"="~%s~","~italic(r)^2~"="~%.3f',
              ifelse(..p.value.. < 0.05,
                     sprintf("%.2e", ..p.value..),
                     sprintf("%.2f", ..p.value..)),
              r.squared)
    )),
    parse = TRUE,
    size = 4, color = "red", label.x = "right", label.y = "top"
  )

png(snakemake@output[[2]], width=12, height=6, units='in', res=200, pointsize=4)
ggarrange(p1, p2,
          labels = c("a", "b"),
          ncol = 2, nrow = 1, font.label=list(size=20))
dev.off()


### number of ROH deserts left in modern dogs:
# merged_wg_modern_mean_deserts <- merged_wg_modern_mean %>% filter(region %in% deserts$region)



####################################### RUN GOfunR for gene ontology analysis  ######################################
if (!require("BiocManager", quietly = TRUE)) {
  options(repos = c(CRAN = "https://cran.rstudio.com"))
  install.packages("BiocManager")
}
#BiocManager::install() #check and update to latest bioconductor
BiocManager::install("TxDb.Cfamiliaris.UCSC.canFam3.refGene")
BiocManager::install("org.Cf.eg.db")

# get first and last 500Kb window of each chromosome to remove:
all_windows_no_first_last <- test %>%
  dplyr::select(region) %>%
  separate(region, into = c("chr", "start_end"), sep = ":") %>%
  separate(start_end, into = c("start", "end"), sep = "-") %>%
  mutate(chr = as.numeric(chr), start = as.numeric(start), end = as.numeric(end)) %>% #make sure these are numeric so I can sort
  arrange(chr, start, end) %>% #sort numerically
  group_by(chr) %>%
  filter(row_number()==1 | row_number()==n()) #choose first and last row for each chr (corresponds to window at beginning and end of chromosome)

all_windows_no_first_last$region <- paste(all_windows_no_first_last$chr,":",all_windows_no_first_last$start,"-",all_windows_no_first_last$end, sep="")

# windows outside of the cnv frequency cutoff to remove
kmer_depth_mean_no_grey_filt_remove <- kmer_depth_mean_no_grey %>% 
  filter(weighted_mean_frequency>=0.13) #cutoff for modern dogs

# extreme coverage windows to remove
final_window_remove <- test %>% 
  filter(Mean<low_wind | Mean>high_wind) %>% dplyr::select(region)  

# Windows to be used as background
final_window_background_keep <- test %>%
  filter(!region %in% final_window_remove$region) %>% #remove extreme coverage windows
  filter(!region %in% all_windows_no_first_last$region) #remove windows at beginning and end of chr

# desert regions windows to be used as foreground (extreme coverage windows already removed)
imputed_modern_regions <- deserts %>% 
  filter(!region %in% kmer_depth_mean_no_grey_filt_remove$region) %>% #remove cnv windows
  filter(!region %in% all_windows_no_first_last$region) #remove windows at beginning and end of chr

# prepare the background and target regions
is_candidate_back <- data.frame(final_window_background_keep$region, is_candidate=c(rep(0,nrow(final_window_background_keep))))
is_candidate_target <- data.frame(imputed_modern_regions$region, is_candidate=c(rep(1,nrow(imputed_modern_regions))))
colnames(is_candidate_back) <- c('regions', 'is_candidate')
colnames(is_candidate_target) <- c('regions', 'is_candidate')
is_candidate_all <- rbind(is_candidate_back, is_candidate_target)

# Running a hypergeometric test with correction for gene length using regions
# I have to set background regions when I use the regions option:
res_hyper1=go_enrich(is_candidate_all, test = 'hyper', n_randsets = 1000, 
                     orgDb='org.Cf.eg.db',
                     txDb = 'TxDb.Cfamiliaris.UCSC.canFam3.refGene', gene_len = TRUE,
                     regions = TRUE)


## first element of go_enrich result has the stats
stats1 = res_hyper1[[1]]

## see which genes are located in the candidate region
input_genes = res_hyper1[[2]]
candidate_genes = input_genes[input_genes[,2]==1, 1]
candidate_genes <- as.data.frame(candidate_genes)

#output table with GO terms and candidate genes:
write.table(stats1, snakemake@output[[3]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(candidate_genes, snakemake@output[[4]], sep = "\t", row.names = FALSE, quote = FALSE)



##############################################################
#Removing the DLA regions to see if the signal still persist
# DLA-chr12	chr12:307171-2872051
# DLA-chr35	chr35:25514060-26406861
# DLA-79	chr18:41142496-41145658

# input DLA regions:
colnames(DLA_regions) <- c('gene', 'region')

# replace chr:
DLA_regions$region <- str_replace_all(DLA_regions$region, 'chr', '')

# Split region column into chr start and end:
DLA_regions <- DLA_regions %>%
  separate(region, into = c("chr", "start_end"), sep = ":") %>%
  separate(start_end, into = c("start", "end"), sep = "-") %>%
  mutate(start = as.numeric(start), end = as.numeric(end))

# correct column names roh deserts (cnv regions and chr start and end regions removed above)
colnames(imputed_modern_regions) <- c("chr", "start", "end", "region")

# Function to check overlap
check_overlap <- function(df2_row, df1) {
  any(df1$chr == df2_row$chr & df1$end >= df2_row$start & df1$start <= df2_row$end)
}

# Filter out ROH deserts overlapping with DLA regions 
imputed_modern_regions_no_DLA <- imputed_modern_regions %>%
  rowwise() %>%
  filter(!check_overlap(cur_data(), DLA_regions)) %>%
  ungroup()

# Filter out background windows overlapping with DLA regions 
final_window_background_keep_no_DLA <- final_window_background_keep %>%
  separate(region, into = c("chr", "start_end"), sep = ":") %>%
  separate(start_end, into = c("start", "end"), sep = "-") %>%
  rowwise() %>%
  filter(!check_overlap(cur_data(), DLA_regions)) %>%
  ungroup()

final_window_background_keep_no_DLA$region <- paste(final_window_background_keep_no_DLA$chr,":",final_window_background_keep_no_DLA$start,"-",final_window_background_keep_no_DLA$end, sep="")

# prepare the background and target regions
is_candidate_back <- data.frame(final_window_background_keep_no_DLA$region, is_candidate=c(rep(0,nrow(final_window_background_keep_no_DLA))))
is_candidate_target <- data.frame(imputed_modern_regions_no_DLA$region, is_candidate=c(rep(1,nrow(imputed_modern_regions_no_DLA))))
colnames(is_candidate_back) <- c('regions', 'is_candidate')
colnames(is_candidate_target) <- c('regions', 'is_candidate')
is_candidate_all <- rbind(is_candidate_back, is_candidate_target)

# Running a hypergeometric test with correction for gene length using regions
# I have to set background regions when I use the regions option:
res_hyper1=go_enrich(is_candidate_all, test = 'hyper', n_randsets = 1000, 
                     orgDb='org.Cf.eg.db',
                     txDb = 'TxDb.Cfamiliaris.UCSC.canFam3.refGene', gene_len = TRUE,
                     regions = TRUE)


## first element of go_enrich result has the stats
stats2 = res_hyper1[[1]]

## see which genes are located in the candidate region
input_genes = res_hyper1[[2]]
candidate_genes = input_genes[input_genes[,2]==1, 1]
candidate_genes <- as.data.frame(candidate_genes)

#output table with GO terms and candidate genes:
write.table(stats2, snakemake@output[[5]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(candidate_genes, snakemake@output[[6]], sep = "\t", row.names = FALSE, quote = FALSE)


