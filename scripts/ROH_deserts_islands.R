#!/usr/bin/env Rscript

# Author:    Katia Bougiouri
# Copyright: Copyright 2025, University of Copenhagen
# Email:     katia.bougiouri@gmail.com
# License:   MIT

#load libraries
library(dplyr)

if (!requireNamespace("windowscanr", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("tavareshugo/windowscanr")
}
library(windowscanr)
library(data.table)
library(tidyverse)
#library(windowscanr)
library(cowplot)
library(grid)
library(ggplotify)
library(patchwork)
library(viridis)
library(scales)
library(ggridges)
library(lme4)
library(broom)
library(ggpubr)

options(scipen=999)

# imputed dogs
hom_sum <-fread(snakemake@input[[1]])
ind_file <- fread(snakemake@input[[2]])

# modern dogs
hom_sum_modern <-fread(snakemake@input[[3]])
ind_file_modern <- read.csv(snakemake@input[[4]], sep="")

#per window coverage estimate file
all_cov_window <- read.delim(snakemake@input[[5]], header=FALSE)


####### Estimating the window cutoffs based on imputed dogs or imputed wolves (from the bams)

#create index column
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

#get mean and std
meann_wind <- mean(joinn$Mean)
stdd_wind <- 2*sd(joinn$Mean)
low_wind <- meann_wind-stdd_wind
high_wind <- meann_wind+stdd_wind

# Extreme coverage window plot
png(snakemake@output[[1]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(joinn, aes(x=Mean, UNAFF_mean))+
  geom_point(data = .%>% filter(Mean>low_wind & Mean<high_wind), alpha=0.6, colour="cyan3")+
  geom_point(data = .%>% filter(Mean<low_wind | Mean>high_wind), alpha=0.6, colour="grey40")+
  geom_smooth (data =  .%>% filter(Mean>low_wind & Mean<high_wind), method='lm', colour="cyan4") +
  geom_vline(xintercept=low_wind, linetype="dashed", 
             color = "red3", size=0.8)+
  geom_vline(xintercept=high_wind, linetype="dashed", 
             color = "red3", size=0.8)+
  theme_bw()+
  ylim(0,1)+
  labs(x="Mean depth per window (500 KB)", y="Mean ROH prevelance per window (500 KB)")
dev.off()



######################################
######### HEAT MAPS ################## 
######################################


####################################### Imputed dogs  ######################################

### only ancient dog heatmap

hom_sum <- hom_sum %>%
  mutate(MB = BP / 1000000,
         KB = BP / 1000,
         index = 1:nrow(.))

# count ROH in running windows of 500 Kb
# UNAFF	is the number of non-cases with a ROH including this SNP
# UNAFF_n is number of SNPs in a window
# UNAFF_mean is mean ROH prevalence in a window
running_roh_imputed <- winScan(x = hom_sum,
                               groups = "CHR",
                               position = "BP",
                               values = "UNAFF",
                               win_size = 500000,
                               win_step = 500000,
                               funs = c("mean"),
                               cores = 8)

# take for df:
running_roh_imputed$region <- paste(running_roh_imputed$CHR,":",running_roh_imputed$win_start,"-",running_roh_imputed$win_end, sep="")

# remove windows without snps
running_roh_imputed %>% 
  mutate(UNAFF_mean = UNAFF_mean/nrow(ind_file)) %>% 
  filter(UNAFF_n > 0) -> running_roh_p_imputed

# prepare colour scale:
fill_cols <- viridis(20, option = "A")
qn <- scales::rescale(quantile(running_roh_p_imputed$UNAFF_mean,
                               probs=seq(0, 1, length.out=length(fill_cols))))

# windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean<low_wind | Mean>high_wind)

# remove windows with extreme depth estimates (may have CNVs):
running_roh_p_imputed_nogrey <- running_roh_p_imputed %>%
  filter(!region %in% final_window_remove$region)

p3 <- ggplot(running_roh_p_imputed_nogrey, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
  geom_tile() +
  geom_tile(data = final_window_remove, aes(x = win_start, y = 0.5), fill = 'grey') +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = seq(0, 125000000, by = 10000000),
                     labels = as.character(seq(0, 125000, 10000)/1000))+
  ylab("Chromosome") +
  scale_fill_gradientn("% of ancient samples with ROH",
                       colors = rev(fill_cols), 
                       limits = c(0,0.55),
                       breaks = c(0.1, 0.3, 0.5),
                       labels = c(10, 30, 50),
                       na.value="grey") +
  facet_grid(CHR~., switch="both") +
  xlab("Position in Mb") +
  theme_minimal(base_family = "Helvetica", 
                base_size = 13)+
  theme(panel.spacing.y=unit(0.1, "lines"),
        panel.grid=element_blank(),
        axis.title.x = element_text(margin=margin(t=5), size=15),
        axis.title.y = element_text(margin=margin(r=5), size=15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size=14),
        strip.text.y.left = element_text(size = 14, angle = 0),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.75,0.13),
        legend.direction = "horizontal",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
  ) +
  ggtitle('Ancient')+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  guides(fill = guide_colourbar(title.position = "bottom" ,
                                barwidth = 12, barheight = 0.8))
p3

ggsave(snakemake@output[[2]], p3, width = 7, height = 6, bg='transparent')


#density plot (add manually to previous plot)
x <- running_roh_p_imputed_nogrey$UNAFF_mean
y <- density(x, n = 2^12)

p4 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) + 
  geom_line() + 
  geom_segment(aes(xend = x, yend = 0, color = x)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits = c(0, 0.55),
                     breaks = c(0.1,0.3, 0.5),
                     labels = c(10, 30, 50)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = c(2,4,6,8,10,12)) +
  scale_color_gradientn(colors = rev(fill_cols), 
                        limits = c(0, 0.55),
                        breaks = c(0.1, 0.3, 0.5),
                        labels = c(10, 30, 50)) +
  theme_minimal(base_family = "Helvetica", 
                base_size = 13)+
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()) + 
  ylab("Density")

p4

ggsave(snakemake@output[[3]], p4, width = 5, height = 2.5)



####################################### Modern dogs  ######################################

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

#prepare colour scale:
fill_cols <- viridis(20, option = "A")
qn <- scales::rescale(quantile(running_roh_p_modern$UNAFF_mean,
                               probs=seq(0, 1, length.out=length(fill_cols))))

#windows below and above threshold to be coloured grey:
final_window_remove <- joinn %>% filter(Mean<low_wind | Mean>high_wind)

#remove windows with extreme depth estimates (may have CNVs):
running_roh_p_modern_nogrey <- running_roh_p_modern %>%
  filter(!region %in% final_window_remove$region)

p5 <- ggplot(running_roh_p_modern_nogrey, aes(x = win_start, y = 0.5, fill = UNAFF_mean)) + 
  geom_tile() +
  geom_tile(data = final_window_remove, aes(x = win_start, y = 0.5), fill = 'grey') +
  scale_y_continuous(expand = c(0,0))+
  scale_x_continuous(expand = c(0,0), 
                     breaks = seq(0, 125000000, by = 10000000),
                     labels = as.character(seq(0, 125000, 10000)/1000))+
  ylab("Chromosome") +
  scale_fill_gradientn("% of present-day samples with ROH",
                       colors = rev(fill_cols), 
                       limits = c(0,0.55),
                       breaks = c(0.1, 0.3, 0.5),
                       labels = c(10, 30, 50),
                       na.value="grey") +
  facet_grid(CHR~., switch="both") +
  xlab("Position in Mb") +
  theme_minimal(base_family = "Helvetica", 
                base_size = 13)+
  theme(panel.spacing.y=unit(0.1, "lines"),
        panel.grid=element_blank(),
        axis.title.x = element_text(margin=margin(t=5), size=15),
        axis.title.y = element_text(margin=margin(r=5), size=15),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size=14),
        strip.text.y.left = element_text(size = 14, angle = 0),
        axis.ticks.x = element_line(linewidth = 0.3),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.5, unit = "cm"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(0.75,0.13),
        legend.direction = "horizontal",
        legend.text = element_text(size=14),
        legend.title = element_text(size=15),
  ) +
  ggtitle('Present-day')+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  guides(fill = guide_colourbar(title.position = "bottom" ,
                                barwidth = 12, barheight = 0.8))
p5

ggsave(snakemake@output[[4]], p5, width = 7, height = 6, bg='transparent')


#density plot (add manually to previous plot)
x <- running_roh_p_modern_nogrey$UNAFF_mean
y <- density(x, n = 2^12)

p6 <- ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) + 
  geom_line() + 
  geom_segment(aes(xend = x, yend = 0, color = x)) +
  scale_x_continuous(expand = c(0, 0), 
                     limits = c(0, 0.55),
                     breaks = c(0.1,0.3, 0.5),
                     labels = c(10, 30, 50)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_gradientn(colors = rev(fill_cols), 
                        limits = c(0, 0.55),
                        breaks = c(0.1, 0.3, 0.5),
                        labels = c(10, 30, 50)) +
  theme_minimal(base_family = "Helvetica", 
                base_size = 13)+
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=22),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank()) + 
  ylab("Density")

p6

ggsave(snakemake@output[[5]], p6, width = 5, height = 2.5)


# Export window files: 
write.table(running_roh_p_imputed_nogrey, snakemake@output[[6]], sep = "\t", row.names = FALSE, quote = FALSE)
write.table(running_roh_p_modern_nogrey, snakemake@output[[7]], sep = "\t", row.names = FALSE, quote = FALSE)



########### Ancient imputed dogs:

#### DESERTS
imputed_low <- running_roh_p_imputed_nogrey %>%
  filter(UNAFF_mean<0.05) %>%
  select(region)

# prepare bed format for bedtools intersect:
imputed_low_bed <- imputed_low %>% 
  separate_wider_delim(region, ":", names=c('chr','region')) %>% 
  separate_wider_delim(region, "-", names=c('start','end'))

imputed_low_bed$chr <- as.numeric(imputed_low_bed$chr)
imputed_low_bed$start <- as.numeric(imputed_low_bed$start)
imputed_low_bed$end <- as.numeric(imputed_low_bed$end)

# sort 
imputed_low_bed_sorted <- imputed_low_bed %>%
  arrange(chr, start)


########### Modern  dogs:

#### DESERTS
modern_low <- running_roh_p_modern_nogrey %>%
  filter(UNAFF_mean<0.05) %>%
  select(region)

# prepare bed format for bedtools intersect:
modern_low_bed <- modern_low %>% 
  separate_wider_delim(region, ":", names=c('chr','region')) %>% 
  separate_wider_delim(region, "-", names=c('start','end'))

modern_low_bed$chr <- as.numeric(modern_low_bed$chr)
modern_low_bed$start <- as.numeric(modern_low_bed$start)
modern_low_bed$end <- as.numeric(modern_low_bed$end)

# sort 
modern_low_bed_sorted <- modern_low_bed %>%
  arrange(chr, start)

# export bed files of windows for gene annotation:
write.table(imputed_low_bed_sorted, snakemake@output[[8]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(modern_low_bed_sorted, snakemake@output[[9]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


##############################################################
# FIND COMMON OVERLAPS BETWEEN ANCIENT AND MODERN for DESERTS:

modern_imputed_overlap_deserts <- as.data.frame(intersect(modern_low$region,imputed_low$region))
colnames(modern_imputed_overlap_deserts) <- 'Common_regions'

#prepare bed format for bedtools intersect:
modern_imputed_overlap_deserts_bed <- modern_imputed_overlap_deserts %>% 
  separate_wider_delim(Common_regions, ":", names=c('chr','region')) %>% 
  separate_wider_delim(region, "-", names=c('start','end'))

modern_imputed_overlap_deserts_bed$chr <- as.numeric(modern_imputed_overlap_deserts_bed$chr)
modern_imputed_overlap_deserts_bed$start <- as.numeric(modern_imputed_overlap_deserts_bed$start)
modern_imputed_overlap_deserts_bed$end <- as.numeric(modern_imputed_overlap_deserts_bed$end)

#sort 
modern_imputed_overlap_deserts_bed_sorted <- modern_imputed_overlap_deserts_bed %>%
  arrange(chr, start)

#export windows for gene annotation:
write.table(modern_imputed_overlap_deserts_bed_sorted, snakemake@output[[10]], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


####################################### merge_plots  ######################################

#supp figure heatmap:
png(snakemake@output[[11]], width=17, height=10, units='in', res=200, pointsize=4)
ggarrange(p3, p5,
          labels = c("a", "b"),
          ncol = 2, nrow = 1, font.label=list(size=20))
dev.off()


# Export pdf file for main text
ggarrange(p3, p5,
          labels = c("A", "B"),
          ncol = 2, nrow = 1, font.label=list(size=20))

ggsave(snakemake@output[[12]], width=17, height=10, units='in', dpi=300)
