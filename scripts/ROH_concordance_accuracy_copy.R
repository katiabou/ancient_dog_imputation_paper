#! /usr/bin/env Rscript

#import libraries
library(stringr)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(viridis)

#import data
args <- commandArgs(trailingOnly = TRUE)

list_of_concordance_phased <- args[1]
list_of_validation <- args[2]
list_of_phased <- args[3]
sample_name <- args[4]
cov_sample  <- args[5]
info_sample  <- args[6]
site_type <- args[7]
ref_fasta_chr_size <- read.delim(args[13], header=FALSE)


aa <- str_split(list_of_concordance_phased, pattern=',')
bb <- str_split(list_of_validation, pattern=',')
ff <- str_split(list_of_phased, pattern=',')


a = read.csv(aa[[1]][1], sep="") 
b = read.csv(aa[[1]][2], sep="") 
c = read.csv(aa[[1]][3], sep="") 
d = read.csv(aa[[1]][4], sep="") 
e = read.csv(aa[[1]][5], sep="") 
f = read.csv(aa[[1]][6], sep="") 

g = read.csv(ff[[1]][1], sep="") 

validation = read.csv(bb[[1]][1], sep="") 

#remove extra line used for plotting (in case no ROH was found) and use GRanges to create approprite format
get_gr <- function(roh){
  test <- roh %>% filter(IID != 0 & IID==sample_name) %>%
  mutate(range = paste(POS1, "-", POS2, sep="")) 
  GRanges(
  seqnames=test$CHR,
  ranges=test$range,
  sample=test$FID,
  cov=test$cov
)
}

get_val <- get_gr(validation)
get_a <- get_gr(a)
get_b <- get_gr(b)
get_c <- get_gr(c)
get_d <- get_gr(d)
get_e <- get_gr(e)
get_f <- get_gr(f)
get_g <- get_gr(g)

#find overlaps between the two samples
hits_a <-findOverlaps(get_val,get_a)
hits_b <-findOverlaps(get_val,get_b)
hits_c <-findOverlaps(get_val,get_c)
hits_d <-findOverlaps(get_val,get_d)
hits_e <-findOverlaps(get_val,get_e)
hits_f <-findOverlaps(get_val,get_f)
hits_g <-findOverlaps(get_val,get_g)


#get overlaps
get_overlap <- function(query, subject,hits){
  overlaps<- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
  overlaps$proportion_subject <- width(overlaps) / width(subject[subjectHits(hits)]) #this is the proportion of the subject that is contained within the query
  overlaps$proportion_query <- width(overlaps) / width(query[queryHits(hits)]) #this is the proportion of the query that is overlapping with the subject
  
  return(overlaps)
}

#ab <- get_overlap(gr, gr1, hits)
a_overlap <- get_overlap(get_val, get_a, hits_a)
b_overlap <- get_overlap(get_val, get_b, hits_b)
c_overlap <- get_overlap(get_val, get_c, hits_c)
d_overlap <- get_overlap(get_val, get_d, hits_d)
e_overlap <- get_overlap(get_val, get_e, hits_e)
f_overlap <- get_overlap(get_val, get_f, hits_f)
g_overlap <- get_overlap(get_val, get_g, hits_g)



############ SEGMENT BASED ESTIMATION ###########

## Now get stats for segment based overlaps:

get_seg_stats <- function(hits, get_phased, get_val, df){
  #ROHs identified in the phased_target that are also present in the val (so are actually true)
  if (length(unique(queryHits(hits))) != 0){
  TP <- length(unique(queryHits(hits)))
  } else {
    TP <- 1e-16
  }
  
  #false positives will be a ROH found in the phased_target that was not present in the val 
  if (length(get_phased[!(get_phased %in% get_phased[subjectHits(hits)])]) != 0){
  FP <- length(get_phased[!(get_phased %in% get_phased[subjectHits(hits)])])
  } else {
    FP <- 1e-16
  }
  #false negatives will be a ROH not found in the phased_target that was present in the val
  if (length(get_val[!(get_val %in% get_val[queryHits(hits)])]) != 0){
  FN <- length(get_val[!(get_val %in% get_val[queryHits(hits)])])
  } else {
    FN <- 1e-16
  }
  
  ## Now get F1 stat for segment based overlaps:
  recall_seg = TP/(TP+FN)
  precision_seg =  TP/(TP+FP)
  F1_seg=2*(precision_seg*recall_seg/(precision_seg+recall_seg))
  
  #False discovery rate:
  FDR <- FP / (TP + FP)
  
  #summarize:
  all_stats_segment <- data.frame(sample=sample_name,
                          #cov=unique(get_phased$cov),
                          cov=unique(df$cov),
                          TP=TP,
                          FP=FP,
                          FN=FN,
                          sensitivity=recall_seg,
                          precision=precision_seg,
                          F1=F1_seg,
                          FDR=FDR)
}

get_seg_stats_a <- get_seg_stats(hits_a, get_a, get_val, a)
get_seg_stats_b <- get_seg_stats(hits_b, get_b, get_val, b)
get_seg_stats_c <- get_seg_stats(hits_c, get_c, get_val, c)
get_seg_stats_d <- get_seg_stats(hits_d, get_d, get_val, d)
get_seg_stats_e <- get_seg_stats(hits_e, get_e, get_val, e)
get_seg_stats_f <- get_seg_stats(hits_f, get_f, get_val, f)
get_seg_stats_g <- get_seg_stats(hits_g, get_g, get_val, g)



all_seg_stats <- rbind(get_seg_stats_a,
                       get_seg_stats_b,
                       get_seg_stats_c,
                       get_seg_stats_d,
                       get_seg_stats_e,
                       get_seg_stats_f,
                       get_seg_stats_g)





############ TOTAL LENGTH BASED ESTIMATION ###########
genome_size <- ref_fasta_chr_size[1,2]

library(mltools)

get_length_stats <- function(overlap, gen_size, get_phased, hits, get_val, df){
  # TP should be the total ROH of the subject overlapping the validation
  TP_overlaps_sum <- sum(width(overlap))
  if (TP <- TP_overlaps_sum / gen_size != 0){
  TP <- TP_overlaps_sum / gen_size
  } else {
    TP <- 1e-16
  }
   
  # FP should be the total ROH of the subject minus the TP:
  FP_overlaps_sum <- sum(width(get_phased)) - sum(width(overlap))
  if (FP <- FP_overlaps_sum / gen_size != 0){
  FP <- FP_overlaps_sum / gen_size
  } else {
    FP <- 1e-16
  }

  # FN should be the total ROH of the query minus the TP:
  FN_overlaps_sum <- sum(width(get_val)) - sum(width(overlap))
  if (FN <- FN_overlaps_sum / gen_size != 0){
  FN <- FN_overlaps_sum / gen_size
  } else {
    FN <- 1e-16
  }
  
  # TN should be the total genome size minus the TP, FP, FN:
  TN_sum <- gen_size - TP_overlaps_sum - FP_overlaps_sum - FN_overlaps_sum
  TN <- TN_sum / gen_size
  
  recall_length = TP/(TP+FN)
  precision_length =  TP/(TP+FP)
  F1_length=2*(precision_length*recall_length/(precision_length+recall_length))
  
  mcc <- mcc(TP=TP, FP=FP, TN=TN, FN=FN)
  
  mcc_norm <- (mcc+1)/2 # IF 0 until <0.5 it's worse than random, 0.5 random and >0.5 is better estimations
  
  #accuracy (worst value: 0; best value: 1):
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  #Specificity 
  specificity <- TN / (TN+FP)
  
  #False discovery rate:
  FDR <- FP / (TP + FP)

  all_stats_length <- data.frame(sample=sample_name,
                          #cov=unique(get_phased$cov),
                          cov=unique(df$cov),
                          TP=TP,
                          FP=FP,
                          FN=FN,
                          TN=TN,
                          sensitivity=recall_length,
                          precision=precision_length,
                          F1=F1_length,
                          FDR=FDR,
                          mcc=mcc,
                          mcc_norm=mcc_norm,
                          accuracy=accuracy,
                          specificity=specificity)
}
  
get_length_stats_a <- get_length_stats(a_overlap, genome_size, get_a, hits_a, get_val, a)
get_length_stats_b <- get_length_stats(b_overlap, genome_size, get_b, hits_b, get_val, b)
get_length_stats_c <- get_length_stats(c_overlap, genome_size, get_c, hits_c, get_val, c)
get_length_stats_d <- get_length_stats(d_overlap, genome_size, get_d, hits_d, get_val, d)
get_length_stats_e <- get_length_stats(e_overlap, genome_size, get_e, hits_e, get_val, e)
get_length_stats_f <- get_length_stats(f_overlap, genome_size, get_f, hits_f, get_val, f)
get_length_stats_g <- get_length_stats(g_overlap, genome_size, get_g, hits_g, get_val, g)


all_length_stats <- rbind(
                       get_length_stats_a,
                       get_length_stats_b,
                       get_length_stats_c,
                       get_length_stats_d,
                       get_length_stats_e,
                       get_length_stats_f,
                       get_length_stats_g)


#export both files:
write.table(all_seg_stats, file=args[8], quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
write.table(all_length_stats, file=args[9], quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)


#merge both df together (segment and length) based on common columns (mcc is only for length)
all_seg_stats$type <- 'Segments'
all_seg_stats$mcc_norm <- NA
all_length_stats$type <- 'Length'

common_cols <- intersect(colnames(all_seg_stats), colnames(all_length_stats))
all_stats <- rbind(
  subset(all_seg_stats, select = common_cols), 
  subset(all_length_stats, select = common_cols)
)

#melt df
all_stats_final <- all_stats %>%
  select(sample, cov, F1, mcc_norm, type) %>%
  gather(metric, values, -sample, -cov, -type)

all_stats_final <- na.omit(all_stats_final)

## plotting

name_HC_imputed <- paste("HC (",cov_sample,'x)',sep="")
all_stats_final$cov <- gsub("HC_imputed",name_HC_imputed, all_stats_final$cov)
all_stats_final$metric <- gsub("mcc_norm", "nMCC", all_stats_final$metric)

name_title_2 <- gsub('_',' ',info_sample)
name_title_final <- paste(name_title_2, ' - ', sample_name, " (", site_type, ")", sep='')


#F1 and nMCC (segment and length together)
png(args[10], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(all_stats_final, aes(x=cov, y=values, colour=type))+
  geom_line(aes(group = interaction(type, metric), linetype=metric), linewidth = 1.2)+
  geom_point(size=3)+
  scale_color_manual(values = c("steelblue","orange"))+
  ylim(0,1)+
  labs(x = 'Coverage', colour='Count type', linetype= 'Metric') +
  ggtitle(name_title_final) +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 30, size = 16, vjust = 0.5),  
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=18),
        panel.grid.minor = element_blank(),
        legend.key.width= unit(1.5, 'cm'))
dev.off()



### Looking into other accuracy metrics for length:

#melt all_stats_final df
all_stats_final <- all_length_stats %>%
  select(sample, cov, type, FDR, specificity, sensitivity) %>%
  gather(metric, values, -sample, -cov, -type)

all_stats_final <- na.omit(all_stats_final)
all_stats_final$cov <- gsub("HC_imputed",name_HC_imputed, all_stats_final$cov)

## plotting FDR, specificity and sensitivity
png(args[11], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(all_stats_final, aes(x=cov, y=values, colour=metric))+
  geom_line(aes(group = interaction(type, metric)), linewidth = 1.2)+
  geom_point(size=3)+
  scale_color_manual(values = c("gold2","darkolivegreen4", "lightpink3"))+
  ylim(0,1)+
  labs(x = 'Coverage', colour='Metric') +
  ggtitle(name_title_final) +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 30, size = 16, vjust = 0.5),  
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=18),
        panel.grid.minor = element_blank(),
        legend.key.width= unit(1.5, 'cm'))
dev.off()


## plotting specificity against sensitivity
all_length_stats$cov <- gsub("HC_imputed",name_HC_imputed, all_length_stats$cov)

png(args[12], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(all_length_stats, aes(x=specificity, y=sensitivity, colour=cov))+
  geom_point(size=3)+
  scale_colour_viridis(discrete = TRUE, option='D') + 
  ylim(0,1)+
  xlim(0,1)+
  labs(x = 'Specificity', y='Sensitivity', colour='Coverage') +
  ggtitle(name_title_final) +
  theme_bw()+
  theme(axis.text.x=element_text(size = 16),  
        axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=18), 
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size=18),
        panel.grid = element_blank(),
        legend.key.width= unit(1.5, 'cm'))
dev.off()