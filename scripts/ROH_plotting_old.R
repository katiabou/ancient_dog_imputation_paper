#! /usr/bin/env Rscript

#import libraries
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
chrom <- args[6]

a <- str_split(list_of_concordance_phased, pattern=',')
b <- str_split(list_of_validation, pattern=',')
f <- str_split(list_of_phased, pattern=',')

a[[1]][7] <- NA
a[[1]][7] <- b[[1]]
a[[1]][8] <- NA
a[[1]][8] <- f[[1]]

d <- c()

for (i in 1:length(a[[1]])){
  tmp = read.csv(a[[1]][i], sep="") 
  d <- rbind(d,tmp)
}

#chr length
chr_length <- read.delim(args[5], header=FALSE) 
size_chr <- chr_length[1,2]

#plot only sample
sub <- subset(d, d$FID==name) 

#fill in HC_genotyped if no ROHs are found
new_sub <- rbind(sub, 
      sub %>% 
        filter(cov!="HC_genotyped") %>%
        mutate(cov = "HC_genotyped",
               POS1=0,
               POS2=0,
               KB=0,
               NSNP=0,
               DENSITY=0,
               PHOM=0,
               PHET=0))


#make new column with name and cov
new_sub$name <- paste(new_sub$FID, new_sub$cov, sep = '_')

#put same cov for imputd and genotyped HC:
new_sub$cov[new_sub$cov=='HC_imputed'] <- 'HC'
new_sub$cov[new_sub$cov=='HC_genotyped'] <- 'HC'

#replace name 
new_sub$name <- gsub(paste(name, "_", sep=""),"", new_sub$name)
new_sub$name <- gsub("x","x imputed", new_sub$name)
new_sub$name <- gsub("HC_genotyped","HC", new_sub$name)
new_sub$name <- gsub("HC_imputed","HC imputed", new_sub$name)

#reorder names
new_sub$name <- factor(new_sub$name, levels=c("0.01x imputed", "0.05x imputed",
                                              "0.1x imputed", "0.5x imputed",
                                              "1x imputed", "2x imputed",
                                              "HC imputed", "HC"))


#plot
png(args[7], width=16, height=4, units='in', res=250, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
options(scipen=10000)
ggplot(data = new_sub)+
  geom_segment(aes(y = name, yend = name, x = POS1, xend = POS2, colour=cov), linewidth = 10) +
  scale_colour_viridis(discrete = TRUE, option='D') + 
  scale_x_continuous(breaks = seq(0, size_chr, 10000000))+
  theme_classic()+
  theme(legend.position = "none")+
  xlab(paste("Genomic position ",chrom," (bp)", sep=""))+
  theme(axis.title.y=element_blank())+
  theme(axis.text = element_text(size = 11))
  
dev.off()

