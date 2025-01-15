#! /usr/bin/env Rscript

#import libraries
library(devtools)
library(MetBrewer)
library(reshape2)
library(dplyr)
library(ggplot2)
library(qdapRegex)
library(stringr)


input1 <- read.table(snakemake@input[[1]], quote="\"", comment.char="")
input2 <- read.table(snakemake@input[[2]], quote="\"", comment.char="")
input3 <- read.table(snakemake@input[[3]], quote="\"", comment.char="")
input4 <- read.table(snakemake@input[[4]], quote="\"", comment.char="")

input1$ref_panel <- 'all_canids'
input2$ref_panel <- 'all_canids'
input3$ref_panel <- 'dogs_only'
input4$ref_panel <- 'dogs_only'

input <- rbind(input1, input2, input3, input4)

#fixing the MAF_bin column
test <- sapply(strsplit(input$V1, "/\\s*"), tail, 1)
t1 <- input$V1
t2 <- ex_between(test,"MAF_", ".vcf.gz")
t3 <- do.call(rbind.data.frame, t2)
colnames(t3) <- 'MAF_bin'
t3$MAF_bin <- str_replace_all(t3$MAF_bin, "_", ".")
t4 <- data.frame(sub("^([^.]*.[^.]*).", "\\1-",t3$MAF_bin)) 
colnames(t4) <- 'MAF_bin'
t4$MAF_bin[t4$MAF_bin == '0.0-001'] <- '0-0.001'
input$V8 <- t4$MAF_bin
input$V9 <- paste(input$V2, input$V4, sep="_")
colnames(input) <- c('file_name', 'Sample', 'Chr','Coverage','INFO','Sites', 'Ref_panel', 'MAF_bins', 'Sample_cov')

input$INFO <- paste("INFO", input$INFO, sep=" ")
input$INFO[input$INFO == 'INFO 0'] <- 'No INFO cutoff'

#transform
input_b <- input %>% select('Sample_cov', 'INFO','Sites', 'Ref_panel', 'MAF_bins')
melt_data2 <- melt(input_b, id = c("MAF_bins","Sample_cov","Ref_panel", "INFO")) 

cols = met.brewer(name="Demuth", n=8, type="discrete") #cb friendly nice

# New facet label names for reference variable
ref.labs <- c("All canid reference panel", "Dog reference panel")
names(ref.labs) <- c("all_canids", "dogs_only")

sample.labs <- c("Pleistocene wolf (CGG32 1x)", "Neolithic European dog (NGDG 0.5x)")
names(sample.labs) <- c("CGG32_1", "NGDG_0.5")

#change bar order
melt_data2$INFO <- factor(melt_data2$INFO,                                    
                              levels = c("No INFO cutoff", "INFO 0.8", "INFO 0.9", "INFO 0.95"))

#re-scale number of sites
melt_data2 <- melt_data2 %>%
  mutate(value_2 = value / 1000)

#barplot
ggplot(data=melt_data2, aes(x=INFO, y=value_2, fill=MAF_bins, label = round(value_2, digits = 0))) +
  geom_bar(stat="identity")+
  geom_text(size = 3.5, position = position_stack(vjust = 0.5))+
  geom_text(
    aes(label = round(after_stat(y), digits = 0), group = INFO), 
    stat = 'summary', fun = sum, vjust = -0.5
  )+
  facet_grid(Sample_cov ~ Ref_panel, labeller=labeller(Ref_panel=ref.labs, Sample_cov=sample.labs))+
  scale_fill_manual(values = cols, name = "MAF bins")+
  #scale_y_continuous(breaks = seq(0, 1600, by = 100))+
  ylab(bquote("Number of sites " (10^3)))+
  theme_bw()+
  theme(axis.text.x=element_text(size = 14, vjust = 0.5)) +
  theme(axis.text.y=element_text(size = 14, vjust = 0.5)) + 
  theme(axis.title=element_text(size=14)) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_text(size=16)) + 
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  theme(strip.background =element_rect(fill="gray28"))+
  theme(strip.text = element_text(colour = 'white', size=14))+
  theme(legend.spacing.y = unit(0.3, 'cm'))  +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.key.size = unit(0.8, "cm"))
ggsave(snakemake@output[[1]], width = 14, height = 10)

