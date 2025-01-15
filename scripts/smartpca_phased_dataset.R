library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(readr)
library("MetBrewer")


##### PCA dogs 
#samples <- read.table("~/Downloads/merged_phased_ref_panel.all_chr_MAF_0.01_recalibrated_INFO_0.8_dogs.ind", quote="\"", comment.char="")
#eigenval_output <- read.table("~/Downloads/merged_phased_ref_panel.all_chr_MAF_0.01_recalibrated_INFO_0.8_dogs_eigenval_output", quote="\"", comment.char="")
#eigenvec_output <- read_table("~/Downloads/merged_phased_ref_panel.all_chr_MAF_0.01_recalibrated_INFO_0.8_dogs_eigenvec_output", col_names = FALSE)

samples <- read.table(snakemake@input[[1]], quote="\"", comment.char="")
eigenval_output <- read.table(snakemake@input[[2]], quote="\"", comment.char="")
eigenvec_output <- read_table(snakemake@input[[3]], col_names = FALSE)

colnames(samples)[1] <- "Sample"

#import modern metadata file
meta <- read.delim(snakemake@input[[7]])
#meta <- read.delim('~/Downloads/Dog_Wolf_aDNA_WG-Modern.tsv')
colnames(meta)[3] <- 'Sample'
colnames(meta)[colnames(meta) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(meta)[colnames(meta) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"
meta$data <- 'reference_panel'

sample_modern_names_meta <- left_join(samples, meta, "Sample")


#import imputed metadata file
info <- read.delim(snakemake@input[[8]])
#info <- read.delim('~/Downloads/Dog_Wolf_aDNA_WG-Master.tsv')
colnames(info)[colnames(info) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(info)[colnames(info) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"

info$data <- 'target'
colnames(info)[1] <- 'Sample'
sample_imputed_names_meta <- left_join(samples, info, "Sample")

#merge the modern and imputed datasets
final <- left_join(sample_modern_names_meta, sample_imputed_names_meta, "Sample")

#fix duplicate columns
final$Wolf_Dog_PCA <- ifelse(is.na(final$Wolf_Dog_PCA.x), final$Wolf_Dog_PCA.y, final$Wolf_Dog_PCA.x)
final$Dog_PCA <- ifelse(is.na(final$Dog_PCA.x), final$Dog_PCA.y, final$Dog_PCA.x)
final$data <- ifelse(is.na(final$data.x), final$data.y, final$data.x)

#chose columns
vvv <- as.data.frame(final[,c(1,ncol(final)-2,ncol(final)-1,ncol(final))])


#make list with % of each PC:
mylist<-c()
for (s in eigenval_output$V1) {
  print(s / sum(eigenval_output$V1))
  a<-(s / sum(eigenval_output$V1))
  mylist <- c(mylist, a)
}

#round numbers:
PC1 <- round(mylist[1]*100, digits=2)
PC2 <- round(mylist[2]*100, digits=2)
PC3 <- round(mylist[3]*100, digits=2)
PC4 <- round(mylist[4]*100, digits=2)

fn <- eigenvec_output[-1,]
colnames(fn) <- c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop")

#choose the number of eigenvectors (let's say 4 so 5 along with ind column)
b <- as.data.frame(fn[,c(1:5)])

#add dog info to eigenvector file
b[,6] <- vvv$Wolf_Dog_PCA
colnames(b)[6] <- "Wolf_Dog_PCA"
b[,7] <- vvv$Dog_PCA
colnames(b)[7] <- "Dog_PCA"
b[,8] <- vvv$data
colnames(b)[8] <- "data"


b$plot <- ifelse(b$data=="imputed", b$Dog_PCA, 'reference_panel')
#b$target <- ifelse(grepl(name, b$Sample), b$Sample[b$Sample==name], b$data)
b$type_group <- paste(b$data,b$Dog_PCA,sep = "_")

#plotting
#this is the best one!
df_layer_1 <- b[b$data=="reference_panel",]
df_layer_2 <- b[b$data!="reference_panel",]


png(snakemake@output[[1]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)

ggplot(b, aes(x=PC1, y=PC2)) +
  geom_point(data = df_layer_1, aes(col=Dog_PCA, shape=data), size = 2, alpha=0.7)+
  geom_point(data = df_layer_2, aes(col=Dog_PCA, shape=data), size = 3, alpha=0.8)+
  scale_shape_manual(values=c(8, 19), labels=c('Imputed','Reference panel'))+ 
  scale_color_manual(values = met.brewer("Lakota"),labels=c('Africa, Near East, India', 'Americas','Arctic','East Asia','Europe','Americas preContact')) + #yes
  labs(x = paste("PC1 (",PC1,"%)", sep = ""), y = paste("PC2 (",PC2,"%)", sep = ""))+
  theme_classic()+
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  labs(shape="Genotype data", colour="Group")+
  ggtitle("Dog PCA")
dev.off()



#################
### PCA wolves

samples <- read.table(snakemake@input[[4]], quote="\"", comment.char="")
eigenval_output <- read.table(snakemake@input[[5]], quote="\"", comment.char="")
eigenvec_output <- read_table(snakemake@input[[6]], col_names = FALSE)

colnames(samples)[1] <- "Sample"


#import modern metadata file
meta <- read.delim(snakemake@input[[7]])
#meta <- read.delim('~/Downloads/Dog_Wolf_aDNA_WG-Modern.tsv')
colnames(meta)[3] <- 'Sample'
colnames(meta)[colnames(meta) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(meta)[colnames(meta) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"
meta$data <- 'reference_panel'

sample_modern_names_meta <- left_join(samples, meta, "Sample")

#import imputed metadata file
info <- read.delim(snakemake@input[[8]])
#info <- read.delim('~/Downloads/Dog_Wolf_aDNA_WG-Master.tsv')
colnames(info)[colnames(info) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(info)[colnames(info) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"

info$data <- 'target'
colnames(info)[1] <- 'Sample'
sample_imputed_names_meta <- left_join(samples, info, "Sample")


#replace missmatching names:
info$Sample[info$Sample=='Tumat2'] <- 'Tumat'

#merge with metadata
sample_imputed_names_meta <- left_join(samples, info, "Sample")


#merge modern and imputed
final <- left_join(sample_modern_names_meta, sample_imputed_names_meta, "Sample")

final$Wolf_Dog_PCA <- ifelse(is.na(final$Wolf_Dog_PCA.x), final$Wolf_Dog_PCA.y, final$Wolf_Dog_PCA.x)
final$Dog_PCA <- ifelse(is.na(final$Dog_PCA.x), final$Dog_PCA.y, final$Dog_PCA.x)
final$data <- ifelse(is.na(final$data.x), final$data.y, final$data.x)

#vvv <- final %>% select(1,ncol(final)-2,ncol(final)-1,ncol(final))
vvv <- as.data.frame(final[,c(1,ncol(final)-2,ncol(final)-1,ncol(final))])


#make list with % of each PC:
mylist<-c()
for (s in eigenval_output$V1) {
  print(s / sum(eigenval_output$V1))
  a<-(s / sum(eigenval_output$V1))
  mylist <- c(mylist, a)
}

#round numbers:
PC1 <- round(mylist[1]*100, digits=2)
PC2 <- round(mylist[2]*100, digits=2)
PC3 <- round(mylist[3]*100, digits=2)
PC4 <- round(mylist[4]*100, digits=2)

fn <- eigenvec_output[-1,]
colnames(fn) <- c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop")

#choose the number of eigenvectors (let's say 4 so 5 along with ind column)
b <- as.data.frame(fn[,c(1:5)])

#add info to eigenvector file
b[,6] <- vvv$Wolf_Dog_PCA
colnames(b)[6] <- "Wolf_Dog_PCA"
b[,7] <- vvv$Dog_PCA
colnames(b)[7] <- "Dog_PCA"
b[,8] <- vvv$data
colnames(b)[8] <- "data"

b$plot <- ifelse(b$data=="imputed", b$Wolf_Dog_PCA, 'reference_panel')
b$type_group <- paste(b$data,b$Wolf_Dog_PCA,sep = "_")

#plotting

#focus on this one
df_layer_1 <- b[b$data=="reference_panel",]
df_layer_2 <- b[b$data!="reference_panel",]

png(snakemake@output[[2]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(b, aes(x=PC1, y=PC2)) +
  geom_point(data = df_layer_1, aes(col=Wolf_Dog_PCA, shape=data), size = 3, alpha=0.7)+
  geom_point(data = df_layer_2, aes(col=Wolf_Dog_PCA, shape=data), size = 3, alpha=0.8)+
  scale_shape_manual(values=c(8, 19), labels=c('Imputed','Reference panel'))+
  scale_color_manual(values = met.brewer("Hokusai3"),labels=c('Eastern Eurasian wolves', 'North American wolves',
                                                              'Pleistocene wolves','West Eurasian wolves')) + #yes
  labs(x = paste("PC1 (",PC1,"%)", sep = ""), y = paste("PC2 (",PC2,"%)", sep = ""))+
  theme_classic()+
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  labs(shape="Genotype data", colour="Group")+
  ggtitle("Wolf PCA")
dev.off()










