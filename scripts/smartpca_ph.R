library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(readr)
library("MetBrewer")
library(viridis)
#library(ggmagnify)

#sample files
name <- as.character(snakemake@params[["sample"]])
name_title <- as.character(snakemake@params[["info_sample"]])
name_title_2 <- gsub('_',' ',name_title)
cov_hc <- as.character(snakemake@params[["cov_sample"]])

samples <- read.table(snakemake@input[[1]], quote="\"", comment.char="")
eigenval_output <- read.table(snakemake@input[[2]], quote="\"", comment.char="")
eigenvec_output <- read_table(snakemake@input[[3]], col_names = FALSE)

#Plotting
colnames(samples)[1] <- "Sample"

#import modern metadata file
meta <- read.delim(snakemake@input[[4]])
colnames(meta)[3] <- 'Sample'
colnames(meta)[colnames(meta) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(meta)[colnames(meta) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"
meta$data <- 'reference_panel'
sample_modern_names_meta <- left_join(samples, meta, "Sample")


#import imputed concordance metadata file
info <- read.delim(snakemake@input[[5]])
info$data <- 'target'
colnames(info)[1] <- 'Sample'
sample_imputed_names_meta <- left_join(samples, info, "Sample")


#merge the modern and imputed datasets
final <- left_join(sample_modern_names_meta, sample_imputed_names_meta, "Sample")

#fix duplicate columns
final$Wolf_Dog_PCA <- ifelse(is.na(final$Wolf_Dog_PCA.x), final$Wolf_Dog_PCA.y, final$Wolf_Dog_PCA.x)
final$Dog_PCA <- ifelse(is.na(final$Dog_PCA.x), final$Dog_PCA.y, final$Dog_PCA.x)

#chose columns
vvv <- final %>% select(1,ncol(final)-1,ncol(final))

#fill in downsampled and HC imputed fields with matching metadata:
t <- "imputed"
ref <- "reference_panel"

vvv$Wolf_Dog_PCA <- ifelse(grepl(name, vvv$Sample), vvv$Wolf_Dog_PCA[vvv$Sample==name], vvv$Wolf_Dog_PCA)
vvv$Dog_PCA <- ifelse(grepl(name, vvv$Sample), vvv$Dog_PCA[vvv$Sample==name], vvv$Dog_PCA)
vvv$target <- ifelse(grepl(name, vvv$Sample), name, 'reference_panel')
vvv$data <- ifelse(grepl(t, vvv$Sample), 'imputed', 'non-imputed')
vvv$data <- ifelse(grepl(ref, vvv$target), 'reference_panel', vvv$data)
vvv$data <- ifelse(vvv$Sample==name & vvv$data=='non-imputed', 'genotyped', vvv$data)


##### PCA stuff

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
b[,9] <- vvv$target
colnames(b)[9] <- "target"

#fix coverage column
b$coverage <- NA

cov1 <- 0.05
b$coverage <- ifelse(grepl(paste("_",cov1,'x',sep=""), b$Sample), cov1, b$coverage)
cov2 <- 0.1
b$coverage <- ifelse(grepl(paste("_",cov2,'x',sep=""), b$Sample), cov2, b$coverage)
cov3 <- 0.2
b$coverage <- ifelse(grepl(paste("_",cov3,'x',sep=""), b$Sample), cov3, b$coverage)
cov4 <- 0.5
b$coverage <- ifelse(grepl(paste("_",cov4,'x',sep=""), b$Sample), cov4, b$coverage)
cov5 <- 1
b$coverage <- ifelse(grepl(paste("_",cov5,'x',sep=""), b$Sample), cov5, b$coverage)
cov6 <- 2
b$coverage <- ifelse(grepl(paste("_",cov6,'x',sep=""), b$Sample), cov6, b$coverage)

#add high coverage value:
b$coverage <- as.character(ifelse(b$target != 'reference_panel' & is.na(b$coverage), 'HC', b$coverage))
#b$coverage <- ifelse(b$target != 'reference_panel' & is.na(b$coverage), cov_hc, b$coverage)

#plotting
df_layer_1 <- b[b$target=="reference_panel",]
df_layer_2 <- b[b$target!="reference_panel",]


gen_data <- c(
  'imputed' = "Imputed",
  'non-imputed' = "Pseudohaploid",
  'genotyped' = "HC genotyped",
  'reference_panel' = "Reference panel"
)

name_title_final <- paste(name_title_2, ' - ', name,sep='')

png(snakemake@output[[1]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)

a <- ggplot(b, aes(x=PC1, y=PC2)) +
  geom_point(data = df_layer_1, aes(shape=data), col='grey64', size = 2.5, alpha=0.8)+
  geom_point(data = df_layer_2, aes(fill=coverage, shape=data), size = 3.5, alpha=0.9)+
  geom_point(data = df_layer_2 %>% filter(coverage=='HC'), aes(fill=coverage, shape=data), size = 3.5, alpha=0.9)+
  scale_shape_manual(values=c(23, 24, 22,8), labels=gen_data)+
  scale_fill_viridis(discrete = TRUE, option='D', labels=c('0.05x','0.1x', '0.2x', '0.5x','1x','2x',paste("HC (",cov_hc,'x)',sep=""))) + 
  labs(x = paste("PC1 (",PC1,"%)", sep = ""), y = paste("PC2 (",PC2,"%)", sep = ""))+
  theme_bw()+
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  theme(axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),  
        legend.text = element_text(size=18),
        legend.title = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size=18),
        panel.border = element_blank(),
        axis.line = element_line(colour = "gray60"))+
  labs(shape="Genotype data", fill="Coverage")+
  guides(fill=guide_legend(override.aes=list(shape=c(21))))+
  ggtitle(name_title_final)

a

dev.off()


