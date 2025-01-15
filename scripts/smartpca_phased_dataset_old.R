library(dplyr)
library(ggplot2)
library(MetBrewer)
library(ggrepel)

#import files
eigenval_output <- read.table(snakemake@input[[1]], quote="\"", comment.char="")
eigenvec_output <- read.table(snakemake@input[[2]], quote="\"", comment.char="")
fam <- read.table(snakemake@input[[3]], quote="\"", comment.char="")
meta <- read.delim(snakemake@input[[4]])
Meso_dogs <- read.delim(snakemake@input[[5]])
meta$data <- 'reference_panel'
canid_group <- as.character(snakemake@params[["canid_group"]])

#prepare metadata for modern
ref_meta <- meta %>% select(Sample, Dog_PCA)
colnames(ref_meta) <- c('Sample', 'Group')
ref_meta$type <- 'modern'

#prepare metadata for imputed
new_meta <- Meso_dogs %>% select(name_haplo_VCF, Meta.Population)
colnames(new_meta) <- c('Sample', 'Group')
new_meta$type <- 'imputed'

#merge modern and ancient meta
all_meta <- rbind(new_meta, ref_meta)

#chose only metadata for our samples
samples <- as.data.frame(fam$V1)
colnames(samples) <- 'Sample'

final_df <- dplyr::inner_join(samples, all_meta, by='Sample')


#merge to the EMU output
eigenvec_output_meta <- cbind(eigenvec_output, final_df)


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

colnames(eigenvec_output_meta)[1:4] <- c("PC1", "PC2", "PC3", "PC4")

#assign colours to specific groups
values = c("African_NearEast_India_Dogs" = "#04a3bd", 
           "Americas_Dogs" = "#f0be3d", 
           "Arctic_Dogs" = "#931e18",
           "East_Asian_Dogs" = "#da7901",
           "European_Dogs"="#247d3f", 
           "preContact_Dogs"="#20235b", 
           "Eastern_Eurasian_Wolves"="maroon2", 
           "North_American_Wolves"='grey', 
           "Pleistocene_Wolves"='lightgreen', 
           "Western_Eurasian_Wolves"='dodgerblue', 
           "Wolves"='purple2', 
           "Undet"='pink')

#specify order in legend
order_d <- c("African_NearEast_India_Dogs", "Americas_Dogs", "Arctic_Dogs",
                                              "East_Asian_Dogs","European_Dogs", "preContact_Dogs",
                                              "Eastern_Eurasian_Wolves","North_American_Wolves","Pleistocene_Wolves",
                                              "Western_Eurasian_Wolves","Wolves","Undet")

#define modern and ancient layers
df_layer_1 <- eigenvec_output_meta[eigenvec_output_meta$type=="modern",]
df_layer_2 <- eigenvec_output_meta[eigenvec_output_meta$type!="modern",]

#colours <- c("#04a3bd", "#f0be3d", "#931e18", "#da7901", "#247d3f", "#20235b", "pink", 'grey', 'lightgreen','dodgerblue','purple2','maroon2')

plot_title <- paste(canid_group, "PCA", sep=" ")

png(snakemake@output[[1]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(eigenvec_output_meta, aes(x=PC1, y=PC2)) +
  geom_point(data = df_layer_1, aes(col=Group, shape=type), size = 1, alpha=0.3)+
  geom_point(data = df_layer_2, aes(col=Group, shape=type), size = 2, alpha=0.7)+
  scale_shape_manual(values=c(19,8), labels=c('Imputed','Reference panel'))+ 
  #geom_label_repel(data = eigenvec_output_meta %>% filter(Group=='African_NearEast_India_Dogs'), 
  #                 aes(x=PC1, y=PC2, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  #scale_color_manual(values = colours,labels=c('Africa, Near East, India', 'Americas','Arctic','East Asia','Europe','Americas preContact', 'Undet')) + 
  scale_color_manual(values = values, breaks=order_d) + 
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
  ggtitle(plot_title)
dev.off()

png(snakemake@output[[2]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(eigenvec_output_meta, aes(x=PC3, y=PC4)) +
  geom_point(data = df_layer_1, aes(col=Group, shape=type), size = 1, alpha=0.3)+
  geom_point(data = df_layer_2, aes(col=Group, shape=type), size = 2, alpha=0.7)+
  scale_shape_manual(values=c(19,8), labels=c('Imputed','Reference panel'))+ 
  #geom_label_repel(data = eigenvec_output_meta %>% filter(Group=='East_Asian_Dogs'), 
  #                 aes(x=PC3, y=PC4, label=Sample),size=1.3, box.padding = 1, max.overlaps = Inf)+
  #scale_color_manual(values = colours,labels=c('Africa, Near East, India', 'Americas','Arctic','East Asia','Europe','Americas preContact', 'Undet')) + #yes
  scale_color_manual(values = values, breaks=order_d) + 
  labs(x = paste("PC3 (",PC3,"%)", sep = ""), y = paste("PC4 (",PC4,"%)", sep = ""))+
  theme_classic()+
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  labs(shape="Genotype data", colour="Group")+
  ggtitle(plot_title)
  dev.off()


png(snakemake@output[[3]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(eigenvec_output_meta, aes(x=PC1, y=PC2)) +
  geom_point(data = df_layer_1, aes(col=Group, shape=type), size = 1, alpha=0.3)+
  geom_point(data = df_layer_2, aes(col=Group, shape=type), size = 2, alpha=0.7)+
  scale_shape_manual(values=c(19,8), labels=c('Imputed','Reference panel'))+ 
  geom_label_repel(data = eigenvec_output_meta %>% filter(type=='imputed'), 
                   aes(x=PC1, y=PC2, label=Sample),size=1.3, box.padding = 1, max.overlaps = Inf)+
  #scale_color_manual(values = colours,labels=c('Africa, Near East, India', 'Americas','Arctic','East Asia','Europe','Americas preContact', 'Undet')) + 
  scale_color_manual(values = values, breaks=order_d) + 
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
  ggtitle(plot_title)
dev.off()

png(snakemake@output[[4]], width=9, height=6, units='in', res=200, pointsize=4)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 2,
  cex.lab  = 2)
ggplot(eigenvec_output_meta, aes(x=PC3, y=PC4)) +
  geom_point(data = df_layer_1, aes(col=Group, shape=type), size = 1, alpha=0.3)+
  geom_point(data = df_layer_2, aes(col=Group, shape=type), size = 2, alpha=0.7)+
  scale_shape_manual(values=c(19,8), labels=c('Imputed','Reference panel'))+ 
  geom_label_repel(data = eigenvec_output_meta %>% filter(type=='imputed'), 
                   aes(x=PC1, y=PC2, label=Sample),size=1.3, box.padding = 1, max.overlaps = Inf)+
  #scale_color_manual(values = colours,labels=c('Africa, Near East, India', 'Americas','Arctic','East Asia','Europe','Americas preContact', 'Undet')) + #yes
  scale_color_manual(values = values, breaks=order_d) + 
  labs(x = paste("PC3 (",PC3,"%)", sep = ""), y = paste("PC4 (",PC4,"%)", sep = ""))+
  theme_classic()+
  geom_vline(xintercept = 0, size=0.1, linetype = "dashed")+
  geom_hline(yintercept = 0, size=0.1, linetype = "dashed")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),  
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5))+
  labs(shape="Genotype data", colour="Group")+
  ggtitle(plot_title)
  dev.off()