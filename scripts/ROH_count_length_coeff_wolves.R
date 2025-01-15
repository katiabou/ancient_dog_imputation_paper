#############################################################
#ROH counts and length on imputed and modern dogs and wolves#
#############################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(readxl)
library(readr)
library("MetBrewer")


#import imputed ROH data and info
#info <- read.delim('~/Desktop/Copenhagen_PhD/files/IMPUTATION_2023/IMPUTATION_ROH_PAPER/files/Dog_Wolf_aDNA_WG-Master.tsv')
#roh <- read.csv('~/Downloads/merged_phased.allchrom_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogwolf.hom', sep="")
info <- read.delim(snakemake@input[[1]])
roh <- read.csv(snakemake@input[[3]], sep="")

colnames(info)[colnames(info) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(info)[colnames(info) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"
colnames(info)[1] <- "Sample"
colnames(roh)[1] <- "Sample"

#fix Tumat (Tumat2) and WolfHead (Wolf_head_IN18-016)
roh$Sample[roh$Sample == 'Tumat'] <- 'Tumat2'
roh$Sample[roh$Sample == 'WolfHead'] <- 'Wolf_head_IN18-016'


#merge two based on specific columns
final <- left_join(roh, info %>% dplyr::select(Sample, Wolf_Dog_PCA, Dog_PCA, Meta.Population, Species, Age_Mean_BP), "Sample")


#import modern ROH data and info
#info_ref <- read.delim('~/Desktop/Copenhagen_PhD/files/IMPUTATION_2023/IMPUTATION_ROH_PAPER/files/Dog_Wolf_aDNA_WG-Modern.tsv')
#roh_ref <- read.csv('~/Downloads/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_1_dogwolf.hom', sep="")
info_ref <- read.delim(snakemake@input[[2]])
roh_ref <- read.csv(snakemake@input[[4]], sep="")
colnames(info_ref)[colnames(info_ref) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] = "Dog_PCA"
colnames(info_ref)[colnames(info_ref) == "Wolf.Dog_PCA"] = "Wolf_Dog_PCA"

#make new column without number in IDs
roh_ref$Sample <- sub("_[^_]+$", "", roh_ref$FID)
colnames(info_ref)[1] <- "Sample"

#fix sample names which lost part of name in previous step:
roh_ref$Sample[roh_ref$Sample == 'Bern'] <- 'Bern_AlpineDachsbracke'
roh_ref$Sample[roh_ref$Sample == 'CatahoulaLeopardDog01_Reseq'] <- 'CatahoulaLeopardDog01'
roh_ref$Sample[roh_ref$Sample == 'MIX'] <- 'MIX_Dachshund01'
roh_ref$Sample[roh_ref$Sample == 'MIX_AmericanCocker'] <- 'MIX_AmericanCocker_Beagle01'
roh_ref$Sample[roh_ref$Sample == 'MIX_KerryBlueTerrier'] <- 'MIX_KerryBlueTerrier_Beagle01'
roh_ref$Sample[roh_ref$Sample == 'MIX_MiniatureSchnauzer'] <- 'MIX_MiniatureSchnauzer_Beagle01'
roh_ref$Sample[roh_ref$Sample == 'VillDog'] <- 'VillDog_Australia01'
roh_ref$Sample[roh_ref$Sample == 'Wolf_WO001'] <- 'Wolf_WO001_895'
roh_ref$Sample[roh_ref$Sample == 'Wolf_WO002'] <- 'Wolf_WO002_732'
roh_ref$Sample[roh_ref$Sample == 'Wolf_WO003'] <- 'Wolf_WO003_636'

#fix Wolf08 metapopulation:
info_ref$Meta.Population[info_ref$Sample == 'Wolf08'] <- 'Western_Eurasian_Wolves'


#remove duplicate samples from metadata (which were reseq)
info_ref_no_duplicate = info_ref[!duplicated(info_ref$Sample),]

#merge two based on specific columns
final_ref <- left_join(roh_ref, info_ref_no_duplicate %>% dplyr::select(Sample, Wolf_Dog_PCA, Dog_PCA, Meta.Population, Species, Age_Mean_BP), "Sample")
final_ref$Age_Mean_BP <- 0


#import genome sizes per chromosome (for Froh estimation)
#sizes_autosomes <- read.delim('~/Downloads/CanFam31_allchr_size.genome', header=FALSE)
sizes_autosomes <- read.delim(snakemake@input[[5]], header=FALSE)
total_genome_size <- sum(sizes_autosomes$V2)

#add type column
final$type <- 'imputed'
final_ref$type <- 'modern'

#extract only common columns from imputed and modern
common_cols <- intersect(colnames(final), colnames(final_ref))
all <- rbind(
  subset(final, select = common_cols),
  subset(final_ref, select = common_cols)
)

##### ######### #### 
##### ALL ROHS #####
##### ######### #### 

# sum total ROH count per sample:
a <- all %>% 
  group_by(Sample) %>%
  mutate(n=n()) %>%
  distinct(Sample, .keep_all = TRUE)

# sum total ROH length per sample:
b <- all %>%
  group_by(Sample) %>%
  summarise(ROH_tol=sum(KB))

ab <- merge(a, b[c(1,2)], by='Sample')

#estimate Froh for all ROH sizes (have to multiply by 1000 since the ROHs are given in KBs):
ab$froh <- (ab$ROH_tol*1000)/total_genome_size

##########  PLOTTING  ########## 

#############
#Wolves only
#############

my_breaks <- c(5000, 10000, 20000, 50000, 100000)
my_labels <- c('5,000', '10,000', '20,000', '50,000', '100,000')

#filter out wolves and unwanted outgroup populations:
wolves <- ab %>% filter(Dog_PCA=='Wolves' & Wolf_Dog_PCA!='Outgroup' & Wolf_Dog_PCA!='Wolves')

wolves$Age_Mean_BP <- as.numeric(wolves$Age_Mean_BP)

#re-scale x axis and age
wolves <- wolves %>%
  mutate(ROH_tol_2 = ROH_tol / 1e+3,
         Age_Mean_KBP = (Age_Mean_BP / 1000)*(-1),
         Age_Mean_KBP_2 = (Age_Mean_BP / 1000))


wolves$Meta.Population[wolves$Sample == 'Wolf08'] <- 'Western_Eurasian_Wolves'


#make modern and ancient layers for plotting:
df_layer_1 <- wolves[wolves$type=="modern",]
df_layer_2 <- wolves[wolves$type!="modern",]

group_names <- c(
  'Pleistocene_Wolves' = 'Pleistocene',
  'Eastern_Eurasian_Wolves' = 'East Eurasia',
  'Western_Eurasian_Wolves' = 'West Eurasia',
  'North_American_Wolves' = 'North America'
)


#ROH count against size
png(snakemake@output[[1]], width=11, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=ROH_tol_2, y=n)) +
  geom_point(data = df_layer_1, size=4, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Age_Mean_BP), size=6, shape=21, alpha=0.7)+
  scale_fill_viridis_c(trans = 'log', breaks = my_breaks, labels = my_labels, option='F') +
  #geom_smooth(method='lm', se=FALSE, color='gray28', size=0.5, alpha=0.8) +
  #geom_label_repel(data = wolves %>% filter(type=='modern' & c(ROH_tol_2>400 | n>350)), 
  #                 aes(x=ROH_tol_2, y=n, label=Sample),size=3.5, box.padding = 3, max.overlaps = Inf)+
  labs(x = "Total ROH length (Mb)", y='Total # ROH')+
  labs(fill = "Sample age (kya)")+
  facet_wrap(.~factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves')), labeller = as_labeller(group_names))+
  theme_bw()+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



#Froh against time per population
cols <- c('royalblue4','darkturquoise','yellow3')
png(snakemake@output[[2]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span=1, colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.7)+
  scale_fill_manual(values=cols, labels = c("East Eurasian", "Pleistocene", "West Eurasian"), name='Wolf population') + 
  scale_x_continuous(breaks=seq(round(min(wolves$Age_Mean_KBP+1)), 0, 10)) +
  #geom_label_repel(data = wolves %>% filter(type!='modern'), 
  #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3.5, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH])))+
  theme_bw()+
  #theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

#Froh against time per population labelled
cols <- c('royalblue4','darkturquoise','yellow3')
png(snakemake@output[[3]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span=1, colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.7)+
  scale_fill_manual(values=cols, labels = c("East Eurasian", "Pleistocene", "West Eurasian"), name='Wolf population') + 
  scale_x_continuous(breaks=seq(round(min(wolves$Age_Mean_KBP+1)), 0, 10)) +
  geom_label_repel(data = wolves %>% filter(type!='modern'), 
                   aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH])))+
  theme_bw()+
  #theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



########################
#only pleistocene wolves

#filter out samples:
wolves_pleist <- ab %>% filter(Wolf_Dog_PCA=='Pleistocene_Wolves')

wolves_pleist$Age_Mean_BP <- as.numeric(wolves_pleist$Age_Mean_BP)

#re-scale x axis and age
wolves_pleist <- wolves_pleist %>%
  mutate(ROH_tol_2 = ROH_tol / 1e+3,
         Age_Mean_KBP = (Age_Mean_BP / 1000)*(-1),
         Age_Mean_KBP_2 = (Age_Mean_BP / 1000))

#make modern and ancient layers for plotting:
df_layer_1 <- wolves_pleist[wolves_pleist$type=="modern",]
df_layer_2 <- wolves_pleist[wolves_pleist$type!="modern",]

group_names <- c(
  'Pleistocene_Wolves' = 'Pleistocene',
  'Eastern_Eurasian_Wolves' = 'East Eurasia',
  'Western_Eurasian_Wolves' = 'West Eurasia',
  'North_American_Wolves' = 'North America'
)

#Froh against time Pleistocene
cols <- 'darkturquoise'
png(snakemake@output[[4]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves_pleist, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span=1, aes(fill=Meta.Population), colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.8)+
  scale_fill_manual(values=cols) + 
  scale_x_continuous(breaks=seq(round(min(wolves_pleist$Age_Mean_KBP+1)), 0, 10)) +
  #geom_label_repel(data = wolves_pleist %>% filter(type!='modern' & froh>0.0002), 
  #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3.5, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH])))+
  theme_bw()+
  theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

#Froh against time Pleistocene labelled
cols <- 'darkturquoise'
png(snakemake@output[[5]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves_pleist, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span=1, aes(fill=Meta.Population), colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.8)+
  scale_fill_manual(values=cols) + 
  scale_x_continuous(breaks=seq(round(min(wolves_pleist$Age_Mean_KBP+1)), 0, 10)) +
  geom_label_repel(data = wolves_pleist %>% filter(type!='modern' & froh>0.0002), 
                   aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH])))+
  theme_bw()+
  theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()


##### ######### #### #####
##### Long ROHS >=1.6Mb #####
##### ######### #### #####

# sum total long ROH count per sample:

a <- all %>% 
  filter(KB>=1600) %>%
  group_by(Sample) %>%
  mutate(n=n()) %>%
  distinct(Sample, .keep_all = TRUE)


b <- all %>%
  filter(KB>=1600) %>%
  group_by(Sample) %>%
  summarise(ROH_tol=sum(KB))

ab <- merge(a, b[c(1,2)], by='Sample')

#estimate Froh for all ROH sizes:
ab$froh <- (ab$ROH_tol*1000)/total_genome_size


##########  PLOTTING  ########## 

#############
#Wolves only
#############

my_breaks <- c(5000, 10000, 20000, 50000, 100000)
my_labels <- c('5,000', '10,000', '20,000', '50,000', '100,000')

#filter out wolves and unwanted outgroup populations:
wolves <- ab %>% filter(Dog_PCA=='Wolves' & Wolf_Dog_PCA!='Outgroup' & Wolf_Dog_PCA!='Wolves')

wolves$Age_Mean_BP <- as.numeric(wolves$Age_Mean_BP)

#re-scale x axis and age
wolves <- wolves %>%
  mutate(ROH_tol_2 = ROH_tol / 1e+3,
         Age_Mean_KBP = (Age_Mean_BP / 1000)*(-1),
         Age_Mean_KBP_2 = (Age_Mean_BP / 1000))


wolves$Meta.Population[wolves$Sample == 'Wolf08'] <- 'Western_Eurasian_Wolves'


#make modern and ancient layers for plotting:
df_layer_1 <- wolves[wolves$type=="modern",]
df_layer_2 <- wolves[wolves$type!="modern",]

group_names <- c(
  'Pleistocene_Wolves' = 'Pleistocene',
  'Eastern_Eurasian_Wolves' = 'East Eurasia',
  'Western_Eurasian_Wolves' = 'West Eurasia',
  'North_American_Wolves' = 'North America'
)


png(snakemake@output[[6]], width=11, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=ROH_tol_2, y=n)) +
  geom_point(data = df_layer_1, size=4, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Age_Mean_BP), size=6, shape=21, alpha=0.7)+
  scale_fill_viridis_c(trans = 'log', breaks = my_breaks, labels = my_labels, option='F') +
  #geom_smooth(method='lm', se=FALSE, color='gray28', size=0.5, alpha=0.8) +
  #geom_label_repel(data = wolves %>% filter(type=='modern' & c(ROH_tol_2>400 | n>350)), 
  #                 aes(x=ROH_tol_2, y=n, label=Sample),size=3.5, box.padding = 3, max.overlaps = Inf)+
  labs(x = "Total ROH length (Mb) (ROH >= 1.6Mb)", y='Total # ROH (ROH >= 1.6Mb)')+
  labs(fill = "Sample age (kya)")+
  facet_wrap(.~factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves')), labeller = as_labeller(group_names))+
  theme_bw()+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()


#Froh against time per population
cols <- c('royalblue4','darkturquoise','yellow3')
png(snakemake@output[[7]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span=1, colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.7)+
  scale_fill_manual(values=cols, labels = c("East Eurasian", "Pleistocene", "West Eurasian"), name='Wolf population') + 
  scale_x_continuous(breaks=seq(round(min(wolves$Age_Mean_KBP+1)), 0, 10)) +
  #geom_label_repel(data = wolves %>% filter(type!='modern'), 
  #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3.5, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH >= 1.6Mb)',sep="")))+
  theme_bw()+
  #theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

#Froh against time per population labelled
cols <- c('royalblue4','darkturquoise','yellow3')
png(snakemake@output[[8]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span=1, colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.7)+
  scale_fill_manual(values=cols, labels = c("East Eurasian", "Pleistocene", "West Eurasian"), name='Wolf population') + 
  scale_x_continuous(breaks=seq(round(min(wolves$Age_Mean_KBP+1)), 0, 10)) +
  geom_label_repel(data = wolves %>% filter(type!='modern'), 
                   aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH >= 1.6Mb)',sep="")))+
  theme_bw()+
  #theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



########################
#only pleistocene wolves

#filter out samples:
wolves_pleist <- ab %>% filter(Wolf_Dog_PCA=='Pleistocene_Wolves')

wolves_pleist$Age_Mean_BP <- as.numeric(wolves_pleist$Age_Mean_BP)

#re-scale x axis and age
wolves_pleist <- wolves_pleist %>%
  mutate(ROH_tol_2 = ROH_tol / 1e+3,
         Age_Mean_KBP = (Age_Mean_BP / 1000)*(-1),
         Age_Mean_KBP_2 = (Age_Mean_BP / 1000))

#make modern and ancient layers for plotting:
df_layer_1 <- wolves_pleist[wolves_pleist$type=="modern",]
df_layer_2 <- wolves_pleist[wolves_pleist$type!="modern",]

group_names <- c(
  'Pleistocene_Wolves' = 'Pleistocene',
  'Eastern_Eurasian_Wolves' = 'East Eurasia',
  'Western_Eurasian_Wolves' = 'West Eurasia',
  'North_American_Wolves' = 'North America'
)

#Froh against time Pleistocene
cols <- 'darkturquoise'
png(snakemake@output[[9]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves_pleist, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span=1, aes(fill=Meta.Population), colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.8)+
  scale_fill_manual(values=cols) + 
  scale_x_continuous(breaks=seq(round(min(wolves_pleist$Age_Mean_KBP+1)), 0, 10)) +
  #geom_label_repel(data = wolves_pleist %>% filter(type!='modern' & froh>0.0002), 
  #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3.5, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH >= 1.6Mb)',sep="")))+
  theme_bw()+
  theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()

#Froh against time Pleistocene labelled
cols <- 'darkturquoise'
png(snakemake@output[[10]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves_pleist, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span=1, aes(fill=Meta.Population), colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.8)+
  scale_fill_manual(values=cols) + 
  scale_x_continuous(breaks=seq(round(min(wolves_pleist$Age_Mean_KBP+1)), 0, 10)) +
  geom_label_repel(data = wolves_pleist %>% filter(type!='modern' & froh>0.0002), 
                   aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH >= 1.6Mb)',sep="")))+
  theme_bw()+
  theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()



##### ######### #### #####
##### Short ROHS <1.6Mb #####
##### ######### #### #####

# sum total long ROH count per sample:

a <- all %>% 
  filter(KB<1600) %>%
  group_by(Sample) %>%
  mutate(n=n()) %>%
  distinct(Sample, .keep_all = TRUE)


b <- all %>%
  filter(KB<1600) %>%
  group_by(Sample) %>%
  summarise(ROH_tol=sum(KB))

ab <- merge(a, b[c(1,2)], by='Sample')

#estimate Froh for all ROH sizes:
ab$froh <- (ab$ROH_tol*1000)/total_genome_size


##########  PLOTTING  ########## 

#############
#Wolves only
#############

my_breaks <- c(5000, 10000, 20000, 50000, 100000)
my_labels <- c('5,000', '10,000', '20,000', '50,000', '100,000')

#filter out wolves and unwanted outgroup populations:
wolves <- ab %>% filter(Dog_PCA=='Wolves' & Wolf_Dog_PCA!='Outgroup' & Wolf_Dog_PCA!='Wolves')

wolves$Age_Mean_BP <- as.numeric(wolves$Age_Mean_BP)

#re-scale x axis and age
wolves <- wolves %>%
  mutate(ROH_tol_2 = ROH_tol / 1e+3,
         Age_Mean_KBP = (Age_Mean_BP / 1000)*(-1),
         Age_Mean_KBP_2 = (Age_Mean_BP / 1000))


wolves$Meta.Population[wolves$Sample == 'Wolf08'] <- 'Western_Eurasian_Wolves'


#make modern and ancient layers for plotting:
df_layer_1 <- wolves[wolves$type=="modern",]
df_layer_2 <- wolves[wolves$type!="modern",]

group_names <- c(
  'Pleistocene_Wolves' = 'Pleistocene',
  'Eastern_Eurasian_Wolves' = 'East Eurasia',
  'Western_Eurasian_Wolves' = 'West Eurasia',
  'North_American_Wolves' = 'North America'
)


png(snakemake@output[[11]], width=11, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=ROH_tol_2, y=n)) +
  geom_point(data = df_layer_1, size=4, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Age_Mean_BP), size=6, shape=21, alpha=0.7)+
  scale_fill_viridis_c(trans = 'log', breaks = my_breaks, labels = my_labels, option='F') +
  #geom_smooth(method='lm', se=FALSE, color='gray28', size=0.5, alpha=0.8) +
  #geom_label_repel(data = wolves %>% filter(type=='modern' & c(ROH_tol_2>400 | n>350)), 
  #                 aes(x=ROH_tol_2, y=n, label=Sample),size=3.5, box.padding = 3, max.overlaps = Inf)+
  labs(x = "Total ROH length (Mb) (ROH < 1.6Mb)", y='Total # ROH (ROH < 1.6Mb)')+
  labs(fill = "Sample age (kya)")+
  facet_wrap(.~factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves')), labeller = as_labeller(group_names))+
  theme_bw()+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        legend.key.size = unit(0.9, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()


#Froh against time per population
cols <- c('royalblue4','darkturquoise','yellow3')
png(snakemake@output[[12]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span=1, colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.7)+
  scale_fill_manual(values=cols, labels = c("East Eurasian", "Pleistocene", "West Eurasian"), name='Wolf population') + 
  scale_x_continuous(breaks=seq(round(min(wolves$Age_Mean_KBP+1)), 0, 10)) +
  #geom_label_repel(data = wolves %>% filter(type!='modern'), 
  #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3.5, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH < 1.6Mb)',sep="")))+
  theme_bw()+
  #theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
#dev.off()

#Froh against time per population labelled
cols <- c('royalblue4','darkturquoise','yellow3')
png(snakemake@output[[13]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE, span=1, colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.7)+
  scale_fill_manual(values=cols, labels = c("East Eurasian", "Pleistocene", "West Eurasian"), name='Wolf population') + 
  scale_x_continuous(breaks=seq(round(min(wolves$Age_Mean_KBP+1)), 0, 10)) +
  geom_label_repel(data = wolves %>% filter(type!='modern'), 
                   aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH < 1.6Mb)',sep="")))+
  theme_bw()+
  #theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
#dev.off()



########################
#only pleistocene wolves

#filter out samples:
wolves_pleist <- ab %>% filter(Wolf_Dog_PCA=='Pleistocene_Wolves')

wolves_pleist$Age_Mean_BP <- as.numeric(wolves_pleist$Age_Mean_BP)

#re-scale x axis and age
wolves_pleist <- wolves_pleist %>%
  mutate(ROH_tol_2 = ROH_tol / 1e+3,
         Age_Mean_KBP = (Age_Mean_BP / 1000)*(-1),
         Age_Mean_KBP_2 = (Age_Mean_BP / 1000))

#make modern and ancient layers for plotting:
df_layer_1 <- wolves_pleist[wolves_pleist$type=="modern",]
df_layer_2 <- wolves_pleist[wolves_pleist$type!="modern",]

group_names <- c(
  'Pleistocene_Wolves' = 'Pleistocene',
  'Eastern_Eurasian_Wolves' = 'East Eurasia',
  'Western_Eurasian_Wolves' = 'West Eurasia',
  'North_American_Wolves' = 'North America'
)

#Froh against time Pleistocene
cols <- 'darkturquoise'
png(snakemake@output[[14]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves_pleist, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span=1, aes(fill=Meta.Population), colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.8)+
  scale_fill_manual(values=cols) + 
  scale_x_continuous(breaks=seq(round(min(wolves_pleist$Age_Mean_KBP+1)), 0, 10)) +
  #geom_label_repel(data = wolves_pleist %>% filter(type!='modern' & froh>0.0002), 
  #                 aes(x=Age_Mean_KBP, y=froh, label=Sample),size=3.5, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH < 1.6Mb)',sep="")))+
  theme_bw()+
  theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
#dev.off()

#Froh against time Pleistocene labelled
cols <- 'darkturquoise'
png(snakemake@output[[15]], width=12, height=8, units='in', res=200, pointsize=4)
par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i", cex.axis = 2, cex.lab  = 2)
options(scipen = 999)
ggplot(wolves_pleist, aes(x=Age_Mean_KBP, y=froh)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, span=1, aes(fill=Meta.Population), colour="black", size=0.5) +
  geom_point(data = df_layer_1, size=2, alpha=0.6, colour='grey')+
  geom_point(data = df_layer_2, aes(fill=Meta.Population), size=4, shape=21, alpha=0.8)+
  scale_fill_manual(values=cols) + 
  scale_x_continuous(breaks=seq(round(min(wolves_pleist$Age_Mean_KBP+1)), 0, 10)) +
  geom_label_repel(data = wolves_pleist %>% filter(type!='modern' & froh>0.0002), 
                   aes(x=Age_Mean_KBP, y=froh, label=Sample),size=2, box.padding = 1, max.overlaps = Inf)+
  labs(x = "Time (kya)", y=expression(paste(italic('F')[ROH],' (ROH < 1.6Mb)',sep="")))+
  theme_bw()+
  theme(legend.position = "none")+
  #facet_grid(factor(Meta.Population, levels=c('Pleistocene_Wolves','Eastern_Eurasian_Wolves','Western_Eurasian_Wolves', 'North_American_Wolves'))~. , labeller = as_labeller(group_names))+
  theme(strip.background =element_rect(fill="gray28"),
        strip.text = element_text(colour = 'white', size=18),
        axis.text.y=element_text(size=16),
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=18),
        axis.title.x=element_text(size=18),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18),
        #legend.position = c(0.15, 0.9),
        #legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
#dev.off()


