library(ggplot2)
library(dplyr)
library(ggplotify)


#### ROH bands dogs ###

# imputed dogs
info <- read.delim(snakemake@input[[1]])
roh <- read.csv(snakemake@input[[3]], sep = "")
# info <- read.delim('~/Desktop/Copenhagen_PhD/files/IMPUTATION_2023/IMPUTATION_ROH_PAPER/files/Dog_Wolf_aDNA_WG-Master.tsv')
# roh <- read.csv('~/Downloads/merged_phased.allchrom_MAF_0.01_INFO_0.8_all_sites_hom_win_het_1_dogwolf.hom', sep="")

colnames(info)[colnames(info) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] <- "Dog_PCA"
colnames(info)[colnames(info) == "Wolf.Dog_PCA"] <- "Wolf_Dog_PCA"
colnames(info)[1] <- "Sample"
colnames(roh)[1] <- "Sample"

# fix Tumat (Tumat2) and WolfHead (Wolf_head_IN18-016)
roh$Sample[roh$Sample == "Tumat"] <- "Tumat2"
roh$Sample[roh$Sample == "WolfHead"] <- "Wolf_head_IN18-016"

sites <- as.character(snakemake@params[["sites"]])

# merge two based on specific columns
final <- left_join(roh, info %>% dplyr::select(Sample, Wolf_Dog_PCA, Dog_PCA, Meta.Population, Species, Age_Mean_BP), "Sample")
final$name_age <- paste(final$Sample, final$Age_Mean_BP, "bp", sep = "_")


##### reference panel selection of dogs:
info_ref <- read.delim(snakemake@input[[2]])
roh_ref <- read.csv(snakemake@input[[4]], sep = "")
# info_ref <- read.delim('~/Desktop/Copenhagen_PhD/files/IMPUTATION_2023/IMPUTATION_ROH_PAPER/files/Dog_Wolf_aDNA_WG-Modern.tsv')
# roh_ref <- read.csv('~/Downloads/ref-panel_allchrom_sample-snp_filltags_filter_all_sites_hom_win_het_1_dogwolf.hom', sep="")

colnames(info_ref)[colnames(info_ref) == "Dog_PCA..European..Arctic.NA..East.Asia..Near.Eastern.Africa."] <- "Dog_PCA"
colnames(info_ref)[colnames(info_ref) == "Wolf.Dog_PCA"] <- "Wolf_Dog_PCA"

# make new column without number in IDs
roh_ref$Sample <- sub("_[^_]+$", "", roh_ref$FID)
colnames(info_ref)[1] <- "Sample"


# fix sample names which lost part of name in previous step:
roh_ref$Sample[roh_ref$Sample == "Bern"] <- "Bern_AlpineDachsbracke"
roh_ref$Sample[roh_ref$Sample == "CatahoulaLeopardDog01_Reseq"] <- "CatahoulaLeopardDog01"
roh_ref$Sample[roh_ref$Sample == "MIX"] <- "MIX_Dachshund01"
roh_ref$Sample[roh_ref$Sample == "MIX_AmericanCocker"] <- "MIX_AmericanCocker_Beagle01"
roh_ref$Sample[roh_ref$Sample == "MIX_KerryBlueTerrier"] <- "MIX_KerryBlueTerrier_Beagle01"
roh_ref$Sample[roh_ref$Sample == "MIX_MiniatureSchnauzer"] <- "MIX_MiniatureSchnauzer_Beagle01"
roh_ref$Sample[roh_ref$Sample == "VillDog"] <- "VillDog_Australia01"
roh_ref$Sample[roh_ref$Sample == "Wolf_WO001"] <- "Wolf_WO001_895"
roh_ref$Sample[roh_ref$Sample == "Wolf_WO002"] <- "Wolf_WO002_732"
roh_ref$Sample[roh_ref$Sample == "Wolf_WO003"] <- "Wolf_WO003_636"

# fix Wolf08 metapopulation:
info_ref$Meta.Population[info_ref$Sample == "Wolf08"] <- "Western_Eurasian_Wolves"


# remove duplicate samples from metadata (which were reseq)
info_ref_no_duplicate <- info_ref[!duplicated(info_ref$Sample), ]

# merge two based on specific columns
final_ref <- left_join(roh_ref, info_ref_no_duplicate %>% dplyr::select(Sample, Wolf_Dog_PCA, Dog_PCA, Meta.Population, Species, Age_Mean_BP), "Sample")
final_ref$Age_Mean_BP <- 0

dog_list <- c(
    "AfghanHound01", "Basenji02", "IndigenousDogNigeria02", "Saluki01_1233", "AlaskanHusky01",
    "AlaskanMalamute02", "Beagle02", "BelgianMalinois01", "BorderCollie01", "BullTerrier01",
    "CockerSpanielAmerican01_11414", "GreenlandDog01", "GermanShepherd01", "GoldenRetriever01", "LabradorRetriever01",
    "Pomeranian01", "Samoyed01_176", "YorkshireTerrier01"
)

sub_ref <- subset(final_ref, FID %in% dog_list)
sub_ref$name_age <- paste(sub_ref$Sample, sub_ref$Age_Mean_BP, "bp", sep = "_")


# merge imputed and modern based on common columns
common_cols <- intersect(colnames(final), colnames(sub_ref))
merged <- rbind(
    final[, common_cols],
    sub_ref[, common_cols]
)

# only selected dogs
merged_dogs <- merged %>% filter(Dog_PCA != "Wolves" & Dog_PCA != "Americas_Dogs" & Dog_PCA != "East_Asian_Dogs")


# make port au choix arctic:
merged_dogs$Dog_PCA[merged_dogs$Dog_PCA == "preContact_Dogs"] <- "Arctic_Dogs"
merged_dogs$Meta.Population[merged_dogs$Meta.Population == "preContact_Dogs"] <- "Arctic_Dogs"
merged_dogs$Meta.Population[merged_dogs$Meta.Population == "American_European_Dogs"] <- "European_Dogs"

merged_dogs$Age_Mean_BP <- as.numeric(merged_dogs$Age_Mean_BP)

# ~~ ROH for some indiividuals --------------------------------------------------
# Figure 1A, ROH for most and least inbred individuals

all_roh <- merged_dogs %>%
    group_by(IID) %>%
    summarise(sum_roh = sum(KB)) %>%
    ungroup() %>%
    arrange(desc(sum_roh))


num_ind <- length(unique(all_roh$IID))


df_1 <- merged_dogs %>%
    mutate(
        POS1 = POS1 / 1e+6,
        POS2 = POS2 / 1e+6,
        MB = KB / 1000
    )


df_1$IID <- as.factor(df_1$IID)

yax_1 <- df_1 %>%
    arrange(Meta.Population, Age_Mean_BP)

rownames(yax_1) <- NULL

yax <- data.frame(name_age = unique(yax_1$name_age)) %>%
    mutate(yax = seq(
        from = 2,
        to = 2 * length(unique(yax_1$name_age)),
        by = 2
    ))

df <- left_join(df_1, yax, by = "name_age")

shade <- df %>%
    group_by(CHR) %>%
    summarise(min = min(POS1), max = max(POS2)) %>%
    mutate(min = case_when(
        CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
            CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
            CHR == 22 | CHR == 24 | CHR == 26 | CHR == 28 | CHR == 30 |
            CHR == 32 | CHR == 34 | CHR == 36 | CHR == 38 ~ 0,
        TRUE ~ min
    )) %>%
    mutate(max = case_when(
        CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
            CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
            CHR == 22 | CHR == 24 | CHR == 26 | CHR == 28 | CHR == 30 |
            CHR == 32 | CHR == 34 | CHR == 36 | CHR == 38 ~ 0,
        TRUE ~ max
    ))

chr_names <- as.character(1:38)
names(chr_names) <- as.character(1:38)

d <- yax$name_age
c <- yax$yax

cols <- c("#ce4441", "#62929a", "#ffbb44")
col_x <- ifelse(grepl("_0_bp", d), "palevioletred3", "black")

title_final <- paste("Dog ROHs - ", sites, sep = "")


png(snakemake@output[[1]], width = 22, height = 12, units = "in", res = 200, pointsize = 4)
df %>%
    ggplot() +
    geom_rect(
        data = shade, aes(xmin = min, xmax = max, ymin = 0, ymax = num_ind * 2 + 1),
        alpha = 0.5, fill = "gray88"
    ) +
    geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", linewidth = 0.4) +
    geom_rect(aes(
        xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9,
        fill = Meta.Population
    ), size = 0, alpha = 1) +
    scale_fill_manual(values = cols, labels = c("Africa, Near East, India", "Arctic", "Europe")) +
    scale_y_continuous(breaks = c, labels = d) +
    facet_grid(~CHR,
        scales = "free_x", space = "free_x", switch = "x",
        labeller = as_labeller(chr_names)
    ) +
    theme_minimal(base_family = "Helvetica", base_size = 12) +
    theme(
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
        strip.text.x = element_text(size = 14, angle = 30),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11, colour = col_x),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 0), size = 18),
        axis.title.y = element_text(margin = margin(r = 0), size = 18),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(clip = "off") +
    labs(fill = "Population:") +
    xlab("Chromosome") +
    ylab("Individuals") +
    ggtitle(title_final)
dev.off()




### plot wolves

sub_ref_wolf <- final_ref %>% filter(Dog_PCA == "Wolves" & Wolf_Dog_PCA != "Outgroup" & Wolf_Dog_PCA != "Wolves")

# will use these wolves:

# West Eurasia:
# Spain: Wolf27
# Iran: Wolf20
# Portugal: Wolf24
# Iberian: Wolf39

# East Eurasia:
# China: Wolf05,
# India: Wolf19
# Tibet: WolfTibetan02,

# North America
# Mexico: Wolf22
# Great plains: GreatPlainsWolf01
# Isle Royale: Wolf40
# Yellowstone: Wolf28


wolf_list <- c(
    "Wolf27", "Wolf20", "Wolf24", "Wolf39", "Wolf05",
    "Wolf19", "WolfTibetan02", "Wolf22", "GreatPlainsWolf01", "Wolf40", "Wolf28"
)

sub_ref <- subset(final_ref, FID %in% wolf_list)
sub_ref$name_age <- paste(sub_ref$Sample, sub_ref$Age_Mean_BP, "bp", sep = "_")


# merge imputed and modern based on common columns
common_cols <- intersect(colnames(final), colnames(sub_ref))
merged <- rbind(
    final[, common_cols],
    sub_ref[, common_cols]
)

# only selected wolves
merged_wolves <- merged %>% filter(Dog_PCA == "Wolves" & Wolf_Dog_PCA != "Outgroup" & Wolf_Dog_PCA != "Wolves")


all_roh <- merged_wolves %>%
    group_by(IID) %>%
    summarise(sum_roh = sum(KB)) %>%
    ungroup() %>%
    arrange(desc(sum_roh))


num_ind <- length(unique(all_roh$IID))

df_1 <- merged_wolves %>%
    mutate(
        POS1 = POS1 / 1e+6,
        POS2 = POS2 / 1e+6,
        MB = KB / 1000
    )

df_1$IID <- as.factor(df_1$IID)


# this is the put the Pleistocene wolves on the top:
yax_1 <- df_1 %>%
    arrange(Meta.Population = factor(Meta.Population, levels = c(
        "Eastern_Eurasian_Wolves",
        "Western_Eurasian_Wolves",
        "North_American_Wolves",
        "Pleistocene_Wolves"
    )), Age_Mean_BP)


rownames(yax_1) <- NULL

yax <- data.frame(name_age = unique(yax_1$name_age)) %>%
    mutate(yax = seq(
        from = 2,
        to = 2 * length(unique(yax_1$name_age)),
        by = 2
    ))

df <- left_join(df_1, yax, by = "name_age")



shade <- df %>%
    group_by(CHR) %>%
    summarise(min = min(POS1), max = max(POS2)) %>%
    mutate(min = case_when(
        CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
            CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
            CHR == 22 | CHR == 24 | CHR == 26 | CHR == 28 | CHR == 30 |
            CHR == 32 | CHR == 34 | CHR == 36 | CHR == 38 ~ 0,
        TRUE ~ min
    )) %>%
    mutate(max = case_when(
        CHR == 2 | CHR == 4 | CHR == 6 | CHR == 8 | CHR == 10 |
            CHR == 12 | CHR == 14 | CHR == 16 | CHR == 18 | CHR == 20 |
            CHR == 22 | CHR == 24 | CHR == 26 | CHR == 28 | CHR == 30 |
            CHR == 32 | CHR == 34 | CHR == 36 | CHR == 38 ~ 0,
        TRUE ~ max
    ))


chr_names <- as.character(1:38)
names(chr_names) <- as.character(1:38)

d <- yax$name_age
c <- yax$yax

col_x <- ifelse(grepl("_0_bp", d), "palevioletred3", "black")
cols <- c("royalblue4", "darksalmon", "darkturquoise", "yellow3")

title_final <- paste("Wolf ROHs - ", sites, sep = "")

png(snakemake@output[[2]], width = 22, height = 12, units = "in", res = 200, pointsize = 4)
df %>%
    ggplot() +
    geom_rect(
        data = shade, aes(xmin = min, xmax = max, ymin = 0, ymax = num_ind * 2 + 1),
        alpha = 0.5, fill = "gray88"
    ) +
    geom_hline(data = yax, aes(yintercept = yax), color = "#d8dee9", size = 0.4) +
    geom_rect(aes(
        xmin = POS1, xmax = POS2, ymin = yax - 0.5, ymax = yax + 0.9,
        fill = Meta.Population
    ), size = 0, alpha = 1) +
    scale_fill_manual(
        values = cols, breaks = c("Eastern_Eurasian_Wolves", "Western_Eurasian_Wolves", "North_American_Wolves", "Pleistocene_Wolves"),
        labels = c("East Eurasia", "West Eurasia", "North America", "Pleistocene")
    ) +
    scale_y_continuous(breaks = c, labels = d) +
    facet_grid(~CHR,
        scales = "free_x", space = "free_x", switch = "x",
        labeller = as_labeller(chr_names)
    ) +
    theme_minimal(base_family = "Helvetica", base_size = 12) + # start with blank canvas
    theme(
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(r = 0.5, l = 0.1, b = 0.1, t = 0.1, unit = "cm"),
        strip.text.x = element_text(size = 14, angle = 30),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 11, colour = col_x),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 0), size = 18),
        axis.title.y = element_text(margin = margin(r = 0), size = 18),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5)
    ) +
    coord_cartesian(clip = "off") +
    labs(fill = "Population:") +
    xlab("Chromosome") +
    ylab("Individuals") +
    ggtitle(title_final)
dev.off()
