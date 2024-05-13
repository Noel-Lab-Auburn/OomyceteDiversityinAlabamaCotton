################Data Preprocessing############

###### Libraries #####
library(phyloseq)
library(vegan)
library(tidyverse)
library(metagenomeSeq)
library(ggpubr)
library(Biostrings)
library(microbiome)
library(ggrepel)

##### Set global options #####

# no scientific notation
options(scipen=10000) 

# color blind pallets used throughout 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ibm.cbb <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000")
tol.cbb <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

# Read in the RDS file
oomy_raw <- readRDS(file = "AmpliconSeq/oomy_ALABAMA_KEMI_2024-02-19.rds")

# dropping two samples that we couldn't figure out which field they came from 
oomy_raw_filt <- oomy_raw %>%
  subset_samples(!Sample %in% c("KAL385", "KAL389")) %>%
  subset_taxa(Phylum == "Oomycota")

### Rarefaction ### 

sam.data <- data.frame(oomy_raw_filt@sam_data)
oOTU.table <- otu_table(oomy_raw_filt) %>%
  as.data.frame() %>%
  as.matrix()

raremax <- min(rowSums(t(oOTU.table)))
rare.fun <- rarecurve(t(oOTU.table), step = 500, sample = raremax, tidy = T)

oom.rare.curve.extract2 <- left_join(sam.data, rare.fun, by = c("Sample" = "Site"))

oom.rare <- ggplot(oom.rare.curve.extract2, aes(x = Sample.y.y, y = Species, group = Sample)) + 
  #geom_point() +
  geom_line() + 
  xlab("Reads") + 
  ylab("Number of OTUs") +
  theme_classic() + 
  #geom_vline(xintercept = median(sample_sums(oomy_raw_filt)), linetype = "dashed") +
  #geom_vline(xintercept = min(sample_sums(oomy_raw_filt)), linetype = "dashed", color = "red") +
  scale_color_manual(values = cbbPalette)

## Things really start leveling out around 10,000 reads per sample.

####### ALPHA DIVERSITY #########
oomy_raw_filt@sam_data$richness <- estimate_richness(oomy_raw_filt, measures=c("Observed"))$Observed

mean(oomy_raw_filt@sam_data$richness)


oomy_raw_filt %>%
  subset_taxa(Genus %in% c("Pythium",
                           "Globisporangium", 
                           "Phytophthora", 
                           "Phytopythium")) %>% 
  estimate_richness(measures=c("Observed")) %>%
  summary()


# separating the otu table, metadata, and the tax data for ease 
meta.data <- data.frame(oomy_raw_filt@sam_data)
tax.data <- data.frame(oomy_raw_filt@tax_table)
otu <- oomy_raw_filt@otu_table %>%
  as("matrix")

# summary stats
sum(sample_sums(oomy_raw_filt))
mean(sample_sums(oomy_raw_filt))
median(sample_sums(oomy_raw_filt))
#335,385 sequences
#16,000 sequences per sample
#116 OTUs within phylum 

# distribution of species counts
count_otus <- tax.data %>%
  group_by(Species) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  arrange(desc(freq)) %>%
  print(n = 38)
# 6 OTUs were irregulare (5.2%)
# 3 OTUs were ultimum (2.6%) 
# 2 OTUs were heterothallicum (1.7%)
# 1 OTU was nicotianae (0.9%)
# 1 OTU was rostratifigens (0.9%)
# 1 OTU was spinosum (0.9%)

# 3 OTUs were genus Lagenidium (2.6%)

#no Py. deliense OTUs
#no G. sylvaticum OTUs
#no G. acanthophoron OTUs
#no Phyto. cucubitacearum OTUs
#no Phyto. helicoides OTUs

#calculating abundance occupancy
otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- rownames_to_column(as.data.frame(cbind(otu_occ, otu_rel)),'OTU')

occ_abund_tax <- left_join(occ_abun, tax.data, by = "OTU")

occ_abund_tax$Other <- ifelse(occ_abund_tax$otu_occ < 0.31 & occ_abund_tax$otu_rel < 0.05, "Other", occ_abund_tax$Genus)

# abundance occupancy by pathogen
# G. irregulare
irregulare <- occ_abund_tax %>%
  subset(Species == "irregulare")

# G. ultimum
ultimum <- occ_abund_tax %>%
  subset(Species == "ultimum")

# G. heterothallicum
heterothallicum <- occ_abund_tax %>%
  subset(Species == "heterothallicum")

# P. nicotianae
nicotianae <- occ_abund_tax %>%
  subset(Species == "nicotianae")

# G. rostratifingens
rostratifingens <- occ_abund_tax %>%
  subset(Species == "rostratifingens")

# G. spinosum
spinosum <- occ_abund_tax %>%
  subset(Species == "spinosum")

# genus Lagenidium
lagenidium <- occ_abund_tax %>%
  subset(Genus == "Lagenidium")

#### Fig 1a - Abundance occupancy #####
fig1a <- ggplot(occ_abund_tax, aes(y = otu_occ, x = otu_rel, color = Other)) +
  geom_point() + 
  theme_classic2() + 
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) + 
  xlab("% Relative Abundance") +
  ylab("% Occupancy") + 
  scale_color_manual(values = c(cbbPalette, ibm.cbb)) +
  geom_text_repel(data = occ_abund_tax[occ_abund_tax$otu_occ > 0.30 | occ_abund_tax$otu_rel > 0.05,], 
                  aes(y = otu_occ, x = otu_rel, color = Other, label = Label), size = 2) +
  theme(legend.text = element_text(face = "italic", size = 8),
legend.title = element_blank())

#### Fig 1b - Relative abundances ####

# creating a list of taxa that I want to show up in the plot as Other
taxa.other.list <- oomy_raw_filt %>%
  subset_taxa(Genus %in% c("Pythium",
                           "Globisporangium", 
                           "Phytophthora", 
                           "Phytopythium")) %>%
  microbiome::transform("compositional") %>%
  psmelt() %>%
  group_by(Lowest_Taxnomic_Rank) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  mutate(Other = ifelse(MeanRelAbund < 0.01, "Other", Lowest_Taxnomic_Rank)) %>%
  filter(Other == "Other") %>%
  pull(Lowest_Taxnomic_Rank)

# relative abundance of Pythium, Globisporangium, Phytophthora, and Phytopythium species
fig1b <- oomy_raw_filt %>%
  subset_taxa(Genus %in% c("Pythium",
                           "Globisporangium", 
                           "Phytophthora", 
                           "Phytopythium")) %>%
  microbiome::transform("compositional") %>%
  psmelt() %>%
  group_by(Sample.y, Location, Label) %>%
  summarise(MeanRelAbund = mean(Abundance)) %>%
  left_join(tax.data, by = "Label") %>%
  mutate(Other = ifelse(Lowest_Taxnomic_Rank %in% taxa.other.list, "Other < 0.01%", Lowest_Taxnomic_Rank)) %>%
  ggplot(aes(x = Sample.y, y = MeanRelAbund, fill = Other)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values= c(cbbPalette, ibm.cbb, tol.cbb)) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Location, scales = 'free') +
  labs(x = "", y = "Relative abundance (%)",
       title = "") + 
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.text = element_text(face = "italic", size = 5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        legend.position = "bottom") 


### Presence absence of otus
PA_table <- data.frame(t(otu_PA)) %>%
  rownames_to_column() %>%
  pivot_longer(cols = OOTU_103:OOTU_937, values_to = "PA", names_to = "OTU") %>%
  rename("rowname" = "Sample") %>%
  left_join(meta.data, by = "Sample") %>%
  left_join(tax.data, by = "OTU") %>%
  subset(Genus %in% c("Pythium",
                      "Globisporangium", 
                      "Phytophthora", 
                      "Phytopythium"))


fig1c <- ggplot(PA_table, aes(y = reorder(Label, PA), x = reorder(Sample.y.y, Latitude), fill = as.factor(PA))) +
  geom_tile(color = "black") + 
  theme_classic() + 
  scale_fill_manual(values = c("white",cbbPalette[[4]]), name = "", labels = c("Absent", "Present")) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 5),
        axis.text.y = element_text(face = "italic", size = 5)) +
  xlab("") + 
  ylab("") 


## Final figure 
fig1ab <- ggarrange(fig1a, fig1b, labels = "auto", nrow = 2)

ggarrange(fig1ab, fig1c, labels = c("", "c"), ncol = 2)
  
