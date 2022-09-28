if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")



library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)


# color blind friendly pallete 
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

fungi.colors <- c("#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c",
                  "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32",
                  "#fdae6b", "#f16913", "#d94801", "#8c2d04",
                  "#dadaeb", "#bcbddc", "#807dba", "#6a51a3", "#4a1486",
                  "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d",
                  "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525")

library(microshades)
hex_values <-c(microshades_palette("micro_orange",3, lightest = FALSE), 
               microshades_palette("micro_blue",3, lightest = FALSE), 
               microshades_palette("micro_purple",3, lightest = FALSE), microshades_palette("micro_gray", 3, lightest = FALSE), microshades_palette("micro_brown",3, lightest = FALSE) ,microshades_palette("micro_green",2, lightest = FALSE)) 





samp_dat <- read.csv("C:/Users/Kemi Olofintila/Documents/Kemi/My Research/Oomycete Research/Oomycete_Combined/Metadata_Combined.csv", na.strings = "na")
rownames(samp_dat) <- samp_dat$Sample
samp_dat <- samp_dat[,-1]
SAMP <- phyloseq::sample_data(samp_dat)



otu <- read.csv("~/Kemi/My Research/Oomycete Research/Oomycete_Combined/OTU Table_Combined.csv", na.strings = "na")
rownames(otu) <- otu$Species
otu <- otu[,-1]
OTU <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)

tax <- read.csv("C:/Users/Kemi Olofintila/Documents/Kemi/My Research/Oomycete Research/Oomycete_Combined/Taxonomy_Combined.csv")
rownames(tax) <- tax$Species
tax <- tax[,-1]
TAX <- phyloseq::tax_table(as.matrix(tax))



#FASTA <- readDNAStringSet("path to fasta file", format="fasta", seek.first.rec=TRUE, use.names=TRUE)


#Combines all data
oomycetes <- phyloseq::phyloseq(OTU, 
                                TAX, 
                                #FASTA, 
                                SAMP)
#Subset Data for Each year
Oomycete.1 = subset_samples(oomycetes, Year== "2021")
Oomycete1 = prune_taxa(taxa_sums(Oomycete.1) > 0, Oomycete.1) #prune the rest of data to include only 2021
#oomycetes@sam_data$ to pass through phyloseq object

Oomycete.2 = subset_samples(oomycetes, Year== "2022")
Oomycete2 = prune_taxa(taxa_sums(Oomycete.2) > 0, Oomycete.2)
#alpha diversity example 
plot_richness(Oomycete2) # canned version
ggsave("OOMYCETE_RICHNESS.png", dpi = 300)

plot_richness(Oomycete2, measures= "Observed") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Alpha_Diversity.png", dpi = 300)

# manual version, much more control
Oomycete2@sam_data$shannon <- estimate_richness(Oomycete2, measures=c("Shannon"))$Shannon #shannon diversity index
Oomycete2@sam_data$richness <- estimate_richness(Oomycete2, measures=c("Observed"))$Observed #observed number of taxa
Oomycete2@sam_data$even <- Oomycete2@sam_data$shannon/log(Oomycete2@sam_data$richness) #plieou's evenness

#includes above info in a new datafile

oomycete.metadata <- Oomycete2@sam_data


ggplot(oomycete.metadata, aes(x = reorder(Location, richness), y = richness)) + 
  stat_summary(fun.y=mean,geom="bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  geom_jitter(size = 2) +
  #stat_compare_means(aes(group = Location, label = ..p.signif..), label.y = 10) + 
  #stat_compare_means(label = "p.signif", method = "t.test") +
  stat_compare_means(method = "anova")+
  ylab("Richness") + 
  xlab("Location") +
  #ggtitle("Fungi") + 
  xlab("")+
  #scale_color_manual(values=cbbPalette) +
  #scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2))) +
  #stat_compare_means(method = "t.test") +
  theme_classic()
ggsave("Richness_Per_Location.png", dpi = 300)

ggplot(oomycete.metadata, aes(x = reorder(Location, even), y = even)) + 
  stat_summary(fun.y=mean,geom="bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  geom_jitter(size = 2) +
  #stat_compare_means(aes(group = Location, label = ..p.signif..), label.y = 10) + 
  #stat_compare_means(label = "p.signif", method = "t.test") +
  stat_compare_means(method = "anova")+
  ylab("Evenness") + 
  xlab("Location") +
  #ggtitle("Fungi") + 
  xlab("")+
  #scale_color_manual(values=cbbPalette) +
  #scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2))) +
  #stat_compare_means(method = "t.test") +
  theme_classic()
ggsave("Evenness_Per_Location.png", dpi = 300)

# beta diversity example with canned phyloseq ordinations. I can teach you how to do it mannually too so you can have more control over the plotting. 
GP.ord <- ordinate(Oomycete1, "MDS", "bray")
p1 = plot_ordination(Oomycete1, GP.ord, type="samples", color = "Location")
print(p1)


global.nmds.data <- p1$data
#Permanova- Permutational Multivariate ANOVA

oomycete.dist.bray = phyloseq::distance(Oomycete1, "bray") #creates distance matrix (0 (similar) 1 (non-similar))
adonis2(oomycete.dist.bray~pH, as(sample_data(Oomycete1), "data.frame"), permutations = 9999) 

ggplot() + 
  geom_point(data = global.nmds.data, aes(x = Axis.1, y = Axis.2, shape = Location, color = Sand), alpha = 0.8, size = 2) +
  theme_bw() +
  ylab("PCoA2") + 
  xlab("PCoA1") +
  #scale_fill_manual(values=cbbPalette) +
  #stat_ellipse(data = global.nmds.data, aes(x = Axis.1, y = Axis.2, group = Location), type = "norm", linetype = 2) +
  #scale_shape_manual(values=c(21, 22, 20)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))
ggsave("Richness_Per_Location.png", dpi = 300)


#Ploting bars
Plot <- plot_bar(Oomycete1, "Clade", fill = "Sand")
#To save plot for easier manipulation
Plot$data
#fircats for stacked bars
ggplot(Plot$data) +
  geom_bar(aes(x=reorder(Clade, -Abundance, na.rm = F),y= Abundance, fill = Sand), stat= "identity", width = 0.9) +
  theme_classic()+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.position="top",
        axis.text.x=element_text(angle =90, vjust = 0.5, hjust = 1),
        axis.text.y=element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14)) +
  coord_flip() +
  xlab("") +
  ylab("Abundance")
ggsave("Species_Abundance.png", dpi = 300)

Plot$data
#fircats for stacked bars without sand reference
p2021 = ggplot(Plot$data) +
  geom_bar(aes(x=reorder(Clade, -Abundance, FUN = sum),y= Abundance), stat= "identity", width = 0.9) +
  theme_classic()+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 14, color = "black"),
        legend.position="top",
        axis.text.x=element_text(angle =90, vjust = 0.5, hjust = 1),
        axis.text.y=element_text(size = 14), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 13), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14)) +
  coord_flip() +
  xlab("") +
  ylab("Abundance")
ggsave("Species_Abundance.png", dpi = 300)

p2022



#LinearModelPlot4CEC
oomycete1 <- data.frame(Oomycete1@sam_data) #2021
oomycete2 <- data.frame(Oomycete2@sam_data) #2022
oomycete.metadata <- Oomycete1@sam_data

oomycete.metadata <- Oomycete2@sam_data


ggplot(oomycete1) +
  geom_point(aes(x = Sand/100, y = richness, color = Location), size = 3.5)+
  geom_smooth(aes(x = Sand/100, y = richness), method='lm', se = F, color= "black")+
  theme_classic()+
  scale_color_manual(values = cbbPalette) +
  xlab("% Sand") +
  ylab("Oomycete Richness") + 
  theme(axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))+
  #scale_y_continuous(lim = c(0, 1), labels = scales::percent)+ 
  scale_x_continuous(lim = c(0, 1), labels = scales::percent)
ggsave("Sand_OII.png", dpi = 300)

#LinearModelPlot4CEC
ggplot(oomycete) +
  geom_point(aes(x = CEC, y = richness, color = Location), size = 3.5)+
  geom_smooth(aes(x = CEC, y = richness), method='lm', se = F, color= "black")+
  theme_classic()+
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) + #(for points and lines, scale_fill for plots)
  xlab("CEC") +
  ylab("Oomycete Richness") + 
  theme(axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 14), 
        legend.title = element_text(size = 15), 
        legend.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))


ggsave("CEC.png", dpi = 300)

#Shapiro's test (to test that input variable (edaphic factors) follows normal distribution)
shapiro.test(oomycete1$Sand) #value is greater than 0.05 so I assume normality

shapiro.test(oomycete$CEC)


#SpearmanCorrelation (CEC)

cor.test(oomycete2$Sand, oomycete2$richness, method=c("spearman"), exact = FALSE)


#PearsonCorrelation (Clay/Sand)
#cor.test(Oomycete$Clay, Oomycete$OomyceteIsolationIncidence, method=c("pearson"), exact = FALSE)

#INDICATOR SPECIES AND TERNARY PLOT
install.packages("indicspecies")
library(indicspecies)
indicator.tissue.fungi.corn <- indicspecies::multipatt(as.data.frame(t(Oomycete2@otu_table)), cluster = Oomycete2@sam_data$Location, func = "IndVal.g", control = how(nperm=9999))
# summary of results
summary(indicator.tissue.fungi.corn, indvalcomp = TRUE)

#Loading Ternary Plot
install.packages('Ternary')
library('Ternary')

#To create blank ternary plot
TernaryPlot()

fungi.ra.soy <- oomycetes %>%
  #subset_samples(Compartment == "Leaf" & Crop  == "Soy") %>%
  #phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE) %>%
  #phyloseq::transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%                                         # Melt to long format
  group_by(Species.1, Location) %>%
  summarize(mean_size = sum(Abundance,  na.rm = TRUE))
#nest() %>%
#mutate(mean.relabund = purrr::map(data,~mean(.$Abundance))) %>%
#mutate(median.relabund = purrr::map(data,~median(.$Abundance))) %>%
#mutate(SE.relabund = purrr::map(data,~sd(.$Abundance)/sqrt(length(.$Abundance)))) %>%
#unnest(c(mean.relabund, SE.relabund, median.relabund))

mean.Ternary <- pivot_wider(fungi.ra.soy, id_cols = c(Species.1), names_from = Location, values_from = c(mean_size))


library(tidyverse)
packages = c('ggtern', 'plotly', 'readr', 'dplyr', 'tidyr')
for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
  library(p, character.only = T)
}

colnames(mean.Ternary) <- c("Species" ,"Central", "North", "South")


ggtern()+
  geom_point(data = mean.Ternary[1:5,],aes(x=North, y=South , z= Central, fill = Species), size = 6, shape = 21, alpha = 0.5) +
  geom_point(data = mean.Ternary[6:10,],aes(x=North, y=South , z= Central, fill = Species), size = 6, shape = 22, alpha = 0.5) +
  geom_point(data = mean.Ternary[11:15,],aes(x=North, y=South , z= Central, fill = Species), size = 6, shape = 23, alpha = 0.5) +
  geom_point(data = mean.Ternary[16:20,],aes(x=North, y=South , z= Central, fill = Species), size = 6, shape = 24, alpha = 0.5) +
  geom_point(data = mean.Ternary[21:26,],aes(x=North, y=South , z= Central, fill = Species), size = 6, shape = 25, alpha = 0.5) +
  scale_fill_manual(values = c(rep(cbbPalette[1:5], 5), cbbPalette[6])) +
  limit_tern(1.1,1.1,1.1)+
  theme_showgrid() +
  theme_bw()+
  guides(fill = guide_legend(override.aes=list(shape=c(21,21,21,21,21,22,22,22,22,22,23,23,23,23,23,24,24,24,24,24,25,25,25,25,25,25), size = 4))) +
  scale_shape_manual(values = c(21,21,21,21,21,22,22,22,22,22,23,23,23,23,23,24,24,24,24,24,25,25,25,25,25,25)) +
  theme(legend.position="bottom",
        legend.title = element_blank())+
  labs( x   = "",
        xarrow  = "North",
        y       = "",
        yarrow  = "South",
        z       = "",
        zarrow  = "Central") +
  theme_showarrows()
ggsave("TernaryPlot.png", dpi = 300)






#ggplot(Oomycete, aes(Sand, OomyceteIsolationIncidence, color = `Location`))  +
geom_point(stat = "Identity")+
  geom_smooth(aes(x = CEC, y = OomyceteIsolationIncidence), method='lm', se = F, color = "black") +
  geom_bar(stat = "Identity")+
  scale_fill_hue(c = 50) +
  xlab("Sample") +
  ylab("Oomycete Isolation Incidence")

cor.test(CEC, OomyceteIsolationIncidence, method=c("spearman"))

#Oomycete Map
library(ggplot2)
library(tidyverse)
devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)
install.packages("remotes")
remotes::install_github("UrbanInstitute/urbnthemes", build_vignettes = TRUE)
library(urbnthemes)


library(tidyverse)
devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)

##Shows map with legend 
alabama <- countydata %>% 
  left_join(counties, by = "county_fips") %>% 
  filter(state_name =="Alabama") 
ggplot() +
  geom_polygon(data = alabama, mapping = aes(long, lat, group = group), fill = "white", color = "black", size = .25) +
  geom_point(data = samp_dat, aes(x=Longitude, y=Latitude, shape = as.factor(Year), fill = County), color = "black", size= 7, alpha = 0.5) +
  theme_void() +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(name="County", values = cbbPalette)  +
  guides(fill=guide_legend(override.aes=list(shape=21)))+
  labs(shape="Year", colour="County")
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
ggsave("Survey2021Map.png", dpi = 300)

#Plot for Seed Pathogenicity
Pathogenicity <- read.csv("~/Kemi/My Research/Oomycete Research/Oomycete_Combined/All Pathogenicity.csv", na.strings = "N/A")
Pathogenicity2 <- na.omit(Pathogenicity)

Pathogenicity %>%
  group_by(Isolate_Code, Species) %>%
  summarise(mean=mean(DSI)) %>%
  ggplot(aes(x=reorder(Species, mean), mean, fill = ""))  + 
  stat_summary(fun.data = mean_se, geom = "bar", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 1)) +
  stat_summary(fun=mean,geom="point", aes(group = Species), position = position_dodge(width = 1)) +
  geom_point(alpha = 0.5, position = position_jitterdodge()) +
  xlab("") +
  ylab("Disease Severity Index") + 
  theme_classic()+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"),
        legend.position="top",
        axis.text.x=element_text(angle =90, vjust = 1, hjust = 1, size = 14),
        axis.text.y=element_text(size = 15),
        axis.title.y=element_text(size = 15))
ggsave("Pathogenicity.png", dpi = 300)

Pathogenicity_2021 = Pathogenicity %>%
  subset(Year == 2021)


#Subsetting to 2021 repeated only
Pathogenicity %>%
  group_by(Isolate_Code, Species) %>%
  subset(Year == 2021) %>%
  summarise(mean=mean(DSI)) %>%
  ggplot(aes(x=reorder(Species, mean), mean, fill = ""))  + 
  stat_summary(fun.data = mean_se, geom = "bar", position = position_dodge(width = 0.9)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 1)) +
  stat_summary(fun=mean,geom="point", aes(group = Species), position = position_dodge(width = 1)) +
  geom_point(alpha = 0.5, position = position_jitterdodge()) +
  xlab("") +
  ylab("Disease Severity Index") + 
  theme_classic()+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"),
        legend.position="top",
        axis.text.x=element_text(angle =90, vjust = 0.3, hjust = 1, size = 14),
        axis.text.y=element_text(size = 15),
        axis.title.y=element_text(size = 15))
ggsave("Pathogenicity.png", dpi = 300)


#To get the supplementary table
Table = Pathogenicity %>%
  subset(Year == 2021) %>%
  group_by(Species) %>%
  summarize(avg = mean(DSI), 
            n = n(),
            sd = sd(DSI), se = sd/sqrt(n))
  Table <- kable(Table, digits = 3, format = "markdown")
  
  
#to convertn species to factor;
library (tibble)
data <-as_tibble (Pathogenicity2)%>%
  mutate (Species = factor (Species))
library(multcomp)
library(lsmeans)


fit <- aov(DSI ~ Species, data)

set.seed(20140123)
Dunnet <- glht(fit, linfct=mcp(Species = "Tukey"))
summary(Dunnet)


TukeyHSD(Pathogenicity2)

#trying to correct p=values using "bonferroni"


#General model AOV
fit1_20C <- lm(DSI ~ Species, data)
lsmeans_fit1_20C <- lsmeans(fit1_20C,"Species")
anova(fit1_20C)
CvsA_fit1_20C <- contrast(lsmeans_fit1_20C, "trt.vs.ctrl", ref=1, adjust = "bon")
Results_fit1_20C <- summary(CvsA_fit1_20C)
library(knitr)
kable(Results_fit1_20C, digits = 3, format = "markdown")

#General model AOV for Individual Isolates 
fit1_P <- lm(DSI ~ Isolate_Code, data)
lsmeans_fit1_P <- lsmeans(fit1_P,"Isolate_Code")
anova(fit1_P)
CvsA_fit1_P <- contrast(lsmeans_fit1_P, "trt.vs.ctrl", ref=253, adjust = "bon") #FOR POSITIVE CONTROL
Results_fit1_P <- summary(CvsA_fit1_P)
#TO FIND OUT REF FOR POSITIVE CONTROL (exported to excel for exact position)
P_REF = kable(lsmeans_fit1_P, digits = 3, format = "markdown")
library(knitr)
P_Contrast <- kable(Results_fit1_P, digits = 3, format = "markdown")

#Contrast for negative control
fit1_N <- lm(DSI ~ Isolate_Code, data)
lsmeans_fit1_N <- lsmeans(fit1_N,"Isolate_Code")
anova(fit1_N)
CvsA_fit1_N <- contrast(lsmeans_fit1_N, "trt.vs.ctrl", ref=1, adjust = "bon") #FOR NEGATIVE CONTROL
Results_fit1_N <- summary(CvsA_fit1_N)
library(knitr)
N_Contrast <- kable(Results_fit1_N, digits = 3, format = "markdown")

#TO EXTRACT Common Values 
#is.atomic(P_Contrast) #To fix error received due to using $ for atomic vector 
#getElement(P_Contrast, 'p.value')
#new.P_Contrast <- P_Contrast[['p.value' == 0.05]]
#seed.rot.soybean.25 <- subset(seed.rot, Strain %in% new.strains.soybean & Treatment == "+Pythium" & Crop == "Soybean")

Contrast_Negative.Control <- read.csv("~/Kemi/My Research/Oomycete Research/Contrast_Negative Control.csv")
Contrast_Positive.Control <- read.csv("~/Kemi/My Research/Oomycete Research/Contrast_Positive Control.csv")


intersect(Contrast_Negative.Control, Contrast_Positive.Control)


#Plotting diagnostic plots for fit1 model
par(mfrow=c(2,2)) # optional layout 
diag_plot_20C <- plot(fit1_20C)# diagnostic plots

#Plotting residuals
par(mfrow=c(1,1)) # optional layout 
hist_20C_res <- hist(fit1_20C$residuals)

#Test with random effect fit2 using aov
fit1.2_20C <- aov(DSI ~ Species + Species:Isolate_Code + Error(Set), data)
summary(fit1.2_20C, adjust="bon")

install.packages("lme4")
library("lme4")
install.packages("lsmeans")
library("lsmeans")

##Model 2 with fixed effect (species) and random effect (set)
fit2_20C <- lmer(DSI ~ Species + (1|Set), data, REML = FALSE)
summary(fit2_20C)

#Model 2 fitting (fitted vs residuals)
fitvsres_20C <- plot(fit2_20C)

## ----lsmeans 20C model 2, fig.align='center', fig.width=8, fig.height=10----
#lsmeans for model 2
lsmeans_fit2_20C <- lsmeans(fit2_20C,"Species")

#Print summary for data including mean, SE, df, CIs, t-ratio and p.value
summary(lsmeans_fit2_20C, infer=c(TRUE,TRUE), adjust="bon")

plot(lsmeans_fit2_20C)+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"),
        axis.text.x=element_text(size = 14),
        axis.text.y=element_text(size = 15),
        axis.title.y=element_text(size = 15))+
  theme_classic()
#Estimate confidence intervals for model 2
#confint(contrast(lsmeans_fit2, "trt.vs.ctrl", ref=3))

## ----lsmeans_model3_20C--------------------------------------------------
#Model 3 with fixed effect (species) and random effect (set) and nested effect (species:isolate)
fit3_20C <- lmer(DSI ~ Species + (1|Set) + (1|Isolate.Code:Species), data, REML = FALSE)
summary(fit3_20C)

#Model 3 fitting (fitted vs residuals)
plot(fit3_20C)

## ----model comparison, fig.align='center', fig.width=8, fig.height=10----
#lsmeans for model 3
lsmeans_fit3_20C <- lsmeans(fit3_20C, "Species")

#Print summary for data including mean, SE, df, CIs, t-ratio and p.value
#summary(lsmeans_fit3_20C, infer=c(TRUE,TRUE), adjust="bon")

plot(lsmeans_fit3_20C)

#Comparing models 2 and 3 (to determine which is better)
anova(fit2_20C, fit3_20C)
#(always go for lowest AIC & BIC)
## ----contrast for model 3------------------------------------------------
#Contrast for model 3
#CvsA_fit3_20C <- contrast(lsmeans_fit3_20C, "trt.vs.ctrl", ref=1)
#Results_fit3_20C <- summary(CvsA_fit3_20C)
#library(knitr)
#kable(Results_fit3_20C, digits = 3, format = "markdown")

CvsA_fit2_20C <- contrast(lsmeans_fit2_20C, "trt.vs.ctrl", ref=1)
Results_fit2_20C <- summary(CvsA_fit2_20C, adjust="bon")
library(knitr)
kable(Results_fit2_20C, digits = 3, format = "markdown")






library(lme4)

library(emmeans)

library(multcomp)


Pathogenicity_2021$Set <- as.factor(Pathogenicity_2021$Set)

lm.DON <- lmer(DSI ~ Species + (1|Isolate_Code), data = Pathogenicity_2021)
car::Anova(lm.DON)

plot(lm.DON)

lsmeans.DON <- emmeans(lm.DON, ~Species) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans.DON, alpha = 0.05, adjust = "bon", reversed = TRUE, details = TRUE) # contrast with Tukey ajustment
Results_lsmeansEC
kable(Results_lsmeansEC, digits = 3, format = "markdown")

col_vector<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
myCol = c("pink1", "violet", "mediumpurple1", "slateblue1", "purple", "purple3",
          "turquoise2", "skyblue", "steelblue", "blue2", "navyblue",
          "orange", "tomato", "coral2", "palevioletred", "violetred", "red2",
          "springgreen2", "yellowgreen", "palegreen4",
          "wheat2", "tan", "tan2", "tan3", "brown",
          "grey70", "grey50", "grey30")

Final.Virulence.Rating.for.Oomycetes.2021 <- read.csv("~/Final Virulence Rating for Oomycetes 2021.csv")

ggplot(Final.Virulence.Rating.for.Oomycetes.2021, aes(x = Rating, fill = Species )) + 
  geom_bar() +
  scale_fill_manual(values = fungi.colors) +
  xlab("Virulence Level") +
  #scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2))) +
  #stat_compare_means(method = "t.test") +
  theme_classic()
ggsave("Richness_Per_Location.png", dpi = 300)

#subset of table taking out non-confident data
Subset.Virulence.Grouping <- read.csv("~/Kemi/My Research/Oomycete Research/Subset Virulence Grouping.csv")
ggplot(Subset.Virulence.Grouping, aes(x = Rating, fill = Species )) + 
  geom_bar() +
  scale_fill_manual(values = fungi.colors) +
  xlab("Virulence Level") +
  #scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2))) +
  #stat_compare_means(method = "t.test") +
  theme_classic()+
  coord_flip() #used to flip plots

ggsave("Virulence.png", dpi = 300)

