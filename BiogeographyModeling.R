
####### Libraries ######
library(phyloseq)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(lme4)
library(emmeans)
library(ggcorrplot)
library(microshades)
library(psych)
library(broom)
library(AER)
library(MASS)

# for plotting the map
devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)
library(urbnthemes)
########################

###### Color palettes used in this code #####
# color blind friendly pallete 
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

fungi.colors <- c("#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c",
                  "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32",
                  "#fdae6b", "#f16913", "#d94801", "#8c2d04",
                  "#dadaeb", "#bcbddc", "#807dba", "#6a51a3", "#4a1486",
                  "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d",
                  "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525")


hex_values <-c(microshades_palette("micro_orange",3, lightest = FALSE), 
               microshades_palette("micro_blue",3, lightest = FALSE), 
               microshades_palette("micro_purple",3, lightest = FALSE), microshades_palette("micro_gray", 3, lightest = FALSE), microshades_palette("micro_brown",3, lightest = FALSE) ,microshades_palette("micro_green",2, lightest = FALSE)) 
################################

##### Generate phyloseq input from culture collection list ##### 
cultures <- read.csv("2023-05-08_IsolateCollection.csv")

##### Species Count Table #####
OTU.table <- cultures %>%
  group_by(Location_Code, Manuscript_Species) %>%
  count() %>%
  pivot_wider(names_from = Location_Code, values_from = n, values_fill = 0) %>%
  as.data.frame()

#write.csv(OTU.table, "OTUtable.csv")

rnames <- paste("Sp",rownames(OTU.table), sep = "")

rownames(OTU.table) <- rnames
otu_table_names <- OTU.table$Manuscript_Species
OTU.table <- OTU.table[,-1]
OTU <- phyloseq::otu_table(OTU.table, taxa_are_rows = TRUE)

##### Meta data ####
samp_dat <- read.csv("Metadata_Combined.csv", na.strings = "na")
rownames(samp_dat) <- samp_dat$Sample
samp_dat <- samp_dat[,-1]
SAMP <- phyloseq::sample_data(samp_dat)

##### Taxonomy Table ####
tax <- read.csv("Taxonomy.csv")
species_tax_table <- tax$Species
rownames(tax) <- tax$Species_num
tax <- tax[,-1]
TAX <- phyloseq::tax_table(as.matrix(tax))

# Sanity check - does the culture list match the taxonomy table
all.equal(sort(species_tax_table), sort(otu_table_names))

#Combines all data into a phyloseq object
oomycetes <- phyloseq::phyloseq(OTU, TAX, SAMP)

##### RDS File #####

# save an RDS file to make it easy and more reproducible to load in the data. 
saveRDS(oomycetes, "oomycete_phyloseq_final.RDS")
oomycetes <- readRDS("oomycete_phyloseq_final.RDS")

#### Figure 1. Map of sampling ####
meta.data <- data.frame(oomycetes@sam_data) # needed to grab the lat and long.
meta.data$County2 <- paste(meta.data$County, "County") 
alabama <- countydata %>% 
  left_join(counties, by = "county_fips") %>% 
  filter(state_name =="Alabama") 

alabama$sampled <- ifelse(alabama$county_name %in% unique(meta.data$County2), "Sampled", "Not Sampled")

ggplot() +
  geom_polygon(data = alabama, mapping = aes(long, lat, group = group, fill = sampled), color = "black", size = .25) +
  geom_point(data = meta.data, aes(x=Longitude, y=Latitude, shape = as.factor(Year)), fill = "black", size= 2, alpha = 0.5) +
  theme_void() +
  scale_shape_manual(values = c(21,22))+
  scale_fill_manual(name="County", values = c("white", cbbPalette[[1]]))  +
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
ggsave("SurveyMap.png", dpi = 300)

##### Figure 2. General Isolation Stats ####
## How many species were recovered in 2021?
oomycetes %>%
  subset_samples(Year == 2021) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# 25 species 

## How many species were recovered in 2022?
oomycetes %>%
  subset_samples(Year == 2022) %>%
  phyloseq::filter_taxa(function(x) sum(x) > 0, TRUE)

# 28 species 

### What was the distribution of species counts in 2021?
cultures %>%
  subset(Year == 2021) %>%
  group_by(Manuscript_Species) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  arrange(desc(freq)) %>%
  print(n = 30)

### What was the distribution of species counts in 2022?
cultures %>%
  subset(Year == 2022) %>%
  group_by(Manuscript_Species) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  arrange(desc(freq)) %>%
  print(n = 30)

### plot of species counts
fig2.species.counts <- cultures %>%
  group_by(Year, Manuscript_Species) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  arrange(desc(freq)) %>%
  ggplot(aes(x = reorder(Manuscript_Species, freq), y = freq, fill = as.factor(Year))) + 
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  coord_flip() + 
  scale_fill_manual(values = cbbPalette[3:4]) + 
  theme_classic() +
  xlab("") +
  ylab("Percent recovery") +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(title="Year")) +
  theme(legend.position = c(0.8, 0.2))
fig2.species.counts

#### plot of clade counts
clade.counts <- cultures %>%
  group_by(Year, Clade) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  arrange(desc(freq)) %>%
  ggplot(aes(x = reorder(Clade, freq), y = freq, fill = as.factor(Year))) + 
  geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
  coord_flip() + 
  scale_fill_manual(values = cbbPalette[3:4]) + 
  theme_classic() +
  xlab("") +
  ylab("Percent recovery") +
  scale_y_continuous(labels = scales::percent)
clade.counts

#### Alpha diversity ####
# Adds the alpha diversity measure to the metadata
oomycetes@sam_data$shannon <- estimate_richness(oomycetes, measures=c("Shannon"))$Shannon #shannon diversity index
oomycetes@sam_data$richness <- estimate_richness(oomycetes, measures=c("Observed"))$Observed #observed number of taxa
oomycetes@sam_data$even <- oomycetes@sam_data$shannon/log(oomycetes@sam_data$richness) #plieou's evenness

##### Correlation plot ####
meta.data <- data.frame(oomycetes@sam_data)

cor.plot.test <- meta.data %>%
  dplyr::select(Sand, Silt, Clay, CEC,	SOM, pH, Days.after.planting, Latitude,	Longitude, even, richness, Ave_Soil_Moisture, Sum_Precipitation, Ave_Soil_Temp) 

cor.plot <- psych::corr.test(cor.plot.test, method = "spearman", use = "complete.obs", adjust = "holm")

cor.r <- cor.plot$r
cor.p <- cor.plot$p

names <- c("% Sand", 
  "% Silt", 
  "% Clay",
  "CEC", 
  "SOM", 
  "pH", 
  "DAP", 
  "Lat.",
  "Long.",
  "Evenness",
  "Richness",
  "Ave. Soil Moist. (cb)",
  "Sum Precip.",
  "Ave. Soil Temp.")

colnames(cor.r) <- names
row.names(cor.r) <- names

corelation.plot <- ggcorrplot(cor.r, p.mat = cor.p, 
           hc.order = FALSE, 
           tl.cex = 5,
           outline.col = "grey", 
           type = "lower",
           insig = "blank", 
           ggtheme = ggplot2::theme_minimal,
           colors = c(cbbPalette[[3]], "white", cbbPalette[[4]]))

##### Location - North, Central, South ####
# for every X - unit change in X we onserved a exp^b times as many species.  
# 1 exp^b
# 2. Times as many for counts 
# 3. report confidence interval - exp^lower to exp^upper
#### Linear model - Richness
lm <- lm(richness ~ Location*as.factor(Year), data = data.frame(oomycetes@sam_data))
car::Anova(lm)
location.glm <- glm(richness ~ Location, data = data.frame(oomycetes@sam_data), family = "poisson")
car::Anova(location.glm)
AIC(lm)
AIC(location.glm) # glm has better fit

summary(location.glm) 
# do we have bad overdisperssion? residual deviance/ degrees freedom should be equal to 1 more or less
resid.deviance <- summary(location.glm)[[4]] # residual deviance
degree.freedom <- summary(location.glm)[[7]] # residual deviance
resid.deviance/degree.freedom # pretty good for overdispersion
dispersiontest(location.glm) # no significant overdispersion
tidy(location.glm, exponentiate = TRUE, conf.int = TRUE)
exp(location.glm$coefficients[[1]]) # for intercept - central location - mean is modeled as 4.25 species
exp(location.glm$coefficients[[2]]) # for Northern location - had 1.8 times as many species than central
exp(location.glm$coefficients[[3]]) # for Southern location - had 1.17 times as many species than central

# we found that northern fields had 1.82 times as many oomycete species than central or southern locations

# Mean separation
lsmeans <- emmeans(location.glm, ~Location) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE, Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC

# Extracting the letters for the bars
sig.diff.letters <- data.frame(Results_lsmeansEC$emmeans$Location, 
                               str_trim(Results_lsmeansEC$emmeans$.group))
colnames(sig.diff.letters) <- c("Location", 
                                "Letters")

# for plotting with letters from significance test
ave_richness2 <- data.frame(oomycetes@sam_data) %>%
  group_by(Location) %>%
  dplyr::summarize(
    ave.rich = mean(richness, na.rm=TRUE),
    n = n(),
    se = sd(richness)/sqrt(n)) %>%
  left_join(sig.diff.letters) 

meta.data$Location <- factor(meta.data$Location, levels = c("South", "Central ", "North"))

# plot with letters 
richness.location <- ggplot(meta.data, aes(x = Location, y = richness)) + 
  stat_summary(fun=mean,geom="bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5) +
  geom_jitter(size = 1, width = 0.1) +
  geom_text(data = ave_richness2, aes(label = Letters, y = ave.rich+(2*se)), vjust = -0.5) +
  ylab("Richness") + 
  xlab("Location") +
  xlab("")+
  theme_classic() 
richness.location

##### Soil properties, SAND and richness ####
pois.mod <- glm(richness ~ Sand, data=data.frame(oomycetes@sam_data), 
                family="poisson")
summary(pois.mod)
# do we have bad overdisperssion? residual deviance/ degrees freedom should be equal to 1 more or less
resid.deviance <- summary(pois.mod)[[4]] # residual deviance
degree.freedom <- summary(pois.mod)[[7]] # residual deviance
resid.deviance/degree.freedom # pretty good for overdispersion
car::Anova(pois.mod) #sand is a significant predictor of oomycete species richness
confint(pois.mod)
exp(pois.mod$coefficients[[2]]) # for modeled coefficient
exp(confint(pois.mod)[[2]]) # for lower CI of sand
exp(confint(pois.mod)[[4]]) # for upper CI of sand
AIC(pois.mod)
dispersiontest(pois.mod) # no significant overdispersion
tidy(pois.mod, exponentiate = TRUE, conf.int = TRUE)

exp(2.450757)
(0.98*50)+11 #60
exp(60)
exp(1.45)
# equation: y = e^2.45-0.017x
# We found that with every 1 % increase in sand, there were 0.98 times as many species isolated (p < 0.001)

sand.richness <- ggplot(oomycetes@sam_data, aes(x = Sand, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod), aes(x = Sand, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) +
  #scale_y_continuous(labels = scales::percent) + 
  xlab("% Sand") +
  ylab("Richness") 
sand.richness

##### CEC ####
pois.mod.cec <- glm(richness ~ CEC, data=data.frame(oomycetes@sam_data), 
                family="poisson")
summary(pois.mod.cec)
resid.deviance <- summary(pois.mod.cec)[[4]] # residual deviance
degree.freedom <- summary(pois.mod.cec)[[7]] # residual deviance
resid.deviance/degree.freedom # pretty good for overdispersion
car::Anova(pois.mod.cec)
confint(pois.mod.cec)
exp(pois.mod.cec$coefficients[[2]]) # for modeled coefficient
exp(confint(pois.mod.cec)[[2]]) # for lower CI of sand
exp(confint(pois.mod.cec)[[4]]) # for upper CI of sand
AIC(pois.mod.cec)
dispersiontest(pois.mod.cec)

tidy(pois.mod.cec, exponentiate = TRUE, conf.int = TRUE)

# We found that with every 1 unit increase in CEC, there were 1.08 (1.03 to 1.12) times more oomycete species isolated (p < 0.001)

cec.richness <- ggplot(oomycetes@sam_data, aes(x = CEC, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod.cec), aes(x = CEC, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("CEC (meq/100 g soil)") +
  ylab("Richness") 
cec.richness

##### Silty #####
pois.mod.silt <- glm(richness ~ Silt, data=data.frame(oomycetes@sam_data), 
                family="poisson")

exp(pois.mod.silt$coefficients[[2]]) # for modeled coefficient
summary(pois.mod.silt)
AIC(pois.mod.silt)
car::Anova(pois.mod.silt)
dispersiontest(pois.mod.silt)
# equation: y = e^1.08 + 0.023x 
silt.rich <- ggplot(oomycetes@sam_data, aes(x = Silt, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod.silt), aes(x = Silt, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("% Silt") +
  ylab("Richness") 
silt.rich

##### Clay ####
pois.mod.clay <- glm(richness ~ Clay, data=data.frame(oomycetes@sam_data), 
                     family="poisson")
summary(pois.mod.clay)
car::Anova(pois.mod.clay)
dispersiontest(pois.mod.clay)

clay.rich <- ggplot(oomycetes@sam_data, aes(x = Clay, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod.clay), aes(x = Clay, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5)  + 
  xlab("% Clay") +
  ylab("Richness") 
clay.rich

##### SOM #####
# regular poisson model was significanlty overdispersed so using a quasipoisson
pois.mod.SOM2 <- glm(richness ~ SOM, data=data.frame(oomycetes@sam_data), 
                     family="quasipoisson")

summary(pois.mod.SOM2)
car::Anova(pois.mod.SOM2)

som.rich <- ggplot(oomycetes@sam_data, aes(x = SOM, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod.SOM2), aes(x = SOM, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("SOM") +
  ylab("Richness") 
som.rich


##### Latitude ####
pois.mod.lat2 <- glm(richness ~ Latitude, data=data.frame(oomycetes@sam_data), 
                    family="quasipoisson")

summary(pois.mod.lat2)
AIC(pois.mod.lat2)
anova(pois.mod.lat2)

# equation
exp(-4.98056) * exp(0.20313*33)

# we found that for every degree increase in latitude we isolated 1.23 times as many oomycete species


lat.rich <- ggplot(oomycetes@sam_data, aes(x = Latitude, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod.lat2), aes(x = Latitude, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("Latitude") +
  ylab("Richness") 
lat.rich


##### Longitude ####
pois.mod.long2 <- glm(richness ~ Longitude, data=data.frame(oomycetes@sam_data), 
                     family="poisson")

summary(pois.mod.long2)
AIC(pois.mod.lat2)
dispersiontest(pois.mod.long2)

# we found that for every degree increase in longitude we isolated 0.69 times as many oomycete species

long.rich <- ggplot(oomycetes@sam_data, aes(x = Longitude, y = richness, color = Latitude)) + 
  
  geom_line(data = fortify(pois.mod.long2), aes(x = Longitude, y = exp(.fitted)), color = "black") +
  #geom_smooth(method = "glm", se = F,
  #method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("Longitude") +
  ylab("Richness") 
long.rich


##### Days after planting - not sig ####
pois.mod.dap <- glm(richness ~ Days.after.planting, data=data.frame(oomycetes@sam_data), 
                      family="poisson")

summary(pois.mod.dap)
car::Anova(pois.mod.dap)
##### Seed treatment - not sig ####
pois.mod.seedtreatment <- glm(richness ~ Seed.treatment, data=data.frame(oomycetes@sam_data), 
                    family="quasipoisson")

summary(pois.mod.seedtreatment)
car::Anova(pois.mod.seedtreatment)
##### Variety - not sig ####
pois.mod.seedvariety <- glm(richness ~ Seed.Variety, data=data.frame(oomycetes@sam_data), 
                              family="quasipoisson")

summary(pois.mod.seedvariety)
car::Anova(pois.mod.seedvariety)
##### Precipitation - not sig ####
pois.mod.precip <- glm(richness ~ Sum_Precipitation, data=data.frame(oomycetes@sam_data), 
                    family="poisson")

summary(pois.mod.precip)
car::Anova(pois.mod.precip)

##### Soil moisture - not sig ####
pois.mod.sm <- glm(richness ~ Ave_Soil_Moisture, data=data.frame(oomycetes@sam_data), 
                       family="poisson")

summary(pois.mod.sm)
car::Anova(pois.mod.sm)
##### Soil temp - not sig ####
pois.mod.st <- glm(richness ~ Ave_Soil_Temp, data=data.frame(oomycetes@sam_data), 
                   family="poisson")

summary(pois.mod.st)
car::Anova(pois.mod.st)



##### Soil temp - not sig ####
pois.mod.year <- glm(richness ~ as.factor(Year), data=data.frame(oomycetes@sam_data), 
                   family="quasipoisson")

summary(pois.mod.year)
car::Anova(pois.mod.year)

####### Figure 3 - plot #####
plot1 <- ggarrange(lat.rich,
                   clay.rich, 
                   silt.rich, 
          cec.richness, 
          sand.richness,
          som.rich, labels = "auto", common.legend = T, ncol = 2, nrow = 3)

#### Specific species isolation successes #### 

#### Irregulare #### 
irregulare.presence <- oomycetes %>%
  subset_taxa(Species=="Globisporangium irregulare") %>%
  psmelt() %>%
  mutate(present = ifelse(Abundance > 0, 1, 0))

log.glm <- glm(present ~ Latitude*as.factor(Year), data=irregulare.presence, 
                     family=binomial(link = "logit"))
summary(log.glm)
car::Anova(log.glm)
exp(0.02479)
# We found that for each 1 degree increase in Latitude Globisporangium irregulare was  

#### Phytophthora nicotianea #### 
phytophthora.presence <- oomycetes %>%
  subset_taxa(Species=="Phytophthora nicotianae") %>%
  psmelt() %>%
  mutate(present = ifelse(Abundance > 0, 1, 0))

log.glm.phyt <- glm(present ~ Latitude, data=phytophthora.presence, 
               family=binomial(link = "logit"))
car::Anova(log.glm.phyt)
summary(log.glm.phyt)
exp(0.8766) #2.4 times as likely with increase in latitude
confint(log.glm.phyt)
exp(0.2702209)
exp(1.662061)
# We found that for each 1 degree increase in Latitude Phytophthora nicotianae was 2.4 (1.31 to 5.27; 95% CL) times as likely to be isolated (P = 0.012). 


#### Beta Diversity #####

# beta diversity example with canned phyloseq ordinations. I can teach you how to do it mannually too so you can have more control over the plotting. 
pcoa.ordination <- ordinate(oomycetes, "MDS", "jaccard")

global.data <- data.frame(pcoa.ordination$vectors)
global.data$Sample <- rownames(global.data)
meta.data$Sample <- rownames(meta.data)
global.data2 <- left_join(meta.data, global.data, by = "Sample")
global.data.axis1.var <- pcoa.ordination$values$Relative_eig[[1]]
global.data.axis2.var <- pcoa.ordination$values$Relative_eig[[2]]

#Permanova - Permutational Multivariate ANOVA

oomycete.dist.jaccard = phyloseq::distance(oomycetes, "jaccard") #creates distance matrix (0 (similar) 1 (non-similar))
adonis2(oomycete.dist.jaccard~as.factor(Year), as(sample_data(oomycetes), "data.frame"), permutations = 9999) 

### plot ###
ggplot() + 
  geom_point(data = global.data2, aes(x = Axis.1, y = Axis.2, shape = as.factor(Year), fill = Location), alpha = 0.8, size = 2) +
  theme_bw() +
  ylab("PCoA2") + 
  xlab("PCoA1") +
  scale_fill_manual(values=cbbPalette) +
  scale_shape_manual(values=c(21, 22)) +
  guides(fill=guide_legend(override.aes=list(shape=21)))


####### Figure 4 - Plotting bar plot with sand colored in #########
plot.bar <- plot_bar(oomycetes, x = "Species", y = "Abundance", fill = "Sand", facet_grid = "Year")

ggplot(plot.bar$data) +
  geom_col(aes(x=reorder(Label, -Abundance), y= Abundance, fill = Sand), width = 0.9) +
  theme_classic()+
  theme(strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
        strip.text.x = element_text(size = 14, color = "black"),
        axis.text.x=element_text(angle =90, vjust = 0.5, hjust = 1)) +
  scale_fill_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 50) +
  xlab("") +
  ylab("Abundance") +
  facet_wrap(~Year)



