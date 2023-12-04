
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
library(psych)
library(broom)
library(AER)
library(MASS)

# for plotting the map
#devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)
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

################################

##### Figure 1. Map of sampling ####
# Load in acrage planted map
cottonacres <- read.csv("CottonAcresPlanted.csv", na.strings = "na")


meta.data$County2 <- paste(meta.data$County, "County") 
alabama <- countydata %>% 
  left_join(counties, by = "county_fips") %>% 
  filter(state_name =="Alabama") 
head(alabama)
alabama2 <- left_join(alabama, cottonacres, by = c("county_name"))

alabama2$sampled <- ifelse(alabama2$county_name %in% unique(meta.data$County2), "Sampled", "Not Sampled")

ggplot() +
  geom_polygon(data = alabama2, mapping = aes(long, lat, group = group, fill = Hectares), color = "black", size = .25) +
  geom_point(data = meta.data, aes(x=Longitude, y=Latitude, shape = as.factor(Year)), size= 2, alpha = 0.5, position=position_jitter(h=0.02, w=0.03)) +
  theme_void() +
  scale_shape_manual(values = c(21,24))+
  scale_fill_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 5000, na.value = "white", breaks = c(2000, 4000, 6000, 8000, 10000, 12000)) +
  labs(shape="Year", fill="Ave. Hectares \n Planted") +
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




##### Read in weather data #####
weather <- read.csv("WeatherData_Oomycetes/combine_weather.csv", na.strings = "na")
weather$Date_time <- as.Date(weather$Date_time, "%m/%d/%y")
ave_weather_tocollection <- weather %>%
  group_by(Sample) %>%
  summarise(sum_precip = sum(Precipitation),
            mean_min_soil_temp = mean(Min_Soil_Temperature),
            mean_ave_soil_moist = mean(Soil_Moisture)
            )

weather_three_days <- weather %>%
  group_by(Sample) %>%
  top_n(-3, wt = Date_time) %>%
  summarise(sum_precip_threeday = sum(Precipitation), 
            mean_soil_temp_threeday = mean(Min_Soil_Temperature),
            soil_moist_threeday = mean(Soil_Moisture)) # cbar higher is dryer

weather_five_days <- weather %>%
  group_by(Sample) %>%
  top_n(-5, wt = Date_time) %>%
  summarise(sum_precip_fiveday = sum(Precipitation), 
            min_soil_temp_fiveday = mean(Min_Soil_Temperature),
            soil_moist_fiveday = mean(Soil_Moisture)) # cbar higher is dryer


##### Generate phyloseq input from culture collection list ##### 
cultures <- read.csv("2023-05-08_IsolateCollection.csv")

### How many isolates were isolated
n.isolates <- cultures %>%
  group_by(Location_Code) %>%
  summarise(cnt = n()) 


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
str(samp_dat)
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
oomycetes@sam_data

##### RDS File #####
# save an RDS file to make it easy and more reproducible to load in the data. 
saveRDS(oomycetes, "oomycete_phyloseq_final.RDS")
oomycetes <- readRDS("oomycete_phyloseq_final.RDS")

# New meta.data with all the weather information on it. 
meta.data <- data.frame(oomycetes@sam_data) # needed to grab the lat and long.


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
  group_by(Year, Label) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>% 
  arrange(desc(freq)) %>%
  ggplot(aes(x = reorder(Label, freq), y = freq, fill = as.factor(Year))) + 
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), color = "black") +
  coord_flip() + 
  scale_fill_manual(values = cbbPalette[3:4]) + 
  theme_classic() +
  xlab("") +
  ylab("Percent recovery") +
  scale_y_continuous(labels = scales::percent) +
  guides(fill=guide_legend(title="Year")) 
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

##### Alpha diversity ####

# Adds the alpha diversity measure to the metadata
oomycetes@sam_data$shannon <- estimate_richness(oomycetes, measures=c("Shannon"))$Shannon #shannon diversity index
oomycetes@sam_data$simpson <- estimate_richness(oomycetes, measures=c("Simpson"))$Simpson #simpson diversity index
oomycetes@sam_data$richness <- estimate_richness(oomycetes, measures=c("Observed"))$Observed #observed number of taxa
oomycetes@sam_data$even <- oomycetes@sam_data$shannon/log(oomycetes@sam_data$richness) #plieou's evenness

### Alpha diversity on 30-year normals ### 
merged.oomycete.30year <- merge_samples(sample_data(oomycetes), "Location_Code") # averages everything by the location 
merged.oomycete.30year$year_30_norm_precip_may_cm <- merged.oomycete.30year$year_30_norm_precip_may*2.54
merged.oomycete.30year$annual_30_year_normal_cm <- merged.oomycete.30year$annual_30_year_normal_in*2.54
merged.oomycete.30year$spring_30_year_precip_cm <- merged.oomycete.30year$spring_30_year_precip_in*2.54
merged.oomycete.30year$richness <- round(merged.oomycete.30year$richness)

#### Correlation plot with the merged dataset 30-year averages #####
meta.data.30year <- data.frame(merged.oomycete.30year)

cor.plot.test <- meta.data.30year %>%
  dplyr::select(Sand, # % sand
                Silt, # % silt
                Clay, # % clay
                CEC,	# % CEC
                SOM, # % Soil organic matter
                pH, # ph of soil
                Days.to.planting, # days to planting after start of year
                Latitude,	
                Longitude, 
                sum_precip, # sum precipitation from planting to collection
                mean_min_soil_temp, # average minimum soil temp from planting to collection
                mean_ave_soil_moist, # average soil moist (higher is dryer) from planting to collection
                richness, # richness, number of species
                shannon,# shannon diversity index
                simpson, # simpson diversity index
                annual_30_year_normal_cm
  ) 

cor.plot <- psych::corr.test(cor.plot.test, method = "spearman", use = "pairwise")

cor.r <- cor.plot$r
cor.p <- cor.plot$p

names <- c("% Sand", 
           "% Silt", 
           "% Clay",
           "CEC", 
           "SOM", 
           "pH", 
           "Days to planting",
           "Lat.",
           "Long.",
           "Sum precipitation",
           "Mean Min Soil Temp.",
           "Mean Soil Moisture (centibar)",
           "Richness",
           "Shannon",
           "Simpson",
           "30Y-Normal Precip."
)

colnames(cor.r) <- names
row.names(cor.r) <- names

year30norm.r <- as.matrix(cor.r[,16])
year30norm.p <- as.matrix(cor.p[,16])

colnames(year30norm.r) <- "30Y-Normal Precip."


corelation.plot.30year <- ggcorrplot(year30norm.r, 
                              tl.cex = 5,
                              outline.col = "grey", 
                              p.mat = year30norm.p,
                              #type = "lower",
                              insig = "blank", 
                              ggtheme = ggplot2::theme_void,
                              colors = c(cbbPalette[[3]], "white", cbbPalette[[4]]))

#### Correlation plot with the full dataset ####

## NOTE I am including the 30 year precipitation data, but this includes the same values for several fields and thus represents a pseudoreplication. Therefore, I manually remove the significant squares so it matches the merged correlation plot done above. 
meta.data <- data.frame(oomycetes@sam_data)

cor.plot.test <- meta.data %>%
  dplyr::select(Sand, # % sand
                Silt, # % silt
                Clay, # % clay
                CEC,	# % CEC
                SOM, # % Soil organic matter
                pH, # ph of soil
                Days.to.planting, # days to planting after start of year
                Latitude,	
                Longitude, 
                sum_precip, # sum precipitation from planting to collection
                mean_min_soil_temp, # average minimum soil temp from planting to collection
                mean_ave_soil_moist, # average soil moist (higher is dryer) from planting to collection
                richness, # richness, number of species
                shannon,# shannon diversity index
                simpson, # simpson diversity index
                annual_30_year_normal_in
                ) 

cor.plot <- psych::corr.test(cor.plot.test, method = "spearman")

cor.r <- cor.plot$r
cor.p <- cor.plot$p

# manual correction to avoid overinterpretation based on peudoreplication
cor.p[16,8] <- 1
cor.p[16,9] <- 1
cor.p[16,12] <- 1

names <- c("% Sand", 
 "% Silt", 
 "% Clay",
 "CEC", 
 "SOM", 
 "pH", 
 "Days to planting",
 "Lat.",
 "Long.",
 "Sum precipitation",
 "Mean Min Soil Temp.",
 "Mean Soil Moisture (centibar)",
 "Richness",
 "Shannon",
 "Simpson",
 "30Y-Normal Precip."
)

colnames(cor.r) <- names
row.names(cor.r) <- names

corelation.plot <- ggcorrplot(cor.r, 
           tl.cex = 5,
           outline.col = "grey", 
           p.mat = cor.p,
           type = "lower",
           insig = "blank", 
           ggtheme = ggplot2::theme_void,
           colors = c(cbbPalette[[3]], "white", cbbPalette[[4]]))

#### Location - North, Central, South ####
# for every X - unit change in X we observed a exp^b times as many species.  
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

#### Thirty year normals #####
pois.mod.30year <- glm(richness ~ annual_30_year_normal_cm, data=data.frame(merged.oomycete.30year), 
                family="poisson")
summary(pois.mod.30year)

# plot with letters 
richness.30year <- ggplot(merged.oomycete.30year, aes(x = annual_30_year_normal_cm, y = richness, color = Latitude)) + 
  geom_point() +
  stat_cor(method = "spearman") +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) +
  ylab("Richness") + 
  xlab("30-year-normal Precipitation (cm)") +
  theme_classic() 
richness.30year
#### Soil properties, SAND ####
pois.mod <- glm(richness ~ Sand, data=data.frame(oomycetes@sam_data), 
                family="poisson")
summary(pois.mod)
# do we have bad overdisperssion? residual deviance/ degrees freedom should be equal to 1 more or less
resid.deviance <- summary(pois.mod)[[4]] # residual deviance
degree.freedom <- summary(pois.mod)[[7]] # residual deviance
resid.deviance/degree.freedom # pretty good for overdispersion
car::Anova(pois.mod) #sand is a significant predictor of oomycete species richness
confint(pois.mod)
exp(pois.mod$coefficients[[1]]) # for modeled coefficient the intercept
exp(10*-pois.mod$coefficients[[2]]) # for modeled coefficient the beta
exp(10*-confint(pois.mod)[[2]]) # for lower CI of sand
exp(10*-confint(pois.mod)[[4]]) # for upper CI of sand
AIC(pois.mod)
dispersiontest(pois.mod) # no significant overdispersion
tidy(pois.mod, exponentiate = TRUE, conf.int = TRUE)

# equation: y = e^11.8-0.983x
# We found that with every 1 % increase in sand, there were 0.98 times as many species isolated (p < 0.001)

sand.richness <- ggplot(oomycetes@sam_data, aes(x = Sand, y = richness, color = Latitude)) + 
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) +
  #scale_y_continuous(labels = scales::percent) + 
  xlab("% Sand") +
  ylab("Richness") 
sand.richness

#### CEC ####
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
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("CEC (meq/100 g soil)") +
  ylab("Richness") 
cec.richness

#### Silty #####
pois.mod.silt <- glm(richness ~ Silt, data=data.frame(oomycetes@sam_data), 
                family="poisson")

exp(pois.mod.silt$coefficients[[2]]) # for modeled coefficient
summary(pois.mod.silt)
AIC(pois.mod.silt)
car::Anova(pois.mod.silt)
exp(10*pois.mod.silt$coefficients[[2]]) # for modeled coefficient - beta
exp(10*confint(pois.mod.silt)[[2]]) # for lower CI of silt for 10% increase
exp(10*confint(pois.mod.silt)[[4]]) # for upper CI of silt for 10% increase

dispersiontest(pois.mod.silt)
# equation: y = e^1.08 + 0.023x 
silt.rich <- ggplot(oomycetes@sam_data, aes(x = Silt, y = richness, color = Latitude)) + 
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("% Silt") +
  ylab("Richness") 
silt.rich

#### Clay ####
pois.mod.clay <- glm(richness ~ Clay, data=data.frame(oomycetes@sam_data), 
                     family="poisson")
summary(pois.mod.clay)
car::Anova(pois.mod.clay)
dispersiontest(pois.mod.clay)
exp(pois.mod.clay$coefficients[[2]]) # for modeled coefficient
exp(confint(pois.mod.clay)[[2]]) # for lower CI of sand
exp(confint(pois.mod.clay)[[4]]) # for upper CI of sand



clay.rich <- ggplot(oomycetes@sam_data, aes(x = Clay, y = richness, color = Latitude)) + 
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5)  + 
  xlab("% Clay") +
  ylab("Richness") 
clay.rich

#### SOM #####
pois.mod.SOM <- glm(richness ~ SOM, data=data.frame(oomycetes@sam_data), 
                     family="poisson")

dispersiontest(pois.mod.SOM)
exp(pois.mod.SOM$coefficients[[2]]) # for modeled coefficient
exp(confint(pois.mod.SOM)[[2]]) # for lower CI of sand
exp(confint(pois.mod.SOM)[[4]]) # for upper CI of sand
summary(pois.mod.SOM)
car::Anova(pois.mod.SOM)

som.rich <- ggplot(oomycetes@sam_data, aes(x = SOM, y = richness, color = Latitude)) + 
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("SOM") +
  ylab("Richness") 
som.rich


#### Latitude ####
pois.mod.lat <- glm(richness ~ Latitude, data=data.frame(oomycetes@sam_data), 
                    family="poisson")

dispersiontest(pois.mod.lat)
exp(pois.mod.lat$coefficients[[2]]) # for modeled coefficient
exp(confint(pois.mod.lat)[[2]]) # for lower CI of sand
exp(confint(pois.mod.lat)[[4]]) # for upper CI of sand


summary(pois.mod.lat)
AIC(pois.mod.lat)
car::Anova(pois.mod.lat)
exp(0.21078)
# equation
exp(-5.24888) * exp(0.21078*33)

# we found that for every degree increase in latitude we isolated 1.23 times as many oomycete species

# Quadratic regression due to the drop in richness at midlatitudes
meta.data$Latitude2 <- meta.data$Latitude^2
quadratic.model <-lm(richness ~ Latitude + I(Latitude^2), data=meta.data)
summary(quadratic.model)
AIC(quadratic.model)
AIC(pois.mod.lat)
lat.rich <- ggplot(oomycetes@sam_data, aes(x = Latitude, y = richness)) + 
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", se = FALSE) +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("Latitude") +
  ylab("Richness") 
lat.rich


#### Longitude ####
pois.mod.long2 <- glm(richness ~ Longitude, data=data.frame(oomycetes@sam_data), 
                     family="poisson")

summary(pois.mod.long2)
AIC(pois.mod.long2)
dispersiontest(pois.mod.long2)
exp(pois.mod.long2$coefficients[[1]]) # for modeled coefficient - intercept
exp(-pois.mod.long2$coefficients[[2]]) # for modeled coefficient - beta
exp(-confint(pois.mod.long2)[[2]]) # for lower CL
exp(-confint(pois.mod.long2)[[4]]) # for upper CL

# equation
exp(-31.7147) * exp(-0.3841*33)
# we found that for every degree increase in longitude we isolated 0.68 times as many oomycete species

long.rich <- ggplot(oomycetes@sam_data, aes(x = Longitude, y = richness, color = Latitude)) + 
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("Longitude") +
  ylab("Richness") 
long.rich

#### Days after planting - not sig ####
pois.mod.dap <- glm(richness ~ Days.after.planting, data=data.frame(oomycetes@sam_data), 
                      family="poisson")

summary(pois.mod.dap)
car::Anova(pois.mod.dap)
#### Seed treatment - not sig ####
pois.mod.seedtreatment <- glm(richness ~ Seed.treatment, data=data.frame(oomycetes@sam_data), 
                    family="poisson")

summary(pois.mod.seedtreatment)
car::Anova(pois.mod.seedtreatment)
#### Variety - not sig ####
meta.data <- data.frame(oomycetes@sam_data)
seed.treatment.meta.data <- meta.data %>%
  filter(Seed.Variety %in% c("DP 1646 B2XF", "DP 2038 B3XF"))
pois.mod.seedvariety <- glm(richness ~ Seed.Variety, data=seed.treatment.meta.data, 
                              family="poisson")
dispersiontest(pois.mod.seedvariety)
summary(pois.mod.seedvariety)
car::Anova(pois.mod.seedvariety)

lsmeans <- emmeans(pois.mod.seedvariety, ~Seed.Variety) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE, Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC

ggplot(seed.treatment.meta.data, aes(x = Seed.Variety, y = richness)) +
  geom_boxplot()


#### Precipitation from planting to collection - not sig ####
pois.mod.precip <- glm(richness ~ sum_precip, data=data.frame(oomycetes@sam_data), 
                    family="poisson")

summary(pois.mod.precip)
car::Anova(pois.mod.precip)

#### Soil moisture from planting to collection - not sig ####
pois.mod.sm <- glm(richness ~ mean_ave_soil_moist, data=data.frame(oomycetes@sam_data), 
                       family="poisson")

summary(pois.mod.sm)
car::Anova(pois.mod.sm)

som.rich <- ggplot(oomycetes@sam_data, aes(x = mean_ave_soil_moist, y = richness, color = Latitude)) + 
  
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  geom_point() +
  theme_classic() + 
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  xlab("Soil Moisture (cbar)") +
  ylab("Richness") 
som.rich


#### Soil temp from planting to collection - not sig ####
pois.mod.st <- glm(richness ~ mean_min_soil_temp, data=data.frame(oomycetes@sam_data), 
                   family="poisson")

summary(pois.mod.st)
car::Anova(pois.mod.st)



#### Year ####
pois.mod.year <- glm(richness ~ as.factor(Year), data=data.frame(oomycetes@sam_data), 
                   family="poisson")

summary(pois.mod.year)
car::Anova(pois.mod.year)

####### Figure 3 - plot #####

plot1 <- ggarrange(long.rich,
                   clay.rich, 
                   silt.rich, 
          cec.richness, 
          sand.richness,
          som.rich,
          labels = c("b", "c", "d", "e", "f", "g"), common.legend = T, ncol = 3, nrow = 2)

ggarrange(corelation.plot, plot1, labels = "a", nrow = 1)


#### Beta Diversity #####

# The strategy for beta diversity analysis was to use Jaccard index on binary data since we don't feel that with culture dependent methods we fully captured the abundance and thus would not be appropraite to use Bray-Curtis. 

##### 30Y normal precipitation values - based on averages across county.####
merged.oomycete.30year.physeq <- merge_samples(oomycetes, "Location_Code") # averages everything by the location 
meta.data.30year
env <- meta.data.30year %>% 
  dplyr::select(Sand, # % sand
                Silt, # % silt
                Clay, # % clay
                CEC,	# % CEC
                SOM, # % Soil organic matter
                pH,
                annual_30_year_normal_in # ph of soil

  )

# pcoa #
ord <- wcmdscale(vegdist(merged.oomycete.30year.physeq@otu_table, distance = "jaccard", binary = T), eig = T)
plot(ord)
ord$eig/sum(ord$eig)*100 # proportion explained each axis
fit <- envfit(ord, env, perm = 9999, na.rm = TRUE) # fitting the environmental data to the first two pcoa axes


####### Permanova - Permutational Multivariate ANOVA ####
set.seed(12573)
dist.jac <- vegdist(t(oomycetes@otu_table), distance = "jaccard", binary = T) # calcluate distance
adonis2(dist.jac~as.factor(Year) + Location + Seed.treatment + Seed.Variety, as(sample_data(oomycetes), "data.frame"), permutations = 9999) # PERMANOVA

####### Figure 4a - Plotting heat map with soil texture#########
meta.data$Location_Code <- rownames(meta.data)
datf <- cultures %>%
  group_by(Location_Code, Label) %>%
  summarise(cnt = n()) %>%
  mutate(freq = cnt / sum(cnt)) %>%
  left_join(meta.data, by = "Location_Code")

datt <- left_join(fig2.species.counts$data, datf, by = "Label") # just for ordering properties

Fig4heatmap <- ggplot(data = datt, aes(x = reorder(Location_Code, Latitude), y = reorder(Label, freq.x), fill = Silt+Clay)) + 
  geom_tile(colour = "black") +
  scale_size_continuous(name = "Count") +
  scale_fill_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 60) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(face = "italic")) +
    xlab("") + 
    ylab("")

####### Figure 4b - PcoA with Env. fit data #####
## Environment fit analysis
env <- meta.data %>% 
  dplyr::select(Sand, # % sand
                Silt, # % silt
                Clay, # % clay
                CEC,	# % CEC
                SOM, # % Soil organic matter
                pH, # ph of soil
                sum_precip,
                mean_min_soil_temp,
                mean_ave_soil_moist
  )

# pcoa #
ord <- wcmdscale(vegdist(t(oomycetes@otu_table), distance = "jaccard", binary = T), eig = T)
ord$eig/sum(ord$eig)*100 # proportion explained each axis
fit <- envfit(ord, env, perm = 9999, na.rm = TRUE) # fitting the environmental data to the first two pcoa axes

df_ord <- as.data.frame(vegan::scores(ord, display = "sites"))
df_arrows <- as.data.frame(scores(fit, "vectors"))
mult <- scale_arrow(df_arrows, df_ord[ , c("Dim1", "Dim2")])
df_arrows <- mult * df_arrows
df_arrows$var <- rownames(df_arrows)
df_arrows$p.val <- fit$vectors$pvals
df_ord$Sample <- rownames(df_ord)
meta.data$Sample <- rownames(meta.data)
global.data2 <- left_join(meta.data, df_ord, by = "Sample")
df_arrows$sig <- ifelse(df_arrows$p.val <= 0.05, "Sig.", "Not Sig.")

### Figure 4b
pcoa.jaccard <- ggplot() + 
  geom_point(data = global.data2, aes(x = Dim1, y = Dim2, color = Location, shape = as.factor(Year))) + 
  geom_segment(data = df_arrows, aes(x = 0, xend = Dim1, y = 0, yend = Dim2, color = sig), 
               arrow = arrow(length = unit(0.1,"cm"))) + 
  geom_text(data=df_arrows, aes(x=Dim1, y=Dim2, label=var), hjust="outward", size = 2) + 
  theme_classic() + 
  scale_color_manual(values = c(cbbPalette)) + 
  xlab("PcoA1 [26.77 %]") +
  ylab("PcoA2 [19.56 %]")

#### Specific species isolation successes #### 

######## Figure 4c ######
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

# Specific species plots
fig4c <- ggplot(phytophthora.presence, aes(x = Latitude, y = present)) +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = binomial(link = "logit")), color = "black") +
  geom_point(position=position_jitter(h=0.00, w=0.06), alpha = 0.6) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1), breaks = c(0,1)) + 
  ylab("Presence/Absence of \n Phytophthora nicotianae")


######## Figure 4d #######
#### Irregulare #### 
irregulare.presence <- oomycetes %>%
  subset_taxa(Species=="Globisporangium irregulare") %>%
  psmelt() %>%
  mutate(present = ifelse(Abundance > 0, 1, 0))

log.glm <- glm(present ~ Latitude, data=irregulare.presence, 
               family=binomial(link = "logit"))
summary(log.glm)
car::Anova(log.glm)

# We found that for each 1 degree increase in Latitude Globisporangium irregulare was  

# Specific species plots
fig4d <- ggplot(irregulare.presence, aes(x = Latitude, y = present)) +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = binomial(link = "logit")), color = "black") +
  geom_point(position=position_jitter(h=0.00, w=0.06), alpha = 0.6) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1), breaks = c(0,1)) + 
  ylab("Presence/Absence of \n Globisporangium irregulare")

######## Weather Data ##### 

meta.data$rain_per_day <- meta.data$sum_precip/meta.data$Days.after.planting

rain <- lm(rain_per_day ~ Location*as.factor(Year) + cnt, data = meta.data)
anova(rain)

lsmeans <- emmeans(rain, ~Location|as.factor(Year)) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE, Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC

ggplot(meta.data, aes(x = as.factor(Year), y = rain_per_day, fill = Location)) + 
  stat_summary(fun=mean,geom="bar", position = "dodge") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, position = "dodge") 
  
ggplot(meta.data, aes(x = annual_30_year_normal_in, y = richness, color = Location)) + 
  geom_point() +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = "poisson"), color = "black") +
  stat_cor()




