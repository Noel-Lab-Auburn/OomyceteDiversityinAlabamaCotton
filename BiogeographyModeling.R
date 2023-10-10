
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

##### Read in weather data #####
weather <- read.csv("WeatherData_Oomycetes/combine_weather.csv", na.strings = "na")
weather$Date_time <- as.Date(weather$Date_time, "%m/%d/%y")
ave_weather_tocollection <- weather %>%
  group_by(Sample) %>%
  summarise(mean_min_air_temp = mean(Min_Air_Temperature),
            mean_max_air_temp = mean(Max_Air_Temperature),
            mean_ave_air_temp = mean(Average_Temperature),
            sum_precip = sum(Precipitation),
            min_soil_temp = min(Min_Soil_Temperature),
            mean_max_soil_temp = mean(Max_Soil_Temperature),
            mean_ave_soil_temp = mean(Soil_Temperature),
            mean_ave_soil_moist = mean(Soil_Moisture),
            min_min_soil_moist = min(Min_Soil_Moisture), # cbar higher is dryer
            #std
            std_min_air_temp = sd(Min_Air_Temperature),
            std_max_air_temp = sd(Max_Air_Temperature),
            std_ave_air_temp = sd(Average_Temperature),
            std_precip = sd(Precipitation),
            std_min_soil_temp = sd(Min_Soil_Temperature),
            std_max_soil_temp = sd(Max_Soil_Temperature),
            std_ave_soil_temp = sd(Soil_Temperature),
            std_ave_soil_moist = sd(Soil_Moisture),
            std_max_soil_moist = sd(Max_Soil_Moisture),
            std_min_soil_moist = sd(Min_Soil_Moisture)
  )

weather_three_days <- weather %>%
  group_by(Sample) %>%
  top_n(-3, wt = Date_time) %>%
  summarise(sum_precip_threeday = sum(Precipitation), 
            min_soil_temp_threeday = min(Min_Soil_Temperature),
            soil_moist_threeday = min(Soil_Moisture)) # cbar higher is dryer

weather_five_days <- weather %>%
  group_by(Sample) %>%
  top_n(-5, wt = Date_time) %>%
  summarise(sum_precip_fiveday = sum(Precipitation), 
            min_soil_temp_fiveday = min(Min_Soil_Temperature),
            soil_moist_fiveday = max(Soil_Moisture)) # cbar higher is dryer


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
samp_dat2 <- left_join(samp_dat, ave_weather_tocollection, by = "Sample")
samp_dat3 <- left_join(samp_dat2, n.isolates, by = c("Sample" = "Location_Code"))
samp_dat4 <- left_join(samp_dat3, weather_three_days, by = "Sample")
samp_dat5 <- left_join(samp_dat4, weather_five_days, by = "Sample")

rownames(samp_dat5) <- samp_dat5$Sample
samp_dat5 <- samp_dat5[,-1]
SAMP <- phyloseq::sample_data(samp_dat5)

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

##### Get the closest NOAA weather station for searching the 30 year normals ##### 
#### From Pat Schloss Tutorial on how to get local past weather from the closeset NOAA weather station
inventory_url <- "https://www.ncei.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"

inventory <- read_table(inventory_url,
                        col_names = c("station", "lat", "lon", "variable", "start", "end"))

samp_dat5$lat_rads <- samp_dat5$Latitude * 2 * pi / 360
samp_dat5$long_rads <- samp_dat5$Longitude * 2 * pi / 360
#my_lat <- 42.33831964441621 * 2 * pi / 360
#my_lon <- -83.88938389977316 * 2 * pi / 360


# Distance, d = 3963.0 * arccos[(sin(lat1) * sin(lat2)) + cos(lat1) * cos(lat2) * cos(long2 – long1)]
# 
# The obtained distance, d, is in miles. If you want your value to be in units of kilometers, multiple d by 1.609344.
# d in kilometers = 1.609344 * d in miles

inventory$lat_r <- inventory$lat *2 *pi/360
inventory$lon_r <- inventory$lon *2 *pi/360

weather.sum2 <- NULL
for (i in 1:nrow(samp_dat5)) {
  
  inventory$d_miles = 1.609344 * 3963 * acos((sin(inventory$lat_r) * sin(samp_dat5$lat_rads[[i]])) + cos(inventory$lat_r) * cos(samp_dat5$lat_rads[[i]]) * cos(samp_dat5$long_rads[[i]] - inventory$lon_r))
  station <- filter(inventory, start < 1991 & end > 2020) %>%
    top_n(n = -1, d_miles) %>%
    distinct(station, d_miles) 
  samp_dat5_i <- cbind.data.frame(samp_dat5[i,], station)
  weather.sum2 <- rbind.data.frame(samp_dat5_i, weather.sum2)
} 

weather.sum2$kilometers <- weather.sum2$d_miles * 1.60934
weather.sum2$kilometers

write.csv(weather.sum2, "weather_expanded.csv")

### Specifically for EVSmith, since the station it chose does not have 30- year norms. the next closest one does.
inventory$d_miles = 1.609344 * 3963 * acos((sin(inventory$lat_r) * sin(samp_dat5$lat_rads[[11]])) + cos(inventory$lat_r) * cos(samp_dat5$lat_rads[[11]]) * cos(samp_dat5$long_rads[[11]] - inventory$lon_r))
station <- filter(inventory, start < 1991 & end > 2020) %>%
  top_n(n = -10, d_miles) %>%
  distinct(station, d_miles) 

#### Figure 1. Map of sampling ####
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

#### Alpha diversity ####

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




# plot with letters 
richness.30year <- ggplot(merged.oomycete.30year, aes(x = year_30_norm_precip_may_cm, y = richness)) + 
  geom_point() +
  geom_smooth(method = "lm", se = F, color = "black") +
  ylab("Richness") + 
  xlab("30-year-normal Precipitation (cm)") +
  theme_classic() +
  stat_cor(label.y = 2)+ 
  stat_regline_equation(label.y = 2.5) 
richness.30year


##### Correlation plot ####
meta.data <- data.frame(oomycetes@sam_data)

cor.plot.test <- meta.data %>%
  dplyr::select(Sand, # % sand
                Silt, # % silt
                Clay, # % clay
                CEC,	# % CEC
                SOM, # % Soil organic matter
                pH, # ph of soil
                Days.after.planting, # days after planting collection
                Days.to.planting, # days to planting after start of year
                Latitude,	
                Longitude, 
                sum_precip, # sum precipitation from planting to collection
                min_soil_temp, # min soil temp from planting to collection
                min_min_soil_moist, # min soil moist (higher is dryer) from planting to collection
                sum_precip_threeday, # sum precipitation three days post planting
                min_soil_temp_threeday, # min soil temp three days post planting
                soil_moist_threeday, # soil moisture three days post planting
                richness, # richness, number of species
                shannon,# shannon diversity index
                simpson, # simpson diversity index
                even, # Pileou's evenness
                cnt # count of number of isolates recovered
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
 "Days after planting", 
 "Days to planting",
 "Lat.",
 "Long.",
 "Sum precipitation",
 "Min. Soil Temperature",
 "Min. Soil Moisture (centibar)",
 "Sum precipitation \n three days post planting",
 "Min. Soil Temperature \n three days post planting",
 "Min. Soil Moisture (centibar) \n three days post planting",
 "Richness",
 "Shannon",
 "Simpson",
 "Evenness",
 "Number of oomycete \n isolates recovered")

colnames(cor.r) <- names
row.names(cor.r) <- names

corelation.plot <- ggcorrplot(cor.r, 
           hc.order = TRUE, 
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

##### Soil properties, SAND ####
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


sand.richness <- ggplot(oomycetes@sam_data, aes(x = Silt + Clay, y = richness, color = Latitude)) + 
  
  #geom_line(data = fortify(pois.mod), aes(x = Sand, y = exp(.fitted)), color = "black") +
  geom_smooth(method = "glm", se = F,
  method.args = list(family = "poisson"), color = "black") +
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

# Quadratic regression due to the drop in richness at midlatitudes
meta.data$Latitude2 <- meta.data$Latitude^2
quadratic.model <-lm(richness ~ Latitude + I(Latitude^2), data=meta.data)
summary(quadratic.model)
AIC(quadratic.model)
AIC(pois.mod.lat2)
lat.rich <- ggplot(oomycetes@sam_data, aes(x = Latitude, y = richness)) + 
  
  geom_line(data = fortify(pois.mod.lat2), aes(x = Latitude, y = exp(.fitted)), color = "black", size = 0.5) +
  #geom_line(data = fortify(quadratic.model), aes(x = Latitude, y = exp(.fitted)), color = "black") +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 0.5, color = "black", se = FALSE) +
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
##### Precipitation from planting to collection - not sig ####
pois.mod.precip <- glm(richness ~ sum_precip, data=data.frame(oomycetes@sam_data), 
                    family="poisson")

summary(pois.mod.precip)
car::Anova(pois.mod.precip)

##### Soil moisture from planting to collection - not sig ####
pois.mod.sm <- glm(richness ~ min_min_soil_moist, data=data.frame(oomycetes@sam_data), 
                       family="poisson")

summary(pois.mod.sm)
car::Anova(pois.mod.sm)
##### Soil temp from planting to collection - not sig ####
pois.mod.st <- glm(richness ~ min_soil_temp, data=data.frame(oomycetes@sam_data), 
                   family="poisson")

summary(pois.mod.st)
car::Anova(pois.mod.st)

##### Precipitation from three day - not sig ####
pois.mod.precip <- glm(richness ~ sum_precip_threeday, data=data.frame(oomycetes@sam_data), 
                       family="poisson")

summary(pois.mod.precip)
car::Anova(pois.mod.precip)

##### Soil moisture from three day - not sig ####
pois.mod.sm <- glm(richness ~ soil_moist_threeday, data=data.frame(oomycetes@sam_data), 
                   family="poisson")

summary(pois.mod.sm)
car::Anova(pois.mod.sm)
##### Soil temp three day - not sig ####
pois.mod.st <- glm(richness ~ min_soil_temp_threeday, data=data.frame(oomycetes@sam_data), 
                   family="poisson")

summary(pois.mod.st)
car::Anova(pois.mod.st)


##### Year ####
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
adonis2(oomycete.dist.jaccard~as.factor(Year)+
          Latitude+
          Longitude+
          Sand+
          Silt+
          Clay+
          CEC+
          Days.after.planting, as(sample_data(oomycetes), "data.frame"), permutations = 9999, by = "onedf") 

### plot ###
ggplot() + 
  geom_point(data = global.data2, aes(x = Axis.1, y = Axis.2, shape = as.factor(Year), color = Latitude), alpha = 0.8, size = 2) +
  theme_bw() +
  ylab("PCoA2") + 
  xlab("PCoA1") +
  scale_color_gradient2(low = cbbPalette[[3]], mid = "grey", high = cbbPalette[[4]], midpoint = 32.5) + 
  guides(fill=guide_legend(override.aes=list(shape=21)))


####### Figure 4 - Plotting heat map with soil texture#########
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

# We found that for each 1 degree increase in Latitude Globisporangium irregulare was  

# Specific species plots
fig4b <- ggplot(irregulare.presence, aes(x = Latitude, y = present)) +
  geom_smooth(method = "glm", se = F,
              method.args = list(family = binomial(link = "logit")), color = "black") +
  geom_point(position=position_jitter(h=0.00, w=0.06), alpha = 0.6) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 1), breaks = c(0,1)) + 
  ylab("Presence/Absence of \n Globisporangium irregulare")

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

ggarrange(Fig4heatmap, label = "a", 
          ggarrange(fig4b, fig4c, nrow = 2, labels = c("b", "c")))

  
  
  
