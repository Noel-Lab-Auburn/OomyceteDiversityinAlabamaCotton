library(ggplot2)
library(tidyr)
library(ggpubr)
library(lme4)
library(emmeans)
library(stringr)
library(dplyr)

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
options(dplyr.print_max = 1e9)
STAND <- read.csv("Pathogenicity.csv", na.strings = "NA")

counted.STAND <- STAND %>%
  group_by(Label) %>%
  dplyr::count() %>%
  mutate(is.count = n/6) %>%
  mutate(keep = ifelse(is.count > 1, "Keep", "No")) %>%
  mutate(keep = ifelse(Label %in% c("G. acanthophoron",
                                    "Py. oligandrum/cedri",
                                    "G. nunn", 
                                    "G. longandrum", 
                                    "Py. periplocum"), "No", keep)) %>%
  dplyr::filter(keep == "Keep") 

filtered.STAND <- subset(STAND, Label %in% counted.STAND$Label) 

lm1 <- lm(DSI_1 ~ Label, data = filtered.STAND)
lm2 <- lmer(DSI_1 ~ Label + (1|Isolate.Code/Trial), data = filtered.STAND)
lm3 <- lmer(DSI_1 ~ Label + (1|Isolate.Code/Trial) + (1|Isolate.Code/Label), data = filtered.STAND)

anova(lm1, lm2, lm3)


AIC(lm3)
summary(lm2)
plot(lm2)
car::Anova(lm2)

lsmeans <- emmeans(lm2, ~Label) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, adjust = "bon", reversed = TRUE, details = TRUE, Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC

# Extracting the letters for the bars
sig.diff.letters <- data.frame(Results_lsmeansEC$emmeans$Label, 
                               str_trim(Results_lsmeansEC$emmeans$.group))
colnames(sig.diff.letters) <- c("Label", 
                                "Letters")

# for plotting with letters from significance test
ave_dsi <- filtered.STAND %>%
  group_by(Label) %>%
  dplyr::summarize(
    ave.dsi = mean(DSI_1, na.rm=TRUE),
    n = n(),
    se = sd(DSI_1)/sqrt(n)) %>%
  left_join(sig.diff.letters) %>%
  arrange(-ave.dsi)

# plot with letters 
dsi.path <- ggplot(ave_dsi, aes(x = reorder(Label, -ave.dsi), y = ave.dsi)) + 
  geom_col()+
  geom_errorbar(aes(x=reorder(Label, -ave.dsi), ymin=ave.dsi+(se), ymax=ave.dsi-(se)), width = 0.5)+
  geom_text(aes(label = Letters, y = ave.dsi+(2*se)), vjust = -0.5) +
  ylab("% Disease Severity Index") + 
  xlab("Location") +
  xlab("")+
  theme_classic() +
  theme(
    strip.background = element_rect(color="white", fill="white", size=1.5, linetype="solid"),
    strip.text.x = element_text(size = 12, color = "black"),
    legend.position="right",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1))
dsi.path

means <- filtered.STAND %>%
  group_by(Label) %>%
  summarise(Mean = mean(DSI_1), 
            n = n(), 
            sd.dev = sd(DSI_1)) %>%
  mutate(std.err = sd.dev/sqrt(n)) %>%
  arrange(-Mean) 

write.csv(means, "Mean_DSI.csv")




