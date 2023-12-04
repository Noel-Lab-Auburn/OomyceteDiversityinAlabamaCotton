library(ggplot2)
library(tidyr)
library(ggpubr)
library(lme4)
library(emmeans)
library(stringr)
library(dplyr)

cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
options(dplyr.print_max = 1e9)
path <- read.csv("Pathogenicity.csv", na.strings = "NA")

# Number of isolates tested
counted <- path %>%
  group_by(Label) %>%
  dplyr::count() %>%
  mutate(is.count = n/6) # divide by six due to three reps across two trials = 6 plates per isolate

# Number of isolates tested
unique(path$IsolateCode)

lm1 <- lm(DSI ~ Label, data = path)
lm2 <- lmer(DSI ~ Label + (1|Set), data = path)
lm3 <- lmer(DSI ~ Label + (1|Set) + (1|Trial), data = path)
lm4 <- lmer(DSI ~ Label + (1|Set) + (1|IsolateCode/Label), data = path)
lm5 <- lmer(DSI ~ Label + (1|Set) + (1|IsolateCode/Label) + (1|Trial), data = path)

# The LRT tells us that the fourth model is the best, which is the one with isolates within set and isolates within species as a random factor
anova(lm5, lm4, lm3, lm2, lm1) # Going with Model 4

# Basically models that include the Trial as a random factor were no better. So dropping that as a random factor did not matter. 
# Lowest AIC is also lm4

anova(lm4, lm2) # just comparing lm2 to lm4 and lm4 is better with including isolate nested within species as a random factor. 

plot(lm4)
car::Anova(lm4)

lsmeans <- emmeans(lm4, ~Label) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, adjust = "bon", reversed = TRUE, details = TRUE, Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC

# Extracting the letters for the bars
sig.diff.letters <- data.frame(Results_lsmeansEC$emmeans$Label, 
                               str_trim(Results_lsmeansEC$emmeans$.group))
colnames(sig.diff.letters) <- c("Label", 
                                "Letters")

# for plotting with letters from significance test
ave_dsi <- path %>%
  group_by(Label) %>%
  dplyr::summarize(
    ave.dsi = mean(DSI, na.rm=TRUE),
    n = n(),
    se = sd(DSI)/sqrt(n)) %>%
  left_join(sig.diff.letters) %>%
  arrange(-ave.dsi)

ave_dsi$Label <- factor(ave_dsi$Label, levels = c("Positive Control",
                                                  "Negative Control",
                                                  "G. attrantheridium",
                                                  "G. heterothallicum",
                                                  "G. irregulare",
                                                  "G. orthogonon",
                                                  "G. perplexum",
                                                  "G. rostratifingens",
                                                  "G. spinosum",
                                                  "G. sp. Lev1523",
                                                  "G. sylvaticum",
                                                  
                                                  "Phy. cucurbitacearum",
                                                  "Phy. helicoides",
                                                  
                                                  "P. nicotianae",
                                                  
                                                  "Py. acanthicum",
                                                  "Py. deliense",
                                                  "Py. dissotocum",
                                                  "Py. torulosum"))

x.axis.labels <- c("Positive Control",
  "Negative Control",
  "G. attrantheridium (n = 4)",
  "G. heterothallicum (n = 2)",
  "G. irregulare (n = 8)",
  "G. orthogonon (n = 3)",
  "G. perplexum (n = 4)",
  "G. rostratifingens (n = 4)",
  "G. spinosum (n = 3)",
  "G. sp. Lev1523 (n = 2)",
  "G. sylvaticum (n = 4)",
  
  "Phy. cucurbitacearum (n = 4)",
  "Phy. helicoides (n = 3)",
  
  "P. nicotianae (n = 5)",
  
  "Py. acanthicum (n = 3)",
  "Py. deliense (n = 6)",
  "Py. dissotocum (n = 2)",
  "Py. torulosum (n = 3)")

# plot with letters 
dsi.path <- ggplot(ave_dsi, aes(x = Label, y = ave.dsi)) + 
  geom_col()+
  geom_errorbar(aes(x=Label, ymin=ave.dsi+(se), ymax=ave.dsi-(se)), width = 0.5)+
  geom_text(aes(label = Letters, y = ave.dsi+(2*se)), vjust = -0.5) +
  scale_x_discrete(labels = x.axis.labels)+
  ylab("% Disease Severity Index") + 
  xlab("Location") +
  xlab("")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1, face = "italic"))
dsi.path

means <- path %>%
  group_by(Label) %>%
  summarise(Mean = mean(DSI), 
            n = n(), 
            sd.dev = sd(DSI)) %>%
  mutate(std.err = sd.dev/sqrt(n)) %>%
  arrange(-Mean) 

write.csv(means, "Mean_DSI.csv")




