#Exposure Related Effects of Zequanox Data and Analysis
#James Luoma's Data from U.S.G.S.

## @knitr load_library

library(tidyverse)
library(plyr)
library(stats)
library(lme4)
library(lmerTest)
library(nlme)
library(gridExtra)
library(lubridate)
library(formattable)
library(ggthemes)
library(forestplot)
library(broom)
library(moments)

set.seed(123)



## @knitr load_data

#Read in data and format
james <- read.csv("James_Data.csv")
james$Species <- factor(james$Species) #Lake Surgeon and Lake Trout
james$Tank <- factor(james$Tank)
james$Treatment <- factor(james$Treatment) #0, 50, 100 mg/L
james$Condition.Factor <- james$Condition.Factor * 100
james$UniqueTank <- paste(james$Tank, james$Species, sep = "-")
james$UniqueTank <- factor(james$UniqueTank)
james$UniqueFish <- paste(james$Tank, james$Fish, sep = "-")



## @knitr exploratory_graphics

ggplot(aes(x = UniqueTank, y = Condition.Factor, color = Treatment, shape = Species), data = james) + geom_jitter() + coord_flip()
ggplot(aes(x = Treatment, y = Condition.Factor, fill = Species), data = james) + geom_violin()
ggplot(aes(x = Condition.Factor, fill = Treatment), data = james) + geom_density(alpha = 0.5)

levels(james$Treatment) <- c("Control", "50 mg/L", "100 mg/L")
ggplot(aes(x = Condition.Factor, fill = Species), data = james) + facet_grid(~ Treatment) + geom_density(alpha = 0.5) + xlab("Condition Factor") + ylab("Density") + theme_bw() + theme(strip.background = element_blank()) + scale_fill_colorblind(labels = c("Lake Trout", "Lake Surgeon")) #Significant differences (USE?)

#Separate by species for graphics
james_lt <- james %>%
  filter(Species == "LAT")
james_ls <- james %>%
  filter(Species == "LST")

lt_plot <- ggplot(aes(x = Condition.Factor, fill = Treatment), data = james_lt) + geom_density(alpha = 0.5) + ylab("Density") + xlab("Condition Factor") + scale_fill_colorblind() #Good
ls_plot <- ggplot(aes(x = Condition.Factor, fill = Treatment), data = james_ls) + geom_density(alpha = 0.5) + ylab("Density") + xlab("Condition Factor") + theme(legend.position = "none") + scale_fill_colorblind() #Still gross

grid.arrange(lt_plot, ls_plot, ncol = 2) #give fish label



## @knitr james_pooled

james_pool_dat <- james %>%
  group_by(UniqueTank, Species, Treatment) %>%
  dplyr::summarise(Pooled = mean(Condition.Factor), .groups = "keep")
james_pooled <- lm(Pooled ~ Treatment + Species, data = james_pool_dat)
summary(james_pooled)
plot(james_pooled)

#95% CI
confint(james_pooled)



## @knitr james_fixed

james_fixed <- lm(Condition.Factor ~ UniqueTank + Treatment + Species,
                  data = james)
summary(james_fixed)
plot(james_fixed)

#95% CI
confint(james_fixed)



## @knitr james_ml

james_unrest <- lme(Condition.Factor ~ Treatment + Species,
                    random = ~ 1 | UniqueTank/UniqueFish,
                    data = james, method = "ML")
summary(james_unrest)
anova(james_unrest)

#95% CI
james_unrest_ci <- intervals(james_unrest)
james_unrest_ci



## @knitr james_reml

james_rest <- lme(Condition.Factor ~ Treatment + Species,
                  random = ~ 1 | UniqueTank/UniqueFish,
                  data = james, method = "REML")
summary(james_rest)
anova(james_rest)

#95% CI
james_rest_ci <- intervals(james_rest)
james_rest_ci
