#Exposure Related Effects of Zequanox Data and Analysis
#James Luoma's Data from USGS

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
james$Condition.Factor <- james$Condition.Factor * 100000
james$UniqueTank <- paste(james$Tank, james$Species, sep = "-")
james$UniqueTank <- factor(james$UniqueTank)
james$UniqueFish <- paste(james$Tank, james$Fish, sep = "-")



## @knitr exploratory_graphics

ggplot(aes(x = UniqueTank, y = Condition.Factor,
           color = Treatment, shape = Species), data = james) +
  geom_jitter() + coord_flip()
ggplot(aes(x = Treatment, y = Condition.Factor, fill = Species),
       data = james) +
  geom_violin()
ggplot(aes(x = Condition.Factor, fill = Treatment), data = james) +
  geom_density(alpha = 0.5)
ggplot(aes(x = Condition.Factor, fill = Species), data = james) +
  facet_grid(~ Treatment) +
  geom_density(alpha = 0.5) +
  xlab("Condition Factor") +
  ylab("Density") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  scale_fill_colorblind(labels = c("Lake Trout", "Lake Sturgeon"))



## @knitr james_density

levels(james$Treatment) <- c("Control", "50 mg/L", "100 mg/L")
james$Species2 <- james$Species
james$Species2 <- factor(james$Species2,
                                 levels = c("LST", "LAT"))

ggplot(aes(x = Condition.Factor, fill = Treatment), data = james) +
  facet_grid(~ Species2,
             labeller = labeller(Species2 = c("LAT" = "Lake Trout",
                                             "LST" = "Lake Sturgeon")),
             scales = "free_x") +
  geom_density(alpha = 0.5) +
  xlab("Condition Factor") +
  ylab("Density") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_fill_colorblind()



## @knitr count_table

levels(james$Treatment) <- c("Control", "50 mg/L", "100 mg/L")
james$Treatment <- factor(james$Treatment,
                          levels = c("Control", "50 mg/L", "100 mg/L"))

count_table <- james %>%
  group_by(Tank, Species, Treatment) %>%
  dplyr::summarise(Count = n(), .groups = "keep")
formattable(count_table,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr summary_stats

summary_table <- james %>%
  group_by(Tank, Species) %>%
  dplyr::summarise(MinWeight = min(Weight), AvgeWeight = mean(Weight),
                   MaxWeight = max(Weight), MinLength = min(Length),
                   AvgLength = mean(Length), MaxLength = max(Length),
                   MinCon = min(Condition.Factor),
                   AvgCon = mean(Condition.Factor),
                   MaxCon = max(Condition.Factor),
                   Count = n(), .groups = "keep")
summary_table <- data.frame(summary_table)
summary_table$Species <- as.character(summary_table$Species)
summary_table$Species <- ifelse(summary_table$Species == "LAT",
                                "Lake Trout", "Lake Sturgeon")
names(summary_table) <- c("Tank", "Species", "Minimum Weight",
                          "Average Weight", "Maximum Weight", "Minimum Length",
                          "Average Lenght", "Maximum Length",
                          "Minimum Condition Factor",
                          "Average Conition Factor",
                          "Maximum Condition Factor", "Count")
summary_table[, 3:12] <- round(summary_table[, 3:12], 4)
formattable(summary_table,
            align = c("l", "c", "c", "c", "c", "c", "c", "c",
                      "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr james_pooled

james_pool_dat <- james %>%
  group_by(UniqueTank, Species, Treatment) %>%
  dplyr::summarise(Pooled = mean(Condition.Factor), .groups = "keep")
james_pooled <- lm(Pooled ~ Treatment + Species, data = james_pool_dat)
#plot(james_pooled)
pool_sum <- summary(james_pooled)$coefficients
pool_ci <- confint(james_pooled)
pool_table <- cbind(pool_sum[, 1:2], pool_ci, pool_sum[, 3:4])
pool_table <- data.frame(pool_table)
pool_table <- round(pool_table, 4)
colnames(pool_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                          "Upper 95% CI", "t-Value", "p-Value")
pool_table2 <- pool_table
pool_table2$`p-Value` <- as.character(pool_table2$`p-Value`)
pool_table2[c(1, 4), 6] <- "<0.0001"
formattable(pool_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr james_fixed

james_fixed <- lm(Condition.Factor ~ Treatment + Species + UniqueTank,
                  data = james)
#plot(james_fixed)
fixed_sum <- summary(james_fixed)$coefficients
fixed_ci <- confint(james_fixed)
fixed_table <- cbind(fixed_sum[1:4, 1:2], fixed_ci[1:4, ], fixed_sum[1:4, 3:4])
fixed_table <- data.frame(fixed_table)
fixed_table <- round(fixed_table, 4)
colnames(fixed_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                          "Upper 95% CI", "t-Value", "p-Value")
fixed_table2 <- fixed_table
fixed_table2$`p-Value` <- as.character(fixed_table2$`p-Value`)
fixed_table2[c(1, 4), 6] <- "<0.0001"
formattable(fixed_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr james_ml

james_unrest <- lme(Condition.Factor ~ Treatment + Species,
                    random = ~ 1 | UniqueTank / UniqueFish,
                    data = james, method = "ML")
unrest_sum <- summary(james_unrest)$tTable
unrest_ci <- intervals(james_unrest)$fixed
unrest_table <- cbind(unrest_sum[, 1:2], unrest_ci[, -2], unrest_sum[, 4:5])
unrest_table <- data.frame(unrest_table)
unrest_table <- round(unrest_table, 4)
colnames(unrest_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                           "Upper 95% CI", "t-Value", "p-Value")
unrest_table2 <- unrest_table
unrest_table2$`p-Value` <- as.character(unrest_table2$`p-Value`)
unrest_table2[c(1, 4), 6] <- "<0.0001"
unrest_table2[2, 6] <- "0.0440"
formattable(unrest_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr james_reml

james_rest <- lme(Condition.Factor ~ Treatment + Species,
                  random = ~ 1 | UniqueTank / UniqueFish,
                  data = james, method = "REML")
rest_sum <- summary(james_rest)$tTable
rest_ci <- intervals(james_rest)$fixed
rest_table <- cbind(rest_sum[, 1:2], rest_ci[, -2], rest_sum[, 4:5])
rest_table <- data.frame(rest_table)
rest_table <- round(rest_table, 4)
colnames(rest_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                            "Upper 95% CI", "t-Value", "p-Value")
rest_table2 <- rest_table
rest_table2$`p-Value` <- as.character(rest_table2$`p-Value`)
rest_table2[c(1, 4), 6] <- "<0.0001"
formattable(rest_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr anova_comparison

james_pooled_anova <- data.frame(anova(james_pooled))
pooled_den_df <- james_pooled_anova[3, 1]
pooled_df <- cbind(james_pooled_anova[1, 1], pooled_den_df,
                james_pooled_anova[1, 4:5])
pooled_df <- data.frame(pooled_df)
names(pooled_df) <- c("Numerator DF", "Denominator DF", "F-Value", "p-Value")
rownames(pooled_df) <- "Pooled Treatment"

james_fixed_anova <- data.frame(anova(james_fixed))
fixed_den_df <- james_fixed_anova[4, 1]
fixed_df <- cbind(james_fixed_anova[1, 1], fixed_den_df,
                  james_fixed_anova[1, 4:5])
fixed_df <- data.frame(fixed_df)
names(fixed_df) <- c("Numerator DF", "Denominator DF", "F-Value", "p-Value")
rownames(fixed_df) <- "Fixed Treatment"

james_unrest_anova <- data.frame(anova(james_unrest))
names(james_unrest_anova) <- c("Numerator DF", "Denominator DF",
                               "F-Value", "p-Value")
unrest_df <- james_unrest_anova[2, ]
rownames(unrest_df) <- "ML Treatment"

james_rest_anova <- data.frame(anova(james_rest))
names(james_rest_anova) <- c("Numerator DF", "Denominator DF",
                               "F-Value", "p-Value")
rest_df <- james_rest_anova[2, ]
rownames(rest_df) <- "REML Treatment"


james_anova <- data.frame(rbind(pooled_df, fixed_df, unrest_df, rest_df))
james_anova <- round(james_anova, 4)
names(james_anova) <- c("Numerator DF", "Denominator DF",
                               "F-Value", "p-Value")
james_anova$`p-Value` <- as.character(james_anova$`p-Value`)
james_anova[2, 4] <- "<0.0001"
formattable(james_anova,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))
