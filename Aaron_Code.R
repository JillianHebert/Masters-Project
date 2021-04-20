#Acoustic Telemetry Evaluation of Carbon Dioxide Data and Analysis
#Aaron Cupp's Data from USGS

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

#Load original data for graphics and format
aaron <- read.csv("Aaron_Data.csv")
aaron$hour <- factor(aaron$hour)
aaron$species <- factor(aaron$species) #Bighead Carp and Grass Carp
aaron$level <- factor(aaron$level)
aaron$level <- mapvalues(aaron$level,
                         from = c("high", "medium", "low"),
                         to = c("High", "Medium", "Low"))
aaron$trial <- ifelse(aaron$trial == "10", "6", aaron$trial) #Keep ordered
aaron$trial <- factor(as.numeric(aaron$trial))
aaron$period <- ifelse(aaron$period == "acclimation",
                       "Acclimation", "Treatment")
aaron$period <- factor(aaron$period)
aaron$trtside <- ifelse(aaron$trtside == "east", "East", "West")
aaron$trtside <- factor(aaron$trtside)

#Load data from adehabitat package for analysis and format
aaron_new <- read.csv("aaron_data_final.csv")
aaron_new$date <- as.Date(ymd(aaron_new$date))
aaron_new$trial <- ifelse(aaron_new$trial == "10", "6", aaron_new$trial)
aaron_new$trial <- as.factor(aaron_new$trial)
aaron_new$period <- ifelse(aaron_new$period == "acclimation",
                           "Acclimation", "Treatment")
aaron_new$period <- factor(aaron_new$period)
aaron_new$level <- mapvalues(aaron_new$level, from = c("high", "medium", "low"),
                             to = c("High", "Medium", "Low"))
aaron_new$level <- as.factor(aaron_new$level)
aaron_new$trtside <- ifelse(aaron_new$trtside == "east", "East", "West")
aaron_new$trtside <- as.factor(aaron_new$trtside)
aaron_new$hour <- as.factor(aaron_new$hour)
aaron_new$species <- as.factor(aaron_new$species)
aaron_new$TagCodeTrial <- as.factor(aaron_new$TagCodeTrial)

#Separate by species and period for graphics
aaron_bhc <- aaron_new %>% filter(species == "BHC")
aaron_bhc$species <- factor(aaron_bhc$species)
aaron_bhc_a <- aaron_bhc %>% filter(period == "Acclimation")
aaron_bhc_a$period <- factor(aaron_bhc_a$period)
aaron_bhc_t <- aaron_bhc %>% filter(period == "Treatment")
aaron_bhc_t$period <- factor(aaron_bhc_t$period)

aaron_grc <- aaron_new %>% filter(species == "GRC")
aaron_grc$species <- factor(aaron_grc$species)
aaron_grc_a <- aaron_grc %>% filter(period == "Acclimation")
aaron_grc_a$period <- factor(aaron_grc_a$period)
aaron_grc_t <- aaron_grc %>% filter(period == "Treatment")
aaron_grc_t$period <- factor(aaron_grc_t$period)


#BHC
#summarise by period
bhc_a_dis <- aaron_bhc_a %>%
  group_by(trial, level, trtside, TagCodeTrial) %>%
  dplyr::summarise(acc_dist = mean(avg_dist), .groups = "keep")

bhc_t_dis <- aaron_bhc_t %>%
  group_by(trial, level, trtside, TagCodeTrial) %>%
  dplyr::summarise(treat_dist = mean(avg_dist), .groups = "keep")

#Join period together for one data frame
bhc_dis <- left_join(bhc_a_dis, bhc_t_dis,
                     by = c("trial", "level", "trtside", "TagCodeTrial"))
bhc_dis <- data.frame(bhc_dis)
bhc_dis$TagCodeTrial <- factor(bhc_dis$TagCodeTrial) #43 total fish

#Pivot to long structure
bhc_long <- bhc_dis %>% pivot_longer(c(acc_dist, treat_dist),
                                     names_to = "period", values_to = "dist")
bhc_long <- data.frame(bhc_long)
bhc_long$period <- ifelse(bhc_long$period == "acc_dist",
                          "Acclimation", "Treatment")
bhc_long$period <- factor(bhc_long$period)
bhc_long$level <- factor(bhc_long$level, levels = c("Low", "Medium", "High"))


#GRC
#summarise by period
grc_a_dis <- aaron_grc_a %>%
  group_by(trial, level, trtside, TagCodeTrial) %>%
  dplyr::summarise(acc_dist = mean(avg_dist), .groups = "keep")

grc_t_dis <- aaron_grc_t %>%
  group_by(trial, level, trtside, TagCodeTrial) %>%
  dplyr::summarise(treat_dist = mean(avg_dist), .groups = "keep")

#Join period together for one data frame
grc_dis <- left_join(grc_a_dis, grc_t_dis,
                     by = c("trial", "level", "trtside", "TagCodeTrial"))
grc_dis <- data.frame(grc_dis)
grc_dis$TagCodeTrial <- factor(grc_dis$TagCodeTrial) #43 total fish

#Pivot to long structure
grc_long <- grc_dis %>% pivot_longer(c(acc_dist, treat_dist),
                                     names_to = "period", values_to = "dist")
grc_long <- data.frame(grc_long)
grc_long$period <- ifelse(grc_long$period == "acc_dist",
                          "Acclimation", "Treatment")
grc_long$period <- factor(grc_long$period)
grc_long$level <- factor(grc_long$level, levels = c("Low", "Medium", "High"))


## @knitr original_exploratory_graphics

ggplot(aes(x = trial, y = distance, fill = level), data = aaron) +
  geom_boxplot()
ggplot(aes(x = period, y = distance, fill = level), data = aaron) +
  geom_violin()
ggplot(aes(x = trial, y = distance, color = period, shape = trtside),
       data = aaron) + geom_jitter()
ggplot(aes(x = distance, fill = species), data = aaron) +
  geom_density(alpha = 0.5) #Almost identical
ggplot(aes(x = distance, fill = period), data = aaron) +
  geom_density(alpha = 0.5)
ggplot(aes(x = trial, y = distance, fill = species), data = aaron) +
  geom_boxplot() #Almost identical



## @knitr original_distance_species

ggplot(aes(x = trial, y = distance, fill = period), data = aaron) +
  geom_boxplot() +
  facet_grid(~ species,
             labeller = labeller(species = c("BHC" = "Big Head Carp",
                                             "GRC" = "Grass Carp"))) +
  labs(fill = "Period") +
  xlab("Trial") +
  ylab("Distance") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))



## @knitr original_density

aaron$level <- factor(aaron$level, levels = c("Low", "Medium", "High"))

ggplot(aes(x = distance, fill = level), data = aaron) +
  geom_density(alpha = 0.5) +
  facet_grid(trtside ~ period) +
  labs(fill = "Level") +
  xlab("Distance") +
  ylab("Density") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_fill_colorblind()



## @knitr habitat_density

#adehabitat package data
aaron_new$level <- factor(aaron_new$level, levels = c("Low", "Medium", "High"))

ggplot(aes(x = avg_dist, fill = level), data = aaron_new) +
  geom_density(alpha = 0.5) +
  facet_grid(trtside ~ period) +
  labs(fill = "Level") +
  xlab("Absolute Distance") +
  ylab("Density") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_fill_colorblind()



## @knitr habitat_average

#adehabitat package data
ggplot(aes(x = avg_dx, y = avg_dy, color = species), data = aaron_new) +
  geom_point() +
  facet_grid(trtside ~ period) +
  xlab("Average X Distance") +
  ylab("Average Y Distance") +
  labs(color = "Species") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_color_colorblind()



## @knitr total_crossing

ggplot(aes(x = total_cross, fill = period), data = aaron_new) +
  geom_histogram(bins = 5, position = "dodge") +
  facet_wrap(~ species,
             labeller = labeller(species = c("BHC" = "Big Head Carp",
                                             "GRC" = "Grass Carp"))) +
  xlab("Total Crossings") +
  ylab("Count") +
  labs(fill = "Period") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))



## @knitr bhc_exploratory

aaron_bhc$level <- factor(aaron_bhc$level,
                          levels = c("Low", "Medium", "High"))
aaron_bhc_a$level <- factor(aaron_bhc_a$level,
                            levels = c("Low", "Medium", "High"))
aaron_bhc_t$level <- factor(aaron_bhc_t$level,
                            levels = c("Low", "Medium", "High"))

ggplot(aes(x = avg_dist, fill = period), data = aaron_bhc) +
  geom_density(alpha = 0.5)
ggplot(aes(x = avg_dist, fill = level), data = aaron_bhc) +
  geom_density(alpha = 0.5)
ggplot(aes(x = avg_dist, fill = level), data = aaron_bhc) +
  geom_density(alpha = 0.5) + facet_grid(trtside~period)
ggplot(aes(x = total_cross), data = aaron_bhc) +
  geom_histogram(bins = 20) + facet_wrap(~period)

ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_bhc) +
  geom_point() + facet_wrap(~trtside)
#By period
ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_bhc_a) +
  geom_point() + facet_wrap(~trtside)
ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_bhc_t) +
  geom_point() + facet_wrap(~trtside)

ggplot(aes(x = avg_dist, fill = level), data = aaron_bhc) +
  geom_density(alpha = 0.5) + facet_grid(trtside ~ period) +
  labs(fill = "Level") + xlab("Average Distance") + ylab("Density") +
  theme_bw() + theme(strip.background = element_blank()) +
  scale_fill_colorblind()
ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_bhc) +
  geom_point() + facet_grid(trtside ~ period) +
  labs(color = "Level") + xlab("Hour") + ylab("Average Distance") +
  theme_bw() + theme(strip.background = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_color_colorblind()
ggplot(aes(x = total_cross, fill = level), data = aaron_bhc) +
  geom_histogram(bins = 20) + facet_grid(trtside ~ period) +
  labs(fill = "Level") + xlab("Total Crosses") + ylab("Count") +
  theme_bw() + theme(strip.background = element_blank()) +
  scale_fill_colorblind()



## @knitr grc_exploratory

aaron_grc$level <- factor(aaron_grc$level,
                          levels = c("Low", "Medium", "High"))
aaron_grc_a$level <- factor(aaron_grc_a$level,
                            levels = c("Low", "Medium", "High"))
aaron_grc_t$level <- factor(aaron_grc_t$level,
                            levels = c("Low", "Medium", "High"))

ggplot(aes(x = avg_dist, fill = period), data = aaron_grc) +
  geom_density(alpha = 0.5)
ggplot(aes(x = avg_dist, fill = level), data = aaron_grc) +
  geom_density(alpha = 0.5)
ggplot(aes(x = avg_dist, fill = level), data = aaron_grc) +
  geom_density(alpha = 0.5) + facet_grid(trtside~period)
ggplot(aes(x = total_cross), data = aaron_grc) +
  geom_histogram(bins = 20) + facet_wrap(trtside~period)

ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_grc) +
  geom_point() + facet_wrap(~trtside)
#By period
ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_grc_a) +
  geom_point() + facet_wrap(~trtside)
ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_grc_t) +
  geom_point() + facet_wrap(~trtside) #MUCH more movement during treatment

ggplot(aes(x = avg_dist, fill = level), data = aaron_grc) +
  geom_density(alpha = 0.5) + facet_grid(trtside ~ period) +
  labs(fill = "Level") + xlab("Average Distance") + ylab("Density") +
  theme_bw() + theme(strip.background = element_blank()) +
  scale_fill_colorblind()
ggplot(aes(x = hour, y = avg_dist, color = level), data = aaron_grc) +
  geom_point() + facet_grid(trtside ~ period) + labs(color = "Level") +
  xlab("Hour") + ylab("Average Distance") + theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  scale_color_colorblind()
ggplot(aes(x = total_cross, fill = level), data = aaron_grc) +
  geom_histogram(bins = 20) + facet_grid(trtside ~ period) +
  labs(fill = "Level") + xlab("Total Crosses") + ylab("Count") +
  theme_bw() + theme(strip.background = element_blank()) +
  scale_fill_colorblind()



## @knitr bhc_avg_table

bhc_table <- aaron_bhc %>%
  group_by(period, trtside, level) %>%
  dplyr::summarise(AverageDist = mean(avg_dist), .groups = "keep")
names(bhc_table) <- c("Period", "Treatment Side", "Level", "Average Distance")
bhc_table[, 4] <- round(bhc_table[, 4], 4)
formattable(bhc_table,
               align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
               list(`Indicator Name` = formatter("span",
               style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr grc_avg_table

grc_table <- aaron_grc %>%
  group_by(period, trtside, level) %>%
  dplyr::summarise(AverageDist = mean(avg_dist), .groups = "keep")
names(grc_table) <- c("Period", "Treatment Side", "Level", "Average Distance")
grc_table[, 4] <- round(grc_table[, 4], 4)
formattable(grc_table,
              align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
              list(`Indicator Name` = formatter("span",
              style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr bhc_pooled

#Pool over random effect
aaron_bhc_pool <- bhc_long %>%
  group_by(period, level, trtside) %>%
  dplyr::summarise(Pool_Dist = mean(dist), .groups = "keep")
aaron_bhc_pooled <- lm(Pool_Dist ~ trtside + level + period + (period * level),
                       data = aaron_bhc_pool)
#plot(aaron_bhc_pooled)
bhc_pool_sum <- summary(aaron_bhc_pooled)$coefficients
bhc_pool_ci <- confint(aaron_bhc_pooled)
bhc_pool_table <- cbind(bhc_pool_sum[, 1:2], bhc_pool_ci, bhc_pool_sum[, 3:4])
bhc_pool_table <- data.frame(bhc_pool_table)
bhc_pool_table <- round(bhc_pool_table, 4)
colnames(bhc_pool_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                          "Upper 95% CI", "t-Value", "p-Value")
formattable(bhc_pool_table,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr grc_pooled

#Pool over random effect
aaron_grc_pool <- grc_long %>%
  group_by(period, level, trtside) %>%
  dplyr::summarise(Pool_Dist = mean(dist), .groups = "keep")
aaron_grc_pooled <- lm(Pool_Dist ~ trtside + level + period + (period * level),
                       data = aaron_grc_pool)
#plot(aaron_grc_pooled)
grc_pool_sum <- summary(aaron_grc_pooled)$coefficients
grc_pool_ci <- confint(aaron_grc_pooled)
grc_pool_table <- cbind(grc_pool_sum[, 1:2], grc_pool_ci, grc_pool_sum[, 3:4])
grc_pool_table <- data.frame(grc_pool_table)
grc_pool_table <- round(grc_pool_table, 4)
colnames(grc_pool_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                              "Upper 95% CI", "t-Value", "p-Value")
formattable(grc_pool_table,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr bhc_fixed

aaron_bhc_fixed <- lm(dist ~ level + period + (period * level) + trial / level,
                      data = bhc_long)
#plot(aaron_bhc_fixed)
bhc_fixed_sum <- summary(aaron_bhc_fixed)$coefficients
bhc_fixed_ci <- confint(aaron_bhc_fixed)
bhc_fixed_table <- cbind(bhc_fixed_sum[c(1:4, 11:12), 1:2],
                         bhc_fixed_ci[c(1:4, 13:14), ],
                         bhc_fixed_sum[c(1:4, 11:12), 3:4])
bhc_fixed_table <- data.frame(bhc_fixed_table)
bhc_fixed_table <- round(bhc_fixed_table, 4)
colnames(bhc_fixed_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                               "Upper 95% CI", "t-Value", "p-Value")
bhc_fixed_table2 <- bhc_fixed_table
bhc_fixed_table2$`p-Value` <- as.character(bhc_fixed_table2$`p-Value`)
bhc_fixed_table2[c(1:3), 6] <- "<0.0001"
bhc_fixed_table2[6, 6] <- "0.0001"
formattable(bhc_fixed_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr grc_fixed

aaron_grc_fixed <- lm(dist ~ level + period + (period * level) + trial / level,
                      data = grc_long)
#plot(aaron_grc_fixed)
grc_fixed_sum <- summary(aaron_grc_fixed)$coefficients
grc_fixed_ci <- confint(aaron_grc_fixed)
grc_fixed_table <- cbind(grc_fixed_sum[c(1:4, 11:12), 1:2],
                         grc_fixed_ci[c(1:4, 13:14), ],
                         grc_fixed_sum[c(1:4, 11:12), 3:4])
grc_fixed_table <- data.frame(grc_fixed_table)
grc_fixed_table <- round(grc_fixed_table, 4)
colnames(grc_fixed_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                               "Upper 95% CI", "t-Value", "p-Value")
grc_fixed_table2 <- grc_fixed_table
grc_fixed_table2$`p-Value` <- as.character(grc_fixed_table2$`p-Value`)
grc_fixed_table2[c(1:3), 6] <- "<0.0001"
grc_fixed_table2[5, 6] <- "0.0002"
formattable(grc_fixed_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))


## @knitr bhc_ml

aaron_bhc_unrest <- lme(dist ~ trtside + level + period + (period * level),
                        random = ~ 1 | trial / TagCodeTrial,
                        data = bhc_long, method = "ML")
bhc_unrest_sum <- summary(aaron_bhc_unrest)$tTable
bhc_unrest_ci <- intervals(aaron_bhc_unrest, which = "fixed")
bhc_unrest_ci <- data.frame(bhc_unrest_ci$fixed)
bhc_unrest_table <- cbind(bhc_unrest_sum[, 1:2], bhc_unrest_ci[, -2],
                          bhc_unrest_sum[, 4:5])
bhc_unrest_table <- data.frame(bhc_unrest_table)
bhc_unrest_table <- round(bhc_unrest_table, 4)
colnames(bhc_unrest_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                          "Upper 95% CI", "t-Value", "p-Value")
bhc_unrest_table2 <- bhc_unrest_table
bhc_unrest_table2$`p-Value` <- as.character(bhc_unrest_table2$`p-Value`)
bhc_unrest_table2[1, 6] <- "<0.0001"
bhc_unrest_table2[7, 6] <- "0.0002"
formattable(bhc_unrest_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr grc_ml

aaron_grc_unrest <- lme(dist ~ trtside + level + period + (period * level),
                        random = ~ 1 | trial / TagCodeTrial,
                        data = grc_long, method = "ML")
grc_unrest_sum <- summary(aaron_grc_unrest)$tTable
grc_unrest_ci <- intervals(aaron_grc_unrest, which = "fixed")
grc_unrest_ci <- data.frame(grc_unrest_ci$fixed)
grc_unrest_table <- cbind(grc_unrest_sum[, 1:2], grc_unrest_ci[, -2],
                          grc_unrest_sum[, 4:5])
grc_unrest_table <- data.frame(grc_unrest_table)
grc_unrest_table <- round(grc_unrest_table, 4)
colnames(grc_unrest_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                                "Upper 95% CI", "t-Value", "p-Value")
grc_unrest_table2 <- grc_unrest_table
grc_unrest_table2$`p-Value` <- as.character(grc_unrest_table2$`p-Value`)
grc_unrest_table2[1, 6] <- "<0.0001"
grc_unrest_table2[6, 6] <- "0.0005"
formattable(grc_unrest_table2,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr bhc_reml

aaron_bhc_rest <- lme(dist ~ trtside + level + period + (period * level),
                      random = ~ 1 | trial / TagCodeTrial,
                      data = bhc_long, method = "REML")
bhc_rest_sum <- summary(aaron_bhc_rest)$tTable
bhc_rest_ci <- intervals(aaron_bhc_rest, which = "fixed")
bhc_rest_ci <- data.frame(bhc_rest_ci$fixed)
bhc_rest_table <- cbind(bhc_rest_sum[, 1:2], bhc_rest_ci[, -2],
                          bhc_rest_sum[, 4:5])
bhc_rest_table <- data.frame(bhc_rest_table)
bhc_rest_table <- round(bhc_rest_table, 4)
colnames(bhc_rest_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                                "Upper 95% CI", "t-Value", "p-Value")
formattable(bhc_rest_table,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr grc_reml

aaron_grc_rest <- lme(dist ~ trtside + level + period + (period * level),
                      random = ~ 1 | trial / TagCodeTrial,
                      data = grc_long, method = "REML")
grc_rest_sum <- summary(aaron_grc_rest)$tTable
grc_rest_ci <- intervals(aaron_grc_rest, which = "fixed")
grc_rest_ci <- data.frame(grc_rest_ci$fixed)
grc_rest_table <- cbind(grc_rest_sum[, 1:2], grc_rest_ci[, -2],
                          grc_rest_sum[, 4:5])
grc_rest_table <- data.frame(grc_rest_table)
grc_rest_table <- round(grc_rest_table, 4)
colnames(grc_rest_table) <- c("Estimate", "Std. Error", "Lower 95% CI",
                                "Upper 95% CI", "t-Value", "p-Value")
formattable(grc_rest_table,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))


## @knitr bhc_anova

bhc_pooled_anova <- data.frame(anova(aaron_bhc_pooled))
bhc_pooled_den_df <- bhc_pooled_anova[5, 1]
bhc_pooled_df <- cbind(bhc_pooled_anova[4, 1], bhc_pooled_den_df,
                   bhc_pooled_anova[4, 4:5])
bhc_pooled_df <- data.frame(bhc_pooled_df)
names(bhc_pooled_df) <- c("Numerator DF", "Denominator DF",
                          "F-Value", "p-Value")
rownames(bhc_pooled_df) <- "Pooled Level:Period"

bhc_fixed_anova <- data.frame(anova(aaron_bhc_fixed))
bhc_fixed_den_df <- bhc_fixed_anova[5, 1]
bhc_fixed_df <- cbind(bhc_fixed_anova[4, 1], bhc_fixed_den_df,
                       bhc_fixed_anova[4, 4:5])
bhc_fixed_df <- data.frame(bhc_fixed_df)
names(bhc_fixed_df) <- c("Numerator DF", "Denominator DF",
                         "F-Value", "p-Value")
rownames(bhc_fixed_df) <- "Fixed Level:Period"

bhc_unrest_anova <- data.frame(anova(aaron_bhc_unrest))
names(bhc_unrest_anova) <- c("Numerator DF", "Denominator DF",
                               "F-Value", "p-Value")
bhc_unrest_df <- bhc_unrest_anova[5, ]
rownames(bhc_unrest_df) <- "ML Level:Period"

bhc_rest_anova <- data.frame(anova(aaron_bhc_rest))
names(bhc_rest_anova) <- c("Numerator DF", "Denominator DF",
                             "F-Value", "p-Value")
bhc_rest_df <- bhc_rest_anova[5, ]
rownames(bhc_rest_df) <- "REML Level:Period"


bhc_anova <- data.frame(rbind(bhc_pooled_df, bhc_fixed_df,
                              bhc_unrest_df, bhc_rest_df))
bhc_anova <- round(bhc_anova, 4)
names(bhc_anova) <- c("Numerator DF", "Denominator DF",
                        "F-Value", "p-Value")
bhc_anova$`p-Value` <- as.character(bhc_anova$`p-Value`)
bhc_anova[-1, 4] <- "<0.0001"
formattable(bhc_anova,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr grc_anova

grc_pooled_anova <- data.frame(anova(aaron_grc_pooled))
grc_pooled_den_df <- grc_pooled_anova[5, 1]
grc_pooled_df <- cbind(grc_pooled_anova[4, 1], grc_pooled_den_df,
                       grc_pooled_anova[4, 4:5])
grc_pooled_df <- data.frame(grc_pooled_df)
names(grc_pooled_df) <- c("Numerator DF", "Denominator DF",
                          "F-Value", "p-Value")
rownames(grc_pooled_df) <- "Pooled Level:Period"

grc_fixed_anova <- data.frame(anova(aaron_grc_fixed))
grc_fixed_den_df <- grc_fixed_anova[5, 1]
grc_fixed_df <- cbind(grc_fixed_anova[4, 1], grc_fixed_den_df,
                      grc_fixed_anova[4, 4:5])
grc_fixed_df <- data.frame(grc_fixed_df)
names(grc_fixed_df) <- c("Numerator DF", "Denominator DF",
                         "F-Value", "p-Value")
rownames(grc_fixed_df) <- "Fixed Level:Period"

grc_unrest_anova <- data.frame(anova(aaron_grc_unrest))
names(grc_unrest_anova) <- c("Numerator DF", "Denominator DF",
                             "F-Value", "p-Value")
grc_unrest_df <- grc_unrest_anova[5, ]
rownames(grc_unrest_df) <- "ML Level:Period"

grc_rest_anova <- data.frame(anova(aaron_grc_rest))
names(grc_rest_anova) <- c("Numerator DF", "Denominator DF",
                           "F-Value", "p-Value")
grc_rest_df <- grc_rest_anova[5, ]
rownames(grc_rest_df) <- "REML Level:Period"


grc_anova <- data.frame(rbind(grc_pooled_df, grc_fixed_df,
                              grc_unrest_df, grc_rest_df))
grc_anova <- round(grc_anova, 4)
names(grc_anova) <- c("Numerator DF", "Denominator DF",
                      "F-Value", "p-Value")
grc_anova$`p-Value` <- as.character(grc_anova$`p-Value`)
grc_anova[-1, 4] <- "<0.0001"
formattable(grc_anova,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))
