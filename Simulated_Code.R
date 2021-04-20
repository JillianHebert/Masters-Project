#Simulated Data

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



## @knitr simulation

sim_data <- function(alpha = 0.05, btwn_sub, wthn_sub,
                     n_sub, n_rep = 20, delta) {

  alpha <- alpha #Significance level
  rep <- 1000
  final <- data.frame(pool_p = double(),
                      fixed_p = double(),
                      un_p = double(),
                      rest_p = double())

  final_ci_table <- data.frame(pool_lower_ci = double(),
                               pool_upper_ci = double(),
                               fixed_lower_ci = double(),
                               fixed_upper_ci = double(),
                               un_lower_ci = double(),
                               un_upper_ci = double(),
                               rest_lower_ci = double(),
                               rest_upper_ci = double())

  for(j in 1:rep) {

    eps <- 7 #location parameter (-inf, inf)
    lam <- 1.5 #scale parameter (must be positive)
    eta <- 16 #shape parameter (must be positive)

    gam_mean <- eps + (lam * eta)
    gam_sd <- lam * sqrt(eta)
    skew <- 2 / sqrt(eta) #always positive

    n_trt <- 2 # number of treatments
    n_sub <- n_sub # number of subjects per treatment** (trials)
    n_rep <- n_rep #number of replicate measurements per subject
    n <- n_trt * n_sub * n_rep # grand total

    trt_means <- c(10, 10 + (delta * gam_sd)) #treatment difference
    sig_sub <- btwn_sub # between subject st dev **
    sig_e <- wthn_sub # error st dev (within subject st dev) **

    # Create data frame columns that specify trt, sub, and rep id
    rep_k <- as.factor(1:(n_trt * n_sub * n_rep))
    sub_j <- as.factor(rep(1:(n_sub * n_trt), each = n_rep))
    trt_i <- as.factor(rep(1:n_trt, each = n_sub * n_rep))

    #Generate response components, errors and treatment effects
    e_k <- scale(rgamma(n = n_trt * n_sub * n_rep, shape = eta) * sig_e)
    e_sub_j <- rep(scale(rgamma(n = n_sub * n_trt, shape = eta)) * sig_sub,
                   each = n_rep)
    mu_trt_i <- rep(trt_means, each = n_sub * n_rep)

    #Construct the response variable by summing components, errors and trt
    sim.df <- cbind.data.frame(trt_i, sub_j, rep_k, mu_trt_i, e_sub_j, e_k)
    sim.df <- sim.df %>%
      mutate(Y = mu_trt_i + e_sub_j + e_k)

    #Create a version of the df that averages across reps within subjects
    sim.df.group <- sim.df %>%
      group_by(sub_j, trt_i) %>%
      dplyr::summarise(Y_pooled = mean(Y), .groups = "keep")

    #Fit the pooled model and store the associated treatment p-value
    pool_mod <- lm(Y_pooled ~ trt_i, data = sim.df.group)
    pool_sum <- summary(pool_mod)
    pool_p <- pool_sum$coefficients[2, 4]
    pool_lower_ci <- as.numeric(confint(pool_mod)[2, ][1])
    pool_upper_ci <- as.numeric(confint(pool_mod)[2, ][2])

    #Fit the fixed model and store the associate treatment p-value
    #(subject as block)
    fixed_mod <- lm(Y ~ trt_i + sub_j, data = sim.df)
    fixed_sum <- summary(fixed_mod)
    fixed_p <- fixed_sum$coefficients[2, 4]
    fixed_lower_ci <- as.numeric(confint(fixed_mod)[2, ][1])
    fixed_upper_ci <- as.numeric(confint(fixed_mod)[2, ][2])

    #Fit the mixed effects model with ML and store p-value
    un_mod <- lme(Y ~ trt_i, random = list(~1 | sub_j),
                  data = sim.df, method = "ML")
    un_sum <- summary(un_mod)
    un_p <- un_sum$tTable[2, 5]
    un_intervals <- intervals(un_mod)
    un_lower_ci <- as.numeric(un_intervals$fixed[2, ][1])
    un_upper_ci <- as.numeric(un_intervals$fixed[2, ][3])

    #Fit the mixed effects model with REML and store p-value
    rest_mod <- lme(Y ~ trt_i, random = list(~ 1 | sub_j),
                    data = sim.df, method = "REML")
    rest_sum <- summary(rest_mod)
    rest_p <- rest_sum$tTable[2, 5]
    rest_intervals <- intervals(rest_mod)
    rest_lower_ci <- as.numeric(rest_intervals$fixed[2, ][1])
    rest_upper_ci <- as.numeric(rest_intervals$fixed[2, ][3])

    #Assemble p-values and CI from all models and bind them with previous sims
    table <- as.data.frame(cbind(pool_p, fixed_p, un_p, rest_p))
    final <- rbind(final, table)

    #Assemble CI from all models and bind them with previous sims
    ci_table <- as.data.frame(cbind(pool_lower_ci, pool_upper_ci,
                                    fixed_lower_ci, fixed_upper_ci,
                                    un_lower_ci, un_upper_ci, rest_lower_ci,
                                    rest_upper_ci))
    final_ci_table <- rbind(final_ci_table, ci_table)

  } #End of simulation for loop


  #Create power table
  power_table <- final %>%
    dplyr::summarise(pool_power = sum(pool_p < 0.05) / rep,
                     fixed_power = sum(fixed_p < 0.05) / rep,
                     un_power = sum(un_p < 0.05) / rep,
                     rest_power = sum(rest_p < 0.05) / rep)

  power_table <- as.data.frame(power_table)


  #Create CI table
  average_ci_table <- final_ci_table %>%
    dplyr::summarise(avg_pool_lower = mean(pool_lower_ci),
                     avg_pool_upper = mean(pool_upper_ci),
                     avg_fixed_lower = mean(fixed_lower_ci),
                     avg_fixed_upper = mean(fixed_upper_ci),
                     avg_un_lower = mean(un_lower_ci),
                     avg_un_upper = mean(un_upper_ci),
                     avg_rest_lower = mean(rest_lower_ci),
                     avg_rest_upper = mean(rest_upper_ci))

  average_ci_table <- as.data.frame(average_ci_table)

  #Combine tables for output
  output <- cbind(power_table, average_ci_table)
  output <- as.data.frame(output)
  output

}



## @knitr no_effect_power

cases <- read.csv("cases.csv")
n_cases <- dim(cases)[1]
results <- as.data.frame(cbind(case = seq(1:n_cases), pool_power = NA,
                               fixed_power = NA, un_power = NA,
                               rest_power = NA, avg_pool_lower = NA,
                               avg_pool_upper = NA, avg_fixed_lower = NA,
                               avg_fixed_upper = NA, avg_un_lower = NA,
                               avg_un_upper = NA, avg_rest_lower = NA,
                               avg_rest_upper = NA))

#Loop through all 9 cases and save results
for(i in 1:n_cases) {
  args <- as.numeric(cases[i, 1:5])
  results[i, 2:13] <- sim_data(alpha = args[1], btwn_sub = args[2],
                               wthn_sub = args[3], n_sub = args[4],
                               n_rep = args[5], delta = 0)
}

#Format output
equal_sim_results <- cbind(seq(1:9), cases[, 1:5], results[, -1])
equal_sim_results <- as.data.frame(equal_sim_results)
names(equal_sim_results) <- c("Case", "Alpha", "Between Subject SD",
                              "Within Subject SD", "Number of Subjects",
                              "Number of Repetitions", "Pooled Power",
                              "Fixed Power", "ML Power", "REML Power",
                              "Pooled Lower CI", "Pooled Upper CI",
                              "Fixed Lower CI", "Fixed Upper CI",
                              "ML Lower CI", "ML Upper CI", "REML Lower CI",
                              "REML Upper CI")

#Create power table
equal_power <- equal_sim_results[, 1:10]
equal_power <- round(equal_power, 4)
formattable(equal_power,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"),
            list(`Indicator Name` = formatter(
              "span", style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr no_effect_graphics

equal_data <- equal_sim_results %>%
  select(`Pooled Upper CI`, `Pooled Lower CI`, `Fixed Upper CI`,
         `Fixed Lower CI`, `ML Upper CI`, `ML Lower CI`, `REML Upper CI`,
         `REML Lower CI`) %>%
  dplyr::mutate(pool_estimate = (`Pooled Upper CI` + `Pooled Lower CI`) / 2,
                fixed_estimate = (`Fixed Upper CI` + `Fixed Lower CI`) / 2,
                ml_estimate = (`ML Upper CI` + `ML Lower CI`) / 2,
                reml_estimate = (`REML Upper CI` + `REML Lower CI`) / 2)
equal_data <- cbind(Case = 1:9, equal_data)
equal_data <- as.data.frame(equal_data)

#Separate by model
pooled_data <- equal_data %>%
  dplyr::summarise(term = "Pooled",
                   pool_est = pool_estimate,
                   pool_low = `Pooled Lower CI`,
                   pool_up = `Pooled Upper CI`,
                   .groups = "keep") %>%
  t()

fixed_data <- equal_data %>%
  dplyr::summarise(term = "Fixed",
                   fixed_est = fixed_estimate,
                   fixed_low = `Fixed Lower CI`,
                   fixed_up = `Fixed Upper CI`,
                   .groups = "keep") %>%
  t()

ml_data <- equal_data %>%
  dplyr::summarise(term = "ML",
                   ml_est = ml_estimate,
                   ml_low = `ML Lower CI`,
                   ml_up = `ML Upper CI`,
                   .groups = "keep") %>%
  t()

reml_data <- equal_data %>%
  dplyr::summarise(term = "REML",
                   reml_est = reml_estimate,
                   reml_low = `REML Lower CI`,
                   reml_up = `REML Upper CI`,
                   .groups = "keep") %>%
  t()

#Create CI graph for each model and each case

plot_list <- list()
for(i in 1:9) {

  case <- as.data.frame(rbind(pooled_data[, i], fixed_data[, i],
                              ml_data[, i], reml_data[, i]))
  names(case) <- c("Model", "Estimate", "Lower CI", "Upper CI")
  case$Estimate <- as.numeric(case$Estimate)
  case$`Lower CI` <- as.numeric(case$`Lower CI`)
  case$`Upper CI` <- as.numeric(case$`Upper CI`)
  case$Model <- ordered(case$Model, levels = c("REML", "ML", "Fixed", "Pooled"))

  graph <- ggplot(aes(x = Model, y = Estimate), data = case) +
    geom_point() +
    scale_y_continuous(breaks = seq(-8, 8, by = 2), limits = c(-8, 8)) +
    labs(title = paste("Case ", i)) +
    geom_linerange(aes(ymin = `Lower CI`, ymax = `Upper CI`)) +
    coord_flip() +
    geom_hline(yintercept = 0, color = "red") +
    theme_bw() +
    theme(strip.background = element_blank())

  plot_list[[i]] <- graph

}

grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]],
             plot_list[[5]], plot_list[[6]], plot_list[[7]], plot_list[[8]],
             plot_list[[9]], ncol = 3, nrow = 3)



## @knitr effect_power

cases <- read.csv("cases.csv")
n_cases <- dim(cases)[1]
results2 <- as.data.frame(cbind(case = seq(1:n_cases), pool_power = NA,
                               fixed_power = NA, un_power = NA,
                               rest_power = NA, avg_pool_lower = NA,
                               avg_pool_upper = NA, avg_fixed_lower = NA,
                               avg_fixed_upper = NA, avg_un_lower = NA,
                               avg_un_upper = NA, avg_rest_lower = NA,
                               avg_rest_upper = NA))

#Loop through all 9 cases
for(i in 1:n_cases) {
  args <- as.numeric(cases[i, 1:5])
  results2[i, 2:13] <- sim_data(alpha = args[1], btwn_sub = args[2],
                               wthn_sub = args[3], n_sub = args[4],
                               n_rep = args[5], delta = 0.5)
}

#Format output
unequal_sim_results <- cbind(seq(1:9), cases[, 1:5], results2[, -1])
unequal_sim_results <- as.data.frame(unequal_sim_results)
names(unequal_sim_results) <- c("Case", "Alpha", "Between Subject SD",
                                "Within Subject SD", "Number of Subjects",
                                "Number of Repetitions", "Pooled Power",
                                "Fixed Power", "ML Power", "REML Power",
                                "Pooled Lower CI", "Pooled Upper CI",
                                "Fixed Lower CI", "Fixed Upper CI",
                                "ML Lower CI", "ML Upper CI", "REML Lower CI",
                                "REML Upper CI")

#Create power table
unequal_power <- unequal_sim_results[, 1:10]
unequal_power <- round(unequal_power, 4)
formattable(unequal_power,
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr effect_graphics

unequal_data <- unequal_sim_results %>%
  select(`Pooled Upper CI`, `Pooled Lower CI`, `Fixed Upper CI`,
         `Fixed Lower CI`, `ML Upper CI`, `ML Lower CI`, `REML Upper CI`,
         `REML Lower CI`) %>%
  dplyr::mutate(pool_estimate = (`Pooled Upper CI` + `Pooled Lower CI`) / 2,
                fixed_estimate = (`Fixed Upper CI` + `Fixed Lower CI`) / 2,
                ml_estimate = (`ML Upper CI` + `ML Lower CI`) / 2,
                reml_estimate = (`REML Upper CI` + `REML Lower CI`) / 2)
unequal_data <- cbind(Case = 1:9, unequal_data)
unequal_data <- as.data.frame(unequal_data)

#Separate by model type
pooled_data2 <- unequal_data %>%
  dplyr::summarise(term = "Pooled",
                   pool_est = pool_estimate,
                   pool_low = `Pooled Lower CI`,
                   pool_up = `Pooled Upper CI`,
                   .groups = "keep") %>%
  t()

fixed_data2 <- unequal_data %>%
  dplyr::summarise(term = "Fixed",
                   fixed_est = fixed_estimate,
                   fixed_low = `Fixed Lower CI`,
                   fixed_up = `Fixed Upper CI`,
                   .groups = "keep") %>%
  t()

ml_data2 <- unequal_data %>%
  dplyr::summarise(term = "ML",
                   ml_est = ml_estimate,
                   ml_low = `ML Lower CI`,
                   ml_up = `ML Upper CI`,
                   .groups = "keep") %>%
  t()

reml_data2 <- unequal_data %>%
  dplyr::summarise(term = "REML",
                   reml_est = reml_estimate,
                   reml_low = `REML Lower CI`,
                   reml_up = `REML Upper CI`,
                   .groups = "keep") %>%
  t()

#Create CI graph for each model and each case
plot_list2 <- list()
for(i in 1:9) {

  case <- as.data.frame(rbind(pooled_data2[, i], fixed_data2[, i],
                              ml_data2[, i], reml_data2[, i]))
  names(case) <- c("Model", "Estimate", "Lower CI", "Upper CI")
  case$Estimate <- as.numeric(case$Estimate)
  case$`Lower CI` <- as.numeric(case$`Lower CI`)
  case$`Upper CI` <- as.numeric(case$`Upper CI`)
  case$Model <- ordered(case$Model, levels = c("REML", "ML", "Fixed", "Pooled"))

#Change center to difference (10, 10 * (0.5 + gam_sd))
  graph <- ggplot(aes(x = Model, y = Estimate), data = case) +
    geom_point() +
    scale_y_continuous(breaks = seq(-4, 10, by = 2), limits = c(-4, 10)) +
    labs(title = paste("Case ", i)) +
    geom_linerange(aes(ymin = `Lower CI`, ymax = `Upper CI`)) +
    coord_flip() +
    geom_hline(yintercept = 3, color = "red") +
    theme_bw() +
    theme(strip.background = element_blank())

  plot_list2[[i]] <- graph

}

grid.arrange(plot_list2[[1]], plot_list2[[2]], plot_list2[[3]], plot_list2[[4]],
             plot_list2[[5]], plot_list2[[6]], plot_list2[[7]], plot_list2[[8]],
             plot_list2[[9]], ncol = 3, nrow = 3)



## @knitr interval_width

pooled_width <- data.frame(t(pooled_data))
pooled_width$pool_low <- as.numeric(pooled_width$pool_low)
pooled_width$pool_up <- as.numeric(pooled_width$pool_up)
pooled_width$width <- pooled_width$pool_up - pooled_width$pool_low

fixed_width <- data.frame(t(fixed_data))
fixed_width$fixed_low <- as.numeric(fixed_width$fixed_low)
fixed_width$fixed_up <- as.numeric(fixed_width$fixed_up)
fixed_width$width <- fixed_width$fixed_up - fixed_width$fixed_low

ml_width <- data.frame(t(ml_data))
ml_width$ml_low <- as.numeric(ml_width$ml_low)
ml_width$ml_up <- as.numeric(ml_width$ml_up)
ml_width$width <- ml_width$ml_up - ml_width$ml_low

reml_width <- data.frame(t(reml_data))
reml_width$reml_low <- as.numeric(reml_width$reml_low)
reml_width$reml_up <- as.numeric(reml_width$reml_up)
reml_width$width <- reml_width$reml_up - reml_width$reml_low

pooled_width2 <- data.frame(t(pooled_data2))
pooled_width2$pool_low <- as.numeric(pooled_width2$pool_low)
pooled_width2$pool_up <- as.numeric(pooled_width2$pool_up)
pooled_width2$width <- pooled_width2$pool_up - pooled_width2$pool_low

fixed_width2 <- data.frame(t(fixed_data2))
fixed_width2$fixed_low <- as.numeric(fixed_width2$fixed_low)
fixed_width2$fixed_up <- as.numeric(fixed_width2$fixed_up)
fixed_width2$width <- fixed_width2$fixed_up - fixed_width2$fixed_low

ml_width2 <- data.frame(t(ml_data2))
ml_width2$ml_low <- as.numeric(ml_width2$ml_low)
ml_width2$ml_up <- as.numeric(ml_width2$ml_up)
ml_width2$width <- ml_width2$ml_up - ml_width2$ml_low

reml_width2 <- data.frame(t(reml_data2))
reml_width2$reml_low <- as.numeric(reml_width2$reml_low)
reml_width2$reml_up <- as.numeric(reml_width2$reml_up)
reml_width2$width <- reml_width2$reml_up - reml_width2$reml_low

final_width2 <- rbind(pooled_width[, c(1, 5)], fixed_width[, c(1, 5)],
                      ml_width[, c(1, 5)], reml_width[, c(1, 5)])
effect_width <- c(pooled_width2[, 5], fixed_width2[, 5],
                      ml_width2[, 5], reml_width2[, 5])
final_width2 <- cbind(final_width2, effect_width)
final_width2$term <- factor(final_width2$term,
                           levels = c("Pooled", "Fixed", "ML", "REML"))
final_width2 <- cbind(seq(1, 9, by = 1), final_width2)
names(final_width2) <- c("Case", "Model", "No Treatment Effect Interval Width",
                         "Treatment Effect Interval Width")
final_width2$Case <- factor(final_width2$Case)
final_width2 <- final_width2[order(final_width2$Case), ]
final_width2[, 3:4] <- round(final_width2[, 3:4], 4)
rownames(final_width2) <- NULL
formattable(final_width2,
            align = c("l", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr avg_widths

no_pool <- mean(pooled_width$width)
no_fixed <- mean(fixed_width$width)
no_ml <- mean(ml_width$width)
no_reml <- mean(reml_width$width)

pool <- mean(pooled_width2$width)
fixed <- mean(fixed_width2$width)
ml <- mean(ml_width2$width)
reml <- mean(reml_width2$width)

mod <- c("Pooled", "Fixed", "ML", "REML")

no_effect <- rbind(no_pool, no_fixed, no_ml, no_reml)
effect <- rbind(pool, fixed, ml, reml)
avg_table <- cbind(mod, no_effect, effect)
rownames(avg_table) <- NULL
avg_table <- data.frame(avg_table)
avg_table$V2 <- round(as.numeric(avg_table$V2), 4)
avg_table$V3 <- round(as.numeric(avg_table$V3), 4)
names(avg_table) <- c("Model", "Average No Effect Width",
                      "Average Effect Width")
formattable(avg_table,
            align = c("l", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))



## @knitr width_percent

data <- cbind(seq(1:9), ml_width$width, reml_width$width,
              ml_width2$width, reml_width2$width)
data <- data.frame(data)
names(data) <- c("Case", "No Effect ML Width", "No Effect REML Width",
                 "Effect ML Width", "Effect REML Width")
data$No_Effect_Difference <- data$`No Effect REML Width` -
  data$`No Effect ML Width`
data$Effect_Difference <- data$`Effect REML Width` -
  data$`Effect ML Width`
no_effect_avg_constant <- (avg_table[3, 2] + avg_table[4, 2]) / 2
data$No_Effect_Percent <- data$No_Effect_Difference /
  no_effect_avg_constant * 100
effect_avg_constant <- (avg_table[3, 3] + avg_table[4, 3]) / 2
data$Effect_Percent <- data$Effect_Difference /
  effect_avg_constant * 100
data[, -1] <- round(data[, -1], 4)
names(data) <- c("Case", "No Effect ML Width", "No Effect REML Width",
                 "Effect ML Width", "Effect REML Width",
                 "No Effect Difference", "Effect Difference",
                 "No Effect Percent", "Effect Percent")
formattable(data[, c(1:5, 8:9)],
            align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "r"),
            list(`Indicator Name` = formatter("span",
            style = ~ style(color = "grey", font.weight = "bold"))))

