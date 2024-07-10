#===============================================================================
#  CODE CHUNKS AND FIGURES FOR PAPER
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

#-- Install devtools if it is not already installed

install.packages("devtools")

#-- Install the IPD package from GitHub

devtools::install_github("awanafiaz/IPD")

#-- Load the IPD library

library(IPD)

#-- Load additional libraries

library(tidyverse)

library(patchwork)

#--- DATA GENERATION -----------------------------------------------------------

#-- Generate example data for linear regression

set.seed(123)

n <- c(10000, 500, 1000)

dat_ols <- simdat(n = n, effect = 3, sigma_Y = 1,

  model = "ols", shift = 1, scale = 2)

#-- Plot example labeled data

dat_ols_labeled <- dat_ols[dat_ols$set == "labeled",]

p1 <- dat_ols_labeled |>

  pivot_longer(Y:f) |>

  mutate(

    name = factor(name, levels = c("f", "Y"),

      labels = c("Predicted ", "Observed ")) |> fct_rev()) |>

  ggplot(aes(x = X1, y = value, group = name, color = name, fill = name)) +

  theme_bw() + coord_fixed(1/3) +

  geom_point(alpha = 0.5) +

  geom_smooth(method = "lm") +

  scale_x_continuous(limits = c(-2.5, 2.5)) +

  scale_y_continuous(limits = c(-7.5, 7.5)) +

  scale_color_manual(values = c("#FFB500", "#00C1D5")) +

  scale_fill_manual(values = c("#FFB500", "#00C1D5")) +

  labs(x = "\nCovariate", y = "Outcome\n", fill = NULL, color = NULL) +

  theme(

    legend.direction = "horizontal",

    legend.position = "inside",

    legend.position.inside = c(0.3, 0.91),

    axis.title = element_text(size = 14),

    axis.text = element_text(size = 12))

p2 <- dat_ols_labeled |> ggplot(aes(x = f, y = Y)) +

  theme_bw() + coord_fixed(1) +

  geom_point(color = "#1B365D", alpha = 0.5) +

  geom_smooth(method = "lm", color = "#1B365D") +

  scale_x_continuous(limits = c(-7.5, 7.5)) +

  scale_y_continuous(limits = c(-7.5, 7.5)) +

  labs(x = "\nPredicted Outcome", y = "Observed Outcome\n") +

  theme(

    axis.title = element_text(size = 14),

    axis.text = element_text(size = 12))

fig1 <- (p1 + theme(plot.margin = unit(c(0,20,0, 0), "pt"))) +

  (p2 + theme(plot.margin = unit(c(0,20,0,20), "pt"))) +

  plot_annotation(tag_levels = "A")

fig1

ggsave("fig1.png", fig1, width = 10, height = 6, units = "in")

#=== ANALYSIS ==================================================================

#--- MODEL FITTING -------------------------------------------------------------

dat_ols_l <- dat_ols[dat_ols$set == "labeled",]

dat_ols_u <- dat_ols[dat_ols$set == "unlabeled",]

#-- Oracle Regression

fit0 <- lm(Y ~ X1, data = dat_ols_u)

#-- Naive Regression

fit1 <- lm(f ~ X1, data = dat_ols_u)

#-- Classic Regression

fit2 <- lm(Y ~ X1, data = dat_ols_l)

#-- Specify the Formula for IPD

formula <- Y - f ~ X1

#-- PostPI Bootstrap Correction

fit3 <- ipd(formula, method = "postpi_boot", model = "ols",

  data = dat_ols, label = "set", nboot = 200)

#-- PostPI Analytic Correction

fit4 <- ipd(formula, method = "postpi_analytic", model = "ols",

  data = dat_ols, label = "set")

#-- PPI

fit5 <- ipd(formula,  method = "ppi", model = "ols",

  data = dat_ols, label = "set")

#-- PPI++

fit6 <- ipd(formula, method = "ppi_plusplus", model = "ols",

  data = dat_ols, label = "set")

#-- POP-Inf

fit7 <- ipd(formula, method = "popinf", model = "ols",

  data = dat_ols, label = "set")

#--- SUMMARIZE RESULTS ---------------------------------------------------------

mthds <- c("Oracle", "Naive", "Classic",

  "PostPI Bootstrap", "PostPI Analytic", "PPI", "PPI++", "POP-Inf")

fig2_dat <- cbind(mthds,

  rbind(

    c(summary(fit0)$coefficients[2, 1:2], confint(fit0)[2,]),

    c(summary(fit1)$coefficients[2, 1:2], confint(fit1)[2,]),

    c(summary(fit2)$coefficients[2, 1:2], confint(fit2)[2,]),

    summary(fit3)$coefficients[2,],

    summary(fit4)$coefficients[2,],

    summary(fit5)$coefficients[2,],

    summary(fit6)$coefficients[2,],

    summary(fit7)$coefficients[2,])) |>

  as_tibble() |>

  mutate(

    mthds = factor(mthds, levels = mthds, labels = mthds),

    cov = if_else(Lower.CI <= 3 & Upper.CI >= 3, 1, 0) |>

      factor(levels = 0:1, labels = c("No", "Yes")))

fig2 <- fig2_dat |>

  ggplot(aes(x = mthds, y = Estimate, ymin = Lower.CI, ymax = Upper.CI,

    color = cov, fill = cov)) +

    theme_minimal() + coord_flip() +

    geom_hline(yintercept = 3, linetype = 2) +

    geom_pointrange(linewidth = 3) +

    geom_point(color = "#1B365D", size = 3) +

    scale_x_discrete(limits = rev) +

    scale_color_manual(values = c("#FFB500", "#00C1D5")) +

    labs(x = NULL, y = "\nEstimated Coefficient",

      color = "Covers the True Coefficient?",

      fill = "Covers the True Coefficient?") +

    theme(

      legend.direction = "horizontal",

      legend.position = "inside",

      legend.position.inside = c(0.21, 0.125),

      axis.title = element_text(size = 14),

      axis.text = element_text(size = 12))

ggsave("fig2.png", fig2, width = 10, height = 4, units = "in")

#--- PRINTING, SUMMARIZING, AND TIDYING ----------------------------------------

df <- dat_ols[which(dat_ols$set != "training"),]

print(fit3)                                      #- Print

summ <- summary(fit3)                            #- Summary

print(summ)                                      #- Print Summary

tidy(fit3)                                       #- Tidy

glance(fit3)                                     #- Glance

augmented_df <- augment(fit3, data = df)         #- Augment

#=== END =======================================================================
