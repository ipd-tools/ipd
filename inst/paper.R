#===============================================================================
#  CODE CHUNKS AND FIGURES FOR PAPER
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(ipd)

library(tidyverse)

library(patchwork)

#--- DATA GENERATION -----------------------------------------------------------

#-- Generate example data for linear regression

set.seed(123)

n <- c(10000, 500, 1000)

dat_ols <- simdat(n = n, effect = 1, sigma_Y = 4,

  model = "ols", shift = 0, scale = 1)

#-- Plot example labeled data

dat_ols_labeled <- dat_ols[dat_ols$set_label == "labeled",]

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

  scale_color_manual(values = c("#AA4AC4", "#00C1D5")) +

  scale_fill_manual(values = c("#AA4AC4", "#00C1D5")) +

  labs(x = "\nCovariate", y = "Outcome\n", fill = NULL, color = NULL) +

  theme(

    legend.direction = "horizontal",

    legend.position = "inside",

    legend.position.inside = c(0.3, 0.91),

    axis.title = element_text(size = 14),

    axis.text = element_text(size = 12))

p2 <- dat_ols_labeled |> ggplot(aes(x = f, y = Y)) +

  theme_bw() + coord_fixed(1) +

  geom_point(color = "#00C1D5", alpha = 0.5) +

  geom_smooth(method = "lm", color = "#00C1D5") +

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

#=== ANALYSIS ==================================================================

#--- MODEL FITTING -------------------------------------------------------------

dat_ols_l <- dat_ols[dat_ols$set_label == "labeled",]

dat_ols_u <- dat_ols[dat_ols$set_label == "unlabeled",]

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

  data = dat_ols, label = "set_label", rel_func = "gam", nboot = 200)

#-- PostPI Analytic Correction

fit4 <- ipd(formula, method = "postpi_analytic", model = "ols",

  data = dat_ols, label = "set_label")

#-- PPI

fit5 <- ipd(formula,  method = "ppi", model = "ols",

  data = dat_ols, label = "set_label")

#-- PPI++

fit6 <- ipd(formula, method = "ppi_plusplus", model = "ols",

  data = dat_ols, label = "set_label")

#-- PSPA

fit7 <- ipd(formula, method = "pspa", model = "ols",

  data = dat_ols, label = "set_label")

#--- SUMMARIZE RESULTS ---------------------------------------------------------

mthds <- c("Oracle", "Naive", "Classic",

  "PostPI Bootstrap", "PostPI Analytic", "PPI", "PPI++", "PSPA")

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

    geom_point(color = "#00C1D5", size = 3) +

    scale_x_discrete(limits = rev) +

    scale_color_manual(values = c("#AA4AC4", "#00C1D5")) +

    labs(x = NULL, y = "\nEstimated Coefficient",

      color = "Covers the True Coefficient?",

      fill = "Covers the True Coefficient?") +

    theme(

      legend.direction = "horizontal",

      legend.position = "inside",

      legend.position.inside = c(0.21, 0.125),

      axis.title = element_text(size = 14),

      axis.text = element_text(size = 12))

fig2

#--- PRINTING, SUMMARIZING, AND TIDYING ----------------------------------------

df <- dat_ols[which(dat_ols$set_label != "training"),]

print(fit3)                                      #- Print

summ <- summary(fit3)                            #- Summary

print(summ)                                      #- Print Summary

tidy(fit3)                                       #- Tidy

glance(fit3)                                     #- Glance

augmented_df <- augment(fit3, data = df)         #- Augment

#--- SIMULATION FUNCTION -------------------------------------------------------

sim_func <- function(nsims = 1000, n_t = 100, n_l = 100, n_u = 1000,

  nboot = 100, dgm, true_beta) {

  n_tot <- n_t + n_l + n_u

  est_oracle <- est_naive <- est_classic <- est_postpi <- est_ppi <-

     est_plusplus <- est_pspa <- rep(NA, nsims)

  err_oracle <- err_naive <- err_classic <- err_postpi <- err_ppi <-

     err_plusplus <- err_pspa <- rep(NA, nsims)

  for(sim in 1:nsims) {

    if (sim %% 100 == 0) cat("sim:", sim, "of", nsims, "\n")

    #--- SIMULATE DATA ---------------------------------------------------------

    #-- GENERATE DATASET

    dat <- simdat(n = c(n_t, n_l, n_u), effect = true_beta,

      sigma_Y = 4, model = "ols", shift = 0, scale = 1)

    #-- LABELED DATA

    dat_l = dat |>

      filter(set_label == "labeled")

    #-- UNLABELED DATA

    dat_u = dat |>

      filter(set_label == "unlabeled")

    #--- FIT VARIOUS METHODS ---------------------------------------------------

    #-- 0. ORACLE ESTIMATE

    fit_oracle <- lm(Y ~ X1, data = dat_u)

    est_oracle[sim] <- fit_oracle$coefficients[2]

    err_oracle[sim] <- summary(fit_oracle)$coefficients[2, 2]

    #-- 1. NAIVE ESTIMATE

    fit_naive <- lm(f ~ X1, data = dat_u)

    est_naive[sim] <- fit_naive$coefficients[2]

    err_naive[sim] <- summary(fit_naive)$coefficients[2, 2]

    #-- 2. CLASSIC ESTIMATE

    fit_classic <- lm(Y ~ X1, data = dat_l)

    est_classic[sim] <- fit_classic$coefficients[2]

    err_classic[sim] <- summary(fit_classic)$coefficients[2, 2]

    #-- 3. POSTPI ESTIMATE

    fit_postpi <- ipd(Y - f ~ X1, data = dat_l, unlabeled_data = dat_u,

      method = "postpi_boot", model = "ols", rel_func = "gam", n_t = n_t,

      nboot = nboot)

    est_postpi[sim] <- fit_postpi$coefficients[2]

    err_postpi[sim] <- fit_postpi$se[2]

    #-- 5. PPI

    fit_ppi <- ipd(Y - f ~ X1, data = dat_l, unlabeled_data = dat_u,

      method = "ppi", model = "ols")

    est_ppi[sim] <- fit_ppi$coefficients[2]

    err_ppi[sim] <- fit_ppi$se[2]

    #-- 6. PPI++

    fit_plusplus <- ipd(Y - f ~ X1, data = dat_l,

      unlabeled_data = dat_u, method = "ppi_plusplus", model = "ols")

    est_plusplus[sim] <- fit_plusplus$coefficients[2]

    err_plusplus[sim] <- fit_plusplus$se[2]

    #-- 7. PSPA

    fit_pspa <- ipd(Y - f ~ X1, data = dat_l,

      unlabeled_data = dat_u, method = "pspa", model = "ols")

    est_pspa[sim] <- fit_pspa$coefficients[2]

    err_pspa[sim] <- fit_pspa$se[2]
  }

  results <- data.frame(

    est_oracle,   err_oracle,
    est_naive,    err_naive,
    est_classic,  err_classic,
    est_postpi,   err_postpi,
    est_ppi,      err_ppi,
    est_plusplus, err_plusplus,
    est_pspa,   err_pspa)

  results <- results |>

    mutate(

      #- 0. Oracle

      upr_oracle  = est_oracle + 1.96 * err_oracle,
      lwr_oracle  = est_oracle - 1.96 * err_oracle,
      cov_oracle  = ((upr_oracle > true_beta) & (lwr_oracle < true_beta)),
      tval_oracle = est_oracle / err_oracle,
      pval_oracle = 2 * pt(-abs(tval_oracle), df = (n_u - 2)),

      #- 1. Naive

      upr_naive  = est_naive + 1.96 * err_naive,
      lwr_naive  = est_naive - 1.96 * err_naive,
      cov_naive  = ((upr_naive > true_beta) & (lwr_naive < true_beta)),
      tval_naive = est_naive / err_naive,
      pval_naive = 2 * pt(-abs(tval_naive), df = (n_u - 2)),

      #- 2. Classic

      upr_classic  = est_classic + 1.96 * err_classic,
      lwr_classic  = est_classic - 1.96 * err_classic,
      cov_classic  = ((upr_classic > true_beta) & (lwr_classic < true_beta)),
      tval_classic = est_classic / err_classic,
      pval_classic = 2 * pt(-abs(tval_classic), df = (n_l - 2)),

      #- 4. PostPI

      upr_postpi  = est_postpi + 1.96 * err_postpi,
      lwr_postpi  = est_postpi - 1.96 * err_postpi,
      cov_postpi  = ((upr_postpi > true_beta) & (lwr_postpi < true_beta)),
      tval_postpi = est_postpi / err_postpi,
      pval_postpi = 2 * pt(-abs(tval_postpi), df = (n_u - 2)),

      #- 5. PPI

      upr_ppi  = est_ppi + 1.96 * err_ppi,
      lwr_ppi  = est_ppi - 1.96 * err_ppi,
      cov_ppi  = ((upr_ppi > true_beta) & (lwr_ppi < true_beta)),
      tval_ppi = est_ppi / err_ppi,
      pval_ppi = 2 * pt(-abs(tval_ppi), df = (n_l + n_u - 2)),

      #- 6. PPI++

      upr_plusplus  = est_plusplus + 1.96 * err_plusplus,
      lwr_plusplus  = est_plusplus - 1.96 * err_plusplus,
      cov_plusplus  = ((upr_plusplus > true_beta) & (lwr_plusplus < true_beta)),
      tval_plusplus = est_plusplus / err_plusplus,
      pval_plusplus = 2 * pt(-abs(tval_plusplus), df = (n_l + n_u - 2)),

      #- 7. PSPA

      upr_pspa  = est_pspa + 1.96 * err_pspa,
      lwr_pspa  = est_pspa - 1.96 * err_pspa,
      cov_pspa  = ((upr_pspa > true_beta) & (lwr_pspa < true_beta)),
      tval_pspa = est_pspa / err_pspa,
      pval_pspa = 2 * pt(-abs(tval_pspa), df = (n_l + n_u - 2))
    )

  return(results)
}

#=== SIMULATION ================================================================

#-- Number of Simulated Replicates

nsims <- 1000

#-- Sample Sizes

n_t <- 100

n_l <- 100

n_u <- 1000

#-- Effect Size

true_beta <- 1

#-- Number of Bootstrap Replicates

nboot <- 100

results <- sim_func(nsims, n_t, n_l, n_u, nboot, dgm, true_beta)

#=== SUMMARIZE RESULTS =========================================================

methods <- c(

  "oracle" = "Oracle", "naive" = "Naive", "classic"  = "Classical",
  "postpi" = "PostPI", "ppi"   = "PPI",   "plusplus" = "PPI++",
  "pspa" = "PSPA")

coverages <- results |>

  summarise(across(starts_with("cov_"), mean))

cov_text <- coverages |>

  pivot_longer(everything()) |>

  separate(name, c("stat", "method")) |>

  mutate(

    method = factor(method, names(methods), methods),

    sampsize = case_when(

      method %in% c("Oracle", "Naive", "PostPI") ~ n_u,

      method == "Classical" ~ n_l,

      method %in% c("PPI", "PPI++", "PSPA") ~ n_l + n_u),

    ipd = if_else(method %in% c("PostPI", "PPI", "PPI++", "PSPA"),

      "IPD Method", "Benchmark Method"),

    beta = -0.5, id = 750,

    label = paste0("Coverage:\n", value*100, "%\nSamples Used:\n",

      prettyNum(sampsize, big.mark=",")))

results_long <- results |>

  mutate(id = 1:nsims) |>

  pivot_longer(-id) |>

  separate(name, c("stat", "method")) |>

  pivot_wider(names_from = stat, values_from = value) |>

  mutate(

    method = factor(method, names(methods), methods),

    ipd = if_else(method %in% c("PostPI", "PPI", "PPI++", "PSPA"),

      "IPD Method", "Benchmark Method"),

    cov = factor(cov, 0:1, c("No", "Yes"))) |>

  arrange(method, est) |>

  mutate(y = rep(1:nsims, 7))

p1 <- results_long |>

  filter(ipd == "Benchmark Method") |>

  ggplot(aes(y = y)) + theme_bw() +

    facet_grid(vars(ipd), vars(method)) +

    geom_linerange(aes(xmin = lwr, xmax = upr, col = cov)) +

    geom_point(aes(x = est, y = y), col = "white", size = 1) +

    geom_vline(xintercept = true_beta, linetype = 2) +

    geom_text(data = cov_text |> filter(ipd == "Benchmark Method"),

      mapping = aes(x = beta, y = id, label = label), size = 5) +

    scale_fill_manual(values = c("#AA4AC4", "#00C1D5")) +

    scale_color_manual(values = c("#AA4AC4", "#00C1D5")) +

    labs(

      x = "\nPoint Estimate and 95% CI",

      y = "Simulation #\n",

      color = "Covers?",

      fill = "Covers?") +

    theme(

      legend.position = "top",

      legend.text = element_text(size = 14),

      legend.title = element_text(face = "bold", size = 16),

      strip.background = element_rect(fill = "#00C1D5"),

      strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

      axis.text = element_text(size = 14),

      axis.title = element_text(face = "bold", size = 16)) +

    guides(

      color = guide_legend(override.aes = list(size = 5)))

p2 <- results_long |>

  filter(ipd != "Benchmark Method") |>

  ggplot(aes(y = y)) + theme_bw() +

    facet_grid(vars(ipd), vars(method)) +

    geom_linerange(aes(xmin = lwr, xmax = upr, col = cov)) +

    geom_point(aes(x = est, y = y), col = "white", size = 1) +

    geom_vline(xintercept = true_beta, linetype = 2) +

    geom_text(data = cov_text |> filter(ipd != "Benchmark Method"),

      mapping = aes(x = beta, y = id, label = label), size = 5) +

    scale_fill_manual(values = c("#AA4AC4", "#00C1D5")) +

    scale_color_manual(values = c("#AA4AC4", "#00C1D5")) +

    labs(

      x = "\nPoint Estimate and 95% CI",

      y = "Simulation #\n",

      color = "Covers?",

      fill = "Covers?") +

    theme(

      legend.position = "top",

      legend.text = element_text(size = 14),

      legend.title = element_text(face = "bold", size = 16),

      strip.background = element_rect(fill = "#00C1D5"),

      strip.text = element_text(face = "bold", color = "#FFFFFF", size = 14),

      axis.text = element_text(size = 14),

      axis.title = element_text(face = "bold", size = 16)) +

    guides(

      color = guide_legend(override.aes = list(size = 5)))

fig <- p1 / p2 + plot_layout(guides = "collect", axis_titles = "collect") &

  theme(legend.position = "top")

fig

#=== END =======================================================================
