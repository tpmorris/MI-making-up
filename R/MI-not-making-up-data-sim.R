# Load required packages
library(tidyverse)
library(mvtnorm)
library(mice)
library(broom)
library(mirai)
library(ragg)
library(rsimsum)
library(here)

# Define a function for simulating a dataset
simulate_data <- function(
    n_obs = 500,
    n_mis = 500,
    Sigma = matrix(data = c(1, 0.5, 0.5, 1), ncol = 2),
    p_Z = 0.5) {
  YX <- mvtnorm::rmvnorm(n = n_obs, sigma = Sigma)
  Y <- YX[, 1]
  X <- YX[, 2]
  length(Y) <- n_obs + n_mis
  length(X) <- n_obs + n_mis
  Z <- rbinom(n = n_obs + n_mis, size = 1, prob = p_Z)
  out <- data.frame(Y = Y, X = X, Z = Z)
  return(out)
}

do_analysis <- function(data, method) {
  # Ensure the packages are loaded
  # -> - useful for parallel processing with {mirai}
  library(broom)
  library(mice)
  # Different methods, different analysis steps
  if (method == 1) {
    fit1 <- lm(Y ~ X, data = data)
    out1 <- tidy(fit1) |>
      dplyr::filter(term == "X") |>
      dplyr::mutate(
        method = "Full",
        nimps = 0,
        target = "x",
        df = glance(fit1)$df.residual
      ) |>
      dplyr::select(estimate, std.error, method, nimps, target)
    fit2 <- lm(Y ~ X + Z, data = data)
    out2 <- tidy(fit2) |>
      dplyr::filter(term == "Z") |>
      dplyr::mutate(
        method = "Full",
        nimps = 0,
        target = "z",
        df = glance(fit2)$df.residual
      ) |>
      dplyr::select(estimate, std.error, method, nimps, target, df)
  }
  if (method == 2) {
    impdata <- mice::mice(
      data = data,
      m = 5,
      method = "norm",
      visitSequence = "monotone",
      maxit = 1,
      print = FALSE
    )
    fit1 <- with(impdata, lm(Y ~ X))
    est1 <- mice::pool(fit1)
    out1 <- tidy(est1) |>
      dplyr::filter(term == "X") |>
      dplyr::mutate(
        method = "MI",
        nimps = 5,
        target = "x",
        df = dfcom
      ) |>
      dplyr::select(estimate, std.error, method, nimps, target, df, riv, fmi)
    fit2 <- with(impdata, lm(Y ~ X + Z))
    est2 <- mice::pool(fit2)
    out2 <- tidy(est2) |>
      dplyr::filter(term == "Z") |>
      dplyr::mutate(
        method = "MI",
        nimps = 5,
        target = "z",
        df = dfcom
      ) |>
      dplyr::select(estimate, std.error, method, nimps, target, df, riv, fmi)
  }
  if (method == 3) {
    impdata <- mice::mice(
      data = data,
      m = 100,
      method = "norm",
      visitSequence = "monotone",
      maxit = 1,
      print = FALSE
    )
    fit1 <- with(impdata, lm(Y ~ X))
    est1 <- mice::pool(fit1)
    out1 <- tidy(est1) |>
      dplyr::filter(term == "X") |>
      dplyr::mutate(
        method = "MI",
        nimps = 100,
        target = "x",
        df = dfcom
      ) |>
      dplyr::select(estimate, std.error, method, nimps, target, df, riv, fmi)
    fit2 <- with(impdata, lm(Y ~ X + Z))
    est2 <- mice::pool(fit2)
    out2 <- tidy(est2) |>
      dplyr::filter(term == "Z") |>
      dplyr::mutate(
        method = "MI",
        nimps = 100,
        target = "z",
        df = dfcom
      ) |>
      dplyr::select(estimate, std.error, method, nimps, target, df, riv, fmi)
  }
  out <- dplyr::bind_rows(out1, out2)
  return(out)
}

# Setup paraller processing, using all available cores
daemons(parallel::detectCores())

# Run B repetitions using purrr::map()
B <- 10000 
time_start <- Sys.time()
out <- map(
  .x = seq(B),
  .f = in_parallel(
    \(i) {
      simdata <- simulate_data()
      dplyr::bind_rows(
        do_analysis(simdata, method = 1),
        do_analysis(simdata, method = 2),
        do_analysis(simdata, method = 3)
      ) |>
        dplyr::mutate(rep_id = i)
    },
    simulate_data = simulate_data,
    do_analysis = do_analysis
  ),
  .progress = TRUE
)
time_stop <- Sys.time()
time_diff <- difftime(time_stop, time_start, units = "mins")
time_diff

# Reset to sequential processing
daemons(0)

# Combine repetitions in a single dataset and add true values of each estimand
outb <- list_rbind(out) |>
  mutate(true = ifelse(target == "x", 0.5, 0.0))

# Of course we use {rsimsum} for this ;)
relprec <- simsum(
  data = outb,
  estvarname = "estimate",
  se = "std.error",
  methodvar = "nimps",
  ref = "0",
  by = "target",
  true = "true"
)
fmi <- simsum(
  data = outb,
  estvarname = "fmi",
  methodvar = "nimps",
  ref = "0",
  by = "target",
  true = 0.0
)

# Combine the results of the simulation study
# and prepare for plotting
finres <- bind_rows(
  tidy(summary(relprec, stat = "relprec")),
  tidy(summary(fmi, stat = "bias"))
) |>
  filter(nimps != "0") |>
  mutate(
    ref = case_when(
      stat == "bias" ~ 0.5,
      stat == "relprec" ~ 0.0,
      .default = NA_real_
    )
  ) |>
  mutate(
    nimps = factor(
      nimps,
      levels = c("5", "100"),
      labels = c("5 Imputations", "100 Imputations")
    )
  ) |>
  mutate(
    stat = factor(
      stat,
      levels = c("relprec", "bias"),
      labels = c("Relative ~ '%' ~ Precision", "FMI")
    )
  ) |>
  mutate(
    target = factor(
      target,
      levels = c("x", "z"),
      labels = c("Estimand: beta[X]", "Estimand: beta[Z]")
    )
  )

# Plot
pp <- finres |>
  ggplot(aes(x = nimps, y = est)) +
  geom_hline(aes(yintercept = ref), linetype = "dashed") +
  geom_errorbar(
    aes(ymin = est - mcse, ymax = est + mcse),
    # 1 MCSE confidence bands
    width = 0
  ) +
  geom_point() +
  facet_grid(target ~ stat, scales = "free_x", labeller = label_parsed) +
  coord_flip() +
  theme_bw(base_size = 12, base_family = "Atkinson Hyperlegible") +
  labs(x = "", y = "")
pp

# Export
# Plot(s):
ggsave(
  filename = here("R/performance.png"),
  plot = pp,
  dev = ragg::agg_png,
  width = 6,
  height = 4,
  dpi = 600
)
ggsave(
  filename = here("R/performance.pdf"),
  plot = pp,
  dev = cairo_pdf,
  width = 6,
  height = 4,
  dpi = 600
)
# Data:
write_rds(x = outb, file = here("R/estimates.rds"))
write_rds(x = finres, file = here("R/perf_final.rds"))
