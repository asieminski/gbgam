# A complex Poisson (canonical log link) model: a daily-mortality
# distributed-lag model (DLM). Deaths are counted over n = 300 days against
# the expected deaths E from reference rates (an offset log(E) in the
# formula, so exp(eta - log(E)) is the standardised mortality ratio). The
# model has a day-of-week factor, PM2.5 as a linear term, a slow trend in
# time and a heat effect over 7 lags as a linear functional smooth of the
# temperature and lag matrices; the first two weeks are dropped by `subset`.
# It covers list data with matrix covariates, unpenalised learners of
# several types (factor treatment contrasts, a numeric linear term, I(), the
# t2() null-space columns of a plain and of a functional smooth),
# initialise_intercept = FALSE with an offset, tau = 2, stopping by the .632+
# bootstrap estimate, Hutch++ EDF with block weights (checked against the
# exact trace on the same path), a free_selection parametric bootstrap with
# the model widening factor and, for the non-canonical sqrt link, a refit
# with update(), exact EDF and AIC stopping.

# Daily mortality with a distributed-lag heat effect. temp and lag are
# n x n_lag matrices (row i holds the temperatures of days i, i - 1, ...,
# i - n_lag + 1, and the lags 0, ..., n_lag - 1); eta is the true linear
# predictor, offset included
simulate_mortality_dlm <- function(n = 300, n_lag = 7, seed = 501) {
  set.seed(seed)
  n_full    <- n + n_lag - 1
  temp_full <- 15 + 8 * sin(2 * pi * seq_len(n_full) / 365) +
    stats::rnorm(n_full, sd = 3)
  temp <- t(vapply(seq_len(n),
                   function(i) temp_full[i + n_lag - 1 - 0:(n_lag - 1)],
                   numeric(n_lag)))
  lag  <- matrix(0:(n_lag - 1), n, n_lag, byrow = TRUE)
  time <- seq_len(n)
  days <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
  dow  <- factor(days[(time - 1) %% 7 + 1], levels = days)
  E    <- stats::runif(n, 15, 40)
  # Heat: the excess over 15 degrees, with lag weights decaying over days
  w_lag <- exp(-(0:(n_lag - 1)) / 2)
  heat  <- 0.04 * as.vector(pmax(temp - 15, 0) %*% (w_lag / sum(w_lag)))
  pm    <- stats::runif(n, 5, 40)
  eta   <- log(E) + 0.1 * sin(2 * pi * time / 150) +
    0.08 * (dow %in% c("Sat", "Sun")) + heat + 0.004 * pm
  y     <- stats::rpois(n, exp(eta))
  list(y = y, temp = temp, lag = lag, time = time, dow = dow, E = E,
       pm = pm, heat = heat, eta = eta)
}

test_that("Poisson DLM: .632+ stopping, Hutch++, bootstrap, sqrt refit", {
  skip_on_cran()
  d   <- simulate_mortality_dlm()
  fml <- y ~ dow + pm + offset(log(E)) + t2(time, k = 8, bs = "cr") +
    t2(temp, lag, k = c(5, 4), bs = "cr")

  # Fixes the Hutch++ probes, the resampling weights and the bootstrap
  # replicates. With B = 4 bootstrap samples about 0.632^4, or a sixth, of
  # the rows are in every sample, far above the 5% that triggers the warning.
  # The Hutch++ block weights are those of Meyer et al. (2021), with more
  # probes in the second sketch than in the first
  set.seed(9)
  expect_warning(
    mod <- gbgam(
      fml, family = stats::poisson(), data = d, subset = time > 14,
      learner_control   = gb_learner_ctrl(step_size = 0.02,
                                          initialise_intercept = FALSE,
                                          tau = 2),
      stop_control      = gb_stop_ctrl(
        "cv", iter_no_improv = 40, max_iter = 1500,
        resampling = gb_resampling("bootstrap", B = 4, estimator = "632plus",
                                   mode = "free_smoothing", chunk = 50)),
      edf_control       = gb_edf_ctrl("hutchpp", z = 60,
                                      cweights = c(1, 2, 3)),
      uq_control = gb_param_boot_ctrl(n_samples = 8,
                                      mode = "free_selection",
                                      correction = "model"),
      verbose = FALSE),
    "not left out by any resampling fit")

  idx <- as.integer(rownames(mod$model))
  n   <- length(idx)
  expect_identical(idx, 15:300)
  expect_true(mod$converged)
  expect_gt(mod$mstop, 0L)

  # Invariants; the canonical link lets Hutch++ be checked against the exact
  # trace on the same path
  expect_complex_invariants(mod, d, truth_eta = d$eta[idx], min_cor = 0.9,
                            hutchpp_rtol = 0.15)
  # The offset alone correlates highly with the true eta, so the signal is
  # checked on the log-SMR scale, eta - log(E)
  expect_gt(stats::cor(mod$linear.predictors - mod$offset,
                       d$eta[idx] - log(d$E[idx])), 0.85)

  # The heat effect is recovered by the linear functional smooth of the
  # temperature and lag matrices
  nd <- training_rows(mod, d)
  pt <- predict(mod, newdata = nd, type = "terms")
  expect_gt(stats::cor(pt[, "t2(temp,lag)"], d$heat[idx]), 0.6)

  # The learners: the factor's treatment contrasts, the numeric term, and
  # the penalised and null-space parts of both smooths
  labels <- vapply(mod$learners, function(lr) lr$var_label, character(1))
  expect_identical(labels, c(
    "(Intercept)", paste0("dow", c("Tue", "Wed", "Thu", "Fri", "Sat", "Sun")),
    "pm", "t2(time)r", "t2(time)n1", "t2(temp,lag)rr", "t2(temp,lag)nr",
    "t2(temp,lag)rn", "t2(temp,lag)n1"))

  # The resampling estimate decides mstop
  cv  <- mod$cv
  run <- sum(!is.na(mod$metrics$k_opts))
  expect_identical(cv$type, "bootstrap")
  expect_identical(cv$estimator, "632plus")
  expect_identical(cv$B, 4L)
  expect_length(cv$risk, run + 1L)
  expect_identical(cv$mstop, mod$mstop)
  expect_identical(mod$mstop, which.min(cv$risk) - 1L)
  # The stopping rule is checked every `chunk` = 50 iterations
  expect_gte(run - mod$mstop, 40L)
  expect_identical(run %% 50L, 0L)
  expect_identical(dim(cv$weights), c(n, 4L))
  expect_true(all(cv$weights == round(cv$weights)))
  expect_true(all(colSums(cv$weights) == n))
  expect_identical(cv$n_never_out, sum(rowSums(cv$weights == 0) == 0))
  expect_gt(cv$n_never_out, 0L)
  expect_equal(cv$oob, rowSums(cv$rep_risk), tolerance = 1e-12)
  # Where the out-of-sample estimate is below the no-information error, the
  # .632+ risk is (1 - weight) apparent + weight oob, with weight in
  # [0.632, 1]
  ok <- cv$oob <= cv$noinfo
  expect_true(all(ok))
  expect_equal(cv$risk[ok], (1 - cv$weight[ok]) * cv$apparent[ok] +
                 cv$weight[ok] * cv$oob[ok], tolerance = 1e-10)
  expect_true(all(cv$weight >= 0.632 - 1e-12 & cv$weight <= 1 + 1e-12))
  # The weight is 0.632 / (1 - 0.368 R) with Efron and Tibshirani's relative
  # overfitting rate R, from the oob error truncated at the no-information
  # error
  oob_t <- pmin(cv$oob, cv$noinfo)
  R     <- ifelse(cv$noinfo > cv$apparent & oob_t > cv$apparent,
                  (oob_t - cv$apparent) / (cv$noinfo - cv$apparent), 0)
  expect_equal(cv$weight, 0.632 / (1 - 0.368 * R), tolerance = 1e-10)
  expect_equal(mod$metrics$ls[seq_len(run + 1L)], cv$risk, tolerance = 1e-12)

  # Boosting starts from eta = log(E), that is SMR = 1, with edf~ = 0: the
  # apparent error at iteration 0 is the mean deviance of y against E. The
  # intercept-only MLE is still stored: the log of the overall SMR
  y <- mod$y
  E <- d$E[idx]
  expect_equal(mod$offset, log(E))
  expect_identical(mod$edfs[1], 0)
  expect_equal(cv$apparent[1], mean(stats::poisson()$dev.resids(y, E, 1)),
               tolerance = 1e-10)
  expect_equal(mod$intercept0, log(sum(y) / sum(E)), tolerance = 1e-6)
  # The no-information error at mstop is the mean deviance over all pairs of
  # observations and fitted means
  dev_pairs <- outer(y, mod$fitted.values,
                     function(a, b) stats::poisson()$dev.resids(a, b, 1))
  expect_equal(cv$noinfo[mod$mstop + 1L], mean(dev_pairs), tolerance = 1e-10)

  # The mean fitted SMR is close to the true one, both weekend coefficients
  # exceed every weekday one (the weekend effect was 0.08 on the log scale)
  # and the PM2.5 effect is positive
  smr <- exp(mod$linear.predictors - log(E))
  expect_lt(abs(mean(smr) - mean(exp(d$eta[idx] - log(E)))), 0.02)
  b <- mod$coefficients
  expect_gt(min(b[c("dowSat", "dowSun")]),
            max(b[c("dowTue", "dowWed", "dowThu", "dowFri")]))
  expect_gt(b[["pm"]], 0)

  # free_selection bootstrap with the model widening factor
  expect_identical(dim(mod$boot_coef), c(8L, length(b)))
  expect_length(mod$widening_factor, 1L)
  expect_gt(mod$widening_factor, 0.5)
  expect_lt(mod$widening_factor, 2)
  expect_equal(mod$Vp, mod$widening_factor^2 * stats::cov(mod$boot_coef),
               tolerance = 1e-10)

  # Standard errors at the training rows (list newdata with matrices) are
  # those of the linear predictor matrix and Vp
  pr <- predict(mod, newdata = nd, se.fit = TRUE)
  X  <- predict(mod, newdata = nd, type = "lpmatrix")
  expect_equal(as.vector(pr$fit), as.vector(predict(mod, newdata = nd)),
               tolerance = 1e-10)
  expect_true(all(is.finite(pr$se.fit) & pr$se.fit > 0))
  expect_equal(as.vector(pr$se.fit),
               unname(sqrt(rowSums((X %*% mod$Vp) * X))), tolerance = 1e-8)

  # summary(): the stopping loss is the resampling estimate
  s <- summary(mod)
  expect_identical(s$loss_name, ".632+ risk")
  expect_true("%d_stop_loss" %in% names(s$learner_table))
  expect_equal(s$best_loss, min(cv$risk), tolerance = 1e-10)
  expect_equal(sum(s$learner_table$sel), mod$mstop)
  # The shares are rounded to 0.1 per learner
  expect_lte(abs(sum(s$learner_table[["%d_stop_loss"]]) - 100),
             0.05 * nrow(s$learner_table))
  expect_match(paste(cli::cli_fmt(print(s)), collapse = " "),
               "best .632+ risk", fixed = TRUE)

  # The non-canonical sqrt link (variance-stabilising for counts), refitted
  # with update(). A log offset makes no sense on the sqrt scale, so
  # sqrt(E) enters as a linear term. Exact EDF and AIC stopping; Hutch++ is
  # not checked for non-canonical links
  expect_message(
    mod_sqrt <- update(
      mod,
      y ~ dow + pm + I(sqrt(E)) + t2(time, k = 8, bs = "cr") +
        t2(temp, lag, k = c(5, 4), bs = "cr"),
      family          = stats::poisson("sqrt"),
      edf_control     = gb_edf_ctrl("exact"),
      stop_control    = gb_stop_ctrl("aic", iter_no_improv = 20),
      learner_control = gb_learner_ctrl(step_size = 0.02),
      uq_control      = gb_uq_ctrl("none")),
    "Non-canonical link")
  expect_identical(mod_sqrt$family$link, "sqrt")
  expect_identical(as.integer(rownames(mod_sqrt$model)), idx)
  expect_identical(mod_sqrt$edf_compute_method, "exact")
  expect_true(all(mod_sqrt$offset == 0))
  expect_null(mod_sqrt$Vp)
  expect_true("I(sqrt(E))" %in%
                vapply(mod_sqrt$learners, function(lr) lr$var_label, ""))
  expect_complex_invariants(mod_sqrt, d, truth_eta = exp(d$eta[idx] / 2),
                            min_cor = 0.9)
  # The intercept starts at the intercept-only MLE, sqrt(mean(y)), and the
  # fit tracks the true mean about as well as the log-link fit
  expect_equal(mod_sqrt$intercept0, sqrt(mean(y)), tolerance = 1e-6)
  expect_equal(mod_sqrt$edfs[1], 1, tolerance = 1e-8)
  mu_true <- exp(d$eta[idx])
  expect_lt(mean((mod_sqrt$fitted.values - mu_true)^2),
            2 * mean((mod$fitted.values - mu_true)^2))
  # sqrt(E) dominates sqrt(mu) and is in the model, so the signal is checked
  # on the SMR scale, mu / E
  expect_gt(stats::cor(mod_sqrt$fitted.values / E, mu_true / E), 0.8)
  # A sanity check: both fits follow E
  expect_gt(stats::cor(mod_sqrt$fitted.values, mod$fitted.values), 0.95)
})
