# =============================================================================
# Coverage simulation tests for the gbgam uncertainty quantification methods
# =============================================================================
#
# PURPOSE
# -------
# Verify that the parametric bootstrap and Bayesian backfitting produce
# confidence / credible intervals with approximately nominal pointwise
# coverage of the TRUE data-generating function.
#
# WHY NOT IN tests/testthat/
# --------------------------
# These tests run \mathcal{O}(n_sim) model fits and take several minutes.
# They are excluded from R CMD check and devtools::check() via .Rbuildignore.
#
# HOW TO RUN
# ----------
# From the package root (with the package loaded/installed):
#
#   pkgload::load_all()
#   testthat::test_file("tests/coverage_tests/test-uq-coverage.R")
#
# or for a quicker smoke test:
#
#   source("tests/coverage_tests/test-uq-coverage.R")
#
# WHAT IS TESTED
# --------------
# For each combination of family x dimension:
#   1. Simulate n_sim independent datasets from the known DGP.
#   2. Fit gbgam with either the bootstrap or Bayesian backfitting.
#   3. Compute pointwise 95 % CIs on the LINK scale.
#   4. Compare against the true linear predictor (link-scale f).
#   5. Assert that mean pointwise coverage is in [lo_tol, hi_tol].
#
# TOLERANCES
# ----------
# With only n_sim = 30 replications the Monte-Carlo error on coverage is
# roughly +/-sqrt(0.95*0.05/30) ~= +/-0.04.  We therefore use a wide
# tolerance of [0.72, 1.00] as a sanity check, and a tighter [0.80, 1.00] for
# the calibrated bootstrap (which targets training-data coverage = 0.95 by
# construction) and for the Gibbs sampler (exact conditionals for Gaussian).
#
# =============================================================================

library(testthat)

# Load the package (works both when sourced and when run via test_file)
if (!isNamespaceLoaded("gbgam")) pkgload::load_all(quiet = TRUE)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Run n_sim coverage experiments and return the per-rep coverage vector.
#
# sim_fn(seed) must return a list with:
#   $formula  - formula object
#   $data     - data.frame
#   $f_link   - true linear predictor at training points (link scale)
#
# `alpha` is the nominal mis-coverage level; CIs are (1-alpha)*100 % Wald.
run_coverage_sim <- function(
    sim_fn,
    family,
    n_sim              = 30L,
    max_iter           = 500L,
    step_size          = 0.05,
    uq_method   = c("bootstrap", "bayes"),
    n_samples          = 80L,
    n_burnin           = 30L,
    alpha              = 0.05,
    seed_offset        = 0L,
    edf_compute_method = "hutchpp",
    stopping_method    = "aic",
    verbose            = FALSE) {

  uq_method <- match.arg(uq_method)
  z_crit  <- qnorm(1 - alpha / 2)
  coverages <- rep(NA_real_, n_sim)

  for (i in seq_len(n_sim)) {
    setup <- sim_fn(seed = seed_offset + i)
    tryCatch({
      # gb_uq_ctrl() passes the settings to the method's own control
      # function, which refuses the settings of the other method
      ic_args <- list(method = uq_method, n_samples = n_samples)
      if (uq_method == "bootstrap") ic_args$alpha    <- alpha
      if (uq_method == "bayes")     ic_args$n_burnin <- n_burnin

      mod <- gbgam(
        setup$formula,
        data              = setup$data,
        family            = family,
        learner_control   = gb_learner_ctrl(step_size = step_size),
        stop_control      = gb_stop_ctrl(stopping_method = stopping_method,
                                         max_iter = max_iter),
        edf_control = gb_edf_ctrl(edf_compute_method = edf_compute_method),
        uq_control  = do.call(gb_uq_ctrl, ic_args),
        verbose     = verbose
      )

      # Link-scale prediction with SE
      pred   <- predict(mod, type = "link", se.fit = TRUE)
      ci_lo  <- pred$fit - z_crit * pred$se.fit
      ci_hi  <- pred$fit + z_crit * pred$se.fit
      f_true <- setup$f_link

      coverages[i] <- mean(f_true >= ci_lo & f_true <= ci_hi, na.rm = TRUE)

    }, error = function(e) {
      warning(sprintf("rep %d failed: %s", i, conditionMessage(e)),
              call. = FALSE)
    })
  }
  coverages
}

# Assert that the mean of non-NA coverages is in [lo, hi].
expect_coverage_in <- function(coverages, lo, hi,
                                label = deparse(substitute(coverages))) {
  ok   <- !is.na(coverages)
  n_ok <- sum(ok)
  if (n_ok == 0) testthat::fail(paste(label, ": all reps failed"))
  mc   <- mean(coverages[ok])
  msg  <- sprintf("%s: mean coverage = %.3f  (expected in [%.2f, %.2f], n_ok = %d)",
                  label, mc, lo, hi, n_ok)
  testthat::expect_gte(mc, lo, label = msg)
  testthat::expect_lte(mc, hi, label = msg)
  invisible(mc)
}


# ---------------------------------------------------------------------------
# DGP factories
# ---------------------------------------------------------------------------

# 1-D smooth: y = sin(2*pi*x) + noise
sim_1d_gaussian <- function(seed, n = 200L) {
  set.seed(seed)
  x <- sort(runif(n))
  f <- sin(2 * pi * x)
  y <- f + rnorm(n, sd = 0.4)
  list(formula = y ~ t2(x, k = 10, bs = "cr"),
       data    = data.frame(x = x, y = y),
       f_link  = f)          # identity link: f_link = f
}

# 1-D Poisson: log(lambda) = 0.5 + sin(2*pi*x)
sim_1d_poisson <- function(seed, n = 200L) {
  set.seed(seed)
  x   <- sort(runif(n))
  eta <- 0.5 + sin(2 * pi * x)
  y   <- rpois(n, lambda = exp(eta))
  list(formula = y ~ t2(x, k = 10, bs = "cr"),
       data    = data.frame(x = x, y = y),
       f_link  = eta)        # log link: f_link = log(mu)
}

# 1-D Binomial (aggregate): logit(p) = sin(2*pi*x)
sim_1d_binomial <- function(seed, n = 200L) {
  set.seed(seed)
  x        <- sort(runif(n))
  eta      <- sin(2 * pi * x)
  n_trials <- sample(10:30, n, replace = TRUE)
  s        <- rbinom(n, size = n_trials, prob = plogis(eta))
  f        <- n_trials - s
  list(formula = cbind(s, f) ~ t2(x, k = 10, bs = "cr"),
       data    = data.frame(x = x, s = s, f = f, n_trials = n_trials),
       f_link  = eta)        # logit link: f_link = logit(p)
}

# 4-D additive Gaussian using mgcv::gamSim example 1
sim_4d_gaussian <- function(seed, n = 400L) {
  set.seed(seed)
  dat <- mgcv::gamSim(1, n = n, dist = "normal", scale = 1, verbose = FALSE)
  # dat$f is the true linear predictor (sum of component smooth functions)
  list(formula = y ~ t2(x0, k = 8, bs = "cr") +
                     t2(x1, k = 8, bs = "cr") +
                     t2(x2, k = 8, bs = "cr") +
                     t2(x3, k = 8, bs = "cr"),
       data   = dat,
       f_link = dat$f)       # identity link
}

# 4-D Poisson using mgcv::gamSim example 1 (Poisson variant)
sim_4d_poisson <- function(seed, n = 300L) {
  set.seed(seed)
  dat <- mgcv::gamSim(1, n = n, dist = "poisson", scale = 0.2, verbose = FALSE)
  # gamSim stores true log(mu) in dat$f for Poisson
  list(formula = y ~ t2(x0, k = 8, bs = "cr") +
                     t2(x1, k = 8, bs = "cr") +
                     t2(x2, k = 8, bs = "cr") +
                     t2(x3, k = 8, bs = "cr"),
       data   = dat,
       f_link = dat$f)
}


# ---------------------------------------------------------------------------
# 1-D Gaussian -- bootstrap
# ---------------------------------------------------------------------------

test_that("Bootstrap: 1-D Gaussian achieves >= 0.80 mean pointwise coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_1d_gaussian,
    family           = gaussian(),
    n_sim            = 30L,
    uq_method = "bootstrap",
    n_samples        = 120L,
    alpha            = 0.05,
    seed_offset      = 1000L
  )
  expect_coverage_in(covs, lo = 0.80, hi = 1.00,
                     label = "Bootstrap 1D Gaussian")
})

# ---------------------------------------------------------------------------
# 1-D Gaussian -- Bayes (Gibbs)
# ---------------------------------------------------------------------------

test_that("Bayes Gibbs: 1-D Gaussian achieves >= 0.78 mean pointwise coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_1d_gaussian,
    family           = gaussian(),
    n_sim            = 30L,
    uq_method = "bayes",
    n_samples        = 300L,
    n_burnin         = 100L,
    alpha            = 0.05,
    seed_offset      = 2000L
  )
  expect_coverage_in(covs, lo = 0.78, hi = 1.00,
                     label = "Bayes Gibbs 1D Gaussian")
})

# ---------------------------------------------------------------------------
# 1-D Poisson -- bootstrap
# ---------------------------------------------------------------------------

test_that("Bootstrap: 1-D Poisson achieves >= 0.72 mean pointwise coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_1d_poisson,
    family           = poisson(),
    n_sim            = 30L,
    uq_method = "bootstrap",
    n_samples        = 120L,
    alpha            = 0.05,
    seed_offset      = 3000L
  )
  expect_coverage_in(covs, lo = 0.72, hi = 1.00,
                     label = "Bootstrap 1D Poisson")
})

# ---------------------------------------------------------------------------
# 1-D Poisson -- Bayes (MH)
# ---------------------------------------------------------------------------

test_that("Bayes MH: 1-D Poisson achieves >= 0.70 mean pointwise coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_1d_poisson,
    family           = poisson(),
    n_sim            = 30L,
    uq_method = "bayes",
    n_samples        = 400L,
    n_burnin         = 150L,
    alpha            = 0.05,
    seed_offset      = 4000L
  )
  expect_coverage_in(covs, lo = 0.70, hi = 1.00,
                     label = "Bayes MH 1D Poisson")
})

# ---------------------------------------------------------------------------
# 1-D Binomial -- bootstrap
# ---------------------------------------------------------------------------

test_that("Bootstrap: 1-D Binomial achieves >= 0.72 mean pointwise coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_1d_binomial,
    family           = binomial(),
    n_sim            = 30L,
    uq_method = "bootstrap",
    n_samples        = 120L,
    alpha            = 0.05,
    seed_offset      = 5000L
  )
  expect_coverage_in(covs, lo = 0.72, hi = 1.00,
                     label = "Bootstrap 1D Binomial")
})

# ---------------------------------------------------------------------------
# 1-D Binomial -- Bayes (MH)
# ---------------------------------------------------------------------------

test_that("Bayes MH: 1-D Binomial achieves >= 0.70 mean pointwise coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_1d_binomial,
    family           = binomial(),
    n_sim            = 30L,
    # AIC stopping needs ~2000-4000 iterations here; the default cap of 500
    # leaves a heavily under-fitted (biased) boosting fit
    max_iter         = 5000L,
    uq_method = "bayes",
    n_samples        = 400L,
    n_burnin         = 150L,
    alpha            = 0.05,
    seed_offset      = 6000L
  )
  expect_coverage_in(covs, lo = 0.70, hi = 1.00,
                     label = "Bayes MH 1D Binomial")
})

# ---------------------------------------------------------------------------
# 4-D Gaussian (gamSim 1) -- bootstrap, reduced n_sim because slower
# ---------------------------------------------------------------------------

test_that("Bootstrap: 4-D Gaussian (gamSim 1) achieves >= 0.78 mean coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_4d_gaussian,
    family           = gaussian(),
    n_sim            = 15L,
    max_iter         = 150L,
    step_size        = 0.02,
    uq_method = "bootstrap",
    n_samples        = 100L,
    alpha            = 0.05,
    seed_offset      = 7000L
  )
  expect_coverage_in(covs, lo = 0.78, hi = 1.00,
                     label = "Bootstrap 4D Gaussian")
})

# ---------------------------------------------------------------------------
# 4-D Gaussian (gamSim 1) -- Bayes (Gibbs)
# ---------------------------------------------------------------------------

test_that("Bayes Gibbs: 4-D Gaussian (gamSim 1) achieves >= 0.75 mean coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_4d_gaussian,
    family           = gaussian(),
    n_sim            = 15L,
    max_iter         = 150L,
    step_size        = 0.02,
    uq_method = "bayes",
    n_samples        = 300L,
    n_burnin         = 100L,
    alpha            = 0.05,
    seed_offset      = 8000L
  )
  expect_coverage_in(covs, lo = 0.75, hi = 1.00,
                     label = "Bayes Gibbs 4D Gaussian")
})

# ---------------------------------------------------------------------------
# 4-D Poisson (gamSim 1) -- bootstrap only (MH for Poisson 4D is very slow)
# ---------------------------------------------------------------------------

test_that("Bootstrap: 4-D Poisson (gamSim 1) achieves >= 0.70 mean coverage", {
  covs <- run_coverage_sim(
    sim_fn           = sim_4d_poisson,
    family           = poisson(),
    n_sim            = 15L,
    max_iter         = 150L,
    step_size        = 0.02,
    uq_method = "bootstrap",
    n_samples        = 100L,
    alpha            = 0.05,
    seed_offset      = 9000L
  )
  expect_coverage_in(covs, lo = 0.70, hi = 1.00,
                     label = "Bootstrap 4D Poisson")
})

# ---------------------------------------------------------------------------
# Calibration check: bootstrap widening should bring coverage towards 0.95
# ---------------------------------------------------------------------------

test_that("Calibrated bootstrap: inflated coverage exceeds un-inflated coverage", {
  # Single dataset, two fits: with and without widening.
  # The inflated version should have higher (or equal) training-data coverage.
  setup <- sim_1d_gaussian(seed = 9999L, n = 200L)

  fit_base <- function(alpha_val) {
    gbgam(
      setup$formula,
      data            = setup$data,
      family          = gaussian(),
      learner_control = gb_learner_ctrl(tau = 0, step_size = 0.05),
      stop_control    = gb_stop_ctrl(stopping_method = "none", max_iter = 80L),
      edf_control     = gb_edf_ctrl(edf_compute_method = "none"),
      uq_control      = gb_uq_ctrl(method = "bootstrap",
                                             n_samples = 100L, alpha = alpha_val)
    )
  }

  mod_uninflated <- fit_base(alpha_val = NULL)   # raw covariance
  mod_inflated   <- fit_base(alpha_val = 0.05)   # widened

  z_crit <- qnorm(0.975)
  pred_u <- predict(mod_uninflated, type = "link", se.fit = TRUE)
  pred_i <- predict(mod_inflated,   type = "link", se.fit = TRUE)

  cov_u <- mean(setup$f_link >= pred_u$fit - z_crit * pred_u$se.fit &
                setup$f_link <= pred_u$fit + z_crit * pred_u$se.fit)
  cov_i <- mean(setup$f_link >= pred_i$fit - z_crit * pred_i$se.fit &
                setup$f_link <= pred_i$fit + z_crit * pred_i$se.fit)

  # Inflation factor must be >= 1
  expect_gte(mod_inflated$widening_factor, 1.0 - 1e-6)

  # Inflated coverage must be >= uninflated coverage (or equal if both saturate)
  expect_gte(cov_i, cov_u - 0.02)  # allow tiny floating-point slack
})


# ---------------------------------------------------------------------------
# Structural sanity checks (cheap, run even on reduced machines)
# ---------------------------------------------------------------------------

test_that("Bootstrap SEs are finite and positive for all families", {
  for (setup_fn in list(sim_1d_gaussian, sim_1d_poisson, sim_1d_binomial)) {
    setup <- setup_fn(seed = 42L, n = 120L)
    fam   <- if (identical(setup_fn, sim_1d_gaussian)) gaussian()
             else if (identical(setup_fn, sim_1d_poisson)) poisson()
             else binomial()

    mod <- gbgam(
      setup$formula,
      data            = setup$data,
      family          = fam,
      learner_control = gb_learner_ctrl(tau = 0, step_size = 0.05),
      stop_control    = gb_stop_ctrl(stopping_method = "none", max_iter = 40L),
      edf_control     = gb_edf_ctrl(edf_compute_method = "none"),
      uq_control      = gb_uq_ctrl(method = "bootstrap",
                                             n_samples = 20L, alpha = NULL)
    )
    pred <- predict(mod, type = "link", se.fit = TRUE)
    expect_true(all(is.finite(pred$se.fit)),
                info = paste("non-finite SEs for family:", fam$family))
    expect_true(all(pred$se.fit > 0),
                info = paste("non-positive SEs for family:", fam$family))
  }
})

test_that("Bayes posterior Vp is symmetric positive-semi-definite (Gaussian Gibbs)", {
  setup <- sim_1d_gaussian(seed = 101L, n = 150L)
  mod <- gbgam(
    setup$formula,
    data            = setup$data,
    family          = gaussian(),
    learner_control = gb_learner_ctrl(tau = 0, step_size = 0.05),
    stop_control    = gb_stop_ctrl(stopping_method = "none", max_iter = 50L),
    edf_control     = gb_edf_ctrl(edf_compute_method = "exact"),
    uq_control      = gb_uq_ctrl(method = "bayes",
                                           n_samples = 200L, n_burnin = 50L)
  )
  V <- vcov(mod)
  # Symmetry
  expect_equal(V, t(V), tolerance = 1e-10)
  # All eigenvalues non-negative
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > -1e-8),
              info = paste("negative eigenvalue:", min(ev)))
})

test_that("Bayes posterior Vp is symmetric PSD for Poisson MH", {
  setup <- sim_1d_poisson(seed = 102L, n = 150L)
  mod <- gbgam(
    setup$formula,
    data            = setup$data,
    family          = poisson(),
    learner_control = gb_learner_ctrl(tau = 0, step_size = 0.05),
    stop_control    = gb_stop_ctrl(stopping_method = "none", max_iter = 50L),
    edf_control     = gb_edf_ctrl(edf_compute_method = "exact"),
    uq_control      = gb_uq_ctrl(method = "bayes",
                                           n_samples = 200L, n_burnin = 100L)
  )
  V  <- vcov(mod)
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > -1e-8),
              info = paste("negative eigenvalue:", min(ev)))
  # Every sampled learner must mix (NA = learner held fixed at zero)
  ar <- mod$mh_accept_rate
  ar <- ar[!is.na(ar)]
  expect_true(all(ar >= 0 & ar <= 1))
  expect_true(all(ar > 0.05),
              info = paste("MH sampler stuck for some learner:",
                           paste(round(ar, 3), collapse = " ")))
})
