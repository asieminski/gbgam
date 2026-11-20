# bootstrap_gbgam(): output contract, modes, corrections, argument checks and
# early returns; simulation and replays are in test-bootstrap-*.R.

# Check the basic output contract of a bootstrap fit
check_boot_contract <- function(mod_b, n_rep, p) {
  testthat::expect_s3_class(mod_b, "gbgam")
  testthat::expect_equal(dim(mod_b$Vp), c(p, p))
  testthat::expect_equal(dim(mod_b$boot_coef), c(n_rep, p))
  # Vp symmetric and positive semi-definite
  testthat::expect_equal(mod_b$Vp, t(mod_b$Vp))
  testthat::expect_true(
    all(eigen(mod_b$Vp, only.values = TRUE)$values >= -1e-10))
  # All entries finite
  testthat::expect_true(all(is.finite(mod_b$boot_coef)))
  testthat::expect_true(all(is.finite(mod_b$Vp)))
}


# ---- Structure and fields ----------------------------------------------------

test_that("bootstrap_gbgam returns correct structure (alpha = NULL)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 100)
  mod <- fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 20L,
                   uq_control = gb_uq_ctrl("none"))

  set.seed(42)
  n_rep <- 30L
  # edf "none" + gaussian: the dispersion warning is expected
  expect_warning(
    mod_b <- bootstrap_gbgam(mod, n_samples = n_rep, alpha = NULL),
    "edf_compute_method"
  )
  p      <- length(mod$coefficients)

  check_boot_contract(mod_b, n_rep, p)

  # Row and column names on Vp match coefficient names
  expect_equal(rownames(mod_b$Vp), names(mod$coefficients))
  expect_equal(colnames(mod_b$Vp), names(mod$coefficients))

  # alpha = NULL => no widening computed
  expect_null(mod_b$widening_factor)
})

test_that("alpha widens the replicates' covariance by the model factor", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 104)
  # An AIC-stopped fit, whose factor is near 1 (no factor-range warning)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl("aic", iter_no_improv = 50L),
                uq_control = gb_uq_ctrl("none"))

  set.seed(42)
  expect_no_warning(
    mod_infl <- bootstrap_gbgam(mod, n_samples = 40L, alpha = 0.05))
  set.seed(42)
  mod_raw  <- bootstrap_gbgam(mod, n_samples = 40L, alpha = NULL)
  expect_identical(mod_infl$boot_coef, mod_raw$boot_coef)

  # alpha = NULL: no factor, and Vp is the replicates' covariance
  expect_null(mod_raw$widening_factor)
  expect_equal(mod_raw$Vp, stats::cov(mod_raw$boot_coef), ignore_attr = TRUE)
  # alpha = 0.05: the factor is the 0.95 quantile over observations and
  # replicates of |mu^(r) - mu~| / sd(mu^(r)), divided by z_0.975, computed
  # here from predict()'s design (exact up to rounding: at most 6e-15
  # measured over 4 seeds)
  X    <- stats::predict(mod, type = "lpmatrix")
  mu_r <- mod_raw$boot_coef %*% t(X)
  dev  <- abs(sweep(mu_r, 2L, fitted(mod), "-"))
  dev  <- sweep(dev, 2L, apply(mu_r, 2L, stats::sd), "/")
  c_ref <- stats::quantile(dev, 0.95, names = FALSE) / stats::qnorm(0.975)
  expect_equal(mod_infl$widening_factor, c_ref, tolerance = 1e-10)
  expect_equal(mod_infl$Vp, c_ref^2 * stats::cov(mod_raw$boot_coef),
               ignore_attr = TRUE)

  # Poisson (log link): the factor is computed from the fitted means, not
  # from the linear predictors (which give 1.06649 against 1.06571)
  pois <- fit_1d(simulate_poisson_1d(n = 100, seed = 3), stats::poisson(),
                 stop_control = gb_stop_ctrl("aic", iter_no_improv = 50L),
                 uq_control = gb_uq_ctrl("none"))
  set.seed(42)
  expect_no_warning(
    pb <- bootstrap_gbgam(pois, n_samples = 40L, alpha = 0.05))
  X    <- stats::predict(pois, type = "lpmatrix")
  mu_r <- exp(pb$boot_coef %*% t(X))
  dev  <- abs(sweep(mu_r, 2L, fitted(pois), "-"))
  dev  <- sweep(dev, 2L, apply(mu_r, 2L, stats::sd), "/")
  c_ref <- stats::quantile(dev, 0.95, names = FALSE) / stats::qnorm(0.975)
  expect_equal(pb$widening_factor, c_ref, tolerance = 1e-10)
})


# ---- Families (output contract only) -----------------------------------------

test_that("bootstrap_gbgam satisfies structural contract for all supported families", {
  skip_on_cran()
  specs <- list(
    gaussian         = list(family = stats::gaussian(),
                            sim_fn = simulate_gaussian_1d,    seed = 107L),
    poisson          = list(family = stats::poisson(),
                            sim_fn = simulate_poisson_1d,     seed = 108L),
    binomial         = list(family = stats::binomial(),
                            sim_fn = simulate_binomial_1d,    seed = 103L),
    Gamma            = list(family = stats::Gamma(link = "log"),
                            sim_fn = simulate_gamma_1d,       seed = 105L),
    inverse.gaussian = list(family = stats::inverse.gaussian(link = "log"),
                            sim_fn = simulate_invgaussian_1d, seed = 106L)
  )

  for (nm in names(specs)) {
    spec  <- specs[[nm]]
    dat   <- spec$sim_fn(n = 80L, seed = spec$seed)
    mod   <- fit_fixed(y ~ t2(x, k = 5, bs = "cr"), dat, family = spec$family,
                       max_iter = 15L,
                       uq_control = gb_uq_ctrl("none"))
    set.seed(13L)
    n_rep <- 20L
    # edf "none": free-dispersion families warn about the dispersion
    if (nm %in% c("gaussian", "Gamma", "inverse.gaussian")) {
      expect_warning(
        mod_b <- bootstrap_gbgam(mod, n_samples = n_rep, alpha = NULL),
        "edf_compute_method"
      )
    } else {
      mod_b <- bootstrap_gbgam(mod, n_samples = n_rep, alpha = NULL)
    }
    check_boot_contract(mod_b, n_rep, length(mod$coefficients))
  }
})


# ---- Mode validation errors --------------------------------------------------

test_that("bootstrap_gbgam free_* modes error when the main model has no EDF tracking", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 100)
  mod <- fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 15L,
                   uq_control = gb_uq_ctrl("none"))

  for (bad_mode in c("free_smoothing", "free_selection", "free_stopping")) {
    expect_error(
      bootstrap_gbgam(mod, n_samples = 5L, mode = bad_mode),
      "requires EDF tracking",
      info = paste("expected error for mode =", bad_mode)
    )
  }
})

test_that("bootstrap_gbgam errors when free_stopping is used with stopping_method = 'none'", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 100)
  # EDF tracking present but no stopping criterion
  mod <- fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 20L,
                   edf = "hutchpp",
                   uq_control = gb_uq_ctrl("none"))

  expect_error(
    bootstrap_gbgam(mod, n_samples = 5L, mode = "free_stopping"),
    "requires the main fit to stop by"
  )
})


# ---- Bootstrap modes (output contract) ---------------------------------------

test_that("bootstrap_gbgam free_smoothing mode satisfies structural contract", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 111)
  mod <- fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 25L,
                   edf = "hutchpp",
                   uq_control = gb_uq_ctrl("none"))

  set.seed(55)
  n_rep <- 20L
  mod_b  <- bootstrap_gbgam(mod, n_samples = n_rep, alpha = NULL,
                              mode = "free_smoothing")
  check_boot_contract(mod_b, n_rep, length(mod$coefficients))
  expect_equal(mod_b$uq_control$control$mode, "free_smoothing")
})

test_that("bootstrap_gbgam free_stopping mode satisfies structural contract", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 112)
  mod <- gbgam(
    y ~ t2(x, k = 6, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.05),
    stop_control    = gb_stop_ctrl("aic", max_iter = 80L),
    edf_control     = gb_edf_ctrl("hutchpp"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  set.seed(77)
  n_rep <- 15L
  # the replicates reach max_iter = 80 before AIC stops them (reported as
  # censored, with a warning)
  expect_warning(
    mod_b <- bootstrap_gbgam(mod, n_samples = n_rep, alpha = NULL,
                             mode = "free_stopping"),
    "reached max_iter")
  check_boot_contract(mod_b, n_rep, length(mod$coefficients))
  expect_equal(mod_b$uq_control$control$mode, "free_stopping")
})

test_that("bias corrections and smooth_intervals() work under na.exclude", {
  dat <- simulate_gaussian_2d(n = 120, seed = 21)
  dat$x[1:4] <- NA
  fit <- function(na_action, ...) {
    gbgam(y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 6, bs = "cr"), data = dat,
          na.action = na_action, stop_control = gb_stop_ctrl(max_iter = 300),
          verbose = FALSE, ...)
  }
  m_ex <- fit(stats::na.exclude)
  m_om <- fit(stats::na.omit)
  # 8 replicates leave some learner factors below 0.9 (Monte Carlo noise)
  boot <- function(m, corr) {
    without_low_widening_warning(
      bootstrap_gbgam(m, n_samples = 8, correction = corr))
  }
  for (corr in c("model", "smooth", "learner")) {
    set.seed(3)
    b_ex <- boot(m_ex, corr)
    set.seed(3)
    b_om <- boot(m_om, corr)
    expect_equal(b_ex$Vp, b_om$Vp, tolerance = 1e-10, label = corr)
    ci_ex <- smooth_intervals(b_ex, alpha = 0.1)
    ci_om <- smooth_intervals(b_om, alpha = 0.1)
    expect_equal(ci_ex$factor, ci_om$factor, tolerance = 1e-10, label = corr)
    # Without newdata the intervals are aligned with the data, as fitted()
    expect_equal(nrow(ci_ex$fit), nrow(dat))
    expect_true(all(is.na(ci_ex$lower[1:4, ])))
    expect_equal(ci_ex$upper[-(1:4), ], ci_om$upper, tolerance = 1e-10,
                 ignore_attr = TRUE)
  }
  # Inside gbgam(), a failing correction would lose the whole fit
  set.seed(3)
  m_in <- fit(stats::na.exclude, uq_control = gb_uq_ctrl(
    "bootstrap", n_samples = 8, correction = "smooth"))
  expect_false(is.null(m_in$Vp))
})

test_that("bootstrap_gbgam() removes Bayesian backfitting results", {
  dat <- simulate_gaussian_1d(n = 80, seed = 22)
  mod <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl(max_iter = 300))
  bb  <- bayesian_backfitting(mod, n_samples = 20, n_burnin = 5)
  expect_false(is.null(bb$coefficients_bayes))
  set.seed(1)
  b   <- bootstrap_gbgam(bb, n_samples = 4)
  expect_null(b$coefficients_bayes)
  expect_null(b$beta_samples)
  expect_null(b$bb_info)
  expect_equal(dim(b$boot_coef), c(4L, length(mod$coefficients)))
  expect_error(smooth_intervals(b, centre = "posterior"), "needs Bayesian")
})


# ---- Argument checks and early returns ---------------------------------------

test_that("bootstrap_gbgam requires n_rep >= 2 and handles < 2 completed", {
  dat <- simulate_poisson_1d(n = 80L, seed = 51L)
  mod <- gbgam(
    y ~ t2(x, k = 5),
    data         = dat,
    family       = poisson(),
    stop_control = gb_stop_ctrl("none", max_iter = 10L),
    edf_control  = gb_edf_ctrl("none"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  expect_error(bootstrap_gbgam(mod, n_samples = 1L), "n_samples")
  w <- testthat::capture_warnings(
    res <- bootstrap_gbgam(mod, n_samples = 5L, max_runtime_mins = 0)
  )
  expect_true(any(grepl("time limit reached after 0 / 5", w)))
  expect_true(any(grepl("fewer than 2 replicates", w)))
  expect_identical(res, mod)
})

test_that("an early stop keeps earlier results, or records the replicates", {
  dat <- simulate_gaussian_1d(n = 80, seed = 23)
  mod <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl(max_iter = 100))
  set.seed(1)
  bb <- bayesian_backfitting(mod, n_samples = 10L, n_burnin = 2L)
  expect_identical(bb$n_samples_done, 10L)
  # Too few replicates: the posterior and its settings stay
  expect_warning(
    expect_warning(
      out <- bootstrap_gbgam(bb, n_samples = 5L, max_runtime_mins = 0),
      "time limit reached after 0 / 5"),
    "returning the model unchanged")
  expect_identical(out, bb)
  # Some replicates: they are used, and counted; the settings are those run
  real_fun <- replicate_function
  local_mocked_bindings(replicate_function = function(...) {
    one_rep <- real_fun(...)
    function(b) {
      r <- one_rep(b)
      if (b == 3L) Sys.sleep(4)  # the 3 s limit is passed after replicate 3
      r
    }
  })
  set.seed(2)
  expect_warning(b <- bootstrap_gbgam(bb, n_samples = 30L,
                                      max_runtime_mins = 0.05),
                 "time limit reached after 3 / 30")
  expect_identical(b$n_samples_done, 3L)
  expect_identical(nrow(b$boot_coef), b$n_samples_done)
  expect_identical(b$uq_control$control$n_samples, 30L)
  expect_null(b$coefficients_bayes)
})

test_that("the bootstrap progress bar reports the replicates run", {
  dat <- simulate_gaussian_1d(n = 60, seed = 24)
  mod <- fit_fixed(y ~ t2(x, k = 5, bs = "cr"), dat, max_iter = 10)
  real_fun <- replicate_function
  local_mocked_bindings(replicate_function = function(...) {
    one_rep <- real_fun(...)
    function(b) {
      r <- one_rep(b)
      if (b == 3L) Sys.sleep(0.8)  # the 0.6 s limit is passed here
      r
    }
  })
  op <- options(cli.dynamic = FALSE, cli.unicode = FALSE,
                cli.progress_show_after = 0)
  on.exit(options(op), add = TRUE)
  msgs <- character(0)
  suppressWarnings(withCallingHandlers(
    b <- bootstrap_gbgam(mod, n_samples = 30L, max_runtime_mins = 0.01,
                         verbose = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }))
  expect_identical(b$n_samples_done, 3L)
  out <- paste(msgs, collapse = "")
  expect_match(out, "Bootstrap +3 replicates in")
  expect_no_match(out, "30 replicates")
})

test_that("bootstrap_gbgam rejects alpha outside (0, 0.5]", {
  dat <- simulate_gaussian_1d(n = 60L, seed = 111L)
  mod <- gbgam(
    y ~ t2(x, k = 5),
    data         = dat,
    stop_control = gb_stop_ctrl("none", max_iter = 5L),
    edf_control  = gb_edf_ctrl("exact"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  expect_error(bootstrap_gbgam(mod, n_samples = 2L, alpha = 0), "alpha")
  expect_error(bootstrap_gbgam(mod, n_samples = 2L, alpha = 1), "alpha")
  expect_error(bootstrap_gbgam(mod, n_samples = 2L, alpha = 0.95),
               "alpha is the miscoverage")
  expect_error(gb_param_boot_ctrl(alpha = 0.95), "alpha is the miscoverage")
  expect_identical(gb_param_boot_ctrl(alpha = 0.5)$alpha, 0.5)
})

test_that("bootstrap with mstop = 0 warns and keeps the fit", {
  set.seed(1)
  d <- data.frame(x = stats::runif(60), y = stats::rnorm(60))
  expect_warning(
    fit <- gbgam(y ~ t2(x, k = 6), data = d,
                 stop_control = gb_stop_ctrl(iter_no_improv = 20),
                 uq_control   = gb_uq_ctrl("bootstrap",
                                                       n_samples = 5),
                 verbose           = FALSE),
    "mstop = 0"
  )
  expect_equal(fit$mstop, 0L)
  expect_null(fit$Vp)
  expect_true(all(is.finite(fit$linear.predictors)))
})

test_that("a widening factor just below 1 does not warn", {
  # Setting of the bootstrap_gbgam() example, where n_samples = 20 gives
  # c = 0.99 (Monte Carlo noise): no warning about a poor model fit
  set.seed(1)
  n  <- 200
  x  <- stats::runif(n)
  df <- data.frame(x = x, y = sin(2 * pi * x) + stats::rnorm(n, sd = 0.4))
  mod <- gbgam(
    y ~ t2(x, k = 8, bs = "cr"),
    data       = df,
    uq_control = gb_uq_ctrl("none"),
    verbose    = FALSE
  )
  w <- testthat::capture_warnings(b <- bootstrap_gbgam(mod, n_samples = 20))
  expect_gt(b$widening_factor, 0.9)
  expect_false(any(grepl("widening factor", w)))
})

test_that("widening factors far from 1 warn, naming the groups", {
  # model: one factor, the two messages
  expect_warning(warn_widening_factors(3.5, 50L, "model"),
                 "widening factor = 3.50 > 3: this typically")
  expect_warning(warn_widening_factors(0.5, 50L, "model"),
                 "widening factor = 0.50 < 0.9: .*n_samples = 50")
  expect_silent(warn_widening_factors(1.2, 50L, "model"))
  # learners: one warning per side, naming the groups beyond it
  fac <- c("(Intercept)" = 1, "t2(x)r" = 4.1, "t2(x)n1" = 1.3,
           "t2(z)r" = 0.4, "t2(z)n1" = 3.52)
  w <- testthat::capture_warnings(warn_widening_factors(fac, 50L, "learner"))
  expect_length(w, 2L)
  expect_match(w[1], paste0("widening factors > 3 for 2 of 5 learners: ",
                            "t2(x)r (4.10), t2(z)n1 (3.52); this typically"),
               fixed = TRUE)
  expect_match(w[2], paste0("widening factors < 0.9 for 1 of 5 learners: ",
                            "t2(z)r (0.40); this may indicate a poor model ",
                            "fit (or too few replicates, n_samples = 50)"),
               fixed = TRUE)
  # within [0.9, 3], or exactly 1 (no bootstrap variation): no warning
  expect_silent(warn_widening_factors(fac[c(1, 3)], 50L, "smooth"))
  # more than 5 groups beyond a threshold: the 5 farthest, then the count
  many <- setNames(c(3.1, 9, 4, 5, 6, 7, 8), paste0("s", 1:7))
  w <- testthat::capture_warnings(warn_widening_factors(many, 50L, "smooth"))
  expect_length(w, 1L)
  expect_match(w, paste0("7 of 7 smooths and parametric coefficients: ",
                         "s2 (9.00), s7 (8.00), s6 (7.00), s5 (6.00), ",
                         "s4 (5.00) and 2 more"), fixed = TRUE)
})

test_that("bootstrap_gbgam() names the learner whose widening factor is large", {
  set.seed(2)
  n  <- 120
  df <- data.frame(x = stats::runif(n), z = stats::runif(n))
  df$y <- sin(2 * pi * df$x) + df$z + stats::rnorm(n, sd = 0.4)
  mod <- gbgam(y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 6, bs = "cr"),
               data = df, uq_control = gb_uq_ctrl("none"), verbose = FALSE,
               stop_control = gb_stop_ctrl("none", max_iter = 30))
  # factor 4 for the second learner only, 1.5 for the others
  calls <- 0L
  local_mocked_bindings(wood_widening_factor = function(...) {
    calls <<- calls + 1L
    if (calls == 2L) 4 else 1.5
  })
  w <- testthat::capture_warnings(
    b <- bootstrap_gbgam(mod, n_samples = 5, correction = "learner"))
  lab <- names(b$widening_factor)[2L]
  expect_identical(unname(b$widening_factor[[2L]]), 4)
  expect_length(grep("widening factors", w), 1L)
  expect_match(w[grep("widening factors", w)],
               sprintf("> 3 for 1 of %d learners: %s (4.00)",
                       length(mod$learners), lab), fixed = TRUE)
})

test_that("replicates without variation give a model widening factor of 1", {
  # Every replicate response equals the fitted means, so no fitted mean
  # varies: the factor is 1 at fit time and in smooth_intervals()
  set.seed(1)
  n  <- 100
  x  <- stats::runif(n)
  df <- data.frame(x = x, y = sin(2 * pi * x) + stats::rnorm(n, sd = 0.4))
  mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = df,
               uq_control = gb_uq_ctrl("none"),
               verbose = FALSE)
  local_mocked_bindings(simulate_response = function(family, mu, ...) mu)
  b <- bootstrap_gbgam(mod, n_samples = 5)
  expect_equal(b$widening_factor, 1)
  expect_equal(unname(smooth_intervals(b)$factor), 1)
})

test_that("free modes need EDF tracking and phi must be positive", {
  skip_on_cran()
  set.seed(14)
  d <- data.frame(x = stats::runif(150))
  f <- sin(2 * pi * d$x)
  d$y <- f + stats::rnorm(150, 0, 0.3)
  mod <- gbgam(y ~ t2(x, bs = "ps", k = 8), data = d,
               stop_control = gb_stop_ctrl("true_f", true_f = f,
                                           iter_no_improv = 100L),
               edf_control = gb_edf_ctrl("none"), verbose = FALSE)
  set.seed(15)
  b <- bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL, phi = 0.09)
  expect_length(b$boot_mstop, 3L)
  for (md in c("free_smoothing", "free_selection", "free_stopping"))
    expect_error(bootstrap_gbgam(mod, n_samples = 3L, mode = md),
                 paste0('mode = "', md, '" requires EDF tracking'),
                 fixed = TRUE)
  expect_error(bootstrap_gbgam(mod, n_samples = 3L, phi = 0), "positive")
  expect_error(gb_uq_ctrl("bootstrap", phi = 0), "positive")
})

test_that("bootstrap_gbgam warns when no residual df are left for phi", {
  dat <- simulate_gaussian_1d(n = 80L, seed = 71L)
  mod <- gbgam(
    y ~ t2(x, k = 5),
    data         = dat,
    weights      = rep(1 / 80, 80),
    stop_control = gb_stop_ctrl("none", max_iter = 10L),
    edf_control  = gb_edf_ctrl("exact"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  expect_lte(sum(mod$prior.weights) - mod$edf.total, 0)
  set.seed(1)
  expect_warning(bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL),
                 "no residual degrees of freedom")
  set.seed(1)
  expect_no_warning(bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL,
                                    phi = 0.1))
})

test_that("bootstrap_gbgam warns about dispersion only for free-dispersion families", {
  dat <- simulate_gaussian_1d(n = 80L, seed = 71L)
  mod <- gbgam(
    y ~ t2(x, k = 5),
    data         = dat,
    stop_control = gb_stop_ctrl("none", max_iter = 10L),
    edf_control  = gb_edf_ctrl("none"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  set.seed(1)
  expect_warning(bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL),
                 "edf_compute_method")

  pdat <- simulate_poisson_1d(n = 80L, seed = 72L)
  pmod <- gbgam(
    y ~ t2(x, k = 5),
    data         = pdat,
    family       = poisson(),
    stop_control = gb_stop_ctrl("none", max_iter = 10L),
    edf_control  = gb_edf_ctrl("none"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  set.seed(1)
  w <- testthat::capture_warnings(
    bootstrap_gbgam(pmod, n_samples = 3L, alpha = NULL)
  )
  expect_false(any(grepl("edf_compute_method", w)))
})


# ---- Modes and corrections ---------------------------------------------------

test_that("each mode sets the replicates' overrides", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  rec <- recorded_paths(mod, mod$mstop)
  set <- function(mode, max_iter = NULL)
    replicate_settings(mode, rec$path, rec$rhos, mod$mstop, mod,
                       post_rhos = rec$post_rhos, max_iter = max_iter)
  s <- lapply(setNames(nm = names(boot_modes)), set)
  expect_identical(s$fixed$edf_compute_method, "none")
  expect_identical(s$fixed$preselected_rhos, rec$rhos)
  for (md in c("free_smoothing", "free_selection", "free_stopping")) {
    expect_identical(s[[md]]$edf_compute_method, mod$edf_compute_method)
    expect_null(s[[md]]$preselected_rhos)
  }
  expect_identical(s$free_smoothing$preselected_learner_indices, rec$path)
  expect_null(s$free_selection$preselected_learner_indices)
  for (md in c("fixed", "free_smoothing", "free_selection")) {
    expect_identical(s[[md]]$stopping_method, "none")
    expect_identical(s[[md]]$max_iter, mod$mstop)
  }
  expect_identical(s$free_stopping$stopping_method, mod$stopping_method)
  expect_identical(s$free_stopping$max_iter, mod$max_iter)
  expect_identical(s$free_stopping$iter_no_improv, mod$iter_no_improv)
  expect_identical(s$free_stopping$iter_no_improv_rel,
                   mod$iter_no_improv_rel)
  for (md in c("fixed", "free_smoothing", "free_selection")) {
    expect_null(s[[md]]$iter_no_improv)
    expect_null(s[[md]]$iter_no_improv_rel)
  }
  expect_identical(set("free_stopping", max_iter = 7L)$max_iter, 7L)
  # An extended fit that had no cap gives its replicates none
  mod$max_iter <- 500L
  mod$uncapped <- TRUE
  expect_null(set("free_stopping")$max_iter)
  expect_identical(set("free_stopping", max_iter = 7L)$max_iter, 7L)
})

test_that("free_stopping stops on a model without iter_no_improv_rel", {
  dat <- simulate_gaussian_1d(n = 60, seed = 1)
  old <- fit_1d(dat, stop_control = gb_stop_ctrl(max_iter = 10))
  old$iter_no_improv_rel <- NULL
  expect_error(bootstrap_gbgam(old, 2L, mode = "free_stopping", phi = 1),
               "older gbgam.*refit it")
})

test_that("free_stopping replicates stop by the inherited patience", {
  dat <- simulate_gaussian_1d(n = 120, seed = 1)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(iter_no_improv = 3,
                                                 iter_no_improv_rel = 0.5))
  rec <- recorded_paths(mod, mod$mstop)
  s   <- replicate_settings("free_stopping", rec$path, rec$rhos,
                            mod$mstop, mod, post_rhos = rec$post_rhos)
  set.seed(3)
  y_new <- as.vector(fitted(mod)) + stats::rnorm(120, sd = 0.25)
  b   <- refit_replicate(mod, y_new, s)
  run <- sum(!is.na(b$metrics$k_opts))
  cat("bootstrap replicate: mstop", b$mstop, "run", run, "\n")
  expect_identical(b$iter_no_improv_rel, 0.5)
  expect_true(b$converged)
  expect_gt(ceiling(0.5 * b$mstop), 3)
  expect_identical(run, as.integer(b$mstop + ceiling(0.5 * b$mstop)))
})

test_that("free_stopping bootstrap works for Gamma / IG with zero weights", {
  set.seed(91)
  n   <- 150L
  dat <- data.frame(x = stats::runif(n))
  mu  <- exp(0.5 + sin(2 * pi * dat$x))
  dat$y_ga <- stats::rgamma(n, shape = 5, rate = 5 / mu)
  dat$y_ig <- mu * stats::rgamma(n, shape = 20, rate = 20)
  w   <- c(rep(0, 5), rep(1, n - 5))
  for (spec in list(list("y_ga", stats::Gamma(link = "log")),
                    list("y_ig", stats::inverse.gaussian(link = "log")))) {
    # (quiet: non-canonical link message)
    mod <- suppressMessages(gbgam(
      stats::as.formula(paste(spec[[1]], "~ t2(x, k = 8)")),
      data              = dat,
      weights           = w,
      family            = spec[[2]],
      learner_control   = gb_learner_ctrl(step_size = 0.3),
      stop_control      = gb_stop_ctrl("aic", iter_no_improv = 10L,
                                       max_iter = 100L),
      edf_control = gb_edf_ctrl("exact"),
      uq_control  = gb_uq_ctrl("none"),
      verbose     = FALSE
    ))
    expect_gt(mod$mstop, 0L)
    set.seed(4)
    # every replicate stops by AIC before max_iter = 100, so none is censored
    expect_no_warning(
      b <- bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL,
                           mode = "free_stopping"))
    expect_equal(dim(b$boot_coef), c(3L, length(coef(mod))))
    expect_true(all(is.finite(b$boot_coef)))
    expect_true(all(b$boot_mstop < 100L))
    expect_false(any(b$boot_censored))
  }
})

test_that("bootstrap replicates that reach max_iter are recorded", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  set.seed(12)
  # replicates cannot stop by AIC within 5 iterations
  expect_warning(
    b <- bootstrap_gbgam(mod, n_samples = 4L, alpha = NULL,
                         mode = "free_stopping", max_iter = 5L),
    "reached max_iter")
  expect_true(all(b$boot_censored))
  expect_true(all(b$boot_mstop <= 5L))
  expect_identical(b$uq_control$control$max_iter, 5L)
  # with room to stop (mstop 165 to 324 measured), none is censored
  set.seed(12)
  expect_no_warning(
    b_free <- bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL,
                              mode = "free_stopping", max_iter = 5000L))
  expect_true(all(b_free$boot_mstop < 5000L))
  expect_false(any(b_free$boot_censored))
  # max_iter is the same as the main fit's max_iter set to it
  mod5 <- mod
  mod5$max_iter <- 5L
  set.seed(12)
  b5 <- suppressWarnings(bootstrap_gbgam(mod5, n_samples = 4L, alpha = NULL,
                                         mode = "free_stopping"))
  expect_identical(b5$boot_coef, b$boot_coef)
  expect_identical(b5$boot_mstop, b$boot_mstop)
  expect_error(bootstrap_gbgam(mod, n_samples = 4L, max_iter = 5L),
               "only used with mode")
  set.seed(12)
  b2 <- bootstrap_gbgam(fit_aic(f_col, sim_collinear()), n_samples = 3L,
                        alpha = NULL)
  expect_false(any(b2$boot_censored))       # replayed paths are never capped
})

test_that("fixed-mode Gaussian replicates have covariance phi A W^-1 A^T", {
  # For a Gaussian identity fit a fixed-mode replicate is linear in its
  # response, coef^(r) = A y^(r); column i of A is the replay on e_i. This
  # checks the simulated responses and the linearity of the replay; the
  # replay itself is checked against the main fit only through A y = coef
  dat <- simulate_gaussian_1d(n = 60, seed = 98)
  w   <- rep(c(0.25, 1, 4), 20)
  mod <- fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 20L,
                   weights = w, uq_control = gb_uq_ctrl("none"))
  n <- nrow(dat)
  A <- vapply(seq_len(n), function(i) {
    replay_on_y(mod, diag(n)[, i], "fixed")$coefficients
  }, numeric(length(mod$coefficients)))
  expect_equal(as.vector(A %*% mod$y), unname(mod$coefficients),
               tolerance = 1e-10)
  phi <- 0.09
  B   <- 200L
  set.seed(7)
  b <- bootstrap_gbgam(mod, n_samples = B, alpha = NULL, phi = phi)
  # the replicate responses are N(mu~, phi / w) draws, in replicate order
  set.seed(7)
  Y <- vapply(seq_len(B), function(r) {
    stats::rnorm(n, fitted(mod), sqrt(phi / w))
  }, numeric(n))
  expect_equal(b$boot_coef, t(A %*% Y), tolerance = 1e-10,
               ignore_attr = TRUE)
  # Vp against the exact covariance, through the summed variance of the
  # fitted values: Monte Carlo relative SD 0.08 at B = 200 (from the
  # eigenvalues of the exact covariance), so the tolerance is 3.7 SD; ratios
  # 0.97-1.02 over 4 seeds; simulating without the weights gives 3.25
  V_exact <- phi * A %*% (t(A) / w)
  X <- stats::predict(mod, type = "lpmatrix")
  ratio <- sum(diag(X %*% b$Vp %*% t(X))) / sum(diag(X %*% V_exact %*% t(X)))
  expect_lt(abs(ratio - 1), 0.3)
})

# Share of replicate deviations t^(r) - t0 outside z * factor * sd at level
# 0.95, over observations and replicates
share_outside <- function(t0, tr, factor) {
  s <- apply(tr, 2L, stats::sd)
  mean(abs(sweep(tr, 2L, t0, "-")) >
         stats::qnorm(0.975) * factor * rep(s, each = nrow(tr)))
}

test_that("the model factor is a quantile of the scaled fitted means", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  set.seed(8)
  b <- suppressWarnings(bootstrap_gbgam(mod, n_samples = 30L))
  expect_identical(b$uq_control$control$correction, "model")
  # smooth_intervals() at the bootstrap's alpha finds the same factor
  expect_identical(smooth_intervals(b)$factor, b$widening_factor)
  expect_equal(b$Vp, b$widening_factor^2 * stats::cov(b$boot_coef))
  # the stretched intervals around the main fit's fitted means hold a share
  # 1 - alpha of the replicates' fitted means, up to ties and interpolation
  X    <- stats::predict(b, type = "lpmatrix")
  mu_r <- b$family$linkinv(b$boot_coef %*% t(X))
  out  <- share_outside(b$family$linkinv(as.vector(X %*% b$coefficients)),
                        mu_r, b$widening_factor)
  expect_lt(abs(out - 0.05), 0.01)
})

test_that("per-smooth and per-learner corrections widen by group", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  set.seed(8)
  b_s <- bootstrap_gbgam(mod, n_samples = 30L, correction = "smooth")
  set.seed(8)
  b_l <- bootstrap_gbgam(mod, n_samples = 30L, correction = "learner")
  expect_named(b_s$widening_factor,
               c("(Intercept)", "t2(x0)", "t2(x1)", "t2(x2)"))
  expect_length(b_l$widening_factor, length(mod$learners))
  boots <- list(smooth = b_s, learner = b_l)
  for (nm in names(boots)) {
    b <- boots[[nm]]
    expect_true(isSymmetric(unname(b$Vp)))
    expect_true(all(eigen(b$Vp, symmetric = TRUE, only.values = TRUE)$values >
                      -1e-10))
    expect_identical(b$uq_control$control$correction, nm)
  }
  # the factor of a group puts the share alpha of replicate deviations
  # outside z * c * sd; the factor is an empirical quantile of the scaled
  # deviations, so the share is 0.05 only up to ties and interpolation
  X  <- stats::predict(b_s, type = "lpmatrix")
  sm <- b_s$smooth[[1L]]
  i  <- sm$first.para:sm$last.para
  out <- share_outside(as.vector(X[, i] %*% b_s$coefficients[i]),
                       b_s$boot_coef[, i] %*% t(X[, i]),
                       b_s$widening_factor[["t2(x0)"]])
  expect_lt(abs(out - 0.05), 0.01)
  # Vp is D C D, with C the replicates' covariance and D the diagonal of the
  # factors of each coefficient's group (the intercept's for column 1); the
  # factors differ by 14%, so one factor for all would not match
  fac <- b_s$widening_factor
  d   <- rep(fac[["(Intercept)"]], length(mod$coefficients))
  for (s in b_s$smooth) d[s$first.para:s$last.para] <- fac[[s$label]]
  V_ref <- stats::cov(b_s$boot_coef) * outer(d, d)
  expect_gt(diff(range(fac)) / min(fac), 0.05)
  expect_equal(b_s$Vp, V_ref, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(b_s$Vp[1L, 1L], fac[["(Intercept)"]]^2 *
                 stats::var(b_s$boot_coef[, 1L]), tolerance = 1e-10)
  # smooth_intervals() at the bootstrap's alpha: se of a + f_j from V_ref
  ci <- smooth_intervals(b_s)
  for (s in b_s$smooth) {
    i  <- c(1L, s$first.para:s$last.para)
    se <- sqrt(rowSums((X[, i] %*% V_ref[i, i]) * X[, i]))
    expect_equal(ci$se[, s$label], se, tolerance = 1e-10, ignore_attr = TRUE)
  }

  # learner correction: the replicate learner coefficients rebuild the
  # replicate coefficients, and a learner's factor has the same property
  lr <- b_l$learners
  for (k in seq_along(lr)) lr[[k]]$beta <- b_l$boot_learner_coef[[k]][3L, ]
  expect_equal(unname(learners_to_mgcv_coef(lr, b_l$coef_names, b_l$ind0var,
                                            b_l$cen_const,
                                            cen_absorb = b_l$cen_absorb)),
               unname(b_l$boot_coef[3L, ]))
  k  <- which(names(b_l$widening_factor) == "t2(x0)r")
  Xk <- b_l$learners[[k]]$design
  out <- share_outside(as.vector(Xk %*% b_l$learners[[k]]$beta),
                       b_l$boot_learner_coef[[k]] %*% t(Xk),
                       b_l$widening_factor[[k]])
  expect_lt(abs(out - 0.05), 0.01)

  # every learner factor set to 2 doubles each replicate learner's deviation
  # from the main fit, so Vp is 4 times the replicates' covariance
  local_mocked_bindings(wood_widening_factor = function(...) 2)
  set.seed(8)
  b_2 <- bootstrap_gbgam(mod, n_samples = 10L, correction = "learner")
  expect_equal(unname(b_2$widening_factor), rep(2, length(mod$learners)))
  expect_equal(b_2$Vp, 4 * stats::cov(b_2$boot_coef), tolerance = 1e-10,
               ignore_attr = TRUE)
})
