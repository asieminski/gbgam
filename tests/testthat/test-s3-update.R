# update.gbgam(): refits with changed arguments and uncertainty
# quantification.


# ---- Refits ------------------------------------------------------------------

test_that("update.gbgam re-fits with a changed stop_control", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 70, seed = 51)
  mod1 <- gbgam(
    y ~ t2(x, k = 6, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(),
    stop_control    = gb_stop_ctrl("none", max_iter = 10),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  mod2 <- update(mod1, stop_control = gb_stop_ctrl("none", max_iter = 20))
  expect_s3_class(mod2, "gbgam")
  expect_equal(mod2$mstop, 20L)
})

test_that("update.gbgam with evaluate=FALSE returns a call", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 52)
  mod <- gbgam(
    y ~ t2(x, k = 6, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(),
    stop_control    = gb_stop_ctrl("none", max_iter = 10),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  cl <- update(mod,
    stop_control = gb_stop_ctrl("none", max_iter = 5),
    evaluate = FALSE)
  expect_true(is.call(cl))
})

test_that("update() refits correctly after NA removal or a subset", {
  set.seed(208)
  xx <- stats::runif(100)
  yy <- sin(2 * pi * xx) + stats::rnorm(100, sd = 0.3)
  wv <- stats::runif(100, 0.5, 2)
  xx[c(3, 7)] <- NA
  m1 <- gbgam(yy ~ t2(xx, k = 6), weights = wv, na.action = stats::na.omit,
              stop_control = sc_ff, edf_control = ec_ff, verbose = FALSE)
  expect_equal(fitted(update(m1)), fitted(m1))
  x2 <- stats::runif(100)
  y2 <- sin(2 * pi * x2) + stats::rnorm(100, sd = 0.3)
  m2 <- gbgam(y2 ~ t2(x2, k = 6), subset = 41:100, stop_control = sc_ff,
              edf_control = ec_ff, verbose = FALSE)
  expect_equal(nobs(update(m2)), 60L)
  expect_equal(fitted(update(m2)), fitted(m2))
  # Without NA removal or subset, the model frame still stands in for data
  m3 <- gbgam(y2 ~ t2(x2, k = 6), stop_control = sc_ff, edf_control = ec_ff,
              verbose = FALSE)
  expect_identical(update(m3, evaluate = FALSE)$data, m3$model)
})


# ---- Uncertainty quantification ----------------------------------------------

test_that("update.gbgam re-runs the bootstrap when refitting a bootstrapped model", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 53)
  mod <- gbgam(
    y ~ t2(x, k = 5, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.05),
    stop_control    = gb_stop_ctrl("none", max_iter = 10),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  set.seed(1)
  mod_b <- suppressWarnings(bootstrap_gbgam(mod, n_samples = 15))
  mod_b2 <- suppressWarnings(
    update(mod_b, stop_control = gb_stop_ctrl("none", max_iter = 12))
  )
  expect_s3_class(mod_b2, "gbgam")
  expect_equal(mod_b2$mstop, 12L)
  expect_false(is.null(mod_b2$Vp))
  expect_equal(dim(mod_b2$Vp), rep(length(mod_b2$coefficients), 2))
  expect_equal(mod_b2$uq_control$method, "bootstrap")
  expect_equal(mod_b2$uq_control$control$n_samples, 15L)
})

test_that("update.gbgam preserves the bootstrap from gbgam(uq_control=)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 54)
  mod_b <- suppressWarnings(gbgam(
    y ~ t2(x, k = 5, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.05),
    stop_control    = gb_stop_ctrl("none", max_iter = 10),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("bootstrap", n_samples = 15L),
    verbose         = FALSE
  ))
  mod_b2 <- suppressWarnings(
    update(mod_b, stop_control = gb_stop_ctrl("none", max_iter = 12))
  )
  expect_s3_class(mod_b2, "gbgam")
  expect_equal(mod_b2$mstop, 12L)
  expect_false(is.null(mod_b2$Vp))
  expect_equal(mod_b2$uq_control$method, "bootstrap")
})

test_that("update.gbgam preserves Bayesian backfitting (Gaussian)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 55)
  mod <- gbgam(
    y ~ t2(x, k = 5, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.05, tau = 0),
    stop_control    = gb_stop_ctrl("none", max_iter = 10),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  set.seed(1)
  mod_bb <- bayesian_backfitting(mod, n_samples = 50L, n_burnin = 20L)
  mod_bb2 <- bayesian_backfitting(
    update(mod_bb, stop_control = gb_stop_ctrl("none", max_iter = 12)),
    n_samples = 50L, n_burnin = 20L
  )
  expect_s3_class(mod_bb2, "gbgam")
  expect_equal(mod_bb2$mstop, 12L)
  expect_false(is.null(mod_bb2$Vp))
})

test_that("update.gbgam uq_control can be overridden to 'none'", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 56)
  mod_b <- suppressWarnings(gbgam(
    y ~ t2(x, k = 5, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.05),
    stop_control    = gb_stop_ctrl("none", max_iter = 10),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("bootstrap", n_samples = 15L),
    verbose         = FALSE
  ))
  mod2 <- update(mod_b, uq_control = gb_uq_ctrl("none"))
  expect_null(mod2$Vp)
  expect_equal(mod2$uq_control$method, "none")
})
