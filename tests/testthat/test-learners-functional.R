# Linear functional terms (matrix covariates, as in distributed-lag models;
# see ?mgcv::linear.functional.terms) as t2() learners.

# ---- 2D functional term: t2(x, lag) ------------------------------------------

test_that("2D functional predictor: model fits without warnings", {
  skip_on_cran()
  dat <- simulate_functional_dlm(n = 200, n_lag = 6, seed = 40)

  expect_no_warning(
    mod <- gbgam(
      y ~ t2(x, lag, k = c(6, 5), bs = c("cr", "cr")),
      family            = stats::poisson(),
      data              = dat,
      stop_control      = gb_stop_ctrl("none", max_iter = 10),
      edf_control       = gb_edf_ctrl("none"),
      verbose           = FALSE
    )
  )

  expect_model_contract(mod, nrow(dat$x))
  expect_true(all(is.na(mod$edf)))
})


test_that("2D functional predictor: coefficient reconstruction reproduces training eta", {
  skip_on_cran()
  dat <- simulate_functional_dlm(n = 200, n_lag = 6, seed = 41)

  mod <- gbgam(
    y ~ t2(x, lag, k = c(6, 5), bs = c("cr", "cr")),
    family            = stats::poisson(),
    data              = dat,
    stop_control      = gb_stop_ctrl("none", max_iter = 15),
    edf_control       = gb_edf_ctrl("none"),
    verbose           = FALSE
  )

  eta_rebuilt <- rebuild_training_eta(mod)
  expect_equal(eta_rebuilt, as.vector(mod$linear.predictors), tolerance = 1e-8)
})


test_that("2D functional predictor: EDF vector has correct length and non-negative entries", {
  skip_on_cran()
  dat <- simulate_functional_dlm(n = 200, n_lag = 6, seed = 42)

  # EDF is tracked (with edf_compute_method = "none" every entry is NA)
  mod <- gbgam(
    y ~ t2(x, lag, k = c(6, 5), bs = c("cr", "cr")),
    family            = stats::poisson(),
    data              = dat,
    stop_control      = gb_stop_ctrl("none", max_iter = 10),
    edf_control       = gb_edf_ctrl("exact"),
    verbose           = FALSE
  )

  # edf must parallel coefficients
  expect_equal(length(mod$edf), length(mod$coefficients))
  expect_true(all(mod$edf >= 0))
  # Dropped (zero-variance) columns get zero EDF; this basis has five
  expect_length(mod$ind0var, 5L)
  expect_true(all(mod$edf[mod$ind0var] == 0))
})


# ---- 3D functional term: t2(x, z, lag) ---------------------------------------

test_that("3D functional predictor: model fits without warnings", {
  skip_on_cran()
  dat <- simulate_functional_dlm(n = 250, n_lag = 5, seed = 50)

  # Add a second matrix covariate by adding Gaussian noise to x
  set.seed(51)
  dat$z <- dat$x + matrix(stats::rnorm(length(dat$x), sd = 0.5),
                           nrow = nrow(dat$x), ncol = ncol(dat$x))

  expect_no_warning(
    mod <- gbgam(
      y ~ t2(x, z, lag, k = c(4, 4, 4), bs = c("cr", "cr", "cr")),
      family            = stats::poisson(),
      data              = dat,
      stop_control      = gb_stop_ctrl("none", max_iter = 10),
      edf_control       = gb_edf_ctrl("none"),
      verbose           = FALSE
    )
  )

  expect_model_contract(mod, nrow(dat$x))
})


# ---- Mixed formula: scalar smooth + functional term --------------------------

test_that("Mixed formula with scalar and functional term fits without warnings", {
  skip_on_cran()
  set.seed(60)
  n     <- 200L
  n_lag <- 5L
  dlm   <- simulate_functional_dlm(n = n, n_lag = n_lag, seed = 60)

  # Add a plain scalar covariate
  dat        <- dlm
  dat$scalar <- stats::rnorm(nrow(dlm$x))

  expect_no_warning(
    mod <- gbgam(
      y ~ t2(scalar, k = 6, bs = "cr") +
          t2(x, lag, k = c(5, 4), bs = c("cr", "cr")),
      family            = stats::poisson(),
      data              = dat,
      stop_control      = gb_stop_ctrl("none", max_iter = 10),
      edf_control       = gb_edf_ctrl("none"),
      verbose           = FALSE
    )
  )

  expect_model_contract(mod, nrow(dat$x))
  eta_rebuilt <- rebuild_training_eta(mod)
  expect_equal(eta_rebuilt, as.vector(mod$linear.predictors), tolerance = 1e-8)
})


# ---- partition_t2_basis_matrix invariant -------------------------------------

test_that("partition_t2_basis_matrix excludes zero-variance columns from all t2 basis partitions", {
  # With matrix covariates, zero-variance columns can lie in penalised
  # partitions too, not only in the null space; all must be dropped
  skip_on_cran()
  dat <- simulate_functional_dlm(n = 200, n_lag = 6, seed = 70)

  gp  <- mgcv::interpret.gam(y ~ t2(x, lag, k = c(6, 5), bs = c("cr", "cr")))
  mf  <- stats::model.frame(gp$fake.formula, dat)
  sms <- construct_smooths(mf, gp$smooth.spec)
  sm  <- sms[[1]]

  weights <- rep(1, nrow(mf))
  w_means <- apply(sm$X, 2, stats::weighted.mean, w = weights)
  cenX    <- sweep(sm$X, 2, w_means, "-")

  # Detect zero-variance columns as the package does internally
  is0var  <- apply(cenX, 2, function(col) stats::var(col) < 1e-8)
  ind0var <- if (any(is0var)) which(is0var) else integer(0)

  expect_gt(length(ind0var), 0L)
  Xs <- partition_t2_basis_matrix(sm, cenX, ind0var)

  # Total retained columns across all partitions must equal ncol(X) - |ind0var|
  total_betas <- sum(sapply(Xs, ncol))
  expect_equal(total_betas, ncol(sm$X) - length(ind0var))
})
