# EDF tracking: exact and Hutch++ paths, and the EDF fields of the model.


# ---- EDF paths ---------------------------------------------------------------

test_that("exact EDF path is non-decreasing over iterations (no intercept refit)", {
  skip_on_cran()
  # The non-decreasing EDF property is guaranteed only for Gaussian with exact
  # EDF computation and no intercept refit.  It does NOT hold in general for
  # other families or stochastic EDF estimators (e.g. Hutch++).
  # Force-fitting the intercept after each iter can paradoxically decrease EDF,
  # so we exclude it from selection without force-fitting it.
  dat <- simulate_gaussian_1d(n = 80, seed = 50)

  mod <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"),
    data              = dat,
    family            = stats::gaussian(),
    learner_control   = gb_learner_ctrl(tau = 0, step_size = 0.05,
                                        excluded_learner_indices = 1L),
    stop_control = gb_stop_ctrl("none", max_iter = 20),
    edf_control  = gb_edf_ctrl("exact"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )

  expect_true("edfs" %in% names(mod))
  expect_length(mod$edfs, mod$mstop + 1L)
  expect_true(all(diff(mod$edfs) >= -1e-8))
})

test_that("Hutch++ and exact EDF give similar trajectories on the same problem", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 70, seed = 53)

  mod_exact <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"),
    data            = dat,
    family          = stats::gaussian(),
    learner_control = gb_learner_ctrl(tau = 1, step_size = 0.05),
    stop_control    = gb_stop_ctrl("none", max_iter = 18),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  # Seed the stochastic Hutch++ estimator for reproducibility
  set.seed(200)
  mod_hpp <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"),
    data            = dat,
    family          = stats::gaussian(),
    learner_control = gb_learner_ctrl(tau = 1, step_size = 0.05),
    stop_control    = gb_stop_ctrl("none", max_iter = 18),
    edf_control     = gb_edf_ctrl("hutchpp"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  expect_true("edfs" %in% names(mod_exact))
  expect_true("edfs" %in% names(mod_hpp))
  expect_equal(length(mod_exact$edfs), length(mod_hpp$edfs))

  expect_gt(stats::cor(mod_exact$edfs, mod_hpp$edfs), 0.90)
  expect_gt(stats::cor(mod_exact$linear.predictors, mod_hpp$linear.predictors), 0.90)
})


# ---- EDF fields of the model -------------------------------------------------

test_that("edf_compute_method = 'none' gives NA EDF fields, not 0", {
  d <- core_data()
  m <- gbgam(ff_core, data = d, verbose = FALSE,
             stop_control = gb_stop_ctrl("none", max_iter = 20),
             edf_control  = gb_edf_ctrl("none"))
  expect_true(is.na(m$edf.total))
  expect_true(all(is.na(m$edf)) && all(is.na(m$edf1)))
  expect_true(is.na(m$smooth[[1]]$edf))
  expect_true(is.na(m$rank))
  expect_true(is.na(m$df.residual))
  # The scale still uses 0 model degrees of freedom
  m0 <- m
  m0$edf.total <- 0
  expect_equal(m$scale, estimate_dispersion(m0))
  ll <- logLik(m)
  expect_true(is.na(attr(ll, "df")))
  expect_equal(as.numeric(ll), -neg2loglik(m0) / 2)
  expect_null(summary(m)$edf)
  expect_output(print(m), "degrees of freedom = not tracked")
})

test_that("per-smooth EDF is attached to smooth[[j]]$edf when edf is computed", {
  skip_on_cran()
  mod <- local_gaussian_mod_edf()
  expect_true(length(mod$smooth) >= 1L)
  for (j in seq_along(mod$smooth)) {
    edf_j <- mod$smooth[[j]]$edf
    expect_true(is.numeric(edf_j) && length(edf_j) == 1L,
                info = paste("smooth", j, "missing scalar edf"))
    expect_true(is.finite(edf_j),
                info = paste("smooth", j, "has non-finite edf:", edf_j))
    expect_gte(edf_j, 0,
               label = paste("smooth", j, "edf must be non-negative"))
  }
})
