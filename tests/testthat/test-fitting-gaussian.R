# Gaussian fits: the model object, the RSS path, equivalent selection
# losses, preselected paths and rhos, intercept options and EDF tracking.


# ---- Gaussian fits -----------------------------------------------------------

test_that("basic Gaussian fit returns a coherent model object", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 30)

  mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, max_iter = 20)

  expect_model_contract(mod, nrow(dat))
  expect_equal(mod$mstop, 20)
  expect_length(mod$metrics$k_opts, 20)
  expect_length(mod$metrics$rhos, 20)
  # stopping_method = "none": no losses are tracked
  expect_true(all(is.na(mod$metrics$ls)))
})

test_that("Gaussian RSS is monotone non-increasing over iterations", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 31, sd = 0.35)

  mod <- fit_fixed(y ~ t2(x, k = 7, bs = "cr"), dat,
                   save_fitted_values = TRUE)

  rss <- deviance_path(dat$y, mod$fitted_values)
  expect_true(all(diff(rss) <= 1e-10))
})

test_that("Gaussian deviance and RSS selection give identical paths and coefficients", {
  # In the Gaussian case, deviance == RSS (up to the constant 2/phi), so
  # selection_loss_type = "deviance" and "rss" must pick the same learner at
  # every iteration and produce identical coefficients.  The same equivalence
  # holds for shrinking predictions before selection vs not.
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 36)

  sc <- gb_stop_ctrl("none", max_iter = 20)
  ec <- gb_edf_ctrl("exact")
  ic <- gb_uq_ctrl("none")

  mod_dev <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"), data = dat,
    learner_control   = gb_learner_ctrl(tau = 0, step_size = 0.05,
                                        selection_loss_type = "deviance",
                                        use_step_size_for_selection = TRUE),
    stop_control = sc, edf_control = ec, uq_control = ic,
    verbose = FALSE
  )

  mod_rss <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"), data = dat,
    learner_control   = gb_learner_ctrl(tau = 0, step_size = 0.05,
                                        selection_loss_type = "rss",
                                        use_step_size_for_selection = TRUE),
    stop_control = sc, edf_control = ec, uq_control = ic,
    verbose = FALSE
  )

  mod_unshrunk <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"), data = dat,
    learner_control   = gb_learner_ctrl(tau = 0, step_size = 0.05,
                                        selection_loss_type = "deviance",
                                        use_step_size_for_selection = FALSE),
    stop_control = sc, edf_control = ec, uq_control = ic,
    verbose = FALSE
  )

  # All three must select the same learner at every iteration
  expect_equal(mod_dev$metrics$k_opts, mod_rss$metrics$k_opts)
  expect_equal(mod_dev$metrics$k_opts, mod_unshrunk$metrics$k_opts)

  # Coefficients must be numerically identical
  expect_equal(mod_dev$coefficients, mod_rss$coefficients,      tolerance = 1e-10)
  expect_equal(mod_dev$coefficients, mod_unshrunk$coefficients, tolerance = 1e-10)

  # EDF paths must be identical (same path => same EDF increments)
  expect_equal(mod_dev$edfs, mod_rss$edfs,      tolerance = 1e-10)
  expect_equal(mod_dev$edfs, mod_unshrunk$edfs, tolerance = 1e-10)
})

test_that("preselected learner indices are respected exactly", {
  skip_on_cran()
  dat <- simulate_gaussian_2d(n = 100, seed = 32)

  learner_obj <- new_learners(
    y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 6, bs = "cr"),
    data = dat,
    weights = rep(1, nrow(dat))
  )

  # Cycle through all learner indices 3 times for a robust, non-degenerate path
  n_learners <- length(learner_obj$learners)
  path <- rep(seq_len(n_learners), 3L)

  mod <- gbgam(
    y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 6, bs = "cr"),
    data              = dat,
    family            = stats::gaussian(),
    learner_control   = gb_learner_ctrl(step_size = 0.05,
                                        preselected_learner_indices = path),
    stop_control = gb_stop_ctrl("none", max_iter = length(path)),
    edf_control  = gb_edf_ctrl("none"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )

  expect_equal(mod$metrics$k_opts, path)
})

test_that("preselected rhos are applied exactly when supplied", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 90, seed = 33)

  learner_obj <- new_learners(
    y ~ t2(x, k = 7, bs = "cr"),
    data = dat,
    weights = rep(1, nrow(dat))
  )

  k_smooth <- which(vapply(learner_obj$learners,
                            function(x) inherits(x, "demmler_reinsch"),
                            logical(1)))[1]
  path <- rep(k_smooth, 10)
  rhos <- seq(1.5, -1, length.out = 10)

  mod <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"),
    data              = dat,
    family            = stats::gaussian(),
    learner_control   = gb_learner_ctrl(tau = 0, step_size = 0.05,
                                        preselected_learner_indices = path,
                                        preselected_rhos = rhos),
    stop_control = gb_stop_ctrl("none", max_iter = 10),
    edf_control  = gb_edf_ctrl("exact"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )

  expect_equal(mod$metrics$k_opts, path)
  expect_equal(mod$metrics$rhos,   rhos, tolerance = 1e-12)
})

test_that("intercept handling options all produce finite valid fits", {
  skip_on_cran()
  set.seed(35)
  dat <- data.frame(
    x = seq(0, 1, length.out = 100),
    y = 3 + stats::rnorm(100, sd = 0.2)
  )

  # Three intercept-handling modes (not equivalent in general; this test only
  # checks that each one gives a finite, valid fit):
  # (1) intercept competes for selection (default)
  # (2) intercept excluded from selection, refitted after each iter to the
  #     re-evaluated gradient (cf. bamlss::boost(always = ))
  # (3) intercept excluded from selection, never force-fitted
  ctrl_list <- list(
    list(),
    list(learners_to_fit_after_each_iter = 1L, excluded_learner_indices = 1L),
    list(excluded_learner_indices = 1L)
  )

  for (ctrl in ctrl_list) {
    mod <- gbgam(
      y ~ t2(x, k = 8, bs = "cr"),
      data              = dat,
      family            = stats::gaussian(),
      learner_control   = do.call(gb_learner_ctrl,
                                  c(list(step_size = 0.05), ctrl)),
      stop_control = gb_stop_ctrl("none", max_iter = 20),
      edf_control  = gb_edf_ctrl("none"),
      uq_control   = gb_uq_ctrl("none"),
      verbose      = FALSE
    )

    expect_model_contract(mod, nrow(dat))
    expect_lt(abs(mean(mod$fitted.values) - mean(dat$y)), 0.5)
  }
})

test_that("exact EDF and no-EDF paths agree when tau is zero and stopping is fixed", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 34)

  # Without EDF tracking ASR is off, so only the exact fit needs tau = 0
  mod_none <- fit_fixed(y ~ t2(x, k = 7, bs = "cr"), dat, max_iter = 18)

  mod_exact <- gbgam(
    y ~ t2(x, k = 7, bs = "cr"),
    data            = dat,
    family          = stats::gaussian(),
    learner_control = gb_learner_ctrl(tau = 0, step_size = 0.05),
    stop_control    = gb_stop_ctrl("none", max_iter = 18),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  expect_equal(mod_none$metrics$k_opts,    mod_exact$metrics$k_opts)
  expect_equal(mod_none$linear.predictors, mod_exact$linear.predictors, tolerance = 1e-10)
  expect_equal(mod_none$coefficients,      mod_exact$coefficients,      tolerance = 1e-10)
})
