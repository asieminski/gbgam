# Resumable fits: extend_gbgam() and rewinding a recorded fit.


# ---- Extending a fit ---------------------------------------------------------

test_that("extending a fit equals fitting the longer path at once", {
  dat <- simulate_poisson_1d(n = 120, seed = 3)
  for (edf in list(gb_edf_ctrl("exact"), gb_edf_ctrl("hutchpp", z = 12))) {
    set.seed(1)
    short <- fit_1d(dat, stats::poisson(), gb_stop_ctrl("none", max_iter = 40),
                    edf = edf)
    set.seed(1)
    long  <- fit_1d(dat, stats::poisson(), gb_stop_ctrl("none", max_iter = 100),
                    edf = edf)
    ext <- extend_gbgam(short, n_iter = 60)

    expect_equal(ext$mstop, long$mstop)
    expect_equal(ext$coefficients, long$coefficients, tolerance = 1e-10)
    expect_equal(ext$linear.predictors, long$linear.predictors,
                 tolerance = 1e-10)
    expect_equal(ext$edfs, long$edfs, tolerance = 1e-10)
    expect_equal(ext$metrics$k_opts, long$metrics$k_opts)
    expect_equal(ext$metrics$rhos, long$metrics$rhos, tolerance = 1e-10)
    expect_equal(ext$deviance, long$deviance, tolerance = 1e-10)
  }
})

test_that("extending an AIC fit chooses the best iteration over the whole path", {
  dat <- simulate_gaussian_1d(n = 120, seed = 1)
  short <- fit_1d(dat, stop_control = gb_stop_ctrl(max_iter = 5,
                                                   iter_no_improv_rel = 0.3))
  long  <- fit_1d(dat, stop_control = gb_stop_ctrl(max_iter = 300,
                                                   iter_no_improv = 1000))
  ext <- extend_gbgam(short, n_iter = 295)
  expect_equal(ext$mstop, long$mstop)
  expect_equal(ext$coefficients, long$coefficients, tolerance = 1e-10)
  expect_equal(ext$metrics$ls, long$metrics$ls, tolerance = 1e-10)
  # The user's patience is kept
  expect_equal(ext$iter_no_improv, short$iter_no_improv)
  expect_identical(ext$iter_no_improv_rel, 0.3)
  expect_true(ext$metrics$stopping_conditions[["mstop_reached"]])
})

test_that("the patience is not applied to an extension", {
  dat <- simulate_gaussian_1d(n = 120, seed = 1)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(iter_no_improv = 5,
                                                 iter_no_improv_rel = 0.2))
  run <- sum(!is.na(mod$metrics$k_opts))
  expect_true(mod$converged)
  ext <- extend_gbgam(mod, n_iter = 40)
  expect_identical(sum(!is.na(ext$metrics$k_opts)), run + 40L)
  expect_true(ext$metrics$stopping_conditions[["mstop_reached"]])
  expect_identical(ext[c("iter_no_improv", "iter_no_improv_rel")],
                   mod[c("iter_no_improv", "iter_no_improv_rel")])
})

test_that("extend_gbgam() keeps saved fitted values and can be chained", {
  dat <- simulate_gaussian_1d(n = 80, seed = 2)
  a <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = dat,
             stop_control = gb_stop_ctrl("none", max_iter = 10),
             edf_control = gb_edf_ctrl("exact"),
             save_fitted_values = TRUE, verbose = FALSE)
  b <- extend_gbgam(extend_gbgam(a, 5), 5)
  expect_equal(ncol(b$fitted_values), 21L)
  expect_equal(b$fitted_values[, 1:11], a$fitted_values)
  expect_equal(b$fitted_values[, 21], as.vector(b$fitted.values),
               tolerance = 1e-12)
})

test_that("extend_gbgam() needs a resumable fit", {
  dat <- simulate_gaussian_1d(n = 60, seed = 3)
  m <- fit_1d(dat, stop_control = gb_stop_ctrl("none", max_iter = 5,
                                                resumable = FALSE))
  expect_null(m$resume)
  expect_error(extend_gbgam(m, 5), "cannot be extended")
  m2 <- fit_1d(dat, stop_control = gb_stop_ctrl("none", max_iter = 5))
  expect_error(extend_gbgam(m2, 0))
})

test_that("extend_gbgam() stops on a model without iter_no_improv_rel", {
  dat <- simulate_gaussian_1d(n = 60, seed = 3)
  old <- fit_1d(dat, stop_control = gb_stop_ctrl(max_iter = 5))
  old$iter_no_improv_rel <- NULL
  expect_error(extend_gbgam(old, 5), "older gbgam.*refit it")
})

test_that("extending a CV fit equals running the longer CV fit", {
  dat <- simulate_poisson_1d(n = 100, seed = 3)
  for (keep in c(FALSE, TRUE)) {
    rs <- gb_resampling("bootstrap", B = 3, chunk = 20L, keep_fits = keep)
    # B = 3 bootstrap samples leave some observations in every sample
    set.seed(5)
    expect_warning(
      short <- fit_1d(dat, stats::poisson(),
                      gb_stop_ctrl("cv", max_iter = 30, resampling = rs),
                      edf = gb_edf_ctrl("hutchpp", z = 12)),
      "not left out")
    set.seed(5)
    expect_warning(
      long <- fit_1d(dat, stats::poisson(),
                     gb_stop_ctrl("cv", max_iter = 70, iter_no_improv = 1000,
                                  resampling = rs),
                     edf = gb_edf_ctrl("hutchpp", z = 12)),
      "not left out")
    ext <- extend_gbgam(short, 40)
    expect_equal(is.null(short$resume$cv$fits), !keep)
    expect_equal(ext$cv$rep_risk, long$cv$rep_risk, tolerance = 1e-10)
    expect_equal(ext$cv$risk, long$cv$risk, tolerance = 1e-10)
    expect_equal(ext$mstop, long$mstop)
    expect_equal(ext$coefficients, long$coefficients, tolerance = 1e-10)
  }
})

test_that("extending keeps one loss entry per iteration", {
  dat <- simulate_gaussian_1d(n = 60, seed = 3)
  m <- fit_1d(dat, stop_control = gb_stop_ctrl("none", max_iter = 10))
  e <- extend_gbgam(m, 10)
  # No loss is tracked, so there is no iteration-0 entry (as in a single fit)
  expect_equal(length(e$metrics$ls), length(e$metrics$k_opts))
  expect_equal(length(e$metrics$ls), 20L)
})

test_that("extending drops uncertainty quantification it cannot recompute", {
  dat <- simulate_gaussian_1d(n = 80, seed = 5)
  set.seed(1)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(max_iter = 20),
                uq_control = gb_uq_ctrl("bootstrap",
                                        n_samples = 3))
  expect_false(is.null(mod$Vp))
  ext <- suppressWarnings(extend_gbgam(mod, 10, max_runtime_mins = 1e-9))
  expect_false(identical(ext$coefficients, mod$coefficients))
  expect_null(ext$Vp)
  expect_null(ext$boot_coef)
  expect_null(ext$widening_factor)
  expect_null(ext$n_samples_done)
  expect_identical(ext$uq_control, gb_uq_ctrl("none"))
})

test_that("print() reports the time limit that ended an extension", {
  dat <- simulate_gaussian_1d(n = 80, seed = 5)
  for (orig_limit in list(NULL, 5)) {
    mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
      max_iter = 10, max_runtime_mins = orig_limit))
    ext <- extend_gbgam(mod, 500, max_runtime_mins = 0)
    expect_true(ext$metrics$stopping_conditions[["time_limit_exceeded"]])
    expect_identical(ext$max_runtime_mins, orig_limit)
    expect_output(print(ext), "time limit of 0 mins reached", fixed = TRUE)
    expect_match(summary(ext)$stop_msgs, "time limit of 0 mins reached",
                 fixed = TRUE, all = FALSE)
    # A later extension that ends at max_iter reports no time limit
    ext2 <- extend_gbgam(ext, 2)
    expect_false(any(grepl("time limit", summary(ext2)$stop_msgs)))
  }

  # The same for a fit stopped by a resampling estimate
  set.seed(9)
  cvm <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", max_iter = 10, resampling = gb_resampling(B = 2, chunk = 10L)))
  ext <- extend_gbgam(cvm, 500, max_runtime_mins = 0)
  expect_output(print(ext), "time limit of 0 mins reached", fixed = TRUE)
})

test_that("a failing uncertainty quantification keeps the extended model", {
  dat <- simulate_gaussian_1d(n = 80, seed = 5)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl("none", max_iter = 10))
  mod$uq_control <- gb_uq_ctrl("bootstrap", n_samples = 3)
  local_mocked_bindings(bootstrap_gbgam = function(...) stop("replicate 2"))
  e <- tryCatch(extend_gbgam(mod, 5),
                gbgam_uq_error = function(e) e)
  expect_s3_class(e, "gbgam_uq_error")
  expect_match(conditionMessage(e), "failed after the main fit: replicate 2")
  expect_s3_class(e$fit, "gbgam")
  expect_equal(e$fit$mstop, 15)
  expect_null(e$fit$Vp)
  expect_identical(e$fit$uq_control, gb_uq_ctrl("none"))
})


# ---- Rewinding ---------------------------------------------------------------

test_that("rewinding a recorded fit reproduces the shorter fit", {
  dat <- simulate_binomial_1d(n = 120, seed = 5)
  short <- fit_1d(dat, stats::binomial(), gb_stop_ctrl("none", max_iter = 12))
  # Recorded run via the cv machinery: rewind to iteration 12
  set.seed(9)
  cvm <- fit_1d(dat, stats::binomial(),
                gb_stop_ctrl("cv", max_iter = 30, iter_no_improv = 1000,
                             resampling = gb_resampling(B = 2, chunk = 30L)))
  st <- fit_state_resume(cvm, Sys.time(), FALSE)
  rew <- fit_state_rewind(st, 12)
  expect_equal(as.vector(rew$linear.predictors),
               as.vector(short$linear.predictors), tolerance = 1e-10)
  expect_equal(
    learners_to_mgcv_coef(rew$learners, rew$coef_names, rew$ind0var,
                          rew$cen_const, rew$cen_absorb),
    short$coefficients, tolerance = 1e-10)
  expect_equal(rew$edf.total, short$edfs[13], tolerance = 1e-10)
})
