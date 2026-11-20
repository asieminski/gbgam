# A step size too large for the data stops the fit with an informative
# error (as mboost's "Infinite residuals"), not a later failure or bad model.


# ---- Divergence of the fit ---------------------------------------------------

# Large Poisson counts (mean ~115): the learners fit the score y - mu at the
# prior weights, so one boosting step is about step_size * mu Fisher-scoring
# steps, which overshoots once it exceeds 2: step_size ~ 2 / 115 = 0.017.
# The limit is lower in practice: 0.014 converges, 0.015 diverges at
# iteration 2
large_count_data <- function(mean_y = 115, n = 300) {
  set.seed(1)
  x <- stats::runif(n)
  z <- stats::runif(n)
  data.frame(x = x, z = z,
             y = stats::rpois(n, mean_y * exp(0.5 * sin(2 * pi * x))))
}
ff_div <- y ~ t2(x, k = 10, bs = "cr") + t2(z, k = 10, bs = "cr")

fit_div <- function(dat, lc = gb_learner_ctrl(), ec = gb_edf_ctrl(),
                    sc = gb_stop_ctrl(iter_no_improv = 10, max_iter = 500)) {
  gbgam(ff_div, data = dat, family = stats::poisson(), learner_control = lc,
        edf_control = ec, stop_control = sc, verbose = FALSE)
}

test_that("the default step size is 0.1", {
  expect_identical(gb_learner_ctrl()$step_size, 0.1)
})

test_that("large-count Poisson diverges at the first update", {
  dat <- large_count_data()
  for (nu in c(0.1, 0.05, 0.02)) {
    expect_error(
      fit_div(dat, lc = gb_learner_ctrl(step_size = nu)),
      paste0("Boosting diverged at iteration 1: the training deviance ",
             "increased.*Decrease step_size in gb_learner_ctrl\\(\\)\\.$"))
  }
  # The divergence is caught before the EDF computation, for every EDF
  # method
  for (ec in list(gb_edf_ctrl("hutchpp"), gb_edf_ctrl("exact"),
                  gb_edf_ctrl("none"))) {
    expect_error(fit_div(dat, ec = ec,
                         sc = gb_stop_ctrl("none", max_iter = 50)),
                 "Boosting diverged at iteration 1")
  }
})

test_that("a step size below the stability limit gives a converging fit", {
  dat <- large_count_data()
  mod <- fit_div(dat, lc = gb_learner_ctrl(step_size = 0.01))
  expect_gt(mod$mstop, 5)
  expect_true(is.finite(AIC(mod)))
  expect_true(all(is.finite(mod$edfs[seq_len(mod$mstop + 1)])))
})

test_that("deviance selection stops instead of returning a garbage model", {
  # The divergence is an error, not a model chosen by AIC
  dat <- large_count_data()
  expect_error(
    fit_div(dat, lc = gb_learner_ctrl(selection_loss_type = "deviance",
                                      use_step_size_for_selection = TRUE)),
    "Boosting diverged at iteration 2")
})

test_that("an overflowing linear predictor gives 'Infinite residuals'", {
  dat <- large_count_data(mean_y = 1e8, n = 100)
  expect_error(fit_div(dat, sc = gb_stop_ctrl("none", max_iter = 5)),
               paste0("Infinite residuals at iteration 1: the step size is ",
                      "too large for this data; decrease step_size"))
})

test_that("check_boost_step accepts descent steps and ignores zero weights", {
  set.seed(3)
  n   <- 50
  mod <- list(y                 = stats::rpois(n, 3),
              prior.weights     = c(rep(1, n - 5), rep(0, 5)),
              family            = stats::poisson(),
              linear.predictors = rep(log(3), n),
              step_size         = 0.1)
  u <- pseudo_residuals(mod$y, mod$linear.predictors, mod$family)
  # A small step along the weighted negative gradient lowers the deviance
  expect_null(check_boost_step(mod, mod$linear.predictors + 0.01 * u, iter = 1))
  # Zero-weight observations neither enter the deviance nor the check
  eta <- mod$linear.predictors
  eta[n] <- 1e3
  expect_equal(train_deviance(mod, eta),
               train_deviance(mod, mod$linear.predictors))
  # An overshooting step is caught, with the iteration in the message
  expect_error(check_boost_step(mod, mod$linear.predictors + 5 * u, iter = 7),
               "Boosting diverged at iteration 7")
  eta[1] <- Inf
  expect_error(check_boost_step(mod, eta, iter = 2),
               "Infinite residuals at iteration 2")
})

test_that("an exact fit at the round-off floor of the deviance is no error", {
  # The deviance falls to ~1e-27, where its changes are round-off of any
  # relative size; these do not count as divergence
  set.seed(1)
  d <- data.frame(x = stats::runif(200))
  d$y <- 1 + 2 * d$x
  mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = d, verbose = FALSE,
               stop_control = gb_stop_ctrl("none", max_iter = 400))
  expect_equal(mod$mstop, 400)
  expect_lt(mod$deviance, 1e-20)
  # With AIC stopping it runs past the floor to its own stopping point
  mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = d, verbose = FALSE,
               stop_control = gb_stop_ctrl(max_iter = 1000))
  expect_gt(mod$mstop, 300)
})

test_that("errors in bootstrap replicates and resampling fits say so", {
  # Resampling weights can push a fit near the stability limit over it; the
  # error then names the replicate
  dat <- simulate_gaussian_1d(n = 60L, seed = 11L)
  mod <- gbgam(y ~ t2(x, k = 6), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl(max_iter = 30))
  local_mocked_bindings(
    refit_replicate = function(...) stop("Boosting diverged at iteration 3"))
  expect_error(bootstrap_gbgam(mod, n_samples = 2, verbose = FALSE),
               paste0("Bootstrap replicate 1 failed \\(the main fit is not ",
                      "affected\\): Boosting diverged at iteration 3"))
  local_mocked_bindings(
    advance_resample_fit = function(...) stop("Boosting diverged at x"))
  expect_error(
    gbgam(y ~ t2(x, k = 6), data = dat, verbose = FALSE,
          stop_control = gb_stop_ctrl(
            "cv", max_iter = 30,
            resampling = gb_resampling(B = 2, chunk = 30L))),
    paste0("Resampling fit 1 failed \\(the main fit is not affected\\): ",
           "Boosting diverged at x"))
})


# ---- Classed condition -------------------------------------------------------

test_that("divergence raises a gbgam_diverged condition", {
  set.seed(1)
  n <- 300
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rpois(n, 115 * exp(0.5 * sin(2 * pi * d$x)))
  e <- tryCatch(
    gbgam(y ~ t2(x, k = 10, bs = "cr"), data = d, family = poisson(),
          verbose = FALSE, stop_control = gb_stop_ctrl("none", max_iter = 5)),
    gbgam_diverged = function(e) e)
  expect_s3_class(e, "gbgam_diverged")
  expect_equal(e$iter, 1)
  expect_null(conditionCall(e))
})

test_that("the divergence message gives the size of the increase", {
  set.seed(3)
  n   <- 50
  mod <- list(y                 = stats::rpois(n, 3),
              prior.weights     = rep(1, n),
              family            = stats::poisson(),
              linear.predictors = rep(log(3), n),
              step_size         = 0.1)
  u       <- pseudo_residuals(mod$y, mod$linear.predictors, mod$family)
  eta_new <- mod$linear.predictors + 5 * u
  incr    <- train_deviance(mod, eta_new) -
    train_deviance(mod, mod$linear.predictors)
  msg <- tryCatch(check_boost_step(mod, eta_new, iter = 7),
                  error = conditionMessage)
  expect_match(msg, paste0("increased by ", format(signif(incr, 3))),
               fixed = TRUE)
  expect_match(msg, "of its value (from ", fixed = TRUE)
})
