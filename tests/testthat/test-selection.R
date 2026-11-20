# Learner selection: selection losses, defaults and excluded learners.


# ---- Selection loss ----------------------------------------------------------

test_that("default selection is the unshrunk rss of the gradient", {
  lc <- gb_learner_ctrl()
  expect_identical(lc$selection_loss_type, "rss")
  expect_false(lc$use_step_size_for_selection)
})

test_that("default selection fits large-count poisson data sensibly", {
  # An unshrunk provisional deviance update overshoots for large counts and
  # selected only the intercept; the default rss criterion must not.
  set.seed(5)
  n <- 300
  d <- data.frame(x = stats::runif(n), z = stats::runif(n))
  d$y <- stats::rpois(n, exp(2 + sin(2 * pi * d$x)))
  m <- quiet_fit(y ~ t2(x, k = 8) + t2(z, k = 8), data = d,
                 family = stats::poisson(),
                 stop_control = gb_stop_ctrl("aic", max_iter = 3000L))
  expect_false(m$metrics$stopping_conditions[["mstop_reached"]])
  ref <- mgcv::gam(y ~ s(x) + s(z), data = d, family = stats::poisson(),
                   method = "REML")
  expect_gt(stats::cor(m$linear.predictors, ref$linear.predictors), 0.99)
})

test_that("poisson selection loss has the same argmin as the deviance", {
  dat <- simulate_poisson_1d(n = 120, seed = 116)
  m <- quiet_fit(y ~ t2(x, k = 8, bs = "cr"), data = dat,
                 family = stats::poisson(),
                 stop_control = gb_stop_ctrl("none", max_iter = 5L),
                 edf_control = gb_edf_ctrl("none"))
  # The fast poisson path is only used by the deviance criterion
  m$selection_loss_type         <- "deviance"
  m$use_step_size_for_selection <- TRUE
  g  <- pseudo_residuals(m$y, m$linear.predictors, m$family)
  Wg <- m$prior.weights * g
  fast <- dev <- numeric(length(m$learners))
  for (k in seq_along(m$learners)) {
    lf <- fit_learner(m$learners[[k]], Wg)
    fast[k] <- learner_selection_loss(m$learners[[k]], lf, m, g)
    nu_sel <- if (isTRUE(m$use_step_size_for_selection)) m$step_size else 1
    nu  <- nu_sel * (m$learners[[k]]$nu_k %||% 1)
    mu  <- m$family$linkinv(m$linear.predictors + nu * lf$ghat)
    dev[k] <- sum(m$family$dev.resids(m$y, mu, m$prior.weights))
  }
  expect_equal(which.min(fast), which.min(dev))
  # Differences from the deviance are a learner-independent constant
  expect_equal(2 * fast - dev, rep((2 * fast - dev)[1], length(dev)),
               tolerance = 1e-10)
})

# Gamma data with means from 0.05 to 20, where step_size = 1 overshoots
gamma_wide_data <- function() {
  set.seed(1)
  n <- 100
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rgamma(n, shape = 5, rate = 5 / (0.05 + 20 * d$x^4))
  d
}

lc_dev <- gb_learner_ctrl(step_size = 1, selection_loss_type = "deviance",
                          use_step_size_for_selection = TRUE)

test_that("deviance selection with a too large step stops as diverged", {
  # The full fit stops at the first increase of the training deviance
  # (divergence check), before any selection loss is non-finite
  expect_error(
    suppressWarnings(quiet_fit(
      y ~ t2(x, k = 6, bs = "cr"), data = gamma_wide_data(),
      family = stats::Gamma(), learner_control = lc_dev,
      stop_control = gb_stop_ctrl(iter_no_improv = 20L, max_iter = 500L))),
    "Boosting diverged"
  )
})

test_that("non-finite selection losses for every learner give a clear error", {
  # A provisional step that leaves the valid range for every candidate (the
  # intercept, whose fit is near zero, is excluded)
  m1 <- suppressWarnings(quiet_fit(
    y ~ t2(x, k = 6, bs = "cr"), data = gamma_wide_data(),
    family = stats::Gamma(), learner_control = lc_dev,
    stop_control    = gb_stop_ctrl("none", max_iter = 1L),
    edf_control     = ec_ff))
  m1$step_size <- 1e3
  expect_error(suppressWarnings(iterate_mod(m1, exclude = 1L)),
               "non-finite selection loss")
})

test_that("rss selection with a too large step reports the deviance", {
  # The selection losses stay finite; the step check reports the invalid
  # linear predictor (the stopping-loss guard is unit-tested in
  # test-stopping-loss.R)
  expect_error(
    suppressWarnings(quiet_fit(
      y ~ t2(x, k = 6, bs = "cr"), data = gamma_wide_data(),
      family = stats::Gamma(),
      learner_control = gb_learner_ctrl(step_size = 1),
      stop_control = gb_stop_ctrl(iter_no_improv = 20L, max_iter = 500L))),
    "training deviance became non-finite at iteration 1"
  )
})

test_that("iterate_mod() reports an empty candidate set", {
  dat <- simulate_gaussian_1d(n = 40L, seed = 206L)
  mod <- quiet_fit(y ~ t2(x, k = 6), data = dat,
                   stop_control = gb_stop_ctrl("none", max_iter = 1L),
                   edf_control = ec_ff)
  expect_error(iterate_mod(mod, exclude = seq_along(mod$learners)),
               "No base learner is left")
})

# ---- Excluded learners -------------------------------------------------------

test_that("excluded_learner_indices prevents those learners from being selected", {
  skip_on_cran()
  dat <- simulate_gaussian_2d(n = 100, seed = 60)

  k_exclude <- 2L  # exclude the first smooth partition

  mod <- gbgam(
    y ~ t2(x, k = 5, bs = "cr") + t2(z, k = 5, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(excluded_learner_indices = k_exclude),
    stop_control    = gb_stop_ctrl("none", max_iter = 30),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  expect_true(!any(mod$metrics$k_opts[seq_len(mod$mstop)] == k_exclude))
})
