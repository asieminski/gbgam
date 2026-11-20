# Replaying preselected selection and rho paths
# (preselected_learner_indices, preselected_rhos).


# ---- Replaying rho paths -----------------------------------------------------

test_that("replaying rho path reproduces a fit with learner step sizes", {
  skip_on_cran()
  set.seed(107)
  n <- 150
  d <- data.frame(x = stats::runif(n), z = stats::runif(n))
  d$y <- sin(2 * pi * d$x) + 0.5 * d$z + stats::rnorm(n, sd = 0.3)
  fml <- y ~ t2(x, k = 10, bs = "cr") + t2(z, k = 6, bs = "cr")
  m <- quiet_fit(fml, data = d,
                 learner_control = gb_learner_ctrl(
                   step_size = 0.1, tau = 1,
                   use_learner_step_size_for_update = TRUE),
                 stop_control = gb_stop_ctrl("none", max_iter = 80L),
                 edf_control = gb_edf_ctrl("exact"))
  # ASR has relaxed at least one learner, so nu_k has changed
  expect_true(length(unique(stats::na.omit(m$metrics$rhos))) > 2L)
  r <- quiet_fit(fml, data = d,
                 learner_control = gb_learner_ctrl(
                   step_size = 0.1, tau = 0,
                   use_learner_step_size_for_update = TRUE,
                   preselected_learner_indices = m$metrics$k_opts,
                   preselected_rhos = ifelse(is.na(m$metrics$rhos), 0,
                                             m$metrics$rhos)),
                 stop_control = gb_stop_ctrl("none", max_iter = 80L),
                 edf_control = gb_edf_ctrl("none"))
  expect_equal(r$linear.predictors, m$linear.predictors, tolerance = 1e-10)
})

test_that("rho replay reproduces the fit when a learner saturates (lambda = 0)", {
  skip_on_cran()
  set.seed(114)
  n <- 120
  d <- data.frame(x = stats::runif(n), z = stats::runif(n))
  d$y <- sin(2 * pi * d$x) + 0.5 * d$z + stats::rnorm(n, sd = 0.3)
  # k = 4 gives a 2-column penalised block, which tau = 3 saturates quickly
  fml <- y ~ t2(x, k = 4, bs = "cr") + t2(z, k = 6, bs = "cr")
  for (use_nu_k in c(TRUE, FALSE)) {
    m <- quiet_fit(fml, data = d,
                   learner_control = gb_learner_ctrl(
                     step_size = 0.1, tau = 3,
                     use_learner_step_size_for_update = use_nu_k),
                   stop_control = gb_stop_ctrl("none", max_iter = 60L),
                   edf_control = gb_edf_ctrl("exact"))
    saturated <- vapply(m$learners, function(lr) isTRUE(lr$lambda_at_zero),
                        logical(1))
    expect_true(any(saturated))
    # The saturated learner keeps sedf = 1: nu_k = 1 - sqrt(1 - 1 / p)
    lr_sat <- m$learners[[which(saturated)[1]]]
    expect_equal(lr_sat$nu_k, 1 - sqrt(1 - 1 / length(lr_sat$beta)),
                 tolerance = 1e-12)
    r <- quiet_fit(fml, data = d,
                   learner_control = gb_learner_ctrl(
                     step_size = 0.1, tau = 0,
                     use_learner_step_size_for_update = use_nu_k,
                     preselected_learner_indices = m$metrics$k_opts,
                     preselected_rhos = ifelse(is.na(m$metrics$rhos), 0,
                                               m$metrics$rhos)),
                   stop_control = gb_stop_ctrl("none", max_iter = 60L),
                   edf_control = gb_edf_ctrl("none"))
    expect_equal(r$linear.predictors, m$linear.predictors, tolerance = 1e-8)
    expect_equal(r$coefficients, m$coefficients, tolerance = 1e-8)
    nu_ks <- function(mod) {
      vapply(mod$learners, function(lr) lr$nu_k %||% NA_real_, numeric(1))
    }
    expect_equal(nu_ks(r), nu_ks(m), tolerance = 1e-12)
  }
})

test_that("post-iteration learner rhos are recorded and can be replayed", {
  set.seed(205)
  n <- 150
  d <- data.frame(x = stats::runif(n), g = factor(sample(1:12, n, TRUE)))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(12)[d$g] +
    stats::rnorm(n, sd = 0.3)
  fml  <- y ~ t2(x, k = 8, bs = "cr") + t2(g, bs = "re")
  k_re <- 4L
  n_it <- 40L
  lc <- list(step_size = 0.3, learners_to_fit_after_each_iter = k_re,
             excluded_learner_indices = k_re)
  mod <- quiet_fit(fml, data = d,
                   learner_control = do.call(gb_learner_ctrl, lc),
                   stop_control = gb_stop_ctrl("none", max_iter = n_it),
                   edf_control = gb_edf_ctrl("exact"))
  expect_identical(mod$learners[[k_re]]$var_label, "t2(g)r")
  pr <- mod$metrics$post_rhos
  expect_true(is.matrix(pr))
  expect_equal(dim(pr), c(length(mod$metrics$rhos), 1L))
  # ASR moved the random-effect rho away from rho0 and the last row is the
  # learner's final rho
  expect_lt(pr[n_it, 1L], mod$learners[[k_re]]$rho0 - 1)
  expect_equal(pr[n_it, 1L], mod$learners[[k_re]]$rho)

  # Fixed-path replay without EDF tracking (as in the 'fixed' bootstrap)
  expect_equal(mod$mstop, n_it)
  replay <- function(post_rhos) replay_on_y(mod, mod$y, "fixed", post_rhos)
  f <- replay(pr)
  expect_equal(f$linear.predictors, mod$linear.predictors, tolerance = 1e-10)
  expect_equal(f$metrics$post_rhos, pr)
  # Without the post-iteration rhos the replay keeps rho0 and differs
  f0 <- replay(NULL)
  expect_gt(max(abs(f0$linear.predictors - mod$linear.predictors)), 1e-3)
  # A wrongly shaped replay matrix is rejected
  expect_error(replay(pr[1:5, , drop = FALSE]), "preselected_post_rhos")

  # Unknown-length runs record one row per iteration; no post learners give
  # zero columns
  m2 <- quiet_fit(fml, data = d,
                  learner_control = do.call(gb_learner_ctrl, lc),
                  stop_control = gb_stop_ctrl("aic", iter_no_improv = 5L),
                  edf_control = gb_edf_ctrl("exact"))
  expect_equal(dim(m2$metrics$post_rhos), c(length(m2$metrics$rhos), 1L))
  m3 <- quiet_fit(fml, data = d, stop_control = sc_ff, edf_control = ec_ff)
  expect_equal(dim(m3$metrics$post_rhos), c(15L, 0L))
})

test_that("preselected path shorter than max_iter gives a clear error", {
  dat <- simulate_gaussian_1d(n = 60, seed = 108)
  expect_error(
    quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = dat,
              learner_control = gb_learner_ctrl(
                preselected_learner_indices = c(2L, 3L, 2L)),
              stop_control = gb_stop_ctrl("none", max_iter = 5L),
              edf_control = gb_edf_ctrl("none")),
    "preselected_learner_indices has length 3"
  )
  expect_error(
    quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = dat,
              learner_control = gb_learner_ctrl(
                preselected_learner_indices = c(2L, 3L, 2L))),
    "max_iter is not set"
  )
})
