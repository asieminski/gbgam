# Loss-based stopping: AIC (with the log-likelihood it uses) and the true_f
# oracle.


# ---- AIC stopping ------------------------------------------------------------

test_that("AIC stopping returns the best iteration on the tracked loss path", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 90, seed = 52, sd = 0.30)

  mod <- gbgam(
    y ~ t2(x, k = 8, bs = "cr"),
    data            = dat,
    family          = stats::gaussian(),
    learner_control = gb_learner_ctrl(tau = 0, step_size = 0.05),
    stop_control    = gb_stop_ctrl("aic", max_iter = 40, iter_no_improv = 8),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  expect_model_contract(mod, nrow(dat))

  # The returned model is within the iteration budget
  expect_gte(mod$mstop, 1L)
  expect_lte(mod$mstop, 40L)

  # The returned iteration achieves the minimum tracked AIC
  non_na_ls <- mod$metrics$ls[!is.na(mod$metrics$ls)]
  expect_true(length(non_na_ls) >= 2L)
  expect_lte(mod$metrics$ls[mod$mstop + 1L], min(non_na_ls) + 1e-10)
})

test_that("after AIC stopping with rollback, linear predictor is consistent with coefficients", {
  # AIC stopping runs to mstop + iter_no_improv and then rolls back to the best
  # iteration.  This test verifies that after rollback both mod$coefficients and
  # mod$linear.predictors refer to the same iteration, i.e. eta == X %*% beta.
  # We use a large max_iter with a generous iter_no_improv so that AIC fires
  # before the budget is exhausted and rollback definitely occurs.
  # Higher noise (sd = 3) ensures overfitting occurs within the budget.
  skip_on_cran()
  dat <- simulate_gaussian_2d(n = 110, seed = 54, sd = 3)

  mod <- gbgam(
    y ~ t2(x, k = 7, bs = "cr") + t2(x, z, k = c(4, 4), bs = c("cr", "cr")),
    data            = dat,
    family          = stats::gaussian(),
    learner_control = gb_learner_ctrl(tau = 1, step_size = 0.05),
    stop_control    = gb_stop_ctrl("aic", max_iter = 200, iter_no_improv = 8),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )

  # Rollback must have fired: AIC stopped before the iteration budget
  expect_lt(mod$mstop, 200L)

  # Coefficients and linear predictor must agree exactly
  eta_rebuilt <- rebuild_training_eta(mod)
  expect_equal(eta_rebuilt, as.vector(mod$linear.predictors), tolerance = 1e-8)
})

test_that("converged reports whether loss-based stopping ended the run", {
  d  <- core_data()
  m1 <- gbgam(ff_core, data = d, verbose = FALSE,
              stop_control = gb_stop_ctrl(iter_no_improv = 20))
  m2 <- gbgam(ff_core, data = d, verbose = FALSE,
              stop_control = gb_stop_ctrl(max_iter = 10))
  expect_true(m1$converged)
  expect_false(m2$converged)
})

test_that("non-integer poisson responses give a finite AIC and stop", {
  set.seed(105)
  n <- 100
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rpois(n, exp(1 + sin(2 * pi * d$x))) + 0.5
  expect_warning(
    m <- quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = d,
                   family = stats::poisson(),
                   learner_control = gb_learner_ctrl(step_size = 0.1),
                   stop_control = gb_stop_ctrl("aic", iter_no_improv = 10,
                                               max_iter = 2000)),
    "not an integer"
  )
  ls <- m$metrics$ls
  expect_true(all(is.finite(ls[!is.na(ls)])))
  expect_true(m$metrics$stopping_conditions[["no_improvements"]])
  # Continuous extension equals dpois() for integer counts
  m_int <- m
  m_int$y <- round(m$y)
  mu <- m_int$family$linkinv(m_int$linear.predictors)
  expect_equal(neg2loglik(m_int),
               -2 * sum(m_int$prior.weights *
                          stats::dpois(m_int$y, mu, log = TRUE)))
  expect_error(check_loss_finite(Inf, 3L), "non-finite")
  expect_null(check_loss_finite(NULL, 3L))
})

test_that("inverse Gaussian log-likelihood matches SuppDists::dinvGauss", {
  skip_if_not_installed("SuppDists")
  set.seed(106)
  n   <- 50
  mu  <- stats::runif(n, 1, 3)
  mod <- list(
    family            = new_family(stats::inverse.gaussian(link = "log")),
    y                 = SuppDists::rinvGauss(n, nu = mu, lambda = 3),
    linear.predictors = log(mu),
    prior.weights     = stats::runif(n, 0.5, 2),
    edf.total         = 3
  )
  # neg2loglik() (used by AIC stopping and logLik()) against SuppDists
  phi <- estimate_dispersion(mod)
  expect_equal(fitted_mu(mod), mu)
  expect_equal(
    neg2loglik(mod),
    -2 * sum(mod$prior.weights *
               SuppDists::dinvGauss(mod$y, nu = mu, lambda = 1 / phi,
                                    log = TRUE)),
    tolerance = 1e-10
  )
})


# ---- Oracle stopping ---------------------------------------------------------

test_that("true_f oracle stopping runs and stops before max_iter", {
  set.seed(1)
  n <- 200
  x <- runif(n)
  f <- sin(2 * pi * x)
  d <- data.frame(y = f + rnorm(n, sd = 0.5), x = x)
  m <- gbgam(y ~ t2(x), data = d, verbose = FALSE,
             stop_control = gb_stop_ctrl("true_f", true_f = f,
                                         max_iter = 2000, iter_no_improv = 20))
  expect_lt(m$mstop, 2000)
  expect_error(
    gbgam(y ~ t2(x), data = d, verbose = FALSE,
          stop_control = gb_stop_ctrl("true_f", true_f = f[-1], max_iter = 50)),
    "one value per observation"
  )
})


# ---- Patience ----------------------------------------------------------------

# First iteration at which the patience rule stops a run along the loss path
# ls (iteration 0 first); a tie moves the best iteration to the later one
first_stop <- function(ls, n_abs, rel) {
  best <- 0L
  for (m in seq_along(ls[-1L])) {
    if (ls[m + 1L] <= ls[best + 1L]) best <- m
    if (m - best >= max(n_abs, ceiling(rel * best))) return(m)
  }
  NA_integer_
}

test_that("the patience grows with the best iteration", {
  set.seed(1)
  n <- 200
  x <- runif(n)
  f <- sin(2 * pi * x)
  d <- data.frame(y = f + rnorm(n, sd = 0.5), x = x)
  fit <- function(...) {
    gbgam(y ~ t2(x), data = d, verbose = FALSE,
          stop_control = gb_stop_ctrl("true_f", true_f = f,
                                      max_iter = 2000, ...))
  }
  m_rel <- fit(iter_no_improv = 5, iter_no_improv_rel = 0.5)
  m_abs <- fit(iter_no_improv = 5, iter_no_improv_rel = 0)
  run_rel <- sum(!is.na(m_rel$metrics$k_opts))
  run_abs <- sum(!is.na(m_abs$metrics$k_opts))
  cat("true_f: mstop", m_rel$mstop, "run", run_rel, "(rel 0.5);",
      "mstop", m_abs$mstop, "run", run_abs, "(rel 0)\n")

  # The relative part is the larger one here, so it sets the stop
  expect_gt(ceiling(0.5 * m_rel$mstop), 5)
  expect_true(m_rel$converged)
  expect_identical(run_rel, as.integer(m_rel$mstop +
                                         ceiling(0.5 * m_rel$mstop)))
  expect_identical(run_rel, first_stop(m_rel$metrics$ls, 5, 0.5))
  # rel = 0 is the rule of iter_no_improv alone
  expect_true(m_abs$converged)
  expect_identical(run_abs, as.integer(m_abs$mstop + 5))
  expect_identical(run_abs, first_stop(m_abs$metrics$ls, 5, 0))
  # Same path up to the earlier stop
  expect_identical(m_abs$metrics$k_opts[seq_len(run_abs)],
                   m_rel$metrics$k_opts[seq_len(run_abs)])
  expect_output(print(m_rel), paste0("no MSE improvement after ",
                                     run_rel - m_rel$mstop, " iterations"))
})

test_that("AIC stopping with the default patience stops by its rule", {
  d <- core_data()
  m <- gbgam(ff_core, data = d, verbose = FALSE)
  expect_identical(m$iter_no_improv, 50)
  expect_identical(m$iter_no_improv_rel, 0.1)
  run <- sum(!is.na(m$metrics$k_opts))
  expect_true(m$converged)
  expect_identical(run, as.integer(m$mstop +
                                     max(50, ceiling(0.1 * m$mstop))))
  expect_identical(run, first_stop(m$metrics$ls, 50, 0.1))
})
