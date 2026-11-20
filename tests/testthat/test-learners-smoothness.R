# Adaptive smoothness relaxation (ASR) of the learners: rho updates,
# lambda set to zero, init_rho() and update_learner_trA_smoothness().


# ---- Adaptive smoothness relaxation ------------------------------------------

test_that("ASR reduces smoothing penalty when the same smooth is repeatedly selected", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 51, sd = 0.15)

  learner_obj <- new_learners(
    y ~ t2(x, k = 9, bs = "cr"),
    data = dat,
    weights = rep(1, nrow(dat))
  )
  k_smooth <- which(vapply(learner_obj$learners,
                            function(x) inherits(x, "demmler_reinsch"),
                            logical(1)))[1]
  # 60 steps of 0.05: the learner's EDF passes its tr(A), so ASR relaxes
  n_it <- 60L
  path <- rep(k_smooth, n_it)
  rho0 <- learner_obj$learners[[k_smooth]]$rho

  fit_tau <- function(tau) {
    gbgam(
      y ~ t2(x, k = 9, bs = "cr"),
      data              = dat,
      family            = stats::gaussian(),
      learner_control   = gb_learner_ctrl(tau = tau, step_size = 0.05,
                                          preselected_learner_indices = path,
                                          learners_to_fit_after_each_iter = 1L,
                                          excluded_learner_indices = 1L),
      stop_control = gb_stop_ctrl("none", max_iter = n_it),
      edf_control  = gb_edf_ctrl("exact"),
      uq_control   = gb_uq_ctrl("none"),
      verbose      = FALSE
    )
  }
  final_rho_tau0 <- fit_tau(0)$learners[[k_smooth]]$rho
  final_rho_tau1 <- fit_tau(1)$learners[[k_smooth]]$rho

  # tau = 0 must not change rho from initialisation
  expect_equal(final_rho_tau0, rho0, tolerance = 1e-10)
  # tau = 1 relaxes the penalty: rho = log(lambda) falls, the fit gets rougher
  # (it falls by 1.8 here)
  expect_lt(final_rho_tau1, rho0 - 1)
})

test_that("lambda-to-zero branch keeps the selection edf at its target", {
  set.seed(5)
  n   <- 100
  lb  <- new_learners(y ~ t2(x, k = 5, bs = "cr"),
                      data = data.frame(x = stats::runif(n), y = stats::rnorm(n)),
                      weights = rep(1, n))
  lr  <- dr_learners(lb)[[1]]
  p   <- length(lr$beta)
  lr$edf <- p + 0.5
  lr2 <- suppressMessages(update_learner_trA_smoothness(
    lr, target_tr2AminAA = 1, tau = 1, verbose = FALSE
  ))
  expect_true(lr2$lambda_at_zero)
  expect_equal(sedf(lr2), 1, tolerance = 1e-8)
  expect_equal(lr2$nu_k, 1 - sqrt(1 - 1 / p), tolerance = 1e-8)
})

test_that("lambda_at_zero flag suppresses repeated messages once rho=-1e16", {
  # A learner from new_learners() with its EDF set above p by hand: a fit
  # cannot reliably push a learner's EDF past p
  dat <- simulate_gaussian_1d(n = 40, seed = 55)
  lb  <- new_learners(y ~ t2(x, k = 5, bs = "cr"), data = dat,
                      weights = rep(1, 40))
  lr <- dr_learners(lb)[[1L]]
  lr$edf <- length(lr$beta) + 10

  # First call: exactly one message, and the flag is set
  msgs <- character(0)
  withCallingHandlers(
    { lr <- update_learner_trA_smoothness(lr,
                                          notify_lambda_zero = TRUE) },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(isTRUE(lr$lambda_at_zero))
  expect_equal(lr$rho, -1e16)
  expect_length(msgs, 1L)

  # Second call: the flag is set, so no further message
  msgs2 <- character(0)
  withCallingHandlers(
    { lr <- update_learner_trA_smoothness(lr,
                                          notify_lambda_zero = TRUE) },
    message = function(m) {
      msgs2 <<- c(msgs2, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_length(msgs2, 0L)
})

test_that("lambda-to-zero message is emitted only when mod$verbose is TRUE", {
  set.seed(5)
  n  <- 100
  lb <- new_learners(y ~ t2(x, k = 5, bs = "cr"),
                     data = data.frame(x = stats::runif(n), y = stats::rnorm(n)),
                     weights = rep(1, n))
  lr <- dr_learners(lb)[[1]]
  lr$edf <- length(lr$beta) + 0.5
  upd <- function(verbose) {
    update_learner(lr, mod = list(tau = 1, target_tr2AminAA = 1,
                                  verbose = verbose),
                   dbeta = 0, deta = 0, dedf = 0)
  }
  expect_silent(lr_q <- upd(FALSE))
  expect_true(lr_q$lambda_at_zero)
  expect_message(lr_v <- upd(TRUE), "Setting lambda to zero")
  expect_identical(lr_q, lr_v)
})


# ---- init_rho() and update_learner_trA_smoothness() --------------------------

test_that("init_rho extends its search interval for very large weights", {
  set.seed(1)
  n   <- 200
  dat <- data.frame(x = stats::runif(n))
  dat$y <- sin(2 * pi * dat$x) + stats::rnorm(n, 0, 0.3)
  f   <- y ~ t2(x, k = 8, bs = "cr")
  l1  <- new_learners(f, data = dat, weights = rep(1, n))$learners[[2]]
  lw  <- new_learners(f, data = dat, weights = rep(1e12, n))$learners[[2]]
  sedf_unscaled <- function(lr) learner_edf(1 / (1 + exp(lr$rho) * lr$D))
  expect_equal(sedf_unscaled(lw), 1, tolerance = 1e-4)
  expect_equal(lw$rho - l1$rho, log(1e12), tolerance = 1e-4)
  expect_error(init_rho(c(1, 2), target_edf = 5, label = "t2(x)r"),
               "t2\\(x\\)r")
})

test_that("update_learner_trA_smoothness reports updates with message()", {
  set.seed(1)
  dat <- data.frame(x = stats::runif(100), y = stats::rnorm(100))
  lr  <- new_learners(y ~ t2(x, k = 8, bs = "cr"), data = dat)$learners[[2]]
  lr$edf <- 3
  expect_silent(update_learner_trA_smoothness(lr, tol = 1e-6))
  expect_message(lr2 <- update_learner_trA_smoothness(lr, tol = 1e-6,
                                                      verbose = TRUE),
                 "updated to tr\\(A\\) = 3.00")
  expect_equal(sum(1 / (1 + lr2$D * exp(lr2$rho))), 3, tolerance = 1e-6)
})

test_that("update_learner_trA_smoothness validates its arguments", {
  set.seed(1)
  dat <- data.frame(x = stats::runif(100), y = stats::rnorm(100),
                    z = stats::rnorm(100))
  lb  <- new_learners(y ~ z + t2(x, k = 8, bs = "cr"), data = dat)
  lz  <- lb$learners[[2]]  # unpenalised linear term
  lr  <- lb$learners[[3]]  # penalised learner of t2(x)
  lz$edf <- lr$edf <- 3
  expect_error(update_learner_trA_smoothness(lz),
               "Must inherit from class 'demmler_reinsch'")
  expect_error(update_learner_trA_smoothness(lr, tol = 0),
               "tol must be positive")
  expect_error(update_learner_trA_smoothness(lr, tol = Inf), "tol")
  expect_error(update_learner_trA_smoothness(lr, tau = NA), "tau")
  expect_error(update_learner_trA_smoothness(lr, tau = -1), "tau")
  expect_error(update_learner_trA_smoothness(lr, target_tr2AminAA = 0),
               "target_tr2AminAA must be positive")
  expect_error(update_learner_trA_smoothness(lr, target_tr2AminAA = NA),
               "target_tr2AminAA")
  expect_error(update_learner_trA_smoothness(lr, verbose = NA), "verbose")
  expect_error(update_learner_trA_smoothness(lr, notify_lambda_zero = "yes"),
               "notify_lambda_zero")
})
