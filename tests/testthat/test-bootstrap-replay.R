# Bootstrap refits: reset_gbgam() and replays of the recorded selection and
# rho paths, which reproduce the fit on the original response.

# Two-smooth Gaussian fit with ASR and exact EDF
fit_asr_model <- function() {
  dat <- simulate_gaussian_2d(n = 150L, seed = 21L)
  gbgam(
    y ~ t2(x, k = 6) + t2(z, k = 6),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.1),
    stop_control    = gb_stop_ctrl("aic", iter_no_improv = 20L),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
}


# ---- reset_gbgam() -----------------------------------------------------------

test_that("reset_gbgam restores the construction-time nu_k (sedf = 1)", {
  mod <- fit_asr_model()
  dr  <- vapply(mod$learners, inherits, logical(1), what = "demmler_reinsch")
  # ASR must have relaxed at least one learner for the test to be meaningful
  expect_true(any(vapply(mod$learners[dr], function(lr) lr$nu_k < 1 &&
                           !isTRUE(all.equal(lr$rho, lr$rho0)), logical(1))))

  r <- reset_gbgam(mod, mod$y, mod$edf_compute_method)
  for (lr in r$learners[dr]) {
    expect_equal(lr$rho, lr$rho0)
    expect_false(lr$lambda_at_zero)
    expect_equal(lr$nu_k, min(1, find_nu_k(lr, target_df = 1)))
    if (lr$nu_k < 1) expect_equal(sedf(lr), 1, tolerance = 1e-8)
  }
})

test_that("reset_gbgam initialises the intercept per initialise_intercept", {
  dat <- simulate_poisson_1d(n = 120L, seed = 31L)
  for (init in c(TRUE, FALSE)) {
    for (edf_m in c("none", "exact")) {
      mod <- gbgam(
        y ~ t2(x, k = 6),
        data              = dat,
        family            = poisson(),
        learner_control   = gb_learner_ctrl(step_size = 0.1,
                                            initialise_intercept = init),
        stop_control = gb_stop_ctrl("none", max_iter = 30L),
        edf_control  = gb_edf_ctrl(edf_m),
        uq_control   = gb_uq_ctrl("none"),
        verbose      = FALSE
      )
      set.seed(1)
      y_new <- stats::rpois(nrow(dat), fitted(mod))
      ahat  <- log(mean(y_new))
      r     <- reset_gbgam(mod, y_new, edf_m)
      a0    <- if (init) ahat else 0
      expect_equal(r$learners[[1L]]$beta, a0, tolerance = 1e-6)
      expect_equal(r$learners[[1L]]$edf,
                   as.numeric(init && edf_m != "none"))
      expect_equal(r$linear.predictors, a0 + mod$offset, tolerance = 1e-6)
      expect_equal(r$edf.total, as.numeric(init))
    }
  }
})


# ---- Replays on the original response ----------------------------------------

test_that("replaying every bootstrap mode on the original y reproduces the fit", {
  mod <- fit_asr_model()
  expect_null(mod$max_iter)  # free_stopping must cope with max_iter = NULL
  for (mode in c("fixed", "free_smoothing", "free_selection",
                 "free_stopping")) {
    f <- replay_on_y(mod, mod$y, mode)
    expect_identical(f$mstop, mod$mstop, label = mode)
    expect_equal(f$metrics$k_opts[seq_len(f$mstop)],
                 mod$metrics$k_opts[seq_len(mod$mstop)], label = mode)
    expect_equal(f$coefficients, mod$coefficients, tolerance = 1e-10,
                 label = mode)
  }
  # bootstrap_gbgam() itself must not drop max_iter = NULL in free_stopping
  set.seed(4)
  b <- bootstrap_gbgam(mod, n_samples = 2L, alpha = NULL, mode = "free_stopping")
  expect_equal(dim(b$boot_coef), c(2L, length(mod$coefficients)))
})

test_that("fixed-mode replay reproduces coef and eta (initialise_intercept = TRUE)", {
  dat <- simulate_poisson_1d(n = 120L, seed = 32L)
  mod <- fit_fixed(y ~ t2(x, k = 6), dat, family = poisson(), max_iter = 30L,
                   step_size = 0.1,
                   uq_control = gb_uq_ctrl("none"))
  f <- replay_on_y(mod, mod$y, "fixed")
  expect_equal(f$coefficients, coef(mod), tolerance = 1e-10)
  expect_equal(f$linear.predictors, mod$linear.predictors, tolerance = 1e-10)
})

test_that("fixed-mode replay reproduces coef and eta (initialise_intercept = FALSE)", {
  # gbgam() itself starts the intercept at zero when
  # initialise_intercept = FALSE, as reset_gbgam() does.
  set.seed(33)
  probe_dat <- data.frame(x = stats::runif(60))
  probe_dat$y <- 50 + stats::rnorm(60)
  probe <- gbgam(
    y ~ t2(x, k = 5),
    data            = probe_dat,
    learner_control = gb_learner_ctrl(initialise_intercept = FALSE),
    stop_control    = gb_stop_ctrl("none", max_iter = 1L),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  expect_gt(abs(mean(fitted(probe)) - 50), 10)

  dat <- simulate_poisson_1d(n = 120L, seed = 32L)
  mod <- gbgam(
    y ~ t2(x, k = 6),
    data              = dat,
    family            = poisson(),
    learner_control   = gb_learner_ctrl(step_size = 0.1,
                                        initialise_intercept = FALSE),
    stop_control = gb_stop_ctrl("none", max_iter = 30L),
    edf_control  = gb_edf_ctrl("none"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  f <- replay_on_y(mod, mod$y, "fixed")
  expect_equal(f$coefficients, coef(mod), tolerance = 1e-10)
  expect_equal(f$linear.predictors, mod$linear.predictors, tolerance = 1e-10)
})

test_that("recorded_paths handles post_rhos (present, NA, absent)", {
  mod <- list(
    metrics = list(k_opts = c(2L, 3L, 2L, 3L), rhos = c(1, NA, 2, 5),
                   post_rhos = cbind(c(NA, NA, NA, NA), c(4, 3, NA, 9))),
    learners = list(list(beta = 0), list(rho = 7), list(beta = 0),
                    list(rho = -1)),
    learners_to_fit_after_each_iter = c(1L, 4L)
  )
  rec <- recorded_paths(mod, m = 3L)
  expect_identical(rec$path, c(2L, 3L, 2L))
  expect_identical(rec$rhos, c(1, 0, 2))
  expect_identical(rec$post_rhos, cbind(c(0, 0, 0), c(4, 3, -1)))

  s_fixed <- replicate_settings("fixed",
                                rec$path, rec$rhos, 3L, mod,
                                post_rhos = rec$post_rhos)
  expect_identical(s_fixed$preselected_post_rhos, rec$post_rhos)
  for (mode in c("free_smoothing", "free_selection", "free_stopping")) {
    s <- replicate_settings(mode, rec$path,
                            rec$rhos, 3L, mod, post_rhos = rec$post_rhos)
    expect_null(s$preselected_post_rhos, label = mode)
  }

  # Models fitted before post_rhos were recorded
  mod$metrics$post_rhos <- NULL
  expect_null(recorded_paths(mod, m = 3L)$post_rhos)
  mod$learners_to_fit_after_each_iter <- NULL
  expect_null(recorded_paths(mod, m = 3L)$post_rhos)
})

test_that("fixed-mode replay reproduces a fit with a post-iteration RE learner", {
  set.seed(101)
  n   <- 200L
  dat <- data.frame(x = stats::runif(n),
                    g = factor(sample(seq_len(15), n, replace = TRUE)))
  re  <- stats::rnorm(15, sd = 0.7)
  dat$y <- sin(2 * pi * dat$x) + re[dat$g] + stats::rnorm(n, sd = 0.3)
  fml <- y ~ t2(x, k = 8, bs = "cr") + t2(g, bs = "re")
  lb  <- new_learners(fml, data = dat, weights = rep(1, n))
  k_re <- which(vapply(lb$learners, function(lr) {
    grepl("t2(g)", lr$var_label, fixed = TRUE)
  }, logical(1)))
  expect_length(k_re, 1L)

  mod <- gbgam(
    fml,
    data              = dat,
    learner_control   = gb_learner_ctrl(
      step_size                       = 0.3,
      learners_to_fit_after_each_iter = k_re,
      excluded_learner_indices        = k_re
    ),
    stop_control = gb_stop_ctrl("none", max_iter = 150L),
    edf_control  = gb_edf_ctrl("exact"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  )
  expect_true(is.matrix(mod$metrics$post_rhos))
  # ASR must have relaxed the RE learner for the test to be meaningful
  lr_re <- mod$learners[[k_re]]
  expect_false(isTRUE(all.equal(lr_re$rho, lr_re$rho0)))
  expect_false(k_re %in% mod$metrics$k_opts[seq_len(mod$mstop)])

  f <- replay_on_y(mod, mod$y, "fixed")
  expect_equal(f$learners[[k_re]]$rho, lr_re$rho)
  expect_equal(f$coefficients, coef(mod), tolerance = 1e-10)
  expect_equal(f$linear.predictors, mod$linear.predictors, tolerance = 1e-10)
})
