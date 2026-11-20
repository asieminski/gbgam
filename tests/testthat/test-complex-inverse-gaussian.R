# Complex inverse Gaussian model with the canonical 1/mu^2 link: response
# times (s) in a perception task. The inverse Gaussian is the first-passage
# time of a drift-diffusion process, and 1/mu^2 (the squared drift) rises
# with stimulus contrast, falls for incongruent trials and with the squared
# horizontal eccentricity of the stimulus (not with its side), changes over
# the trial number (a practice trend plus a rise and fall) and is slightly
# higher for left-handed participants. Each trial has a non-integer
# reliability weight.
# Covers: factor, logical, numeric and I(z^2) learners next to a P-spline
# t2() smooth (penalised range space and null-space column), all of them
# boosted, and a numeric learner with no true effect, not selected with this
# seed; the intercept refitted after every iteration and excluded from
# selection; ASR
# with tau = 1.5; true_f stopping; Hutch++ edf~ (asserted against the exact
# trace, as the link is canonical); a public-API replay of the selection
# and rho paths; the free_smoothing bootstrap without widening, with
# smooth_intervals() unwidened and widened afterwards; Bayesian backfitting
# with kept and free linear learners (also on noisier response times, where
# Metropolis-Hastings candidates can leave the range of the link); and the
# initialise_intercept check of the 1/mu^2 link.

# Response-time data: eta = 1/mu^2 in about (1.07, 2.99), dispersion `phi`
sim_response_times <- function(n = 200, seed = 601, phi = 0.08) {
  set.seed(seed)
  d <- data.frame(
    cond     = factor(sample(c("congruent", "incongruent"), n, TRUE)),
    contrast = runif(n, 0.05, 1),
    trial    = runif(n, 1, 300),
    hand     = runif(n) < 0.3,     # left-handed
    rel      = runif(n, 0.5, 1.5), # reliability weight of the trial
    ecc      = runif(n, -1, 1)     # horizontal eccentricity of the stimulus
  )
  d$eta <- 1 + 1.5 * sqrt(d$contrast) - 0.4 * (d$cond == "incongruent") -
    0.5 * d$ecc^2 + 0.3 * sin(pi * d$trial / 300) + 0.4 * d$trial / 300 +
    0.1 * d$hand
  d$rt  <- simulate_response(stats::inverse.gaussian(), 1 / sqrt(d$eta),
                             d$rel, d$rel, phi = phi)
  d
}

rt_formula <- rt ~ cond + hand + contrast + ecc + I(ecc^2) +
  t2(trial, k = 6, bs = "ps")

# The fit takes several seconds, so the tests share it
rt_cache <- new.env()
rt_fit <- function() {
  if (is.null(rt_cache$mod)) {
    d <- sim_response_times(seed = 603)
    set.seed(11)
    rt_cache$mod <- gbgam(
      rt_formula, family = stats::inverse.gaussian(), data = d,
      weights = rel,
      learner_control = gb_learner_ctrl(
        step_size = 0.2, tau = 1.5,
        learners_to_fit_after_each_iter = 1L,
        excluded_learner_indices = 1L
      ),
      stop_control = gb_stop_ctrl("true_f", true_f = d$eta,
                                  iter_no_improv = 30, max_iter = 2000),
      edf_control = gb_edf_ctrl("hutchpp", z = 60),
      uq_control = gb_param_boot_ctrl(n_samples = 8,
                                      mode = "free_smoothing",
                                      alpha = NULL),
      verbose = FALSE
    )
    rt_cache$data <- d
  }
  rt_cache
}

test_that("the response-time model passes the invariants and recovers the truth", {
  skip_on_cran()
  fit <- rt_fit()
  mod <- fit$mod
  d   <- fit$data

  labels <- vapply(mod$learners, `[[`, "", "var_label")
  expect_identical(labels, c("(Intercept)", "condincongruent", "handTRUE",
                             "contrast", "ecc", "I(ecc^2)", "t2(trial)r",
                             "t2(trial)n1"))
  expect_gt(mod$mstop, 0L)
  expect_lt(mod$mstop, 2000L)

  # Hutch++ edf~ within 15% of the exact trace on the same path
  expect_complex_invariants(mod, d, truth_eta = d$eta, min_cor = 0.9,
                            hutchpp_rtol = 0.15)

  # The intercept is never selected; it is refitted after every iteration
  # (its post-iteration rho column is NA: it is unpenalised)
  k_path <- mod$metrics$k_opts[seq_len(mod$mstop)]
  expect_false(any(k_path == 1L))
  expect_true(all(k_path %in% 2:8))
  # Every unpenalised learner type and the null-space column is boosted
  expect_true(all(c(2:4, 6:8) %in% k_path))
  post_rhos <- mod$metrics$post_rhos
  expect_identical(ncol(post_rhos), 1L)
  expect_true(all(is.na(post_rhos)))
  # The refit keeps the intercept's score sum(w (y - mu)) near zero (about
  # 10 times further from zero without it)
  w <- mod$prior.weights
  r <- mod$y - mod$fitted.values
  expect_lt(abs(sum(w * r)) / sum(w * abs(r)), 0.02)

  # ASR only lowers the smoothing parameter of t2(trial)r
  rho_trial <- mod$metrics$rhos[seq_len(mod$mstop)][k_path == 7L]
  expect_gt(length(rho_trial), 1L)
  expect_true(all(diff(rho_trial) <= 0))
  expect_lt(rho_trial[length(rho_trial)], rho_trial[1])

  # Dispersion near the simulated 0.08
  expect_equal(mod$sig2, mod$scale)
  expect_gt(mod$scale, 0.06)
  expect_lt(mod$scale, 0.11)

  # Effects: incongruent trials lower the drift (truth -0.4, shrunk by
  # boosting); the drift rises with contrast and falls with the squared
  # eccentricity, on both sides; the trial smooth follows
  # 0.3 sin(pi trial / 300) + 0.4 trial / 300
  b_cond <- unname(mod$coefficients["condincongruent"])
  expect_gt(b_cond, -0.55)
  expect_lt(b_cond, -0.15)
  expect_gt(unname(mod$coefficients["contrast"]), 0)
  expect_lt(unname(mod$coefficients["I(ecc^2)"]), 0)
  grid <- data.frame(cond = factor("congruent", levels(d$cond)),
                     hand = FALSE, contrast = 0.5,
                     ecc = seq(-1, 1, length.out = 21), trial = 150)
  p_ecc <- as.vector(predict(mod, newdata = grid))
  expect_true(all(diff(p_ecc[1:11]) > 0))
  expect_true(all(diff(p_ecc[11:21]) < 0))
  f_trial <- predict(mod, newdata = d, type = "terms")[, "t2(trial)"]
  expect_gt(cor(f_trial, 0.3 * sin(pi * d$trial / 300) + 0.4 * d$trial / 300),
            0.9)
})

test_that("replaying the recorded paths through gb_learner_ctrl() reproduces the response-time fit", {
  skip_on_cran()
  fit <- rt_fit()
  mod <- fit$mod
  d   <- fit$data
  M   <- mod$mstop
  rec <- recorded_paths(mod, M)
  expect_length(rec$rhos, M)
  expect_false(anyNA(rec$rhos))

  # tau is left out: it would be ignored with a warning without EDF
  mod_pre <- gbgam(
    rt_formula, family = stats::inverse.gaussian(), data = d, weights = rel,
    learner_control = gb_learner_ctrl(
      step_size = 0.2,
      learners_to_fit_after_each_iter = 1L,
      excluded_learner_indices = 1L,
      preselected_learner_indices = rec$path,
      preselected_rhos = rec$rhos
    ),
    stop_control = gb_stop_ctrl("none", max_iter = M),
    edf_control = gb_edf_ctrl("none"),
    verbose = FALSE
  )
  expect_identical(mod_pre$mstop, M)
  expect_equal(mod_pre$coefficients, mod$coefficients, tolerance = 1e-8)
  expect_equal(mod_pre$linear.predictors, mod$linear.predictors,
               tolerance = 1e-8)
  expect_true(is.na(mod_pre$edf.total))
})

test_that("the unwidened free_smoothing bootstrap gives smooth intervals, widened on request", {
  skip_on_cran()
  fit <- rt_fit()
  mod <- fit$mod
  d   <- fit$data

  expect_identical(mod$uq_control$method, "bootstrap")
  expect_null(mod$widening_factor)
  expect_identical(dim(mod$boot_coef),
                   c(8L, length(mod$coefficients)))
  expect_true(all(is.finite(mod$boot_coef)))

  # Unwidened intervals: centred on the intercept plus the trial smooth
  ci <- smooth_intervals(mod, correction = "none")
  expect_null(ci$factor)
  expect_identical(colnames(ci$fit), "t2(trial)")
  expect_identical(nrow(ci$fit), nrow(d))
  pt <- predict(mod, newdata = d, type = "terms")
  expect_equal(as.vector(ci$fit),
               unname(pt[, "t2(trial)"] + mod$coefficients["(Intercept)"]),
               tolerance = 1e-8)
  expect_true(all(ci$se > 0))
  expect_true(all(ci$lower < ci$fit & ci$fit < ci$upper))
  # The default correction is "none" for a bootstrap with alpha = NULL
  expect_equal(smooth_intervals(mod)$se, ci$se)

  # Widened by the model-wide factor found from the replicates
  ci_m <- smooth_intervals(mod, correction = "model")
  expect_length(ci_m$factor, 1L)
  expect_gt(ci_m$factor, 0)
  expect_equal(ci_m$fit, ci$fit)
  expect_equal(ci_m$se, ci$se * ci_m$factor, tolerance = 1e-8)
})

test_that("Bayesian backfitting with kept linear learners gives posterior intervals for the response-time model", {
  skip_on_cran()
  mod <- rt_fit()$mod
  set.seed(4)
  mb  <- bayesian_backfitting(mod, n_samples = 200, n_burnin = 50,
                              linear = "kept")
  expect_identical(mb$uq_control$method, "bayes")
  expect_null(mb$boot_coef)
  expect_length(mb$coefficients_bayes, length(mod$coefficients))
  expect_true(all(is.finite(mb$coefficients_bayes)))
  expect_true(all(is.finite(mb$Vp)))

  # Learners boosting never selected are held at zero
  never <- setdiff(2:8, mod$metrics$k_opts[seq_len(mod$mstop)])
  expect_gt(length(never), 0L)
  for (k in never) {
    expect_true(all(mb$beta_samples[[k]] == 0))
    expect_true(is.na(mb$mh_accept_rate[k]))
  }

  si <- smooth_intervals(mb, centre = "posterior")
  expect_identical(si$centre, "posterior")
  expect_null(si$factor)
  expect_true(all(si$se > 0))
  expect_true(all(si$lower < si$fit & si$fit < si$upper))
  # The posterior mean and the boosting estimate agree on the trial smooth
  ci <- smooth_intervals(mb)
  expect_gt(cor(as.vector(si$fit), as.vector(ci$fit)), 0.95)
  expect_true(all(abs(si$fit - ci$fit) < si$se))
})

test_that("Bayesian backfitting with free linear learners samples the response-time model", {
  skip_on_cran()
  mod <- rt_fit()$mod
  set.seed(1)
  mb <- bayesian_backfitting(mod, n_samples = 300, n_burnin = 50,
                             linear = "free")
  expect_s3_class(mb, "gbgam")
  expect_length(mb$coefficients_bayes, length(mod$coefficients))
  expect_true(all(is.finite(mb$coefficients_bayes)))
  expect_gt(mb$bb_info$linear_accept, 0)
  expect_lt(mb$bb_info$linear_accept, 1)
})

test_that("Bayesian backfitting with free linear learners rejects candidates outside the 1/mu^2 range", {
  skip_on_cran()
  # Noisier response times (phi = 0.2): the posterior of the linear
  # learners is wide enough for candidates with eta <= 0 in some rows
  d   <- sim_response_times(phi = 0.2)
  mod <- gbgam(
    rt_formula, family = stats::inverse.gaussian(), data = d, weights = rel,
    learner_control = gb_learner_ctrl(step_size = 0.2),
    stop_control = gb_stop_ctrl("aic", iter_no_improv = 30),
    edf_control = gb_edf_ctrl("exact"),
    verbose = FALSE
  )
  # The fit itself lies well inside the range of the link
  expect_gt(min(mod$linear.predictors), 0.5)
  expect_gt(cor(as.vector(mod$linear.predictors), d$eta), 0.8)

  set.seed(1)
  expect_no_error(
    bayesian_backfitting(mod, n_samples = 300, n_burnin = 50,
                         linear = "kept")
  )
  # Such a candidate has posterior density 0: Metropolis-Hastings rejects
  # it and keeps the current state, without warnings from evaluating the
  # model at the candidate
  set.seed(1)
  mb    <- NULL
  warns <- character()
  expect_no_error(withCallingHandlers(
    mb <- bayesian_backfitting(mod, n_samples = 300, n_burnin = 50,
                               linear = "free"),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ))
  expect_identical(warns, character())
  expect_length(mb$coefficients_bayes, length(mod$coefficients))
  expect_true(all(is.finite(mb$coefficients_bayes)))
  expect_gt(mb$bb_info$linear_accept, 0)
  expect_lt(mb$bb_info$linear_accept, 1)
})

test_that("initialise_intercept = FALSE is rejected for the 1/mu^2 link without a positive offset", {
  skip_on_cran()
  set.seed(5)
  d <- sim_response_times(n = 60)
  expect_error(
    gbgam(rt_formula, family = stats::inverse.gaussian(), data = d,
          weights = rel,
          learner_control = gb_learner_ctrl(initialise_intercept = FALSE),
          stop_control = gb_stop_ctrl("none", max_iter = 5),
          edf_control = gb_edf_ctrl("none"),
          verbose = FALSE),
    "initialise_intercept"
  )
})
