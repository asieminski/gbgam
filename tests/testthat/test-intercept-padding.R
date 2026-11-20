# Intercept-padded learners (gb_learner_ctrl(intercept_padding = TRUE)):
# with 1^T W V = 0, [1 | V] fits the intercept and V to the same gradient.

# Fit helper: fixed number of iterations, no uncertainty quantification,
# messages silenced
fit_pad <- function(fml, dat, family = stats::gaussian(), weights = NULL,
                    lc = list(), edf = "none", max_iter = 30L,
                    step_size = 0.1) {
  suppressMessages(gbgam(
    fml, data = dat, family = family, weights = weights,
    learner_control   = do.call(gb_learner_ctrl,
                                c(list(step_size = step_size), lc)),
    stop_control = gb_stop_ctrl("none", max_iter = max_iter),
    edf_control  = gb_edf_ctrl(edf),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE
  ))
}

# Fitted values A u of learner `lr` for the gradient u (prior weights w), by
# penalised weighted least squares (does not rely on V^T W V = I)
manual_hat_fit <- function(lr, u, w) {
  X <- as.matrix(lr$design)
  if (inherits(lr, "demmler_reinsch")) {
    P <- diag(exp(lr$rho) * lr$D, ncol(X))
  } else {
    P <- matrix(0, ncol(X), ncol(X))
  }
  as.vector(X %*% solve(crossprod(X, w * X) + P, crossprod(X, w * u)))
}

# Manual intercept-padded boosting for canonical links (u = y - mu): the
# intercept and learner k are both fitted to the same u. If `path` is NULL the
# learner is selected by the padded provisional deviance, as in gbgam().
manual_padded_path <- function(mod, n_iter, path = NULL, nu = 0.1) {
  fam  <- mod$family
  w    <- mod$prior.weights
  y    <- mod$y
  eta  <- mod$intercept0 * isTRUE(mod$initialise_intercept) + mod$offset
  lrs  <- mod$learners
  ks   <- integer(n_iter)
  # Provisional selection step: shrunk by nu only if requested
  nu_sel <- if (isTRUE(mod$use_step_size_for_selection)) nu else 1
  for (m in seq_len(n_iter)) {
    u     <- y - fam$linkinv(eta)
    g_int <- rep(sum(w * u) / sum(w), length(u))
    if (is.null(path)) {
      cand <- seq_along(lrs)[-1L]
      dev  <- vapply(cand, function(k) {
        nu_k <- lrs[[k]]$nu_k %||% 1
        eh   <- eta + nu_sel * (nu_k * manual_hat_fit(lrs[[k]], u, w) + g_int)
        sum(fam$dev.resids(y, fam$linkinv(eh), w))
      }, numeric(1))
      ks[m] <- cand[which.min(dev)]
    } else {
      ks[m] <- path[m]
    }
    eta <- eta + nu * (manual_hat_fit(lrs[[ks[m]]], u, w) + g_int)
  }
  list(eta = eta, k = ks)
}


# ---- (a) Gaussian: padding is a no-op with an initialised intercept ----------

test_that("gaussian: padding == post-iteration refit == no refit", {
  dat <- simulate_gaussian_2d(n = 120L, seed = 71L)
  set.seed(71)
  w   <- stats::runif(nrow(dat), 0.5, 2)
  fml <- y ~ t2(x, k = 6) + t2(z, k = 6)
  pad <- fit_pad(fml, dat, weights = w, edf = "exact", max_iter = 60L,
                 lc = list(intercept_padding = TRUE))
  aft <- fit_pad(fml, dat, weights = w, edf = "exact", max_iter = 60L,
                 lc = list(learners_to_fit_after_each_iter = 1L,
                           excluded_learner_indices = 1L))
  non <- fit_pad(fml, dat, weights = w, edf = "exact", max_iter = 60L,
                 lc = list(excluded_learner_indices = 1L))
  for (other in list(aft, non)) {
    expect_identical(pad$metrics$k_opts, other$metrics$k_opts)
    expect_equal(pad$linear.predictors, other$linear.predictors,
                 tolerance = 1e-10)
    expect_equal(pad$edfs, other$edfs, tolerance = 1e-10)
    expect_equal(pad$metrics$rhos, other$metrics$rhos, tolerance = 1e-10)
  }
  # The padded intercept fit is (numerically) zero in every iteration
  expect_equal(pad$learners[[1L]]$beta, non$learners[[1L]]$beta,
               tolerance = 1e-10)
})


# ---- (b) Gaussian, no initial intercept: same-gradient update ----------------

test_that("gaussian, no initial intercept: matches manual same-gradient fits", {
  dat <- simulate_gaussian_2d(n = 100L, seed = 72L)
  dat$y <- dat$y + 3
  set.seed(72)
  w    <- stats::runif(nrow(dat), 0.5, 2)
  fml  <- y ~ t2(x, k = 6) + t2(z, k = 6)
  path <- c(2L, 3L, 5L, 3L, 4L)
  mod  <- fit_pad(fml, dat, weights = w, max_iter = 5L,
                  lc = list(intercept_padding = TRUE, tau = 0,
                            initialise_intercept = FALSE,
                            preselected_learner_indices = path))
  man  <- manual_padded_path(mod, n_iter = 5L, path = path)
  expect_equal(mod$linear.predictors, man$eta, tolerance = 1e-10)
  # The intercept moved towards the mean of y (it is not initialised): by
  # more than half of one step from 0
  expect_gt(mod$learners[[1L]]$beta, 0.5 * 0.1 * mean(dat$y))
  expect_equal(as.vector(stats::predict(mod, newdata = dat, type = "link")),
               mod$linear.predictors, tolerance = 1e-10)
})

test_that("single-learner model: padding equals the explicit [1 | V] learner", {
  dat <- simulate_gaussian_1d(n = 80L, seed = 73L)
  dat$y <- dat$y + 2
  set.seed(73)
  w    <- stats::runif(nrow(dat), 0.5, 2)
  fml  <- y ~ t2(x, k = 8)
  probe <- fit_pad(fml, dat, weights = w, max_iter = 1L,
                   lc = list(intercept_padding = TRUE, tau = 0,
                             initialise_intercept = FALSE))
  labs <- vapply(probe$learners, `[[`, character(1), "var_label")
  k_r  <- which(labs == "t2(x)r")
  n_it <- 20L
  mod  <- fit_pad(fml, dat, weights = w, max_iter = n_it,
                  lc = list(intercept_padding = TRUE, tau = 0,
                            initialise_intercept = FALSE,
                            preselected_learner_indices = rep(k_r, n_it)))
  # Explicit padded hat matrix [1 | V] with an unpenalised intercept column
  lr   <- mod$learners[[k_r]]
  Xp   <- cbind(1, lr$design)
  P    <- diag(c(0, exp(lr$rho) * lr$D))
  Apad <- Xp %*% solve(crossprod(Xp, w * Xp) + P, t(Xp * w))
  eta  <- rep(0, nrow(dat))
  for (m in seq_len(n_it)) eta <- eta + 0.1 * as.vector(Apad %*% (dat$y - eta))
  expect_equal(mod$linear.predictors, eta, tolerance = 1e-10)
})


# ---- (c) Poisson and binomial: update, selection, refit after ----------------

test_that("poisson/binomial: padded update equals manual same-gradient path", {
  dat_p <- simulate_poisson_1d(n = 100L, seed = 74L)
  set.seed(74)
  dat_p$z <- stats::runif(nrow(dat_p))
  dat_p$o <- 0.3 * stats::rnorm(nrow(dat_p))
  dat_b <- simulate_binomial_prop(n = 100L, seed = 75L)
  dat_b$z <- stats::runif(nrow(dat_b))
  cases <- list(
    list(dat = dat_p, family = stats::poisson(), w = NULL,
         fml = y ~ t2(x, k = 6) + t2(z, k = 5) + offset(o)),
    list(dat = dat_b, family = stats::binomial(), w = dat_b$w,
         fml = y ~ t2(x, k = 6) + t2(z, k = 5)))
  for (cs in cases) {
    # Fixed path
    path <- c(3L, 2L, 3L, 5L, 3L, 4L)
    lc   <- list(tau = 0, preselected_learner_indices = path)
    pad  <- fit_pad(cs$fml, cs$dat, cs$family, cs$w, max_iter = 6L,
                    lc = c(lc, intercept_padding = TRUE))
    man  <- manual_padded_path(pad, n_iter = 6L, path = path)
    expect_equal(pad$linear.predictors, man$eta, tolerance = 1e-10)
    aft  <- fit_pad(cs$fml, cs$dat, cs$family, cs$w, max_iter = 6L,
                    lc = c(lc, learners_to_fit_after_each_iter = 1L))
    # The post-update refit uses a re-evaluated gradient, so it differs (by
    # a term of second order in the step size when the intercept is
    # initialised at the MLE)
    expect_gt(max(abs(pad$linear.predictors - aft$linear.predictors)), 1e-8)

    # Data-driven selection by the padded provisional deviance
    sel <- fit_pad(cs$fml, cs$dat, cs$family, cs$w, max_iter = 8L,
                   lc = list(tau = 0, intercept_padding = TRUE,
                             selection_loss_type = "deviance",
                             use_step_size_for_selection = TRUE))
    man <- manual_padded_path(sel, n_iter = 8L)
    expect_identical(sel$metrics$k_opts, man$k)
    expect_false(1L %in% sel$metrics$k_opts)
    expect_equal(sel$linear.predictors, man$eta, tolerance = 1e-10)
    expect_equal(as.vector(stats::predict(sel, newdata = cs$dat, type = "link")),
                 sel$linear.predictors, tolerance = 1e-10)
  }
})

test_that("rss selection is unaffected by the padded intercept", {
  dat <- simulate_poisson_1d(n = 100L, seed = 76L)
  set.seed(76)
  dat$z <- stats::runif(nrow(dat))
  mod <- fit_pad(y ~ t2(x, k = 6) + t2(z, k = 5), dat, stats::poisson(),
                 max_iter = 10L,
                 lc = list(tau = 0, intercept_padding = TRUE,
                           selection_loss_type = "rss"))
  # Replay: at each iteration the selected learner minimises the unpadded
  # weighted RSS of the gradient
  w      <- mod$prior.weights
  eta    <- mod$intercept0 + mod$offset
  nu_sel <- if (isTRUE(mod$use_step_size_for_selection)) 0.1 else 1
  for (m in 1:10) {
    u   <- mod$y - exp(eta)
    rss <- vapply(2:length(mod$learners), function(k) {
      lr <- mod$learners[[k]]
      sum(w * (u - nu_sel * (lr$nu_k %||% 1) * manual_hat_fit(lr, u, w))^2)
    }, numeric(1))
    k   <- 1L + which.min(rss)
    expect_identical(mod$metrics$k_opts[m], k)
    eta <- eta + 0.1 * (manual_hat_fit(mod$learners[[k]], u, w) +
                          sum(w * u) / sum(w))
  }
  expect_equal(mod$linear.predictors, eta, tolerance = 1e-10)
})


# ---- (d) EDF: exact trace equals the FD trace of dmu/dy ----------------------

test_that("padding: exact dmu/dy and EDF match finite differences", {
  skip_on_cran()
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  set.seed(77)
  n <- 25
  d <- data.frame(x = sort(stats::runif(n)), z = stats::runif(n),
                  o = 0.5 * stats::rnorm(n), w = stats::runif(n, 0.5, 2))
  # Non-integer counts keep y +- eps valid (the resulting warning is expected)
  d$y <- stats::rpois(n, exp(0.3 + 0.8 * sin(2 * pi * d$x) + d$o)) + 0.2
  fml <- y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 5, bs = "cr") + offset(o)
  eps <- 1e-4
  n_it <- 30L
  fit_mu <- function(dd, lc, edf = "none") {
    suppressWarnings(fit_pad(fml, dd, stats::poisson(), d$w, lc = lc,
                             edf = edf, max_iter = n_it))
  }
  for (init in c(TRUE, FALSE)) {
    lc  <- list(tau = 0, intercept_padding = TRUE, initialise_intercept = init)
    ref <- fit_mu(d, lc, edf = "exact")
    lc1 <- c(lc, list(preselected_learner_indices = ref$metrics$k_opts))
    expect_identical(fit_mu(d, lc1)$linear.predictors, ref$linear.predictors)
    js  <- if (init) seq_len(n) else 1:6
    fd  <- vapply(js, function(j) {
      dp <- d
      dp$y[j] <- dp$y[j] + eps
      dm <- d
      dm$y[j] <- dm$y[j] - eps
      (fit_mu(dp, lc1)$fitted.values - fit_mu(dm, lc1)$fitted.values) /
        (2 * eps)
    }, numeric(n))
    expect_equal(unname(ref$derivs$dmudy[, js]), unname(fd), tolerance = 1e-6)
    # The tracked EDF is the trace of the tracked Jacobian ...
    expect_equal(ref$edfs[ref$mstop + 1L], sum(diag(ref$derivs$dmudy)),
                 tolerance = 1e-10)
    # ... which equals the finite-difference trace
    if (init) expect_equal(ref$edfs[ref$mstop + 1L], sum(diag(fd)),
                           tolerance = 1e-6)
    # The learner EDFs add up to the model EDF
    expect_equal(sum(vapply(ref$learners, `[[`, numeric(1), "edf")),
                 ref$edfs[ref$mstop + 1L], tolerance = 1e-10)
    # The intercept is credited with dedf_int, the trace of its own increment
    # of dmu/dy, tr(diag(h'(eta)) 1 d(Delta_int)/dy), in each iteration:
    # compare the last iteration with a central FD of that increment
    prev <- suppressWarnings(fit_pad(fml, d, stats::poisson(), d$w, lc = lc,
                                     edf = "exact", max_iter = n_it - 1L))
    lc0  <- c(lc, list(preselected_learner_indices = ref$metrics$k_opts[-n_it]))
    expect_identical(
      suppressWarnings(fit_pad(fml, d, stats::poisson(), d$w, lc = lc0,
                               max_iter = n_it - 1L))$linear.predictors,
      prev$linear.predictors)
    int_step <- function(dd) {  # the intercept's increment in iteration n_it
      b1 <- suppressWarnings(fit_pad(fml, dd, stats::poisson(), d$w, lc = lc1,
                                     max_iter = n_it))$learners[[1L]]$beta
      b0 <- suppressWarnings(fit_pad(fml, dd, stats::poisson(), d$w, lc = lc0,
                                     max_iter = n_it - 1L))$learners[[1L]]$beta
      b1 - b0
    }
    d_int <- vapply(seq_len(n), function(j) {
      dp <- d
      dp$y[j] <- dp$y[j] + eps
      dm <- d
      dm$y[j] <- dm$y[j] - eps
      (int_step(dp) - int_step(dm)) / (2 * eps)
    }, numeric(1))
    share_fd <- sum(exp(ref$linear.predictors) * d_int)
    expect_equal(ref$learners[[1L]]$edf - prev$learners[[1L]]$edf,
                 share_fd, tolerance = 1e-6)
    # ... and the selected learner with the rest of the iteration's change
    k_last <- ref$metrics$k_opts[n_it]
    expect_equal(ref$learners[[k_last]]$edf - prev$learners[[k_last]]$edf,
                 ref$edfs[n_it + 1L] - ref$edfs[n_it] - share_fd,
                 tolerance = 1e-6)
  }
})

test_that("padding: gaussian EDF split gives init-invariant ASR (tau = 1)", {
  dat <- simulate_gaussian_2d(n = 120L, seed = 81L)
  dat$y <- dat$y + 3
  set.seed(81)
  w   <- stats::runif(nrow(dat), 0.5, 2)
  fml <- y ~ t2(x, k = 8) + t2(z, k = 6)
  n_it <- 120L
  fits <- lapply(c(TRUE, FALSE), function(init) {
    fit_pad(fml, dat, weights = w, edf = "exact", max_iter = n_it,
            lc = list(intercept_padding = TRUE, tau = 1,
                      initialise_intercept = init))
  })
  a <- fits[[1L]]
  b <- fits[[2L]]
  # The smooth parts evolve identically: same path, rhos, coefficients and
  # learner EDFs; only the intercept differs
  expect_identical(a$metrics$k_opts, b$metrics$k_opts)
  expect_equal(a$metrics$rhos, b$metrics$rhos, tolerance = 1e-10)
  smooth_par <- function(m) {
    unlist(lapply(m$learners[-1L], function(lr) c(lr$beta, lr$edf)))
  }
  expect_equal(smooth_par(a), smooth_par(b), tolerance = 1e-10)
  # Intercept EDF: 1 when initialised, 1 - (1 - nu)^m when boosted from 0,
  # which is the trace of the intercept's d(beta)/dy
  expect_equal(a$learners[[1L]]$edf, 1, tolerance = 1e-10)
  expect_equal(b$learners[[1L]]$edf, 1 - 0.9^n_it, tolerance = 1e-10)
  # A post-iteration intercept refit coincides with padding (gaussian)
  aft <- fit_pad(fml, dat, weights = w, edf = "exact", max_iter = n_it,
                 lc = list(learners_to_fit_after_each_iter = 1L,
                           excluded_learner_indices = 1L, tau = 1,
                           initialise_intercept = FALSE))
  expect_identical(aft$metrics$k_opts, b$metrics$k_opts)
  expect_equal(aft$linear.predictors, b$linear.predictors, tolerance = 1e-10)
  expect_equal(aft$edfs, b$edfs, tolerance = 1e-10)
  expect_equal(vapply(aft$learners, `[[`, numeric(1), "edf"),
               vapply(b$learners, `[[`, numeric(1), "edf"), tolerance = 1e-10)
  # Hutch++: the intercept share is estimated from the sketch, the learner
  # EDFs still add up to the model EDF, and padding still coincides with the
  # post-iteration refit (same probes)
  set.seed(3)
  h <- fit_pad(fml, dat, weights = w, edf = "hutchpp", max_iter = n_it,
               lc = list(intercept_padding = TRUE, tau = 1,
                         initialise_intercept = FALSE))
  expect_equal(sum(vapply(h$learners, `[[`, numeric(1), "edf")),
               h$edfs[n_it + 1L], tolerance = 1e-10)
  expect_gt(h$learners[[1L]]$edf, 0.5)
  set.seed(3)
  h_aft <- fit_pad(fml, dat, weights = w, edf = "hutchpp", max_iter = n_it,
                   lc = list(learners_to_fit_after_each_iter = 1L,
                             excluded_learner_indices = 1L, tau = 1,
                             initialise_intercept = FALSE))
  expect_identical(h_aft$metrics$k_opts, h$metrics$k_opts)
  expect_equal(h_aft$metrics$rhos, h$metrics$rhos, tolerance = 1e-8)
  expect_equal(h_aft$linear.predictors, h$linear.predictors,
               tolerance = 1e-8)
})


# ---- (e) Hutch++ sketch ------------------------------------------------------

test_that("padding: hutchpp sketch equals the exact Jacobian times the probes", {
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  dat <- simulate_poisson_1d(n = 80L, seed = 78L)
  set.seed(78)
  dat$z <- stats::runif(nrow(dat))
  fml <- y ~ t2(x, k = 8) + t2(z, k = 5)
  fit_ec <- function(edf) {
    set.seed(1)
    fit_pad(fml, dat, stats::poisson(), edf = edf, max_iter = 25L,
            lc = list(tau = 0, intercept_padding = TRUE))
  }
  ex <- fit_ec("exact")
  hp <- fit_ec("hutchpp")
  expect_identical(ex$metrics$k_opts, hp$metrics$k_opts)
  expect_equal(ex$derivs$dmudy %*% hp$derivs$GGG, hp$derivs$dmudy,
               tolerance = 1e-10)
  expect_true(all(is.finite(hp$edfs)))

  # Default AIC stopping with ASR runs end to end
  set.seed(2)
  mod <- suppressMessages(gbgam(
    fml, data = dat, family = stats::poisson(),
    learner_control   = gb_learner_ctrl(step_size = 0.1,
                                        intercept_padding = TRUE),
    stop_control = gb_stop_ctrl("aic", iter_no_improv = 10L),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE))
  expect_model_contract(mod, nrow(dat))
  expect_false(1L %in% mod$metrics$k_opts)
})


# ---- (f) Validation ----------------------------------------------------------

test_that("intercept_padding is validated", {
  expect_false(gb_learner_ctrl()$intercept_padding)
  expect_true(gb_learner_ctrl(intercept_padding = TRUE)$intercept_padding)
  expect_error(gb_learner_ctrl(intercept_padding = NA), "intercept_padding")
  expect_error(gb_learner_ctrl(intercept_padding = "yes"), "intercept_padding")
  expect_error(gb_learner_ctrl(intercept_padding = TRUE,
                               learners_to_fit_after_each_iter = 1L),
               "learners_to_fit_after_each_iter")
  expect_error(gb_learner_ctrl(intercept_padding = TRUE,
                               preselected_learner_indices = c(2L, 1L)),
               "preselected_learner_indices")
  # Learner 1 must be the intercept, and there must be something to pad
  expect_error(check_intercept_padding(list(new_learner("x"),
                                            new_learner("z"))),
               "\"\\(Intercept\\)\"")
  expect_error(check_intercept_padding(list(new_learner("(Intercept)"))),
               "at least one learner")
  dat <- simulate_gaussian_1d(n = 40L, seed = 79L)
  expect_error(fit_pad(y ~ 1, dat, lc = list(intercept_padding = TRUE)),
               "at least one learner")
  # Excluding every learner besides the intercept leaves nothing to select
  expect_error(fit_pad(y ~ t2(x, k = 6), dat,
                       lc = list(intercept_padding = TRUE,
                                 excluded_learner_indices = 2:3)),
               "No base learner is left to select.*intercept_padding")
  expect_error(fit_pad(y ~ t2(x, k = 6), dat,
                       lc = list(excluded_learner_indices = 1:3)),
               "No base learner is left to select")
  # A preselected path does not need candidates
  expect_silent(fit_pad(y ~ t2(x, k = 6), dat, max_iter = 2L,
                        lc = list(intercept_padding = TRUE,
                                  excluded_learner_indices = 2:3,
                                  preselected_learner_indices = c(2L, 3L))))
})


# ---- (g) Bootstrap: the flag is carried by reset_gbgam() refits --------------

test_that("padding: fixed and free_smoothing bootstrap replays reproduce the fit", {
  dat <- simulate_poisson_1d(n = 100L, seed = 80L)
  set.seed(80)
  dat$z <- stats::runif(nrow(dat))
  mod <- suppressMessages(gbgam(
    y ~ t2(x, k = 8) + t2(z, k = 5), data = dat, family = stats::poisson(),
    learner_control   = gb_learner_ctrl(step_size = 0.1,
                                        intercept_padding = TRUE),
    stop_control = gb_stop_ctrl("aic", iter_no_improv = 15L),
    edf_control  = gb_edf_ctrl("exact"),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE))
  expect_true(mod$intercept_padding)
  for (mode in c("fixed", "free_smoothing")) {
    f <- suppressMessages(replay_on_y(mod, mod$y, mode))
    expect_true(f$intercept_padding)
    expect_equal(f$coefficients, coef(mod), tolerance = 1e-8, label = mode)
    expect_equal(f$linear.predictors, mod$linear.predictors,
                 tolerance = 1e-8, label = mode)
  }
  set.seed(5)
  b <- suppressMessages(bootstrap_gbgam(mod, n_samples = 3L, alpha = NULL,
                                        mode = "fixed"))
  expect_equal(dim(b$boot_coef), c(3L, length(coef(mod))))
})
