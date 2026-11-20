# Data simulators, expectations and fit wrappers shared by the test files
# (testthat sources helper-*.R before the tests).


# ---- Expectations and model reconstruction -----------------------------------

# Evaluate expr without the warning about widening factors below 0.9, which
# bootstraps with few replicates give through Monte Carlo noise
without_low_widening_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("widening factors? (= [0-9.]+ )?<", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}

# Fields, classes and finiteness every fitted gbgam model must have
expect_model_contract <- function(mod, n) {
  required_fields <- c(
    "call", "family", "model", "terms", "pterms",
    "weights", "offset", "linear.predictors",
    "learners", "coef_names", "coefficients",
    "cen_const", "ind0var", "metrics",
    "elapsed_secs", "mstop"
  )

  testthat::expect_s3_class(mod, "gbgam")
  testthat::expect_true(inherits(mod, "gam"))
  testthat::expect_true(inherits(mod, "glm"))
  testthat::expect_true(inherits(mod, "lm"))

  testthat::expect_true(all(required_fields %in% names(mod)))
  testthat::expect_length(mod$linear.predictors, n)
  testthat::expect_length(mod$fitted.values, n)
  testthat::expect_true(all(is.finite(mod$linear.predictors)))
  testthat::expect_true(all(is.finite(mod$fitted.values)))
  testthat::expect_true(all(is.finite(mod$coefficients)))
  testthat::expect_true(length(mod$coefficients) == length(mod$coef_names))
  testthat::expect_true(is.list(mod$learners))
  testthat::expect_true(length(mod$learners) >= 1L)
}

# Full model matrix [parametric | smooths] with mgcv's column names
create_concat_design <- function(mod) {
  pX <- stats::model.matrix(mod$pterms, mod$model)
  smooth_X_list <- lapply(mod$smooth, function(sm) {
    X <- sm$X
    # Centred smooths predict from the centred basis (see set_smooth_centring)
    if (inherits(sm, "gbgam_centred")) X <- sweep(X, 2, sm$gbgam_cen)
    colnames(X) <- paste0(sm$label, ".", seq_len(ncol(X)))
    X
  })
  if (length(smooth_X_list) == 0L) return(pX)
  cbind(pX, do.call(cbind, smooth_X_list))
}

# Training linear predictor rebuilt from the mgcv-scale coefficients
rebuild_training_eta <- function(mod) {
  X    <- create_concat_design(mod)
  beta <- mod$coefficients[colnames(X)]
  as.vector(mod$offset + X %*% beta)
}

# Rows of `data` that the fit used (after subset and na.action): `data`
# must have the default row names, which model.frame() keeps. For a list
# with matrix elements (linear functional smooths) the matrices lose rows
# and the vectors elements.
training_rows <- function(mod, data) {
  idx <- as.integer(rownames(mod$model))
  if (is.data.frame(data)) return(data[idx, , drop = FALSE])
  lapply(data, function(v) {
    if (is.matrix(v)) v[idx, , drop = FALSE] else v[idx]
  })
}

# Family-agnostic invariants of a fitted model, for the complex-model tests.
# `data` is the data the model was fitted to (all rows, before subset and
# na.action). Checks, with `tol` for the exact identities:
# * the model contract (expect_model_contract());
# * predict() at the training rows, with the offset argument added back
#   (mgcv::predict.gam() leaves it out for newdata), equals the linear
#   predictor, on the link and the response scale; predict() without newdata
#   is padded with NA only at the rows excluded by na.exclude;
# * the mgcv-scale coefficients rebuild the linear predictor;
# * the terms, their constant and the offset add up to the linear predictor;
# * with EDF tracking, the learners' edf~_k add up to the edf~ tracked at
#   mstop and to edf.total;
# * a fixed-mode replay of the recorded selection, rho and post-iteration
#   rho paths on the observed y reproduces the fit;
# * fitted() and every residual type follow na.action;
# * with `truth_eta` (the true linear predictor at the training rows), its
#   correlation with the fit is at least `min_cor`;
# * with `hutchpp_rtol` (canonical links and Hutch++ EDF only), the
#   relative error of the Hutch++ edf~ against the exact trace on the same
#   path is at most `hutchpp_rtol`.
expect_complex_invariants <- function(mod, data, truth_eta = NULL,
                                      min_cor = 0.9, tol = 1e-8,
                                      hutchpp_rtol = NULL) {
  n   <- nrow(mod$model)
  eta <- as.vector(mod$linear.predictors)
  expect_model_contract(mod, n)

  # Predictions at the training rows
  nd      <- training_rows(mod, data)
  arg_off <- mod$model[["(offset)"]] %||% 0
  testthat::expect_equal(as.vector(predict(mod, newdata = nd)) + arg_off, eta,
                         tolerance = tol, label = "predict(newdata) + offset")
  testthat::expect_equal(
    as.vector(predict(mod, newdata = nd, type = "response")),
    as.vector(mod$family$linkinv(eta - arg_off)),
    tolerance = tol, label = "predict(newdata, type = \"response\")")

  # Without newdata: NA padding only under na.exclude. As in
  # mgcv::predict.gam(), the offset argument is then included only when the
  # formula has an offset() term as well
  is_excl <- inherits(mod$na.action, "exclude")
  n_full  <- n + if (is_excl) length(mod$na.action) else 0L
  p0      <- as.vector(predict(mod))
  testthat::expect_length(p0, n_full)
  kept <- if (is_excl) -mod$na.action else seq_len(n)
  if (is_excl) testthat::expect_true(all(is.na(p0[mod$na.action])))
  fml_off <- !is.null(attr(mod$smoothers, "offset"))
  testthat::expect_equal(p0[kept] + if (fml_off) 0 else arg_off, eta,
                         tolerance = tol, label = "predict() at the kept rows")

  # Coefficients and terms
  testthat::expect_equal(rebuild_training_eta(mod), eta, tolerance = tol,
                         label = "rebuild_training_eta()")
  pt <- predict(mod, newdata = nd, type = "terms")
  testthat::expect_equal(
    unname(rowSums(pt) + attr(pt, "constant") + mod$offset), eta,
    tolerance = tol, label = "terms + constant + offset")

  # Learner edf~_k add up to the tracked edf~ at mstop
  if (mod$edf_compute_method != "none") {
    acc <- sum(vapply(mod$learners, function(lr) lr$edf, numeric(1)))
    testthat::expect_equal(acc, mod$edfs[mod$mstop + 1L], tolerance = tol,
                           label = "sum of learner edf~_k")
    testthat::expect_equal(mod$edf.total, acc, tolerance = tol)
  } else {
    testthat::expect_true(is.na(mod$edf.total))
  }

  # Fixed-mode replay on the observed y
  if (mod$mstop > 0L) {
    rep_fit <- replay_on_y(mod, mod$y, "fixed")
    testthat::expect_equal(as.vector(rep_fit$linear.predictors), eta,
                           tolerance = tol, label = "fixed-mode replay eta")
    testthat::expect_equal(rep_fit$coefficients, mod$coefficients,
                           tolerance = tol, label = "fixed-mode replay coef")
  }

  # fitted() and residuals() follow na.action
  fv <- fitted(mod)
  testthat::expect_length(fv, n_full)
  testthat::expect_equal(as.vector(fv[kept]), as.vector(mod$fitted.values),
                         tolerance = tol)
  for (type in c("deviance", "response", "working", "pearson",
                 "scaled.pearson")) {
    r <- residuals(mod, type = type)
    testthat::expect_length(r, n_full)
    testthat::expect_true(all(is.finite(r[kept])), label = type)
    if (is_excl) testthat::expect_true(all(is.na(r[mod$na.action])))
  }
  testthat::expect_equal(as.vector(residuals(mod, type = "response")[kept]),
                         as.vector(mod$y - mod$fitted.values), tolerance = tol)

  if (!is.null(truth_eta)) {
    testthat::expect_length(truth_eta, n)
    testthat::expect_gt(stats::cor(eta, truth_eta), min_cor)
  }

  # Hutch++ edf~ against the exact trace on the same path. Only for
  # canonical links: for the others the Jacobian has a Frobenius norm much
  # larger than its trace, and NA-Hutch++ can be far off
  if (!is.null(hutchpp_rtol)) {
    testthat::expect_identical(mod$edf_compute_method, "hutchpp")
    testthat::expect_true(is_canonical_link(mod$family))
    edf_ex <- exact_edf_on_path(mod)
    testthat::expect_lte(abs(mod$edf.total - edf_ex) / edf_ex, hutchpp_rtol)
  }
  invisible(mod)
}

# edf~ at mstop by the exact trace, on the recorded selection, rho and
# post-iteration rho paths of `mod` (as a fixed-mode bootstrap replicate on
# the observed y, but with exact EDF tracking)
exact_edf_on_path <- function(mod) {
  m   <- mod$mstop
  rec <- recorded_paths(mod, m)
  s   <- replicate_settings("fixed", rec$path, rec$rhos, m, mod,
                            post_rhos = rec$post_rhos)
  s$edf_compute_method <- "exact"
  ex <- refit_replicate(mod, mod$y, s)
  testthat::expect_equal(as.vector(ex$linear.predictors),
                         as.vector(mod$linear.predictors), tolerance = 1e-8,
                         label = "exact-EDF replay eta")
  ex$edf.total
}

# Largest absolute error relative to the largest reference entry
max_rel_err <- function(a, b) {
  max(abs(a - b)) / max(abs(b))
}

weighted_colmeans <- function(X, w) {
  colSums(X * w) / sum(w)
}

deviance_path <- function(y, fitted_mat, family = stats::gaussian(), w = rep(1, length(y))) {
  apply(fitted_mat, 2L, function(mu) {
    sum(family$dev.resids(y, mu, wt = w))
  })
}


# ---- Data simulators ---------------------------------------------------------

simulate_gaussian_1d <- function(n = 120L, seed = 1L, sd = 0.25) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  f <- sin(2 * pi * x)
  y <- f + stats::rnorm(n, sd = sd)
  data.frame(y = y, x = x, f = f)
}

simulate_gaussian_2d <- function(n = 150L, seed = 2L, sd = 0.2) {
  set.seed(seed)
  x <- stats::runif(n)
  z <- stats::runif(n)
  f <- sin(2 * pi * x) + cos(2 * pi * z) + 0.75 * sin(2 * pi * x * z)
  y <- f + stats::rnorm(n, sd = sd)
  data.frame(y = y, x = x, z = z, f = f)
}

simulate_poisson_1d <- function(n = 120L, seed = 3L) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  eta <- 0.4 + sin(2 * pi * x)
  mu <- exp(eta)
  y <- stats::rpois(n, lambda = mu)
  data.frame(y = y, x = x, eta = eta, mu = mu)
}

simulate_binomial_1d <- function(n = 120L, seed = 5L) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  eta <- sin(2 * pi * x)
  prob <- stats::plogis(eta)
  y <- stats::rbinom(n, size = 1L, prob = prob)
  data.frame(y = y, x = x, prob = prob)
}

simulate_gamma_1d <- function(n = 120L, seed = 6L) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  eta <- 0.5 + sin(2 * pi * x)
  mu  <- exp(eta)          # log link is default for Gamma
  y   <- stats::rgamma(n, shape = 2, rate = 2 / mu)
  data.frame(y = y, x = x, mu = mu)
}

simulate_invgaussian_1d <- function(n = 120L, seed = 7L) {
  testthat::skip_if_not_installed("SuppDists")
  set.seed(seed)
  x   <- sort(stats::runif(n))
  eta <- 0.5 + 0.5 * sin(2 * pi * x)
  mu  <- exp(eta)                # log link for numerical stability
  # SuppDists parameterisation: Var(y) = nu^3 / lambda. lambda = mu gives
  # Var(y) = mu^2, i.e. a non-constant dispersion phi = 1 / mu (not phi = 1)
  y   <- SuppDists::rinvGauss(n, nu = mu, lambda = mu)
  data.frame(y = y, x = x, mu = mu)
}

simulate_binomial_formats <- function(n = 80L, seed = 4L) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  f <- sin(2 * pi * x)
  prob <- stats::plogis(f)
  n_trials <- sample(5:20, n, replace = TRUE)
  successes <- stats::rbinom(n, size = n_trials, prob = prob)
  failures <- n_trials - successes

  agg <- data.frame(
    x = x,
    f = f,
    n_trials = n_trials,
    successes = successes,
    failures = failures,
    prop = successes / n_trials
  )

  id <- rep(seq_len(n), times = n_trials)
  long <- data.frame(
    id = id,
    x = x[id],
    y = unlist(Map(function(s, m) c(rep(1, s), rep(0, m - s)), successes, n_trials)),
    stringsAsFactors = FALSE
  )
  long$y_fac <- factor(long$y, levels = 0:1, labels = c("fail", "success"))

  list(agg = agg, long = long)
}

# Distributed-lag data: a complete lag matrix of a longer input series, the
# lag index matrix and Poisson counts with a constant mean
simulate_functional_dlm <- function(n = 200L, n_lag = 6L, seed = 40L) {
  set.seed(seed)
  # Longer series so every row has n_lag complete lags
  x_full <- stats::rnorm(n + n_lag - 1L)
  lag_fn <- function(x, n.lag) {
    nv <- length(x)
    X  <- matrix(NA_real_, nv, n.lag)
    for (i in seq_len(n.lag)) X[i:nv, i] <- x[seq_len(nv - i + 1L)]
    X
  }
  x_mat <- lag_fn(x_full, n_lag)
  # Keep only complete rows (no NAs)
  keep  <- complete.cases(x_mat)
  x_mat <- x_mat[keep, , drop = FALSE]
  # lag index matrix: each row is 0, 1, ..., n_lag-1
  lag_mat <- matrix(seq_len(n_lag) - 1L, nrow(x_mat), n_lag, byrow = TRUE)
  # A constant mean: the tests check the construction, not the fit
  eta <- rep(0.3, nrow(x_mat))
  y   <- stats::rpois(nrow(x_mat), lambda = exp(eta))
  list(y = y, x = x_mat, lag = lag_mat)
}

# Proportional binomial data: y = successes / n_trials, w = n_trials.
simulate_binomial_prop <- function(n = 30L, seed = 5L) {
  set.seed(seed)
  x        <- sort(stats::runif(n))
  n_trials <- sample(5L:15L, n, replace = TRUE)
  prob     <- stats::plogis(sin(2 * pi * x))
  succ     <- stats::rbinom(n, n_trials, prob)
  data.frame(y = succ / n_trials, x = x, w = n_trials)
}

# ---- Fit wrappers and learner helpers ----------------------------------------

# A short fixed-length fit, without EDF tracking
sc_ff <- gb_stop_ctrl("none", max_iter = 15L)

ec_ff <- gb_edf_ctrl("none")

# One-dimensional Gaussian data and a formula for quick fits
core_data <- function(n = 150, seed = 1) {
  set.seed(seed)
  d   <- data.frame(x = stats::runif(n))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(n, sd = 0.3)
  d
}

ff_core <- y ~ t2(x, k = 6, bs = "cr")

# Two collinear covariates (x1, x2) and one independent covariate (x0)
sim_collinear <- function(n = 250, seed = 21) {
  set.seed(seed)
  x0 <- stats::runif(n)
  x1 <- stats::runif(n)
  x2 <- 0.9 * x1 + 0.1 * stats::runif(n)
  data.frame(x0 = x0, x1 = x1, x2 = x2,
             y = sin(2 * pi * x0) + x1 + stats::rnorm(n, 0, 0.3))
}

# AIC-stopped fit with exact EDF
fit_aic <- function(formula, data, family = stats::gaussian(), ...) {
  gbgam(formula, data = data, family = family,
        learner_control = gb_learner_ctrl(step_size = 0.1),
        edf_control  = gb_edf_ctrl("exact"),
        stop_control = gb_stop_ctrl("aic", iter_no_improv = 200L),
        verbose = FALSE, ...)
}

f_col <- y ~ t2(x0, bs = "ps", k = 10) + t2(x1, bs = "ps", k = 10) +
  t2(x2, bs = "ps", k = 10)

# A short fixed-length fit (no stopping rule), without progress output or
# messages
fit_fixed <- function(formula, data, family = stats::gaussian(),
                      max_iter = 25, edf = "none", step_size = 0.05, ...) {
  suppressMessages(gbgam(
    formula,
    data            = data,
    family          = family,
    learner_control = gb_learner_ctrl(step_size = step_size),
    stop_control    = gb_stop_ctrl("none", max_iter = max_iter),
    edf_control     = gb_edf_ctrl(edf),
    verbose         = FALSE,
    ...
  ))
}

# One small model per family, fitted afresh in each test block
local_gaussian_mod <- function() {
  fit_fixed(y ~ t2(x, k = 8, bs = "cr"),
            simulate_gaussian_1d(n = 80, seed = 99),
            uq_control = gb_uq_ctrl("none"))
}

local_gaussian_mod_edf <- function() {
  fit_fixed(y ~ t2(x, k = 8, bs = "cr"),
            simulate_gaussian_1d(n = 80, seed = 99), edf = "exact",
            uq_control = gb_uq_ctrl("none"))
}

local_poisson_mod <- function() {
  fit_fixed(y ~ t2(x, k = 8, bs = "cr"),
            simulate_poisson_1d(n = 80, seed = 99),
            family = stats::poisson(),
            uq_control = gb_uq_ctrl("none"))
}

# One-dimensional fit with a given stopping rule
fit_1d <- function(dat, family = stats::gaussian(), stop_control,
                   edf = gb_edf_ctrl("exact"), ...) {
  gbgam(y ~ t2(x, k = 10, bs = "cr"), data = dat, family = family,
        learner_control = gb_learner_ctrl(step_size = 0.1),
        stop_control = stop_control, edf_control = edf, verbose = FALSE, ...)
}

# gbgam() without progress output or messages
quiet_fit <- function(...) {
  suppressMessages(gbgam(..., verbose = FALSE))
}

# The penalised (Demmler-Reinsch) learners of a learner bundle
dr_learners <- function(lb) {
  Filter(\(l) inherits(l, "demmler_reinsch"), lb$learners)
}

# Selection EDF of a penalised learner (1 by construction)
sedf <- function(lr) {
  s <- if (isTRUE(lr$lambda_at_zero)) rep(1, length(lr$D)) else
    1 / (1 + exp(lr$rho) * lr$D)
  2 * lr$nu_k * sum(s) - lr$nu_k^2 * sum(s^2)
}

# Replay one bootstrap replicate on a given response, exactly as one_rep()
# inside bootstrap_gbgam() does; post_rhos overrides the recorded post-
# iteration rhos
replay_on_y <- function(mod, y, mode, post_rhos = rec$post_rhos) {
  m   <- mod$mstop
  rec <- recorded_paths(mod, m)
  s   <- replicate_settings(mode, rec$path,
                            rec$rhos, m, mod, post_rhos = post_rhos)
  refit_replicate(mod, y, s)
}
