# Centred smooths (class "gbgam_centred") in mgcv's prediction machinery:
# term-wise effects and SEs are centred, full predictions are unchanged.

ct_sim_data <- function(n = 160L, seed = 61L) {
  set.seed(seed)
  d <- data.frame(
    x   = stats::runif(n),
    z   = stats::runif(n),
    lin = stats::rnorm(n),
    zz  = stats::rnorm(n),
    g   = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    e   = stats::runif(n, 0.5, 2),
    w   = stats::runif(n, 0.5, 1.5)
  )
  eta <- sin(2 * pi * d$x) + 0.5 * d$lin + cos(2 * pi * d$z) +
    0.3 * d$zz * d$x + 0.2 * as.numeric(d$g)
  d$y  <- eta + stats::rnorm(n, sd = 0.3)
  d$yp <- stats::rpois(n, exp(0.3 * eta + log(d$e)))
  d
}

ct_fit <- function(formula, data, family = stats::gaussian(), max_iter = 60L,
                   edf = "none", ...) {
  gbgam(
    formula,
    data         = data,
    family       = family,
    stop_control = gb_stop_ctrl("none", max_iter = max_iter),
    edf_control  = gb_edf_ctrl(edf),
    verbose      = FALSE,
    ...
  )
}

# Coefficients of the previous mapping, which absorbed the centring shift of
# every smooth into the intercept and predicted from the uncentred basis
ct_old_coef <- function(mod, coef = mod$coefficients) {
  centred <- !mod$cen_absorb
  coef[1L] <- coef[1L] - sum((mod$cen_const * coef)[centred])
  coef
}

# Uncentred lpmatrix: add the centring constants back to centred smooths
ct_uncentre_lpmatrix <- function(mod, X) {
  for (sm in mod$smooth) {
    if (inherits(sm, "gbgam_centred")) {
      idx      <- sm$first.para:sm$last.para
      X[, idx] <- sweep(X[, idx, drop = FALSE], 2, sm$gbgam_cen, "+")
    }
  }
  X
}

ct_centred_smooths <- function(mod) {
  Filter(function(sm) inherits(sm, "gbgam_centred"), mod$smooth)
}

ct_check_predictions <- function(mod, data) {
  lp <- as.vector(mod$linear.predictors)
  expect_equal(as.vector(predict(mod)), lp, tolerance = 1e-10)
  expect_equal(as.vector(predict(mod, newdata = data)), lp, tolerance = 1e-10)
  # Old coefficient mapping on the uncentred basis gives the same predictor
  X     <- predict(mod, newdata = data, type = "lpmatrix")
  X_old <- ct_uncentre_lpmatrix(mod, X)
  off   <- attr(X, "model.offset") %||% 0
  expect_equal(as.vector(X_old %*% ct_old_coef(mod)) + off, lp,
               tolerance = 1e-10)
  # Term-wise decomposition reproduces the link predictor
  pt <- predict(mod, type = "terms")
  expect_equal(unname(rowSums(pt) + attr(pt, "constant") + mod$offset), lp,
               tolerance = 1e-10)
  invisible(pt)
}


# ---- Construction ------------------------------------------------------------

test_that("smooths without by-variable are centred; by-smooths are not", {
  d <- ct_sim_data()
  m <- ct_fit(y ~ g + t2(x) + t2(x, by = g) + t2(z, by = zz), d,
              max_iter = 5L)
  is_cen <- vapply(m$smooth, inherits, logical(1L), "gbgam_centred")
  n_cons <- vapply(m$smooth, attr, integer(1L), "nCons")
  expect_equal(is_cen, c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(n_cons, as.integer(is_cen))
  sm1 <- m$smooth[[1L]]
  expect_equal(sm1$gbgam_cen,
               unname(m$cen_const[sm1$first.para:sm1$last.para]))
  # Only the parametric and by-smooth centring shifts go into the intercept
  expect_length(m$cen_absorb, length(m$coefficients))
  expect_false(any(m$cen_absorb[sm1$first.para:sm1$last.para]))
  expect_true(all(m$cen_absorb[-(sm1$first.para:sm1$last.para)]))
})

test_that("centred smooths reference the gbgam namespace", {
  # The namespace is serialised by reference, so readRDS() of a saved model
  # loads gbgam and registers Predict.matrix.gbgam_centred()
  d  <- ct_sim_data()
  m  <- ct_fit(y ~ t2(x) + t2(x, by = g), d, max_iter = 5L)
  ns <- m$smooth[[1L]]$gbgam_ns
  expect_true(isNamespace(ns))
  expect_identical(getNamespaceName(ns), c(name = "gbgam"))
  expect_null(m$smooth[[2L]]$gbgam_ns)

  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(m, f)
  m2 <- readRDS(f)
  expect_identical(m2$smooth[[1L]]$gbgam_ns, ns)
  expect_equal(as.vector(predict(m2, newdata = d)),
               as.vector(m$linear.predictors), tolerance = 1e-10)
})

test_that("functional (matrix-argument) smooths keep the intercept shift", {
  dat <- simulate_functional_dlm(n = 120, n_lag = 4, seed = 62)
  m   <- gbgam(y ~ t2(x, lag, k = c(5, 4), bs = c("cr", "cr")),
               family       = stats::poisson(),
               data         = dat,
               stop_control = gb_stop_ctrl("none", max_iter = 10),
               edf_control  = gb_edf_ctrl("none"),
               verbose      = FALSE)
  expect_false(inherits(m$smooth[[1L]], "gbgam_centred"))
  expect_identical(attr(m$smooth[[1L]], "nCons"), 0L)
  expect_true(all(m$cen_absorb))
  expect_equal(as.vector(predict(m, newdata = dat)),
               as.vector(m$linear.predictors), tolerance = 1e-10)
})


# ---- Predictions unchanged ---------------------------------------------------

test_that("weighted gaussian, 1-D and 2-D t2: same predictions, centred terms", {
  d    <- ct_sim_data()
  m    <- ct_fit(y ~ lin + t2(x) + t2(x, z), d, weights = d$w)
  pt   <- ct_check_predictions(m, d)
  cols <- vapply(ct_centred_smooths(m), `[[`, "", "label")
  expect_equal(cols, c("t2(x)", "t2(x,z)"))
  expect_equal(unname(colSums(pt[, cols] * d$w) / sum(d$w)), c(0, 0),
               tolerance = 1e-10)
  # The centred smooth columns of the lpmatrix have weighted mean zero
  X <- predict(m, type = "lpmatrix")
  expect_equal(unname(colSums(X * d$w) / sum(d$w)), m$cmX, tolerance = 1e-10)
  expect_true(all(m$cmX[-seq_len(m$nsdf)] == 0))
})

test_that("poisson with offset, 1-D and 2-D t2: predictions unchanged", {
  d  <- ct_sim_data()
  m  <- ct_fit(yp ~ t2(x) + t2(x, z, k = c(4, 4)) + offset(log(e)), d,
               family = stats::poisson())
  pt <- ct_check_predictions(m, d)
  expect_equal(unname(colMeans(pt)), c(0, 0), tolerance = 1e-10)
})

test_that("factor-by and numeric-by smooths: predictions unchanged", {
  d <- ct_sim_data()
  m_f <- ct_fit(y ~ g + t2(x, by = g) + t2(z), d, weights = d$w)
  pt  <- ct_check_predictions(m_f, d)
  expect_equal(sum(pt[, "t2(z)"] * d$w) / sum(d$w), 0, tolerance = 1e-10)
  # Factor-by columns are zero outside their level (uncentred, as before)
  expect_true(all(pt[d$g != "a", "t2(x):ga"] == 0))

  m_n <- ct_fit(y ~ t2(x) + t2(x, by = zz), d)
  pt  <- ct_check_predictions(m_n, d)
  expect_equal(mean(pt[, "t2(x)"]), 0, tolerance = 1e-10)
})

test_that("rank-deficient 're' smooth: predictions unchanged, terms centred", {
  d  <- ct_sim_data()
  m  <- ct_fit(yp ~ t2(x) + t2(g, bs = "re"), d, family = stats::poisson())
  pt <- ct_check_predictions(m, d)
  expect_equal(unname(colMeans(pt)), c(0, 0), tolerance = 1e-10)
  i  <- c(which(d$g == "c")[1L], which(d$g == "a")[1L])
  nd <- data.frame(x = d$x[i], g = factor(c("c", "a"), levels = levels(d$g)))
  expect_equal(as.vector(predict(m, newdata = nd)),
               as.vector(m$linear.predictors[i]), tolerance = 1e-10)
})

test_that("derivative bases are not shifted by the centring constants", {
  d  <- ct_sim_data()
  m  <- ct_fit(y ~ t2(x, bs = "cr"), d, max_iter = 5L)
  sm <- m$smooth[[1L]]
  expect_false(smooth_has_deriv(sm))
  sm$margin[[1L]]$deriv <- 1
  expect_true(smooth_has_deriv(sm))
  raw <- sm
  class(raw) <- setdiff(class(raw), "gbgam_centred")
  nd <- data.frame(x = c(0.1, 0.5, 0.9))
  expect_equal(mgcv::PredictMat(sm, nd), mgcv::PredictMat(raw, nd))
  sm$margin[[1L]]$deriv <- NULL
  raw$margin[[1L]]$deriv <- NULL
  expect_equal(mgcv::PredictMat(sm, nd),
               sweep(mgcv::PredictMat(raw, nd), 2, sm$gbgam_cen),
               ignore_attr = TRUE)
})


# ---- Standard errors after bootstrap -----------------------------------------

test_that("bootstrap: same full-predictor SEs, term SEs drop the constant", {
  d <- ct_sim_data()
  m <- ct_fit(yp ~ t2(x) + t2(x, z, k = c(4, 4)) + t2(z, by = g) +
                offset(log(e)), d, family = stats::poisson(), max_iter = 40L)
  set.seed(63)
  b <- bootstrap_gbgam(m, n_samples = 12)

  # Previous coefficient mapping of the same bootstrap draws
  c2     <- b$widening_factor^2
  old_bc <- t(apply(b$boot_coef, 1L, ct_old_coef, mod = b))
  Vp_old <- c2 * stats::cov(old_bc)
  X      <- predict(b, newdata = d, type = "lpmatrix")
  X_old  <- ct_uncentre_lpmatrix(b, X)
  se_old <- unname(sqrt(rowSums((X_old %*% Vp_old) * X_old)))
  se_new <- predict(b, newdata = d, se.fit = TRUE)$se.fit
  expect_equal(unname(as.vector(se_new)), se_old, tolerance = 1e-8)

  # Term-wise: on average over the training data, the old variance exceeds
  # the new one by exactly the variance of the constant c^T beta
  pt <- predict(b, type = "terms", se.fit = TRUE)
  for (sm in ct_centred_smooths(b)) {
    idx     <- sm$first.para:sm$last.para
    Xs_old  <- X_old[, idx, drop = FALSE]
    v_old   <- rowSums((Xs_old %*% b$Vp[idx, idx]) * Xs_old)
    v_const <- drop(t(sm$gbgam_cen) %*% b$Vp[idx, idx] %*% sm$gbgam_cen)
    expect_equal(mean(pt$se.fit[, sm$label]^2), mean(v_old) - v_const,
                 tolerance = 1e-8)
    expect_gt(v_const, 0)
  }
  # By-smooth SEs are unchanged
  for (lab in c("t2(z):ga", "t2(z):gb")) {
    sm  <- b$smooth[[which(vapply(b$smooth, `[[`, "", "label") == lab)]]
    idx <- sm$first.para:sm$last.para
    Xs  <- X[, idx, drop = FALSE]
    expect_equal(unname(pt$se.fit[, lab]),
                 unname(sqrt(rowSums((Xs %*% b$Vp[idx, idx]) * Xs))),
                 tolerance = 1e-8)
  }

  # plot() with seWithMean and on one page
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(b, seWithMean = TRUE))
  expect_no_error(plot(b, pages = 1, se = TRUE))
  expect_no_error(predict(b, type = "iterms", se.fit = TRUE))
})


# ---- Bayesian backfitting ----------------------------------------------------

test_that("bayesian_backfitting maps posterior samples onto the centred basis", {
  d <- ct_sim_data()
  m <- ct_fit(y ~ lin + t2(x) + t2(x, by = g), d, weights = d$w,
              edf = "exact", max_iter = 40L)
  set.seed(64)
  bb <- bayesian_backfitting(m, n_samples = 20L, n_burnin = 5L)
  # Posterior mean predictor from the learners' designs and mean betas
  eta_bar <- Reduce(`+`, Map(function(lr, bs) {
    as.vector(lr$design %*% rowMeans(bs))
  }, bb$learners, bb$beta_samples))
  X <- predict(bb, type = "lpmatrix")
  expect_equal(as.vector(X %*% bb$coefficients_bayes), eta_bar,
               tolerance = 1e-10)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(bb, seWithMean = TRUE, pages = 1))
})
