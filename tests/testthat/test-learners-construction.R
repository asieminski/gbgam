# Construction of the learners by new_learners(): centring, the selection
# EDF calibration, the Demmler-Reinsch form, aliased blocks and sp values.


# ---- Learner bundle ----------------------------------------------------------

test_that("new_learners returns a coherent learner bundle", {
  set.seed(20)
  n <- 90
  dat <- data.frame(
    y = stats::rnorm(n),
    x = stats::runif(n),
    z = stats::runif(n),
    lin = stats::rnorm(n)
  )
  w <- sample(1:4, n, replace = TRUE)

  obj <- new_learners(
    y ~ lin + t2(x, k = 6, bs = "cr") + t2(x, z, k = c(4, 4), bs = c("cr", "cr")),
    data = dat,
    weights = w
  )

  expect_true(is.list(obj))
  expect_true(all(c("learners", "coef_names", "ind0var", "cen_const") %in% names(obj)))
  expect_type(obj$coef_names, "character")
  expect_true(is.list(obj$learners))
  expect_true(length(obj$learners) >= 3L)
  expect_equal(length(obj$coef_names), length(obj$cen_const))
})

test_that("non-intercept learner designs are weighted-centred", {
  set.seed(21)
  n <- 100
  dat <- data.frame(
    y = stats::rnorm(n),
    x = stats::runif(n),
    z = stats::runif(n),
    lin = stats::rnorm(n)
  )
  w <- sample(1:5, n, replace = TRUE)

  obj <- new_learners(
    y ~ lin + t2(x, k = 6, bs = "cr") + t2(x, z, k = c(4, 4), bs = c("cr", "cr")),
    data = dat,
    weights = w
  )

  for (j in seq_along(obj$learners)) {
    lr <- obj$learners[[j]]
    X <- lr$design

    if (lr$var_label == "(Intercept)") {
      expect_true(all(as.vector(X) == 1))
      next
    }

    mu <- weighted_colmeans(as.matrix(X), w)
    expect_true(max(abs(mu)) < 1e-10)
  }
})

test_that("demmler_reinsch learners satisfy unbiased-selection EDF calibration", {
  set.seed(22)
  n <- 120
  dat <- data.frame(
    y = stats::rnorm(n),
    x = stats::runif(n),
    z = stats::runif(n)
  )

  obj <- new_learners(
    y ~ t2(x, k = 8, bs = "cr") + t2(x, z, k = c(5, 5), bs = c("cr", "cr")),
    data = dat,
    weights = rep(1, n)
  )

  expect_true(length(dr_learners(obj)) >= 1L)

  for (lr in dr_learners(obj)) {
    expect_equal(sedf(lr), 1, tolerance = 1e-4)
    expect_true(lr$nu_k > 0)
    expect_true(lr$nu_k <= 1 + 1e-12)
  }
})

test_that("learner-to-mgcv coefficient reconstruction reproduces training linear predictor", {
  skip_on_cran()
  dat <- simulate_gaussian_2d(n = 120, seed = 23)

  mod <- fit_fixed(
    y ~ t2(x, k = 8, bs = "cr") + t2(x, z, k = c(5, 5), bs = c("cr", "cr")),
    dat
  )

  expect_model_contract(mod, nrow(dat))

  eta_rebuilt <- rebuild_training_eta(mod)
  expect_equal(eta_rebuilt, as.vector(mod$linear.predictors), tolerance = 1e-8)
})

test_that("ind0var bookkeeping satisfies structural invariants", {
  set.seed(24)
  n <- 80
  dat <- data.frame(
    y = stats::rnorm(n),
    x = stats::runif(n),
    z = stats::runif(n)
  )

  obj <- new_learners(
    y ~ t2(x, k = 5, bs = "cr") + t2(z, k = 5, bs = "cr"),
    data = dat,
    weights = rep(1, n)
  )

  # The constant column of each t2(cr) basis has zero variance after
  # centring: one dropped column per smooth
  expect_true(is.integer(obj$ind0var))
  expect_length(obj$ind0var, 2L)
  expect_true(all(obj$ind0var >= 1L))
  expect_true(all(obj$ind0var <= length(obj$coef_names)))
  expect_false(is.unsorted(obj$ind0var, strictly = TRUE))
})


# ---- Demmler-Reinsch form ----------------------------------------------------

test_that("demmler_reinsch $D is always a numeric vector (not a matrix)", {
  set.seed(31)
  n <- 80
  dat <- data.frame(y = rnorm(n), x = runif(n), z = runif(n))
  w   <- sample(1:3, n, replace = TRUE)

  obj <- new_learners(
    y ~ t2(x, k = 7, bs = "cr") + t2(x, z, k = c(4, 4), bs = c("cr", "cr")),
    data = dat, weights = w
  )

  expect_true(length(dr_learners(obj)) >= 1L)

  for (lr in dr_learners(obj)) {
    expect_true(is.numeric(lr$D), info = "D must be a numeric vector")
    expect_true(is.vector(lr$D),  info = "D must not be a matrix")
    expect_true(all(lr$D >= 0),   info = "D entries (penalty eigenvalues) must be non-negative")
  }
})

test_that("demmler_reinsch design is column-orthonormal w.r.t. prior weights", {
  set.seed(32)
  n <- 100
  dat <- data.frame(y = rnorm(n), x = runif(n), z = runif(n))
  w   <- runif(n, 0.5, 2)

  obj <- new_learners(
    y ~ t2(x, k = 8, bs = "cr") + t2(x, z, k = c(4, 4), bs = c("cr", "cr")),
    data = dat, weights = w
  )

  for (lr in dr_learners(obj)) {
    X    <- lr$design
    XtWX <- crossprod(sqrt(w) * X)
    expect_equal(XtWX, diag(ncol(X)), tolerance = 1e-8,
                 info = paste("X^T W X != I for learner", lr$var_label))
  }
})


# ---- Rank-deficient blocks ---------------------------------------------------

test_that("DR learners satisfy V^T W V = I and D > 0 for aliased blocks", {
  set.seed(21)
  n   <- 200
  dat <- data.frame(x = stats::runif(n), y = stats::rnorm(n),
                    g = factor(sample(letters[1:6], n, replace = TRUE)))
  w_zero <- ifelse(dat$x > 0.6, 0, 1)
  cases <- list(
    list(f = y ~ t2(x) + t2(g, bs = "re"),        w = rep(1, n)),
    list(f = y ~ t2(x, g, bs = c("cr", "re")),    w = stats::runif(n, 0.5, 2)),
    list(f = y ~ t2(x, k = 10, bs = "cr"),        w = w_zero)
  )
  for (cs in cases) {
    lb <- new_learners(cs$f, data = dat, weights = cs$w)
    for (lr in dr_learners(lb)) {
      V <- lr$design
      expect_lt(max(abs(crossprod(V, cs$w * V) - diag(ncol(V)))), 1e-8)
      expect_true(all(lr$D > 0), info = lr$var_label)
      expect_equal(sedf(lr), 1, tolerance = 1e-4)
    }
  }
})

dr_reference_svd <- function(Z, w, tol = 1e-7) {
  sv  <- svd(sqrt(w) * Z)
  r   <- sum(sv$d > tol * sv$d[1L])
  ord <- rev(seq_len(r))
  list(A = sv$v[, ord, drop = FALSE] %*% diag(1 / sv$d[ord], r),
       D = 1 / sv$d[ord]^2)
}
dr_hat <- function(Z, w, dr, lambda) {
  V <- Z %*% dr$A
  V %*% (t(V * w) / (1 + lambda * dr$D))
}

test_that("QR and SVD Demmler-Reinsch bases agree on full-rank blocks", {
  set.seed(25)
  for (p in c(3L, 12L, 40L)) {
    n  <- 300
    Z  <- matrix(stats::rnorm(n * p), n, p)
    w  <- stats::runif(n, 0.5, 2)
    Z  <- sweep(Z, 2, colSums(w * Z) / sum(w))
    dr <- demmler_reinsch_basis(Z, w)
    ref <- dr_reference_svd(Z, w)
    expect_equal(dr$D, ref$D, tolerance = 1e-10)
    expect_false(is.unsorted(rev(dr$D)))
    V <- Z %*% dr$A
    expect_lt(max(abs(crossprod(V, w * V) - diag(p))), 1e-10)
    expect_lt(max(abs(colSums(w * V))), 1e-8)
    for (lambda in c(0.01, 1, 100)) {
      expect_lt(max(abs(dr_hat(Z, w, dr, lambda) - dr_hat(Z, w, ref, lambda))),
                1e-10)
    }
  }
})

test_that("rank-deficient blocks use the truncated SVD basis", {
  set.seed(26)
  n <- 200
  w <- stats::runif(n, 0.5, 2)
  # Exactly aliased: the last column is a combination of the first two
  Z <- matrix(stats::rnorm(n * 4), n, 4)
  Z <- cbind(Z, Z[, 1] - 2 * Z[, 2])
  dr  <- demmler_reinsch_basis(Z, w)
  ref <- dr_reference_svd(Z, w)
  expect_length(dr$D, 4L)
  expect_equal(dr$D, ref$D, tolerance = 1e-10)
  expect_lt(max(abs(dr_hat(Z, w, dr, 0.5) - dr_hat(Z, w, ref, 0.5))), 1e-10)
  # Full QR rank but sigma_min / sigma_max below tol: orthogonal columns,
  # one of them tiny, so only the singular value threshold drops it
  Q <- qr.Q(qr(matrix(stats::rnorm(n * 4), n, 4)))
  Z <- sweep(Q, 2, c(1, 0.5, 0.2, 1e-9), `*`) / sqrt(w)
  expect_identical(qr(sqrt(w) * Z, tol = 1e-7)$rank, 4L)
  dr  <- demmler_reinsch_basis(Z, w)
  ref <- dr_reference_svd(Z, w)
  expect_length(dr$D, 3L)
  expect_equal(dr$D, ref$D, tolerance = 1e-10)
  expect_lt(max(abs(dr_hat(Z, w, dr, 0.5) - dr_hat(Z, w, ref, 0.5))), 1e-10)
})

test_that("'re' learner smoother equals the exact ridge smoother", {
  set.seed(22)
  n   <- 150
  dat <- data.frame(y = stats::rnorm(n),
                    g = factor(sample(letters[1:5], n, replace = TRUE)))
  lb  <- new_learners(y ~ t2(g, bs = "re"), data = dat, weights = rep(1, n))
  lr  <- dr_learners(lb)[[1]]
  # Centred dummy columns are aliased with the constant: rank 4 of 5
  expect_equal(ncol(lr$design), nlevels(dat$g) - 1L)

  sm     <- lb$smooths[[1]]
  Zc     <- sweep(sm$X, 2, colMeans(sm$X))
  lambda <- exp(lr$rho)
  exact  <- Zc %*% solve(crossprod(Zc) + lambda * diag(ncol(Zc)), t(Zc))
  V      <- lr$design
  dr     <- V %*% (t(V) / (1 + lambda * lr$D))
  expect_lt(max(abs(dr - exact)), 1e-10)
})

test_that("learners_to_mgcv_coef reproduces the learner fit for aliased blocks", {
  set.seed(23)
  n   <- 150
  dat <- data.frame(y = stats::rnorm(n), x = stats::runif(n),
                    g = factor(sample(letters[1:5], n, replace = TRUE)))
  lb  <- new_learners(y ~ t2(x) + t2(g, bs = "re"), data = dat,
                      weights = rep(1, n))
  lrs <- lapply(lb$learners, function(lr) {
    lr$beta <- stats::rnorm(length(lr$beta))
    lr
  })
  eta_learners <- Reduce(`+`, lapply(lrs, \(lr) as.vector(lr$design %*% lr$beta)))
  coefs <- learners_to_mgcv_coef(lrs, lb$coef_names, lb$ind0var, lb$cen_const)
  X     <- cbind(1, do.call(cbind, lapply(lb$smooths, \(sm) sm$X)))
  expect_equal(as.vector(X %*% coefs), eta_learners, tolerance = 1e-10)
})

test_that("a rank-1 penalised block gives an informative error", {
  set.seed(24)
  n   <- 60
  dat <- data.frame(y = stats::rnorm(n),
                    g = factor(sample(c("u", "v"), n, replace = TRUE)))
  expect_error(
    new_learners(y ~ t2(g, bs = "re"), data = dat, weights = rep(1, n)),
    "has rank 1"
  )
})

test_that("rank-deficient 're' learners spread EDF over their original columns", {
  set.seed(12)
  n   <- 120L
  dat <- data.frame(x = stats::runif(n),
                    g = factor(rep(letters[1:5], each = 24L)))
  dat$y <- sin(2 * pi * dat$x) + (as.numeric(dat$g) - 3) / 4 +
    stats::rnorm(n, sd = 0.3)
  expect_warning(
    mod <- gbgam(y ~ t2(x, k = 6, bs = "cr") + t2(g, bs = "re"), data = dat,
                 learner_control = gb_learner_ctrl(step_size = 0.1),
                 stop_control = gb_stop_ctrl("none", max_iter = 50L),
                 edf_control = gb_edf_ctrl("exact"), verbose = FALSE),
    NA
  )
  # The centred 're' block is aliased: fewer DR coefficients than columns
  lr_re <- Filter(function(lr) !is.null(lr$A) && nrow(lr$A) > length(lr$beta),
                  mod$learners)
  expect_length(lr_re, 1L)
  expect_length(mod$edf, length(mod$coefficients))
  edf_sum <- sum(vapply(mod$learners, `[[`, numeric(1), "edf"))
  expect_equal(mod$edf.total, edf_sum, tolerance = 1e-12)
  expect_equal(sum(mod$edf), edf_sum, tolerance = 1e-12)
  re_sm <- mod$smooth[[length(mod$smooth)]]
  expect_equal(sum(mod$edf[re_sm$first.para:re_sm$last.para]),
               lr_re[[1]]$edf, tolerance = 1e-12)
})


# ---- Smoothing parameters, labels and weights --------------------------------

test_that("initial nu_k is stored as nu_k0", {
  dat <- simulate_gaussian_2d(n = 100)
  lb  <- new_learners(y ~ t2(x, z, sp = -c(4, 3, 2)), data = dat,
                      weights = rep(1, 100))
  for (lr in dr_learners(lb)) {
    expect_identical(lr$nu_k0, lr$nu_k)
    expect_lt(lr$nu_k, 1)
  }
})

test_that("user sp values are matched to the right t2 penalty", {
  set.seed(31)
  n   <- 200
  dat <- data.frame(x = stats::runif(n), z = stats::runif(n), y = stats::rnorm(n))
  # rr has one column (fitted as unpenalised), so nr/rn must still get
  # the 2nd and 3rd sp values
  lb  <- new_learners(y ~ t2(x, z, k = c(3, 3), bs = "cr", sp = -c(1.5, 1.2, 1.9)),
                      data = dat, weights = rep(1, n))
  # sp sets the EDF 2 tr(A) - tr(A^2) at nu_k = 1
  sedf_unscaled <- function(lr) learner_edf(1 / (1 + exp(lr$rho) * lr$D))
  edfs <- vapply(dr_learners(lb), sedf_unscaled, 1)
  labs <- vapply(dr_learners(lb), \(l) l$var_label, "")
  expect_equal(unname(edfs[labs == "t2(x,z)nr"]), 1.2, tolerance = 1e-4)
  expect_equal(unname(edfs[labs == "t2(x,z)rn"]), 1.9, tolerance = 1e-4)

  lb4 <- new_learners(y ~ t2(x, z, k = c(4, 4), bs = "cr", sp = -c(3, 2.4, 3.8)),
                      data = dat, weights = rep(1, n))
  edfs4 <- vapply(dr_learners(lb4), sedf_unscaled, 1)
  expect_equal(unname(edfs4), c(3, 2.4, 3.8), tolerance = 1e-4)
})

test_that("sp <= -p gives an informative error", {
  dat <- simulate_gaussian_1d(n = 60)
  expect_error(
    new_learners(y ~ t2(x, k = 5, sp = -5), data = dat, weights = rep(1, 60)),
    "requests sedf"
  )
})

test_that("nullspace learners of a t2 term get distinct labels", {
  dat <- simulate_gaussian_2d(n = 80)
  lb  <- new_learners(y ~ t2(x) + t2(x, z), data = dat, weights = rep(1, 80))
  labels <- vapply(lb$learners, \(l) l$var_label, "")
  expect_false(anyDuplicated(labels) > 0L)
  expect_true(all(c("t2(x)n1", "t2(x,z)n1", "t2(x,z)n2", "t2(x,z)n3") %in% labels))
})

test_that("new_learners() drops weights of rows removed for NA", {
  dat <- simulate_gaussian_1d(n = 60)
  set.seed(41)
  w   <- stats::runif(60, 0.5, 2)
  dat$x[3] <- NA
  lb_na  <- new_learners(y ~ t2(x), data = dat, weights = w)
  lb_ref <- new_learners(y ~ t2(x), data = dat[-3, ], weights = w[-3])
  expect_equal(lb_na$cen_const, lb_ref$cen_const)
  expect_equal(lb_na$learners[[2]]$design, lb_ref$learners[[2]]$design,
               ignore_attr = TRUE)
})

test_that("new_learners() rejects a negative weight", {
  dat <- simulate_gaussian_1d(n = 40, seed = 12)
  w_bad <- rep(1, nrow(dat))
  w_bad[20] <- -1
  expect_error(
    new_learners(y ~ t2(x, k = 6, bs = "cr"), data = dat, weights = w_bad),
    "Assertion on 'weights' failed: Element 20 is not >= 0"
  )
})
