# bayesian_backfitting() against a reference Bayesian backfitting written
# from Hastie and Tibshirani (2000) in the paper's notation (HT below; page
# and equation numbers refer to the article). Test only: the reference uses
# dense n x n smoother matrices.
#
# gbgam API used (adapt here if it changes):
#   gbgam(), bayesian_backfitting(n_samples, n_burnin, linear),
#   smooth_intervals(mod, alpha, centre = "posterior", newdata):
#   $fit and $upper, predict(type = "lpmatrix"),
#   mod$smooth[[j]]$first.para / $last.para,
#   mod$prior.weights, mod$family,
#   mod$bb_info$phi, mod$bb_info$lambda (term, learner, status, rho),
#   mod$learners[[k]]: $label, $var_label, $design, $D, class
#   "demmler_reinsch" for penalised learners.


# ---- Reference, part A: the paper's algorithm (knows nothing of gbgam) -----

# Smoother S = B (B^T B + lambda Omega)^{-1} B^T, HT (42), p. 211; a zero
# eigenvalue of lambda Omega is a flat prior in that direction, as for the
# linear part of a cubic smoothing spline (p. 198, below (7)).
ht_smoother <- function(B_j, lambda_Omega_j) {
  B_j %*% solve(crossprod(B_j) + lambda_Omega_j, t(B_j))
}

# Symmetric square root S^{1/2} of (8), p. 198, and (18), p. 200. HT use the
# series (39)-(40), Algorithm A.1, p. 210; the eigen-decomposition gives the
# same matrix.
ht_sqrt <- function(S) {
  if (max(abs(S - t(S))) > 1e-8 * max(abs(S))) stop("S is not symmetric")
  e <- eigen(S, symmetric = TRUE)
  if (min(e$values) < -1e-10 || max(e$values) > 1 + 1e-10)
    stop("eigenvalues of S outside [0, 1]")   # a smoother, p. 199
  # eigenvalues 0 up to rounding (S has rank ncol(B_j) < n) set to 0, so
  # their square roots (~1e-8) put no noise outside the span of B_j
  values <- ifelse(abs(e$values) < 1e-10, 0, pmin(e$values, 1))
  e$vectors %*% (sqrt(values) * t(e$vectors))
}

# Standard backfitting (16), p. 199, for the starting values: "the fitted
# curves from a standard additive model fit" (p. 199).
ht_backfit <- function(y, S, tol = 1e-10, max_sweeps = 1e4L) {
  n <- length(y)
  f <- replicate(length(S), numeric(n), simplify = FALSE)
  for (sweep in seq_len(max_sweeps)) {
    change <- 0
    for (j in seq_along(S)) {
      f_new  <- as.vector(S[[j]] %*% (y - mean(y) - Reduce(`+`, f[-j],
                                                           numeric(n))))
      change <- max(change, abs(f_new - f[[j]]))
      f[[j]] <- f_new
    }
    if (change < tol) return(list(alpha = mean(y), f = f))
  }
  stop("standard backfitting did not converge")
}

# Algorithm 3.1 (Bayesian backfitting), p. 199, with the conditional (17),
# p. 200, as the update (18): f_j <- S_j r_j + sigma S_j^{1/2} z_j.
# f_0 = 1 alpha, S_0 = 1 1^T / n: alpha ~ N(ave[y - sum_j f_j], sigma^2 / n)
# (p. 200, item 3). Fits centred after smoothing and generation (p. 200),
# for the constraint (15), p. 199.
ht_bayesian_backfitting <- function(y, S, sigma2, n_iter, n_burn) {
  n      <- length(y)
  p      <- length(S)
  sigma  <- sqrt(sigma2)
  S_half <- lapply(S, ht_sqrt)
  start  <- ht_backfit(y, S)                 # initial f_j^0 (Algorithm 3.1)
  alpha  <- start$alpha
  f      <- start$f
  alpha_draws <- numeric(n_iter)
  f_draws     <- replicate(p, matrix(NA_real_, n, n_iter), simplify = FALSE)
  for (t in seq_len(n_burn + n_iter)) {      # "Do for t = 1, 2, 3, ..."
    # j = 0: the intercept
    r_0   <- y - Reduce(`+`, f, numeric(n))
    alpha <- mean(r_0) + sigma * stats::rnorm(1) / sqrt(n)
    # j = 1, ..., p: r_j^t = y - sum_{k<j} f_k^t - sum_{k>j} f_k^{t-1}
    for (j in seq_len(p)) {
      r_j    <- y - alpha - Reduce(`+`, f[-j], numeric(n))
      z_j    <- stats::rnorm(n)
      f_j    <- as.vector(S[[j]] %*% r_j + sigma * S_half[[j]] %*% z_j)
      f[[j]] <- f_j - mean(f_j)                         # centring, p. 200
    }
    if (t > n_burn) {
      alpha_draws[t - n_burn] <- alpha
      for (j in seq_len(p)) f_draws[[j]][, t - n_burn] <- f[[j]]
    }
  }
  list(alpha = alpha_draws, f = f_draws)
}

# sigma^2 from the standard backfitting fit: the REML estimate at fixed
# lambda (the dispersion gbgam documents, see ?bayesian_backfitting), the
# residual sum of squares plus sum_j theta_j^T lambda Omega_j theta_j over
# n - M_p, with theta_j the coefficients of f_j on the basis B_j of S_j
# (attributes "B" and "pen" from gbgam_as_ht()) and M_p the intercept and
# the unpenalised columns
ht_sigma2 <- function(y, S) {
  fit <- ht_backfit(y, S)
  r   <- y - fit$alpha - Reduce(`+`, fit$f)
  pen <- sum(vapply(seq_along(S), function(j) {
    B     <- attr(S[[j]], "B")
    pj    <- attr(S[[j]], "pen")
    theta <- solve(crossprod(B) + diag(pj, ncol(B)),
                   crossprod(B, r + fit$f[[j]]))
    sum(pj * theta^2)
  }, numeric(1)))
  M_p <- 1 + sum(vapply(S, function(S_j) sum(attr(S_j, "pen") == 0),
                        numeric(1)))
  (sum(r^2) + pen) / (length(y) - M_p)
}

# Draws at new x by the basis expansion (43), p. 211, f(x) = b(x)^T theta,
# with theta the exact least-squares fit of the draws on the data basis
# (aliased columns get 0).
ht_at_new_x <- function(draws, B_data, B_new) {
  theta <- qr.coef(qr(B_data), draws)
  theta[is.na(theta)] <- 0
  if (max(abs(B_data %*% theta - draws)) > 1e-8 * max(abs(draws)))
    stop("draws not in the span of the basis")
  B_new %*% theta
}


# ---- Reference, part B: gbgam's prior as the paper's smoothers ---------------
#
#   phi                          -> sigma^2, fixed (p. 198), ht_sigma2()
#   "(Intercept)", flat prior    -> alpha = f_0, S_0 = 1 1^T / n
#   learners of term t2(x_j)     -> columns of B_j, smoother S_j
#   penalised learner: beta_k ~ N(0, phi (lambda_k D_k)^{-1}),
#     lambda_k = exp(rho_k)      -> block lambda_k diag(D_k) of lambda Omega_j
#   linear learner (column x, scaled by gbgam to unit norm, D_k = 1)
#                                -> block lambda_k ||x||^2 on the unscaled x
#   status "zero" (all linear learners when linear = "free")
#                                -> 0: flat, as a smoothing spline's linear
#                                   part in HT
#   status "inf"                 -> column left out
# gbgam's learner columns are centred, so HT's centring changes nothing.
gbgam_as_ht <- function(mod) {
  stopifnot(all(mod$prior.weights == 1), mod$family$family == "gaussian",
            mod$family$link == "identity")
  lam    <- mod$bb_info$lambda
  labels <- vapply(mod$learners, function(lr) lr$label %||% "", "")
  terms  <- unique(lam$term)
  S <- lapply(terms, function(term) {
    idx  <- which(labels == term)
    rows <- lam[lam$term == term, ]
    stopifnot(identical(rows$learner,
                        vapply(mod$learners[idx], `[[`, "", "var_label")))
    keep <- rows$status != "inf"
    B_j  <- do.call(cbind, lapply(mod$learners[idx][keep], function(lr)
      as.matrix(lr$design)))
    pen  <- unlist(Map(function(lr, status, rho) {
      d <- if (inherits(lr, "demmler_reinsch")) lr$D else
        sum(as.matrix(lr$design)^2)
      (if (status == "zero") 0 else exp(rho)) * d
    }, mod$learners[idx][keep], rows$status[keep], rows$rho[keep]))
    structure(ht_smoother(B_j, diag(pen, ncol(B_j))), B = B_j, pen = pen)
  })
  names(S) <- terms
  S
}


# ---- The comparison ----------------------------------------------------------

test_that("bayesian_backfitting() intervals match Hastie-Tibshirani's", {
  skip_on_cran()
  n_iter <- 3000L
  n_burn <- 300L
  set.seed(1)
  n   <- 200
  dat <- data.frame(x1 = stats::runif(n), x2 = stats::runif(n))
  dat$y <- 1 + sin(2 * pi * dat$x1) + 8 * (dat$x2 - 0.5)^2 +
    stats::rnorm(n, sd = 0.4)
  fit <- gbgam(y ~ t2(x1) + t2(x2), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl("aic", resumable = FALSE))
  grid <- seq(0.01, 0.99, length.out = 50)
  newd <- list(data.frame(x1 = grid, x2 = 0.5),
               data.frame(x1 = 0.5, x2 = grid))
  X    <- stats::predict(fit, type = "lpmatrix")
  # smooth_intervals()'s alpha (miscoverage), not HT's intercept alpha
  alpha_ci <- 0.05
  z        <- stats::qnorm(1 - alpha_ci / 2)

  # Monte Carlo error per point, relative to the half-width z sd, for two
  # independent chains with effective size >= n_iter / 2 each (about 0.9
  # n_iter was measured for both samplers): centre sqrt(2 / ess) / z, half-
  # width sqrt(2 / (2 ess)). The grid points are strongly correlated, so 4
  # of these bound the maximum over the grid.
  ess    <- n_iter / 2
  tol_c  <- 4 * sqrt(2 / ess) / z
  tol_h  <- 4 * sqrt(1 / ess)

  for (linear in c("free", "kept")) {
    mod <- bayesian_backfitting(fit, n_samples = n_iter, n_burnin = n_burn,
                                linear = linear)
    S      <- gbgam_as_ht(mod)
    sigma2 <- ht_sigma2(dat$y, S)
    cat(sprintf("linear = %s: phi relative error %.2g\n", linear,
                abs(mod$bb_info$phi / sigma2 - 1)))
    # backfitting converged to 1e-10; measured relative error 1.6e-13
    expect_lt(abs(mod$bb_info$phi / sigma2 - 1), 1e-8)
    ref <- ht_bayesian_backfitting(dat$y, S, sigma2, n_iter, n_burn)
    for (j in seq_along(S)) {
      term <- names(S)[j]
      cols <- c(1, fit$smooth[[j]]$first.para:fit$smooth[[j]]$last.para)
      Xg   <- stats::predict(fit, newdata = newd[[j]],
                             type = "lpmatrix")[, cols]
      draws <- ht_at_new_x(outer(rep(1, n), ref$alpha) + ref$f[[j]],
                           X[, cols], Xg)          # alpha + f_j on the grid
      ref_c <- rowMeans(draws)
      ref_h <- z * apply(draws, 1, stats::sd)
      si    <- smooth_intervals(mod, alpha = alpha_ci, centre = "posterior",
                                newdata = newd[[j]])
      gb_c  <- si$fit[, term]
      gb_h  <- si$upper[, term] - gb_c
      d_c   <- max(abs(gb_c - ref_c) / ref_h)
      d_h   <- max(abs(gb_h - ref_h) / ref_h)
      cat(sprintf(paste("linear = %s, %s: centre %.3f (tol %.3f),",
                        "half-width %.3f (tol %.3f) of the half-width\n"),
                  linear, term, d_c, tol_c, d_h, tol_h))
      expect_lt(d_c, tol_c, label = paste(linear, term, "centre"))
      expect_lt(d_h, tol_h, label = paste(linear, term, "half-width"))
    }
  }
})
