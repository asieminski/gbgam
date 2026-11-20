# bayesian_backfitting(): samplers, priors, posterior mean, argument checks
# and failures.

fit_bb_model <- function(formula, data, family, max_iter = 150L, ...) {
  gbgam(
    formula,
    data         = data,
    family       = family,
    stop_control = gb_stop_ctrl("none", max_iter = max_iter),
    edf_control  = gb_edf_ctrl("exact"),
    verbose      = FALSE,
    ...
  )
}

mean_link_se <- function(mod) {
  mean(stats::predict(mod, type = "link", se.fit = TRUE)$se.fit)
}

posterior_mean_eta <- function(mod) {
  X <- stats::predict(mod, type = "lpmatrix")
  as.vector(X %*% mod$coefficients_bayes) + mod$offset
}

# The prior the posterior used, on the learner (Demmler-Reinsch) scale, read
# from the learners and bb$bb_info$lambda: columns X of the learners in the
# posterior (lambda = Inf leaves a learner out), and the prior precision
# times phi, cDlam (0 for the intercept and learners at lambda = 0)
posterior_prior <- function(bb) {
  lam   <- bb$bb_info$lambda
  lrs   <- bb$learners
  w     <- bb$prior.weights
  stat  <- vapply(lrs, function(lr) {
    if (identical(lr$var_label, "(Intercept)")) "zero"
    else lam$status[lam$learner == lr$var_label]
  }, character(1))
  keep  <- stat != "inf"
  cDlam <- unlist(lapply(which(keep), function(k) {
    lr <- lrs[[k]]
    if (stat[k] == "zero") return(rep(0, length(lr$beta)))
    rho <- lam$rho[lam$learner == lr$var_label]
    if (inherits(lr, "demmler_reinsch")) exp(rho) * lr$D
    else exp(rho) * sum(w * lr$design^2)  # unit W-norm column scaled back
  }))
  is_int <- vapply(lrs[keep], function(lr)
    identical(lr$var_label, "(Intercept)"), logical(1))
  ncols  <- vapply(lrs[keep], function(lr) length(lr$beta), integer(1))
  list(X = do.call(cbind, lapply(lrs[keep], function(lr)
         as.matrix(lr$design))),
       cDlam = cDlam, keep = keep, is_int = rep(is_int, ncols))
}

# Exact posterior of a Gaussian identity model: mode mu, REML dispersion
# phi = (RSS + mu^T S mu) / (n_+ - M_p) with M_p the number of flat columns
# (see ?bayesian_backfitting) and covariance Sigma
exact_gaussian_posterior <- function(bb) {
  pr  <- posterior_prior(bb)
  X   <- pr$X
  w   <- bb$prior.weights
  y   <- bb$y - bb$offset
  P   <- crossprod(X * sqrt(w))
  diag(P) <- diag(P) + pr$cDlam
  mu  <- as.vector(solve(P, crossprod(X, w * y)))
  r   <- y - as.vector(X %*% mu)
  M_p <- sum(pr$cDlam == 0)
  phi <- (sum(w * r^2) + sum(pr$cDlam * mu^2)) / (sum(w > 0) - M_p)
  c(pr, list(mu = mu, phi = phi, M_p = M_p, Sigma = phi * solve(P)))
}

# Random-walk Metropolis on the log posterior log_post: proposals
# b + step * L z with z standard normal; the draws after n_burn, one per
# column
rw_metropolis <- function(log_post, start, L, step, n_iter, n_burn) {
  b     <- start
  lp_b  <- log_post(b)
  draws <- matrix(NA_real_, length(b), n_iter)
  for (t in seq_len(n_burn + n_iter)) {
    cand    <- b + step * as.vector(L %*% stats::rnorm(length(b)))
    lp_cand <- log_post(cand)
    if (log(stats::runif(1)) < lp_cand - lp_b) {
      b    <- cand
      lp_b <- lp_cand
    }
    if (t > n_burn) draws[, t - n_burn] <- b
  }
  draws
}

# Largest error of the posterior means and sds of draws (one row per
# coefficient) against reference draws, in reference sds
draws_error <- function(draws, ref) {
  ref_sd <- apply(ref, 1, stats::sd)
  c(mean = max(abs(rowMeans(draws) - rowMeans(ref)) / ref_sd),
    sd   = max(abs(apply(draws, 1, stats::sd) / ref_sd - 1)))
}

# Number of sweeps in which the draws (one column per sweep) changed
n_moves <- function(draws) {
  last <- ncol(draws)
  sum(colSums(draws[, -1L, drop = FALSE] != draws[, -last, drop = FALSE]) > 0)
}


# ---- Samplers ----------------------------------------------------------------

test_that("MH sampler mixes for Poisson and gives SEs comparable to mgcv", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 200, seed = 3)
  mod <- fit_bb_model(y ~ t2(x, k = 8, bs = "cr"), dat, stats::poisson(),
                      max_iter = 300L)
  set.seed(1)
  bb <- bayesian_backfitting(mod, n_samples = 200L, n_burnin = 50L)

  ar <- bb$mh_accept_rate
  expect_true(all(ar[!is.na(ar)] > 0.3),
              info = paste("accept rates:", paste(round(ar, 3), collapse = " ")))
  # The proposal is not the exact conditional, so some draws are rejected
  # (largest rate here 0.95)
  expect_lt(max(ar, na.rm = TRUE), 0.999)

  se_bb   <- mean_link_se(bb)
  ref     <- mgcv::gam(y ~ t2(x, k = 8, bs = "cr"), data = dat,
                       family = stats::poisson(), method = "REML")
  se_mgcv <- mean(stats::predict(ref, type = "link", se.fit = TRUE)$se.fit)
  # A chain that does not move gives se_bb / se_mgcv of about 0.25
  expect_gt(se_bb / se_mgcv, 0.4)
  expect_lt(se_bb / se_mgcv, 2.5)
})

test_that("MH sampler mixes for binomial (aggregated, weighted) data", {
  skip_on_cran()
  dat <- simulate_binomial_prop(n = 120, seed = 5)
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::binomial(),
                      max_iter = 300L, weights = dat$w)
  set.seed(2)
  bb <- bayesian_backfitting(mod, n_samples = 150L, n_burnin = 50L)
  ar <- bb$mh_accept_rate
  expect_true(all(ar[!is.na(ar)] > 0.3),
              info = paste("accept rates:", paste(round(ar, 3), collapse = " ")))
  se <- stats::predict(bb, type = "link", se.fit = TRUE)$se.fit
  expect_true(all(is.finite(se) & se > 0))
})

test_that("MH on a Gaussian identity model reproduces the Gibbs posterior", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 150, seed = 1, sd = 0.5)
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::gaussian(),
                      max_iter = 200L)
  expect_true(mod$scale < 0.5)  # phi != 1, so the 1 / phi prior scaling matters

  # For identity-link Gaussian the IWLS proposal is the exact full conditional,
  # so MH must accept (almost) always and target the same posterior as Gibbs.
  # A copy with another family name is routed to the MH sampler.
  # It has no saturated likelihood for the dispersion, so it gets the
  # Gibbs chain's phi.
  mod_mh <- mod
  mod_mh$family$family <- "gaussian_mh"
  set.seed(3)
  b_gb <- bayesian_backfitting(mod, n_samples = 1500L, n_burnin = 50L)
  set.seed(3)
  b_mh <- bayesian_backfitting(mod_mh, n_samples = 1500L, n_burnin = 50L,
                               phi = b_gb$bb_info$phi)

  expect_false(is_gaussian_identity(mod_mh$family))
  ar <- b_mh$mh_accept_rate
  expect_true(all(ar[!is.na(ar)] > 0.99))
  sd_mh <- sqrt(diag(b_mh$Vp))
  sd_gb <- sqrt(diag(b_gb$Vp))
  keep  <- sd_gb > 1e-8
  expect_true(all(abs(sd_mh[keep] / sd_gb[keep] - 1) < 0.15),
              info = paste(round(sd_mh[keep] / sd_gb[keep], 3), collapse = " "))
})

test_that("Gaussian with log link is routed to MH on the link scale", {
  skip_on_cran()
  set.seed(2)
  n   <- 200
  x   <- sort(stats::runif(n))
  eta <- 2 + sin(2 * pi * x)
  dat <- data.frame(x = x, y = exp(eta) + stats::rnorm(n, sd = 0.5))
  # The gradient of the gaussian(log) loss scales with mu^2 (up to ~400
  # here), so plain gradient steps are only stable for step_size < 2 / 400.
  # With 0.01 the boosting fit diverges at iteration 34 and the default 0.1
  # at once (gbgam() stops with an error);
  # 0.002 converges to the mgcv fit (same RSS) within 1500 iterations.
  mod <- suppressMessages(fit_bb_model(
    y ~ t2(x, k = 6, bs = "cr"), dat, stats::gaussian(link = "log"),
    max_iter        = 1500L,
    learner_control = gb_learner_ctrl(step_size = 0.002)
  ))
  set.seed(4)
  bb <- bayesian_backfitting(mod, n_samples = 150L, n_burnin = 50L)

  # Routed to the MH sampler: Gibbs is for the identity link only
  expect_false(is_gaussian_identity(mod$family))
  ar <- bb$mh_accept_rate
  expect_true(all(ar[!is.na(ar)] > 0.3))
  # Posterior is centred on the link scale (Gibbs gave eta in [-7, 10])
  expect_lt(max(abs(posterior_mean_eta(bb) - mod$linear.predictors)), 0.3)

  ref     <- mgcv::gam(y ~ t2(x, k = 6, bs = "cr"), data = dat,
                       family = stats::gaussian(link = "log"), method = "REML")
  se_mgcv <- mean(stats::predict(ref, type = "link", se.fit = TRUE)$se.fit)
  expect_gt(mean_link_se(bb) / se_mgcv, 0.4)
  expect_lt(mean_link_se(bb) / se_mgcv, 2.5)
})

test_that("Gibbs draws, phi and posterior mean are exact (Gaussian)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 120, seed = 9, sd = 0.3)
  set.seed(2)
  weights <- list(unit = rep(1, 120), varying = stats::runif(120, 0.2, 3))
  for (nm in names(weights)) {
    # 40 iterations, so that the prior matters (cDlam from 0.2 to 15)
    mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::gaussian(),
                        max_iter = 40L, weights = weights[[nm]])
    set.seed(10)
    bb <- bayesian_backfitting(mod, n_samples = 3000L, n_burnin = 100L)
    ex <- exact_gaussian_posterior(bb)
    sd_ex <- sqrt(diag(ex$Sigma))
    draws <- do.call(rbind, bb$beta_samples[ex$keep])
    # the posterior mean on the mgcv scale, as a linear predictor
    eta_ex  <- bb$offset + as.vector(ex$X %*% ex$mu)
    sd_eta  <- sqrt(rowSums((ex$X %*% ex$Sigma) * ex$X))
    err_eta <- max(abs(posterior_mean_eta(bb) - eta_ex) / sd_eta)
    err     <- c(phi  = abs(bb$bb_info$phi / ex$phi - 1),
                 mean = max(abs(rowMeans(draws) - ex$mu) / sd_ex),
                 sd   = max(abs(apply(draws, 1, stats::sd) / sd_ex - 1)),
                 eta  = err_eta)
    cat(sprintf("%s weights: %s\n", nm,
                paste(names(err), signif(err, 2), collapse = ", ")))
    # phi: measured relative error 0 (the penalty left out of D_p is off
    # by 0.22)
    expect_lt(err[["phi"]], 1e-8, label = paste(nm, "phi"))
    expect_identical(ex$M_p, 2L)  # intercept and the freed linear learner
    # 3000 draws, effective size >= 2600 measured: Monte Carlo standard
    # errors 0.020 (mean) and 0.014 (sd), in sds; measured 0.033, 0.035 and
    # 0.033 (eta), at most 0.033, 0.048 and 0.038 over 6 chain seeds. A
    # Gibbs mean without the weights is off by 4 sds, the prior left out by
    # 12 sds, a prior precision divided by phi by 5.6 sds
    expect_lt(err[["mean"]], 0.1, label = paste(nm, "mean"))
    expect_lt(err[["sd"]], 0.07, label = paste(nm, "sd"))
    expect_lt(err[["eta"]], 0.1, label = paste(nm, "eta"))
  }
})

test_that("a given phi is used, and Gibbs draws scale with sqrt(phi)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 120, seed = 9, sd = 0.3)
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::gaussian(),
                      max_iter = 40L)
  set.seed(10)
  bb_est <- bayesian_backfitting(mod, n_samples = 200L, n_burnin = 10L)
  phi_est <- bb_est$bb_info$phi
  set.seed(10)
  bb_fix <- bayesian_backfitting(mod, n_samples = 200L, n_burnin = 10L,
                                 phi = 4 * phi_est)
  ex <- exact_gaussian_posterior(bb_est)
  # Gibbs draws start at the mode and are affine in the normal draws, so
  # with the same random numbers they are the mode plus sqrt(phi) times
  # the same deviations
  dev_est <- do.call(rbind, bb_est$beta_samples[ex$keep]) - ex$mu
  dev_fix <- do.call(rbind, bb_fix$beta_samples[ex$keep]) - ex$mu
  err <- max(abs(dev_fix - 2 * dev_est)) / max(abs(dev_est))
  cat(sprintf("phi given: relative error of the scaled draws %.2g\n", err))
  expect_identical(bb_fix$bb_info$phi, 4 * phi_est)
  # measured 1e-15; a draw not scaled by sqrt(phi) is off by 0.5
  expect_lt(err, 1e-8)
})

test_that("phi is ignored for Poisson and checked", {
  dat <- simulate_poisson_1d(n = 100, seed = 3)
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::poisson(),
                      max_iter = 30L)
  bb <- bayesian_backfitting(mod, n_samples = 5L, n_burnin = 0L, phi = 3)
  expect_identical(bb$bb_info$phi, 1)
  expect_error(bayesian_backfitting(mod, phi = 0), "phi must be positive")
  expect_error(bayesian_backfitting(mod, phi = -1), "phi")
  expect_error(bayesian_backfitting(mod, phi = Inf), "phi")
})

test_that("MH draws match a random-walk Metropolis reference (Poisson)", {
  skip_on_cran()
  # 40 iterations leave the smooth shrunk: log(lambda) = 4.2 at Stage 1, and
  # the prior moves the posterior mode by up to 4 posterior sds
  dat <- simulate_poisson_1d(n = 200, seed = 3)
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::poisson(),
                      max_iter = 40L)
  set.seed(4)
  bb <- bayesian_backfitting(mod, n_samples = 4000L, n_burnin = 200L)

  pr    <- posterior_prior(bb)
  X     <- pr$X
  cDlam <- pr$cDlam
  y     <- dat$y
  # log likelihood plus log prior, up to a constant (phi = 1)
  log_post <- function(b) {
    eta <- bb$offset + as.vector(X %*% b)
    sum(y * eta - exp(eta)) - 0.5 * sum(cDlam * b^2)
  }
  b <- unlist(lapply(bb$learners[pr$keep], `[[`, "beta"))
  for (i in 1:50) {                        # Newton steps to the mode
    mu <- exp(bb$offset + as.vector(X %*% b))
    H  <- crossprod(X * mu, X)              # minus the Hessian of log_post
    diag(H) <- diag(H) + cDlam
    b  <- b + as.vector(solve(H, crossprod(X, y - mu) - cDlam * b))
  }
  set.seed(5)
  ref <- rw_metropolis(log_post, b, t(chol(solve(H))), step = 0.8,
                       n_iter = 6e4L, n_burn = 5e3L)
  err <- draws_error(do.call(rbind, bb$beta_samples[pr$keep]), ref)
  ar  <- bb$mh_accept_rate
  cat(sprintf("Poisson: mean %.3f, sd %.3f; accept rates %s\n", err[["mean"]],
              err[["sd"]], paste(round(ar, 3), collapse = " ")))
  # Effective sizes >= 2400 (sampler) and >= 2900 (reference) measured:
  # Monte Carlo standard errors of the difference 0.028 (mean) and 0.019
  # (sd), in sds; measured 0.033 and 0.046. Without the proposal densities
  # in the acceptance ratio the sds are off by 0.32
  expect_lt(err[["mean"]], 0.12)
  expect_lt(err[["sd"]], 0.1)

  expect_lt(max(ar, na.rm = TRUE), 0.999)
  # A penalised learner's draws change only when its own block is accepted;
  # the first kept draw's step is not seen, hence the 1
  for (k in which(vapply(bb$learners, inherits, TRUE, "demmler_reinsch"))) {
    gap <- ar[k] * ncol(bb$beta_samples[[k]]) - n_moves(bb$beta_samples[[k]])
    expect_true(-1e-8 < gap && gap < 1 + 1e-8,
                label = paste("accepted minus moved draws:", gap))
  }
})

test_that("MH matches the exact posterior of a skewed Poisson model", {
  skip_on_cran()
  # 12 counts, 8 of them 0: the posterior is far from normal, so the
  # proposal changes with the state of the chain
  set.seed(5)
  n     <- 12
  dat   <- data.frame(x = stats::rnorm(n))
  dat$y <- stats::rpois(n, exp(-0.7 + 0.8 * dat$x))
  mod   <- fit_bb_model(y ~ x, dat, stats::poisson(), max_iter = 300L)
  set.seed(1)
  bb <- bayesian_backfitting(mod, n_samples = 1e4L, n_burnin = 200L)
  expect_identical(bb$bb_info$lambda$status, "zero")

  # Flat prior on both coefficients (linear = "free"): the posterior is the
  # normalised likelihood, on a grid around the glm fit (moments as on a
  # grid 1.2 times as wide to 6 digits)
  g  <- stats::glm(y ~ x, family = stats::poisson(), data = dat)
  se <- sqrt(diag(stats::vcov(g)))
  a  <- seq(coef(g)[[1]] - 14 * se[[1]], coef(g)[[1]] + 6 * se[[1]],
            length.out = 400)
  b  <- seq(coef(g)[[2]] - 10 * se[[2]], coef(g)[[2]] + 10 * se[[2]],
            length.out = 400)
  ll <- matrix(0, length(a), length(b))   # log likelihood at (a[i], b[j])
  for (i in seq_len(n)) {
    eta <- outer(a, b * dat$x[i], `+`)
    ll  <- ll + dat$y[i] * eta - exp(eta)
  }
  p  <- exp(ll - max(ll))
  p  <- p / sum(p)
  edge <- max(p[c(1, nrow(p)), ], p[, c(1, ncol(p))])
  expect_lt(edge, 1e-11)       # the grid holds the posterior (1.3e-12)
  pa <- rowSums(p)
  pb <- colSums(p)
  m  <- c(sum(pa * a), sum(pb * b))
  s  <- sqrt(c(sum(pa * (a - m[1])^2), sum(pb * (b - m[2])^2)))
  err <- c(mean = max(abs(bb$coefficients_bayes - m) / s),
           sd   = max(abs(sqrt(diag(bb$Vp)) / s - 1)))
  cat(sprintf("skewed Poisson: mean %.3f, sd %.3f\n", err[["mean"]],
              err[["sd"]]))
  # Effective size >= 1200 measured: Monte Carlo standard errors 0.03
  # (mean) and about 0.03 (sd), in sds; measured 0.005 and 0.017 (at most
  # 0.021 and 0.041 over 6 seeds). A reverse proposal built at the current
  # state is off by 1.1 and 0.9 here (at least 0.37 and 0.23 over 6 seeds);
  # the boosted coefficients are 0.40 sds from the posterior mean
  expect_lt(err[["mean"]], 0.15)
  expect_lt(err[["sd"]], 0.15)
})

test_that("MH draws and phi match exact references for Gamma (log link)", {
  skip_on_cran()
  dat <- simulate_gamma_1d(n = 200, seed = 6)
  mod <- suppressMessages(fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat,
                                       stats::Gamma(link = "log"),
                                       max_iter = 40L))
  set.seed(4)
  bb <- bayesian_backfitting(mod, n_samples = 4000L, n_burnin = 200L)

  pr    <- posterior_prior(bb)
  X     <- pr$X
  cDlam <- pr$cDlam
  y     <- dat$y
  # The mode minimises deviance / 2 + b^T S b / 2 (phi cancels); Newton
  # steps, with the Gamma deviance's -log(y / mu) + (y - mu) / mu per row
  b <- unlist(lapply(bb$learners[pr$keep], `[[`, "beta"))
  for (i in 1:100) {
    mu <- exp(bb$offset + as.vector(X %*% b))
    H  <- crossprod(X * (y / mu), X)
    diag(H) <- diag(H) + cDlam
    b  <- b + as.vector(solve(H, crossprod(X, y / mu - 1) - cDlam * b))
  }
  mu <- exp(bb$offset + as.vector(X %*% b))
  # REML phi at fixed lambda (Wood, 2011, eq. 4): minimise
  # D_p / (2 phi) - l_s(phi) - M_p / 2 log(2 pi phi), with D_p the penalised
  # deviance and l_s the saturated Gamma log-likelihood
  D_p  <- 2 * sum(-log(y / mu) + (y - mu) / mu) + sum(cDlam * b^2)
  M_p  <- sum(cDlam == 0)
  ls   <- mgcv::fix.family.ls(stats::Gamma())$ls
  nlr  <- function(lphi) {
    phi <- exp(lphi)
    D_p / (2 * phi) - ls(y, rep(1, length(y)), length(y), phi)[1] -
      M_p / 2 * log(2 * pi * phi)
  }
  phi  <- exp(stats::optimize(nlr, c(-10, 5), tol = 1e-12)$minimum)

  # log likelihood -deviance / (2 phi) plus log prior, up to a constant
  log_post <- function(b) {
    eta <- bb$offset + as.vector(X %*% b)
    (sum(-eta - y * exp(-eta)) - 0.5 * sum(cDlam * b^2)) / phi
  }
  set.seed(5)
  ref <- rw_metropolis(log_post, b, t(chol(phi * solve(H))), step = 0.8,
                       n_iter = 6e4L, n_burn = 5e3L)
  err <- c(phi = abs(bb$bb_info$phi / phi - 1),
           draws_error(do.call(rbind, bb$beta_samples[pr$keep]), ref))
  cat(sprintf("Gamma: phi %.2g, mean %.3f, sd %.3f\n", err[["phi"]],
              err[["mean"]], err[["sd"]]))
  # phi at the mode of the fitting sweeps: measured relative error 3.5e-10
  # (D_p / (n - M_p), the Gaussian closed form, is off by 0.094)
  expect_lt(err[["phi"]], 1e-6)
  # Effective sizes >= 2000 (sampler) and >= 2900 (reference) measured:
  # Monte Carlo standard errors of the difference 0.029 (mean) and 0.020
  # (sd), in sds; measured 0.048 and 0.055. phi = 1 makes the sds off by
  # 0.35
  expect_lt(err[["mean"]], 0.12)
  expect_lt(err[["sd"]], 0.1)
})

test_that("the burn-in sweeps are dropped", {
  mods <- list(
    gaussian = fit_bb_model(y ~ t2(x, k = 6, bs = "cr"),
                            simulate_gaussian_1d(n = 100, seed = 11),
                            stats::gaussian(), max_iter = 100L),
    poisson  = fit_bb_model(y ~ t2(x, k = 6, bs = "cr"),
                            simulate_poisson_1d(n = 100, seed = 11),
                            stats::poisson(), max_iter = 100L))
  for (nm in names(mods)) {
    # Stage 1 and the posterior mode use no random numbers, so at the same
    # seed both chains run the same sweeps
    set.seed(7)
    burnt <- bayesian_backfitting(mods[[nm]], n_samples = 10L, n_burnin = 5L)
    set.seed(7)
    all15 <- bayesian_backfitting(mods[[nm]], n_samples = 15L, n_burnin = 0L)
    expect_identical(burnt$beta_samples,
                     lapply(all15$beta_samples,
                            function(b) b[, 6:15, drop = FALSE]),
                     info = nm)
    # without burn-in the first kept draw is already a draw, not the
    # posterior mode the chain starts from: it depends on the seed
    set.seed(8)
    other <- bayesian_backfitting(mods[[nm]], n_samples = 2L, n_burnin = 0L)
    first <- function(bb) lapply(bb$beta_samples, function(b) b[, 1L])
    expect_false(identical(first(other), first(all15)), info = nm)
  }
})


# ---- Priors and the posterior mean -------------------------------------------

test_that("binomial intercept is not shrunk towards zero", {
  skip_on_cran()
  set.seed(6)
  n   <- 300
  x   <- sort(stats::runif(n))
  dat <- data.frame(x = x,
                    y = stats::rbinom(n, 1, stats::plogis(2 + sin(2 * pi * x))))
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::binomial(),
                      max_iter = 300L)
  set.seed(5)
  bb <- bayesian_backfitting(mod, n_samples = 200L, n_burnin = 50L)
  b0_boost <- unname(mod$coefficients["(Intercept)"])
  b0_bayes <- unname(bb$coefficients_bayes["(Intercept)"])
  sd0      <- sqrt(bb$Vp["(Intercept)", "(Intercept)"])
  expect_gt(b0_boost, 1)
  expect_lt(abs(b0_bayes - b0_boost), 4 * sd0)
})

test_that("learners excluded by the user stay at zero", {
  skip_on_cran()
  dat   <- simulate_gaussian_2d(n = 150, seed = 2)
  dat$u <- stats::rnorm(150)
  mod   <- fit_bb_model(
    y ~ u + t2(x, k = 6, bs = "cr"), dat, stats::gaussian(),
    max_iter        = 100L,
    learner_control = gb_learner_ctrl(excluded_learner_indices = 2L)
  )
  expect_identical(mod$learners[[2]]$beta, 0)

  set.seed(6)
  bb <- bayesian_backfitting(mod, n_samples = 100L, n_burnin = 10L)
  expect_true(all(bb$mh_accept_rate == 1))  # Gibbs contract
  expect_true(all(apply(bb$beta_samples[[2]], 1, stats::var) == 0))
  expect_identical(unname(bb$Vp["u", "u"]), 0)
  expect_identical(unname(bb$coefficients_bayes["u"]), 0)
  expect_gt(bb$Vp["(Intercept)", "(Intercept)"], 0)
  lam <- bb$bb_info$lambda
  expect_identical(lam$status[lam$learner == "u"], "inf")

  # Same for the MH sampler
  dat$yc <- stats::rpois(150, exp(0.5 + 0.5 * sin(2 * pi * dat$x)))
  mod_p  <- fit_bb_model(
    yc ~ u + t2(x, k = 6, bs = "cr"), dat, stats::poisson(),
    max_iter        = 100L,
    learner_control = gb_learner_ctrl(excluded_learner_indices = 2L)
  )
  set.seed(7)
  bb_p <- bayesian_backfitting(mod_p, n_samples = 50L, n_burnin = 10L)
  expect_true(is.na(bb_p$mh_accept_rate[2]))
  expect_false(anyNA(bb_p$mh_accept_rate[-2]))
  expect_identical(unname(bb_p$Vp["u", "u"]), 0)
})

test_that("an excluded intercept fitted after each iteration stays active", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 120, seed = 4)
  mod <- fit_bb_model(
    y ~ t2(x, k = 6, bs = "cr"), dat, stats::gaussian(),
    max_iter        = 60L,
    learner_control = gb_learner_ctrl(learners_to_fit_after_each_iter = 1L,
                                      excluded_learner_indices        = 1L)
  )
  set.seed(8)
  bb <- bayesian_backfitting(mod, n_samples = 50L, n_burnin = 5L)
  expect_identical(bb$mh_accept_rate[1], 1)
  expect_gt(bb$Vp["(Intercept)", "(Intercept)"], 0)
})

test_that("an intercept held at zero by boosting stays at zero", {
  skip_on_cran()
  ctrl <- gb_learner_ctrl(initialise_intercept     = FALSE,
                          excluded_learner_indices = 1L)
  fams <- list(stats::gaussian(), stats::poisson())     # Gibbs and MH
  for (fam in fams) {
    dat <- if (fam$family == "gaussian") {
      simulate_gaussian_1d(n = 100, seed = 8)
    } else {
      simulate_poisson_1d(n = 100, seed = 8)
    }
    mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, fam,
                        max_iter = 60L, learner_control = ctrl)
    expect_identical(mod$learners[[1]]$beta, 0)
    m <- model_as_smoothers(mod)
    expect_false(m$has_int)
    # M_p counts only the flat columns of the smooths
    sp <- map_smoothing_parameters(m)$sp
    cDlam <- smoother_prior(m$smoothers[[1]], sp[[1]])$cDlam
    expect_identical(penalty_null_dim(m, sp), sum(cDlam == 0))

    set.seed(9)
    bb <- bayesian_backfitting(mod, n_samples = 30L, n_burnin = 5L)
    expect_true(all(bb$beta_samples[[1]] == 0))
    expect_true(all(bb$Vp["(Intercept)", ] == 0))
    expect_true(all(bb$Vp[, "(Intercept)"] == 0))
    expect_identical(unname(bb$coefficients_bayes["(Intercept)"]),
                     unname(mod$coefficients["(Intercept)"]))
    expect_identical(bb$mh_accept_rate[1],
                     if (fam$family == "gaussian") 1 else NA_real_)
    expect_gt(bb$Vp[2, 2], 0)
  }
})

test_that("coefficients_bayes is the posterior mean on the mgcv scale", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 11)
  mod <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr"), dat, stats::gaussian(),
                      max_iter = 100L)
  set.seed(12)
  bb <- bayesian_backfitting(mod, n_samples = 50L, n_burnin = 5L,
                             linear = "kept")
  expect_identical(names(bb$coefficients_bayes), names(mod$coefficients))
  expect_true(all(is.finite(bb$coefficients_bayes)))
  # Intervals from predict() remain centred on the boosting estimate
  expect_identical(bb$coefficients, mod$coefficients)
  # With the linear learners at their Stage-1 lambda the posterior mean is
  # close to the boosted fit (freed, it moves along the linear directions,
  # see ?bayesian_backfitting)
  eta_bayes <- posterior_mean_eta(bb)
  expect_lt(max(abs(eta_bayes - mod$linear.predictors)), 0.5)
})


# ---- Argument checks, dispersion and failures --------------------------------

test_that("bayesian_backfitting requires n_samples >= 2, not EDF tracking", {
  dat <- simulate_gaussian_1d(n = 60, seed = 1)
  mod <- gbgam(y ~ t2(x, k = 5, bs = "cr"), data = dat,
               stop_control = gb_stop_ctrl("none", max_iter = 10L),
               edf_control  = gb_edf_ctrl("exact"), verbose = FALSE)
  expect_error(bayesian_backfitting(mod, n_samples = 1L, n_burnin = 0L),
               "n_samples")

  # Stage 1 uses the boosted fit only, so EDF tracking is not needed
  mod_none <- gbgam(y ~ t2(x, k = 5, bs = "cr"), data = dat,
                    stop_control = gb_stop_ctrl("none", max_iter = 10L),
                    edf_control  = gb_edf_ctrl("none"), verbose = FALSE)
  set.seed(1)
  out <- bayesian_backfitting(mod_none, n_samples = 5L, n_burnin = 0L)
  expect_true(all(is.finite(out$Vp)))
})

test_that("time limit with fewer than 2 samples returns the model unchanged", {
  dat  <- simulate_poisson_1d(n = 60, seed = 1)
  mods <- lapply(list(stats::poisson(), stats::gaussian()), function(fam) {
    gbgam(y ~ t2(x, k = 5, bs = "cr"), data = dat, family = fam,
          stop_control = gb_stop_ctrl("none", max_iter = 10L),
          edf_control  = gb_edf_ctrl("exact"), verbose = FALSE)
  })
  for (fam_mod in mods) {
    expect_warning(
      out <- bayesian_backfitting(fam_mod, n_samples = 5L, n_burnin = 0L,
                                  max_runtime_mins = 0),
      "fewer than 2 posterior samples"
    )
    expect_identical(out, fam_mod)
  }
})

test_that("rows with zero weight leave the dispersion of the posterior unchanged", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 120, seed = 4, sd = 0.4)
  fit <- function(d, w) fit_bb_model(y ~ t2(x, k = 8, bs = "cr"), d,
                                     stats::gaussian(), max_iter = 60L,
                                     weights = w)
  set.seed(1)
  bb1 <- bayesian_backfitting(fit(dat, rep(1, 120)), n_samples = 5L,
                              n_burnin = 0L)
  set.seed(1)
  bb2 <- bayesian_backfitting(fit(rbind(dat, dat), rep(1:0, each = 120)),
                              n_samples = 5L, n_burnin = 0L)
  expect_equal(bb2$bb_info$phi, bb1$bb_info$phi, tolerance = 1e-8)
})

test_that("dispersion_at_mode() is the REML estimate, not the Pearson one", {
  dat   <- simulate_gaussian_1d(n = 100, seed = 3, sd = 0.3)
  dat$o <- seq(-0.2, 0.2, length.out = 100)
  set.seed(1)
  w     <- c(stats::runif(97, 0.3, 2), 0, 0, 0)
  mod   <- fit_bb_model(y ~ t2(x, k = 6, bs = "cr") + offset(o), dat,
                        stats::gaussian(), max_iter = 30L, weights = w)
  m     <- model_as_smoothers(mod)
  s1    <- map_smoothing_parameters(m)
  sp    <- s1$sp
  sm    <- m$smoothers[[1]]
  sp[[1]]$status[sm$learners$kind == "null"] <- "zero"   # freed linear
  sp[[1]]$rho[sm$learners$kind == "null"]    <- NA_real_
  st    <- penalised_fit(m, sp, s1$state, s1$tol, s1$max_inner,
                         joint_linear = TRUE)$state
  cDlam <- smoother_prior(sm, sp[[1]])$cDlam_all
  b     <- st$beta[[1]]
  pos   <- m$w > 0
  rss   <- sum(m$w * (m$y - st$eta)^2)
  M_p   <- 1L + sum(cDlam == 0)
  phi   <- (rss + sum(cDlam * b^2)) / (sum(pos) - M_p)
  expect_identical(penalty_null_dim(m, sp), 2L)
  expect_equal(dispersion_at_mode(m, sp, st), phi, tolerance = 1e-12)
  # the Pearson estimate with edf = 1 + tr F differs
  tw      <- smoother_working_model(m, 1, st, "test")
  edf     <- 1 + sum(smoother_shrinkage(sm, tw, sp[[1]], "test")$edf)
  pearson <- rss / (sum(pos) - edf)
  expect_gt(abs(pearson / phi - 1), 1e-3)
  # a given phi is returned unchanged, without the df check
  m1 <- m
  m1$w[-1] <- 0
  expect_identical(dispersion_at_mode(m1, sp, st, phi = 2), 2)
})

test_that("dispersion_at_mode() is 1 for Poisson and binomial", {
  for (fam in list(stats::poisson(), stats::binomial())) {
    dat <- if (fam$family == "poisson") simulate_poisson_1d(n = 80, seed = 2)
           else simulate_binomial_1d(n = 80, seed = 2)
    mod <- suppressMessages(fit_bb_model(y ~ t2(x, k = 5, bs = "cr"), dat,
                                         fam, max_iter = 10L))
    m   <- model_as_smoothers(mod)
    s1  <- map_smoothing_parameters(m)
    expect_identical(dispersion_at_mode(m, s1$sp, s1$state), 1)
  }
})

test_that(paste("dispersion_at_mode() stops with class bb_failed without",
                "residual df"), {
  dat <- simulate_gaussian_1d(n = 60, seed = 1)
  mod <- fit_bb_model(y ~ t2(x, k = 5, bs = "cr"), dat, stats::gaussian(),
                      max_iter = 10L)
  m   <- model_as_smoothers(mod)
  s1  <- map_smoothing_parameters(m)
  m$w[-1] <- 0                                  # one row with weight
  err <- tryCatch(dispersion_at_mode(m, s1$sp, s1$state), error = function(e) e)
  expect_s3_class(err, "bb_failed")
  expect_identical(err$status, "no_residual_df")
})

test_that("numerical failures are signalled as map_approx_failed or bb_failed", {
  not_pd <- function() checked_chol(matrix(-1), "a block")
  nonfin <- function()
    signal_numerical("b2b_nonfinite", "a block: non-finite mu")
  for (failed in list(map_failed_condition, backfitting_failed_condition)) {
    cls <- class(failed("", ""))[1L]
    e1  <- tryCatch(reclass_numerical_failure(not_pd(), failed),
                    error = function(e) e)
    e2  <- tryCatch(reclass_numerical_failure(nonfin(), failed),
                    error = function(e) e)
    expect_s3_class(e1, cls)
    expect_s3_class(e2, cls)
    expect_identical(c(e1$status, e2$status), c("not_pd", "nonfinite"))
    expect_match(conditionMessage(e1), "^bayesian_backfitting: .*a block")
  }
})

test_that("allocation failures are told apart from numerical failures", {
  # Messages of R's allocators (R_alloc, R_Calloc); checked_chol() and
  # matrix_root() pass these on instead of reporting not_pd
  alloc <- c("cannot allocate vector of size 8.0 Gb",
             "'R_Calloc' could not allocate memory (1 of 8 bytes)",
             "vector memory limit of 16.0 Gb reached, see mem.maxVSize()")
  for (msg in alloc) expect_true(is_alloc_error(simpleError(msg)))
  not_pd <- tryCatch(chol(matrix(-1)), error = function(e) e)
  expect_false(is_alloc_error(not_pd))
})

test_that("a parametric column with zero weighted norm stays out of the posterior", {
  skip_on_cran()
  set.seed(3)
  n   <- 150
  dat <- data.frame(x = stats::runif(n),
                    g = factor(sample(c("a", "b", "c"), n, TRUE)))
  dat$y <- sin(2 * pi * dat$x) + (dat$g == "b") + stats::rnorm(n, 0, 0.3)
  w   <- ifelse(dat$g == "c", 0, 1)
  mod <- fit_bb_model(y ~ g + t2(x, k = 8, bs = "cr"), dat, stats::gaussian(),
                      max_iter = 60L, weights = w)
  for (linear in c("free", "kept")) {
    set.seed(1)
    bb  <- bayesian_backfitting(mod, n_samples = 5L, n_burnin = 0L,
                                linear = linear)
    lam <- bb$bb_info$lambda
    expect_identical(lam$status[lam$learner == "gc"], "inf")
    k   <- which(vapply(mod$learners, `[[`, "", "var_label") == "gc")
    expect_true(all(bb$beta_samples[[k]] == 0))
    expect_true(all(is.finite(bb$Vp)))
  }
})

test_that("bayesian_backfitting() removes earlier bootstrap results", {
  dat <- simulate_gaussian_1d(n = 80, seed = 22)
  mod <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl(max_iter = 300))
  set.seed(1)
  b  <- bootstrap_gbgam(mod, n_samples = 4)
  expect_false(is.null(b$boot_coef))
  bb <- bayesian_backfitting(b, n_samples = 20, n_burnin = 5)
  expect_false(any(c("boot_coef", "boot_mstop", "widening_factor") %in%
                     names(bb)))
  expect_false(is.null(bb$coefficients_bayes))
  # Unless the time limit leaves no posterior: then the bootstrap stays
  expect_warning(
    same <- bayesian_backfitting(b, n_samples = 5L, n_burnin = 0L,
                                 max_runtime_mins = 0),
    "unchanged, without posterior results")
  expect_identical(same, b)
})

test_that("a time limit after some draws keeps them and counts them", {
  dat <- simulate_gaussian_1d(n = 80, seed = 24)
  mod <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl(max_iter = 100))
  # The sampler as if the time limit stopped it after 3 kept draws
  real_b2b <- boost_to_backfit
  local_mocked_bindings(boost_to_backfit = function(...) {
    post <- real_b2b(...)
    keep <- 3L
    post$beta_samples <- lapply(post$beta_samples,
                                function(b) b[, seq_len(keep), drop = FALSE])
    post$n_done <- list(...)$n_burnin + keep
    post
  })
  set.seed(3)
  expect_warning(bb <- bayesian_backfitting(mod, n_samples = 10L,
                                            n_burnin = 2L),
                 "time limit reached after 3 / 10")
  expect_identical(bb$n_samples_done, 3L)
  expect_identical(bb$uq_control$control$n_samples, 10L)
  expect_true(all(is.finite(bb$Vp)))
})

test_that("the time limit counts Stage 1, not only the sampling sweeps", {
  dat <- simulate_gaussian_1d(n = 60, seed = 25)
  mod <- gbgam(y ~ t2(x, k = 5, bs = "cr"), data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl("none", max_iter = 20L))
  # Stage 1 as if it took longer than the whole limit of 0.5 s
  real_map <- map_smoothing_parameters
  local_mocked_bindings(map_smoothing_parameters = function(...) {
    Sys.sleep(0.8)
    real_map(...)
  })
  expect_warning(
    out <- bayesian_backfitting(mod, n_samples = 1e5L, n_burnin = 0L,
                                max_runtime_mins = 0.5 / 60),
    "fewer than 2 posterior samples"
  )
  expect_identical(out, mod)
})

test_that("in_link_range() follows the family's valideta and validmu", {
  ig <- stats::inverse.gaussian()            # 1/mu^2: eta > 0
  expect_true(in_link_range(ig, c(0.5, 2)))
  expect_false(in_link_range(ig, c(0.5, -0.1)))
  ga <- stats::Gamma()                       # inverse: eta != 0, mu > 0
  expect_true(in_link_range(ga, c(0.5, 2)))
  expect_false(in_link_range(ga, c(0.5, -2)))
  expect_true(in_link_range(stats::binomial(), c(-50, 0, 50)))
})
