# Simulation of the bootstrap responses: prior weights, zero weights,
# binomial trials and the dispersion.


# ---- Simulated responses -----------------------------------------------------

test_that("Poisson replicates have variance mu / w for prior weight w", {
  fam <- stats::poisson()
  n   <- 20000L
  mu  <- rep(2, n)
  w   <- rep(c(1, 3, 0.5), length.out = n)
  set.seed(1)
  y   <- simulate_response(fam, mu, prior_weights = w, user_weights = w,
                           phi = 1)
  for (wj in unique(w)) {
    yj <- y[w == wj]
    expect_equal(mean(yj), 2, tolerance = 0.05)
    expect_equal(stats::var(yj), 2 / wj, tolerance = 0.1)
  }
  # Unit weights: unchanged draws
  set.seed(2)
  y1 <- simulate_response(fam, mu, prior_weights = rep(1, n),
                          user_weights = rep(1, n), phi = 1)
  set.seed(2)
  expect_identical(y1, as.numeric(stats::rpois(n, mu)))
})

test_that("Gaussian, Gamma and IG replicates have variance V(mu) phi / w", {
  # 1e5 draws per weight. Relative Monte Carlo SE of the mean at most 0.005
  # and of the variance at most 0.005 (Gaussian), 0.009 (Gamma) and 0.02
  # (inverse Gaussian, heavy-tailed); over 4 seeds the errors stay within
  # 2.3 SE (largest: 0.024, inverse Gaussian variance). Tolerances are 4 SE
  # or more; ignoring the weights (Gaussian) or phi (Gamma shape w phi, inverse
  # Gaussian without phi) changes a variance by a factor of 2 or more.
  n   <- 2e5L
  w   <- rep(c(0.5, 4), length.out = n)
  phi <- 0.5
  specs <- list(
    list(fam = stats::gaussian(),              mu = 1, V = 1, tol = 0.02),
    list(fam = stats::Gamma("log"),            mu = 2, V = 4, tol = 0.04),
    list(fam = stats::inverse.gaussian("log"), mu = 2, V = 8, tol = 0.08))
  for (sp in specs) {
    set.seed(1)
    # user weights unlike the prior weights: the variance follows the
    # prior weights only
    y <- simulate_response(sp$fam, rep(sp$mu, n), prior_weights = w,
                           user_weights = rep(3, n), phi = phi)
    for (wj in c(0.5, 4)) {
      yj  <- y[w == wj]
      lab <- paste(sp$fam$family, "w =", wj)
      expect_equal(mean(yj), sp$mu, tolerance = 0.02, label = lab)
      expect_equal(stats::var(yj), sp$V * phi / wj, tolerance = sp$tol,
                   label = lab)
    }
  }
})

test_that("inverse Gaussian draws follow the SuppDists distribution", {
  skip_if_not_installed("SuppDists")
  # Kolmogorov-Smirnov test against SuppDists::pinvGauss(), the distribution
  # that stats::inverse.gaussian()$simulate draws from with lambda = w / phi;
  # mu / lambda up to 1e8, where the textbook smaller root cancels to zero
  n     <- 20000L
  specs <- list(c(1, 1), c(2, 0.5), c(0.3, 20), c(10, 0.01), c(1, 1e-8))
  for (sp in specs) {
    set.seed(1)
    x   <- r_inverse_gaussian(n, rep(sp[1], n), rep(sp[2], n))
    lab <- sprintf("mu = %g, lambda = %g", sp[1], sp[2])
    expect_true(all(is.finite(x) & x > 0), label = lab)
    p <- suppressWarnings(stats::ks.test(x, SuppDists::pinvGauss,
                                         nu = sp[1], lambda = sp[2])$p.value)
    expect_gt(p, 0.01, label = lab)
  }
})

test_that("weights = 3 and triplicated rows give the same Poisson SEs", {
  dat  <- simulate_poisson_1d(n = 100L, seed = 81L)
  ctrl <- list(
    family          = stats::poisson(),
    learner_control = gb_learner_ctrl(step_size = 0.1),
    stop_control    = gb_stop_ctrl("none", max_iter = 40L),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  dat3 <- dat[rep(seq_len(nrow(dat)), 3L), ]
  m_w  <- do.call(gbgam, c(list(y ~ t2(x, k = 6), data = dat,
                                weights = rep(3, nrow(dat))), ctrl))
  m_3  <- do.call(gbgam, c(list(y ~ t2(x, k = 6), data = dat3), ctrl))
  expect_equal(coef(m_w), coef(m_3), tolerance = 1e-8)

  set.seed(5)
  b_w <- bootstrap_gbgam(m_w, n_samples = 100L, alpha = NULL)
  set.seed(5)
  b_3 <- bootstrap_gbgam(m_3, n_samples = 100L, alpha = NULL)
  x_new <- data.frame(x = seq(0.05, 0.95, length.out = 25))
  se_w  <- predict(b_w, x_new, se.fit = TRUE)$se.fit
  se_3  <- predict(b_3, x_new, se.fit = TRUE)$se.fit
  ratio <- mean(se_w) / mean(se_3)
  expect_gt(ratio, 0.85)
  expect_lt(ratio, 1.15)
})

test_that("zero-weight rows are simulated at the fitted mean", {
  mu <- seq(0.5, 2, length.out = 10)
  w  <- c(0, 0, rep(1, 8))
  fams <- list(stats::gaussian(), stats::poisson(), stats::binomial(),
               stats::Gamma(link = "log"),
               stats::inverse.gaussian(link = "log"))
  for (fam in fams) {
    mu_f <- if (fam$family == "binomial") mu / 3 else mu
    set.seed(3)
    y <- simulate_response(fam, mu_f, prior_weights = w, user_weights = w,
                           phi = 0.5)
    expect_identical(y[1:2], mu_f[1:2], label = fam$family)
    expect_true(all(is.finite(y)), label = fam$family)
    if (fam$family %in% c("Gamma", "inverse.gaussian")) {
      expect_true(all(y > 0), label = fam$family)
    }
  }
})

test_that("binomial proportion + weights simulates Binomial(n_trials)", {
  fam <- binomial()
  mu  <- rep(0.3, 200)
  nt  <- rep(c(10, 25), 100)
  set.seed(2)
  y_prop <- simulate_response(fam, mu, prior_weights = nt, user_weights = nt,
                              phi = 1)
  expect_false(all(y_prop %in% c(0, 1)))
  expect_equal(y_prop * nt, round(y_prop * nt))
  # Same draws as the cbind specification (prior = n_trials, user = 1)
  set.seed(2)
  y_cbind <- simulate_response(fam, mu, prior_weights = nt,
                               user_weights = rep(1, 200), phi = 1)
  expect_identical(y_prop, y_cbind)
  # Binary data with unit weights stay Bernoulli
  set.seed(2)
  y_bin <- simulate_response(fam, mu, prior_weights = rep(1, 200),
                             user_weights = rep(1, 200), phi = 1)
  expect_true(all(y_bin %in% c(0, 1)))
  # Fractional weights fall back to the prior/user ratio (Bernoulli here)
  set.seed(2)
  y_frac <- simulate_response(fam, mu, prior_weights = rep(0.5, 200),
                              user_weights = rep(0.5, 200), phi = 1)
  expect_true(all(y_frac %in% c(0, 1)))
})

test_that("cbind and proportion + weights give comparable bootstrap SEs", {
  agg  <- simulate_binomial_formats(n = 80L, seed = 41L)$agg
  ctrl <- list(
    learner_control = gb_learner_ctrl(step_size = 0.1),
    stop_control    = gb_stop_ctrl("none", max_iter = 40L),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
  m_cb <- do.call(gbgam, c(list(cbind(successes, failures) ~ t2(x, k = 6),
                                data = agg, family = binomial()), ctrl))
  m_pr <- do.call(gbgam, c(list(prop ~ t2(x, k = 6), data = agg,
                                weights = agg$n_trials,
                                family = binomial()), ctrl))
  expect_equal(m_pr$linear.predictors, m_cb$linear.predictors,
               tolerance = 1e-8)
  set.seed(3); b_cb <- bootstrap_gbgam(m_cb, n_samples = 20L, alpha = NULL)
  set.seed(3); b_pr <- bootstrap_gbgam(m_pr, n_samples = 20L, alpha = NULL)
  se_cb <- sqrt(diag(b_cb$Vp))
  se_pr <- sqrt(diag(b_pr$Vp))
  used  <- se_cb > 1e-10  # coefficients of never-selected learners have SE 0
  expect_true(any(used))
  ratio <- stats::median(se_pr[used] / se_cb[used])
  expect_gt(ratio, 0.8)
  expect_lt(ratio, 1.25)
})

test_that("phi sets the dispersion of the bootstrap replicates", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  set.seed(9)
  tiny <- bootstrap_gbgam(mod, n_samples = 10L, phi = 1e-10, alpha = NULL)
  set.seed(9)
  fitted_phi <- bootstrap_gbgam(mod, n_samples = 10L, alpha = NULL)
  expect_lt(max(diag(tiny$Vp)), 1e-6 * max(diag(fitted_phi$Vp)))
  expect_identical(tiny$uq_control$control$phi, 1e-10)
})
