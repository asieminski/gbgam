# Stopping by resampling estimates of the prediction error
# (gb_stop_ctrl("cv"), gb_resampling()).


# ---- Resampling estimates ----------------------------------------------------

test_that("the risk of a resampling fit is the weighted out-of-sample loss", {
  dat <- simulate_poisson_1d(n = 100, seed = 3)
  max_iter <- 25
  set.seed(4)
  mod <- fit_1d(dat, stats::poisson(),
                gb_stop_ctrl("cv", max_iter = max_iter, iter_no_improv = 1000,
                             resampling = gb_resampling(B = 4, chunk = 10L)))
  W <- mod$cv$weights
  expect_equal(dim(W), c(100L, 4L))
  expect_true(all(rowSums(W == 0) == 1))
  for (b in c(1L, 3L)) {
    # The resampling fit is gbgam() with the training weights as weights
    fb <- fit_1d(dat, stats::poisson(),
                 gb_stop_ctrl("none", max_iter = max_iter), weights = W[, b])
    out <- W[, b] == 0
    mu  <- fb$fitted.values
    dev  <- stats::poisson()$dev.resids(dat$y[out], mu[out], 1)
    # row 1 of rep_risk is iteration 0
    expect_equal(mod$cv$rep_risk[max_iter + 1, b], sum(dev) / nrow(dat),
                 tolerance = 1e-8)
  }
  # k-fold: the estimate is the sum of the fold risks
  expect_equal(mod$cv$oob, rowSums(mod$cv$rep_risk))
  expect_equal(mod$cv$risk, mod$cv$oob)
  expect_equal(mod$mstop, which.min(mod$cv$risk) - 1)
  expect_equal(mod$metrics$ls, mod$cv$risk)
})

test_that("with prior weights, the risk weights the loss and divides by them", {
  dat <- simulate_poisson_1d(n = 100, seed = 3)
  set.seed(1)
  pw <- stats::runif(100, 0.5, 2)
  set.seed(4)
  mod <- fit_1d(dat, stats::poisson(),
                gb_stop_ctrl("cv", max_iter = 25, iter_no_improv = 1000,
                             resampling = gb_resampling(B = 4, chunk = 10L)),
                weights = pw)
  W <- mod$cv$weights
  for (b in c(1L, 3L)) {
    # the resampling fit trains with the prior times the training weights
    fb  <- fit_1d(dat, stats::poisson(), gb_stop_ctrl("none", max_iter = 25),
                  weights = pw * W[, b])
    out <- W[, b] == 0
    dev <- stats::poisson()$dev.resids(dat$y[out], fb$fitted.values[out],
                                       pw[out])
    # 2e-16 measured (relative); the unweighted loss is off by 16-36%
    expect_equal(mod$cv$rep_risk[26, b], sum(dev) / sum(pw), tolerance = 1e-8)
  }
})

test_that("leave-one-out bootstrap averages per observation", {
  W <- cbind(c(0, 2, 1, 0), c(1, 0, 2, 0), c(2, 1, 1, 0))
  ev <- oob_evaluation_weights(W, prior_weights = c(1, 1, 2, 1))
  # obs 1 out once, obs 2 out once, obs 3 never, obs 4 out three times
  expect_equal(ev$n_never_out, 1L)
  expect_equal(ev$evaluation_weights[4, ], rep(1 / 3, 3) / 3)
  expect_equal(ev$evaluation_weights[1, ], c(1, 0, 0) / 3)
  expect_equal(ev$evaluation_weights[3, ], c(0, 0, 0))
  # Losses l_ib: Err1 = sum_i (mean of i's OOB losses) / sum of weights
  L <- matrix(c(1, 5, 9, 2,  3, 6, 9, 4,  7, 8, 9, 6), 4)
  err1_manual <- (L[1, 1] + L[2, 2] + mean(L[4, ])) / 3
  expect_equal(sum(ev$evaluation_weights * L), err1_manual)
})

test_that("the no-information error matches the double sum of deviances", {
  set.seed(1)
  n  <- 30
  w  <- runif(n, 0.5, 2)
  # Efron and Tibshirani (1997, eq. 26), with prior weights
  brute <- function(fam, y, mu) {
    D <- outer(seq_len(n), seq_len(n),
               function(i, j) fam$dev.resids(y[i], mu[j], 1))
    sum(outer(w, w) * D) / sum(w)^2
  }
  cases <- list(
    gaussian = list(fam = stats::gaussian(), y = rnorm(n),
                    mu = rnorm(n)),
    poisson  = list(fam = stats::poisson(), y = c(0, rpois(n - 1, 3)),
                    mu = runif(n, 1, 5)),
    binomial = list(fam = stats::binomial(), y = c(0, 1, runif(n - 2)),
                    mu = runif(n, 0.1, 0.9)),
    Gamma    = list(fam = stats::Gamma("log"), y = rgamma(n, 2),
                    mu = runif(n, 0.5, 3)),
    inverse.gaussian = list(fam = stats::inverse.gaussian(),
                            y = runif(n, 0.5, 2), mu = runif(n, 0.5, 2)),
    quasipoisson = list(fam = stats::quasipoisson(), y = rpois(n, 3),
                        mu = runif(n, 1, 5)),
    negative.binomial = list(fam = MASS::negative.binomial(2),
                             y = rpois(n, 3), mu = runif(n, 1, 5)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    expect_equal(no_information_deviance(cs$y, cs$mu, w, cs$fam),
                 brute(cs$fam, cs$y, cs$mu), tolerance = 1e-10, label = nm)
  }
})

test_that("the no-information error is the mean over all permutations", {
  # All permutations of 1:n, one per row
  perms <- function(n) {
    if (n == 1L) return(matrix(1L))
    do.call(rbind, lapply(seq_len(n), function(i) {
      rest <- seq_len(n)[-i]
      cbind(i, matrix(rest[perms(n - 1L)], ncol = n - 1L))
    }))
  }
  P <- perms(6L)
  expect_identical(dim(P), c(720L, 6L))
  set.seed(4)
  cases <- list(
    gaussian = list(fam = stats::gaussian(), y = rnorm(6), mu = rnorm(6)),
    poisson  = list(fam = stats::poisson(), y = c(0, 0, 1, 3, 5, 2),
                    mu = runif(6, 0.5, 4)),
    binomial = list(fam = stats::binomial(), y = c(0, 1, 1, 0, 0.5, 0.25),
                    mu = runif(6, 0.1, 0.9)),
    Gamma    = list(fam = stats::Gamma("log"), y = rgamma(6, 2),
                    mu = runif(6, 0.5, 3)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    # Pair the responses with the fitted means in every possible order
    perm_mean <- mean(apply(P, 1L, function(p)
      mean(cs$fam$dev.resids(cs$y, cs$mu[p], 1))))
    expect_equal(no_information_deviance(cs$y, cs$mu, rep(1, 6), cs$fam),
                 perm_mean, tolerance = 1e-12, label = nm)
  }
})

test_that("the no-information error treats prior weights as frequencies", {
  set.seed(5)
  n  <- 8
  w  <- c(0, 1, 2, 3, 1, 0, 4, 2)
  y  <- rpois(n, 2)
  mu <- runif(n, 0.5, 4)
  fam <- stats::poisson()
  # Repeating row i w_i times gives the same value; zero weights drop out
  expect_equal(no_information_deviance(y, mu, w, fam),
               no_information_deviance(rep(y, w), rep(mu, w),
                                       rep(1, sum(w)), fam),
               tolerance = 1e-12)
  keep <- w > 0
  expect_equal(no_information_deviance(y, mu, w, fam),
               no_information_deviance(y[keep], mu[keep], w[keep], fam),
               tolerance = 1e-12)
})

test_that("the no-information error is exact at the edge cases", {
  pairwise <- function(fam, y, mu, w) {
    D <- outer(seq_along(y), seq_along(mu),
               function(i, j) fam$dev.resids(y[i], mu[j], 1))
    sum(outer(w, w) * D) / sum(w)^2
  }
  check <- function(fam, y, mu, w = rep(1, length(y)), label) {
    g <- no_information_deviance(y, mu, w, fam)
    expect_true(is.finite(g), label = label)
    expect_equal(g, pairwise(fam, y, mu, w), tolerance = 1e-10,
                 label = label)
    g
  }
  set.seed(6)
  mu <- runif(5, 0.5, 3)
  # ybar on the boundary of the mean space: all-zero counts, binomial
  # responses all 0 or all 1 (with trial weights)
  g <- check(stats::poisson(), rep(0, 5), mu, label = "poisson, all 0")
  expect_equal(g, 2 * mean(mu))
  p <- runif(5, 0.1, 0.9)
  check(stats::binomial(), rep(0, 5), p, w = c(3, 1, 5, 2, 4),
        label = "binomial, all 0")
  check(stats::binomial(), rep(1, 5), p, w = c(3, 1, 5, 2, 4),
        label = "binomial, all 1")
  # Constant response: only the d(ybar, mu_j) part remains
  fam <- stats::Gamma("log")
  g   <- check(fam, rep(2, 5), mu, label = "constant y")
  expect_equal(g, mean(fam$dev.resids(rep(2, 5), mu, 1)))
  # Fitted means all at ybar (intercept-only fit): the null deviance per
  # unit weight
  y <- rgamma(5, 2)
  g <- check(fam, y, rep(mean(y), 5), label = "mu = ybar")
  expect_equal(g, sum(fam$dev.resids(y, mean(y), 1)) / 5)
  # One observation: its own deviance
  check(fam, 1.5, 0.7, label = "n = 1")
  expect_equal(no_information_deviance(1.5, 0.7, 1, fam),
               fam$dev.resids(1.5, 0.7, 1))
  # Responses over six orders of magnitude
  check(fam, 10^(-3:2), runif(6, 0.01, 50), label = "wide range")
})

test_that("bootstrap stopping keeps each fitted mean's exposure in gamma", {
  set.seed(7)
  n   <- 100
  dat <- data.frame(x = runif(n), expo = runif(n, 0.5, 5))
  dat$y <- rpois(n, dat$expo * exp(0.5 + sin(2 * pi * dat$x)))
  mod <- suppressWarnings(gbgam(
    y ~ t2(x, k = 6, bs = "cr") + offset(log(expo)), family = poisson(),
    data = dat, verbose = FALSE,
    edf_control  = gb_edf_ctrl("none"),
    stop_control = gb_stop_ctrl(
      "cv", max_iter = 30,
      resampling = gb_resampling("bootstrap", B = 4, chunk = 30L))))
  # gamma at mstop from all pairs (y_i, mu_j), with mu_j = fitted(mod)[j]
  # at row j's own exposure
  D <- outer(dat$y, fitted(mod), function(yi, mj)
    poisson()$dev.resids(yi, mj, 1))
  expect_equal(mod$cv$noinfo[mod$mstop + 1L], mean(D), tolerance = 1e-10)
})

test_that(".632 and .632+ combine the errors as in Efron and Tibshirani", {
  err   <- c(1.0, 0.8, 0.6)
  err1  <- c(1.2, 1.1, 1.3)
  gamma <- c(1.5, 1.5, 1.25)
  expect_equal(err_632(err, err1), 0.368 * err + 0.632 * err1)
  ep    <- err_632plus(err, err1, gamma)
  err1p <- pmin(err1, gamma)
  R     <- (err1p - err) / (gamma - err)
  expect_equal(ep$relative_overfitting, R)
  # Efron and Tibshirani (1997, eq. 32); the last iteration has
  # Err(1) > gamma, where Err(.632) keeps the untruncated Err(1)
  expect_equal(ep$risk,
               err_632(err, err1) + (err1p - err) * 0.368 * 0.632 * R /
                 (1 - 0.368 * R))
  expect_equal(ep$risk[3], 0.632 * err1[3] + 0.368 * gamma[3])
  # Where Err(1) <= gamma, eq. 32 is the weighted form of eq. 29 (as in
  # pec and peperr)
  w <- 0.632 / (1 - 0.368 * R)
  expect_equal(ep$weight, w)
  expect_equal(ep$risk[1:2], ((1 - w) * err + w * err1)[1:2])
  # No overfitting (err1 <= err): plain .632 weight
  ep0 <- err_632plus(1, 0.9, 2)
  expect_equal(ep0$weight, 0.632)
  expect_equal(ep0$risk, err_632(1, 0.9))
  # The estimator chosen in gb_resampling()
  expect_identical(resampling_estimate("oob", err, err1, NULL)$risk, err1)
  expect_identical(resampling_estimate("632", err, err1, NULL)$risk,
                   err_632(err, err1))
  expect_identical(resampling_estimate("632plus", err, err1, gamma), ep)
})

test_that("bootstrap stopping gives the .632+ estimate and its components", {
  dat <- simulate_gaussian_1d(n = 120, seed = 1)
  set.seed(3)
  expect_warning(
    mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
      "cv", iter_no_improv = 40,
      resampling = gb_resampling("bootstrap", B = 6, chunk = 40L))),
    "not left out by any resampling fit")
  cv <- mod$cv
  m  <- mod$mstop
  expect_equal(cv$estimator, "632plus")
  expect_equal(length(cv$risk), length(cv$apparent))
  # Each piece at mstop, from the fits themselves (all exact up to rounding)
  y   <- dat$y
  err <- mean((y - fitted(mod))^2)  # apparent error
  expect_equal(cv$apparent[m + 1], err, tolerance = 1e-10)
  expect_equal(cv$apparent[m + 1], mod$deviance / 120, tolerance = 1e-10)
  gam <- mean(outer(y, fitted(mod), "-")^2)  # no-information rate
  expect_equal(cv$noinfo[m + 1], gam, tolerance = 1e-10)
  # leave-one-out bootstrap error: each observation's mean loss over the
  # fits that leave it out, averaged over the observations left out
  W    <- cv$weights
  out  <- W == 0
  loss <- vapply(seq_len(ncol(W)), function(b) {
    fb <- fit_1d(dat, stop_control = gb_stop_ctrl("none", max_iter = m),
                 weights = W[, b])
    (y - fitted(fb))^2
  }, numeric(120))
  left <- rowSums(out) > 0
  err1 <- mean(rowSums(loss * out)[left] / rowSums(out)[left])
  expect_equal(cv$oob[m + 1], err1, tolerance = 1e-10)
  # the .632+ rule of Efron and Tibshirani (1997), with Err(1) below gamma
  expect_lt(err1, gam)
  R <- (err1 - err) / (gam - err)
  w <- 0.632 / (1 - 0.368 * R)
  expect_equal(cv$risk[m + 1], (1 - w) * err + w * err1, tolerance = 1e-10)
  expect_true(all(colSums(cv$weights) == 120))
  expect_equal(mod$mstop, cv$mstop)
  expect_true(length(cv$risk) - 1 - mod$mstop >= 40)
})

test_that("custom weights and evaluation weights are used as given", {
  dat <- simulate_gaussian_1d(n = 60, seed = 2)
  W <- cbind(rep(c(1, 0), 30), rep(c(0, 1), 30))
  V <- 2 * (W == 0) / 60
  set.seed(1)
  a <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", max_iter = 30, resampling = gb_resampling(weights = W)))
  set.seed(1)
  b <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", max_iter = 30,
    resampling = gb_resampling(weights = W, eval_weights = V)))
  # evaluation weights twice the default ones double the risk
  expect_equal(b$cv$risk, 2 * a$cv$risk)
  expect_equal(a$cv$type, "custom")
  expect_equal(a$cv$weights, W)
})

test_that("user-supplied folds are respected", {
  dat <- simulate_gaussian_1d(n = 60, seed = 2)
  folds <- rep(1:3, each = 20)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", max_iter = 20, resampling = gb_resampling(folds = folds)))
  expect_equal(mod$cv$weights, outer(folds, 1:3, `!=`) * 1)
})

test_that("a fit that leaves no observation out adds nothing to the estimate", {
  dat <- simulate_binomial_1d(n = 60, seed = 4)
  W   <- cbind(rep(0:1, 30), rep(1:0, 30), 1)
  expect_warning(
    mod <- fit_1d(dat, stats::binomial(), gb_stop_ctrl(
      "cv", max_iter = 30, resampling = gb_resampling(weights = W,
                                                      chunk = 30L))),
    "1 of 3 resampling fits leave no observation out")
  expect_equal(unname(mod$cv$rep_risk[, 3L]), rep(0, 31))
  expect_true(all(is.finite(mod$cv$risk)))
})


# ---- Resampling fits ---------------------------------------------------------

test_that("CV stopping is reproducible and independent of the chunk size", {
  dat <- simulate_poisson_1d(n = 100, seed = 3)
  fit <- function(chunk) {
    set.seed(11)
    fit_1d(dat, stats::poisson(), gb_stop_ctrl(
      "cv", max_iter = 60, iter_no_improv = 1000,
      resampling = gb_resampling(B = 3, chunk = chunk)),
      edf = gb_edf_ctrl("hutchpp", z = 12))
  }
  a <- fit(15L)
  b <- fit(60L)
  expect_equal(a$cv$rep_risk, b$cv$rep_risk, tolerance = 1e-10)
  expect_equal(a$coefficients, b$coefficients, tolerance = 1e-10)
})

test_that("replicates can replay the main fit's selection and smoothing", {
  dat <- simulate_gaussian_1d(n = 80, seed = 1)
  # 300 iterations, so that ASR relaxes (157 of the main fit's smoothing
  # parameters are set): replaying them or not then makes a difference
  fit <- function(mode) {
    set.seed(2)
    fit_1d(dat, stop_control = gb_stop_ctrl(
      "cv", max_iter = 300, iter_no_improv = 1000,
      resampling = gb_resampling(B = 3, mode = mode, keep_fits = TRUE)))
  }
  mod <- fit("fixed")
  expect_equal(nrow(mod$cv$rep_risk), 301L)
  expect_true(all(is.finite(mod$cv$risk)))
  expect_identical(mod$resampling$mode, "fixed")
  expect_identical(resample_fit_replay("fixed", mod),
                   list(selection = TRUE, smoothing = TRUE,
                        edf_compute_method = "none"))
  expect_true(any(!is.na(mod$metrics$rhos)))
  mod2 <- fit("free_smoothing")
  expect_true(all(is.finite(mod2$cv$risk)))
  expect_identical(resample_fit_replay("free_smoothing", mod2),
                   list(selection = TRUE, smoothing = FALSE,
                        edf_compute_method = "exact"))
  # fixed: every resampling fit replays the main fit's learners and
  # smoothing parameters; free_smoothing: the learners only
  for (m in list(fixed = mod, free_smoothing = mod2)) {
    expect_length(m$resume$cv$fits, 3L)
    for (st in m$resume$cv$fits) {
      expect_identical(st$k_opts, m$metrics$k_opts)
      expect_identical(identical(st$rhos, m$metrics$rhos),
                       m$resampling$mode == "fixed")
    }
  }
})

test_that("cv stopping stops at the time limit", {
  dat <- simulate_gaussian_1d(n = 80, seed = 1)
  set.seed(2)
  # the limit is passed by the first chunk of 10 iterations, of 400 allowed
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", max_iter = 400, max_runtime_mins = 1e-9,
    resampling = gb_resampling(B = 3, chunk = 10L)))
  expect_true(mod$metrics$stopping_conditions[["time_limit_exceeded"]])
  expect_false(mod$metrics$stopping_conditions[["mstop_reached"]])
  expect_identical(nrow(mod$cv$rep_risk), 11L)
})

test_that("a factor level without training weight does not break a replicate", {
  set.seed(1)
  n <- 120
  dat <- data.frame(x = runif(n),
                    f = factor(c(rep("a", 58), rep("b", 58), rep("c", 4))))
  dat$y <- sin(2 * pi * dat$x) + 3 * (dat$f == "c") + rnorm(n, sd = 0.3)
  folds <- c(rep(1:3, length.out = 116), rep(1, 4))
  for (md in c("free_smoothing", "free_selection", "fixed")) {
    mod <- gbgam(y ~ f + t2(x, k = 8, bs = "cr"), data = dat,
                 stop_control = gb_stop_ctrl(
                   "cv", max_iter = 60,
                   resampling = gb_resampling(folds = folds, mode = md)),
                 verbose = FALSE)
    expect_true(all(is.finite(mod$cv$risk)))
  }
  # The zero-norm learner has a zero fit and a zero EDF increment
  lr <- new_unpenalised_learner(c(1, 1, 0), "z", weights = c(0, 0, 1))
  expect_equal(fit_learner(lr, c(1, 2, 3))$b, 0)
  expect_equal(as.vector(delta_dgamma_dy(lr, diag(3))),
               rep(0, 3))
})

test_that("a fixed sp valid for the main fit does not stop the replicates", {
  # the replicates clamp it (clamp_sp = TRUE) where their weights would give
  # a selection EDF below 1
  dat <- simulate_gaussian_1d(n = 90, seed = 4)
  base <- gbgam(y ~ t2(x, k = 10, bs = "cr"), data = dat,
                stop_control = gb_stop_ctrl("none", max_iter = 2),
                verbose = FALSE)
  rho0 <- dr_learners(base)[[1]]$rho0
  sp <- exp(rho0 - 0.2)
  ff <- y ~ t2(x, k = 10, bs = "cr", sp = sp)
  expect_silent(gbgam(ff, data = dat, verbose = FALSE,
                      stop_control = gb_stop_ctrl("none", max_iter = 5)))
  mod <- gbgam(ff, data = dat, verbose = FALSE,
               stop_control = gb_stop_ctrl(
                 "cv", max_iter = 30, resampling = gb_resampling(B = 3)))
  expect_true(all(is.finite(mod$cv$risk)))
})

test_that("rebuilt resample fits redraw their probes; the caller's RNG is kept", {
  dat <- simulate_poisson_1d(n = 80, seed = 3)
  W   <- cbind(rep(0:1, 40), rep(1:0, 40))
  set.seed(3)
  a <- fit_1d(dat, stats::poisson(), gb_stop_ctrl(
    "cv", max_iter = 30,
    resampling = gb_resampling(weights = W, keep_fits = TRUE)),
    edf = gb_edf_ctrl("hutchpp", z = 12))
  plan <- plan_resampling(a, inbag_weights = a$cv$weights,
                          rng_state = a$resume$cv$rng_state)
  # Rebuilt in a session with another RNG kind and state
  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[1], old[2], old[3]))
  set.seed(5)
  kind_before  <- RNGkind()
  state_before <- .Random.seed
  fits <- with_rng_state(plan$rng_state, build_resample_fits(plan))
  expect_identical(RNGkind(), kind_before)
  expect_identical(.Random.seed, state_before)
  # The same probes as the stored fits, which drew them consecutively
  for (b in 1:2) {
    expect_identical(fits[[b]]$mod$derivs$GGG,
                     a$resume$cv$fits[[b]]$mod$derivs$GGG)
  }
  expect_false(identical(fits[[1]]$mod$derivs$GGG,
                         fits[[2]]$mod$derivs$GGG))
})

test_that("an extension rebuilds the resample fits under the session's RNG kind", {
  dat <- simulate_poisson_1d(n = 80, seed = 3)
  old <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old[1], old[2], old[3]))
  fit <- function(max_iter) {
    set.seed(8)
    fit_1d(dat, stats::poisson(), gb_stop_ctrl(
      "cv", max_iter = max_iter, iter_no_improv = 1000,
      resampling = gb_resampling(B = 3, chunk = 10L)),
      edf = gb_edf_ctrl("hutchpp", z = 12))
  }
  long <- fit(40)
  ext  <- extend_gbgam(fit(20), n_iter = 20)
  expect_equal(ext$cv$rep_risk, long$cv$rep_risk, tolerance = 1e-10)
  expect_equal(ext$cv$risk, long$cv$risk, tolerance = 1e-10)
})

test_that("a model without a recorded RNG state cannot be extended", {
  dat <- simulate_gaussian_1d(n = 60, seed = 2)
  set.seed(1)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", max_iter = 10, resampling = gb_resampling(B = 3)))
  mod$resume$cv$rng_state <- NULL
  expect_error(extend_gbgam(mod, n_iter = 5), "older gbgam")
})


# ---- Argument checks and use of a cv model -----------------------------------

test_that("gb_resampling() validates its arguments", {
  expect_s3_class(gb_resampling(), "gb_resampling")
  expect_equal(gb_resampling()$B, 10L)
  expect_equal(gb_resampling("bootstrap")$B, 25L)
  expect_equal(gb_resampling("bootstrap")$estimator, "632plus")
  expect_equal(gb_resampling("subsample")$estimator, "oob")
  expect_error(gb_resampling(B = 1), "at least 2")
  expect_error(gb_resampling(estimator = "632"), "only available")
  expect_silent(gb_resampling("bootstrap", estimator = "632"))
  expect_error(gb_resampling(estimator = "mse"), "estimator")
  # The deviance is the only loss
  expect_error(gb_resampling(loss = "mse"), "unused argument")
  expect_identical(gb_resampling()$mode, "free_selection")
  for (md in c("free_selection", "free_smoothing", "fixed"))
    expect_identical(gb_resampling(mode = md)$mode, md)
  expect_error(gb_resampling(mode = "free_stopping"), "should be one of")
  expect_error(gb_resampling(mode = character()), "length 1")
  # Arguments that are not used are errors
  expect_error(gb_resampling(B = 3, folds = rep(1:2, 5)), "B is not used")
  expect_error(gb_resampling(B = 3, weights = matrix(1, 4, 2)),
               "B is not used")
  expect_error(gb_resampling("bootstrap", weights = matrix(1, 4, 2)),
               "type is not used")
  expect_error(gb_resampling(fraction = 0.6), "only used with type")
  expect_identical(gb_resampling("subsample")$fraction, 0.5)
  expect_null(gb_resampling()$fraction)
  expect_identical(gb_resampling(weights = matrix(1, 4, 2))$type, "custom")
  expect_error(gb_resampling(eval_weights = diag(2)), "together with weights")
  expect_error(gb_resampling(weights = matrix(-1, 3, 2)), "non-negative")
  expect_error(gb_resampling("subsample", fraction = 1), "strictly")
  expect_error(gb_stop_ctrl(resampling = gb_resampling()), "only used")
  expect_error(gb_stop_ctrl("cv", resampling = list()), "gb_resampling")
  expect_s3_class(gb_stop_ctrl("cv")$resampling, "gb_resampling")
})

test_that("cv input errors name the problem", {
  dat <- simulate_gaussian_1d(n = 60, seed = 2)
  expect_error(fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", resampling = gb_resampling(folds = rep(1:2, 10)))),
    "one label per observation")
  expect_error(fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", resampling = gb_resampling(weights = matrix(1, 10, 2)))),
    "one row per observation")
  expect_error(fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", resampling = gb_resampling(weights = matrix(1, 60, 2)))),
    "No observation is left out")
})

test_that("k-fold cross-validation needs at most one fold per observation", {
  dat <- simulate_binomial_1d(n = 20, seed = 4)
  expect_error(fit_1d(dat, stats::binomial(), gb_stop_ctrl(
    "cv", max_iter = 10, resampling = gb_resampling(B = 21))),
    "must not exceed the number of observations used in the fit \\(20\\)")
})

test_that("cv stopping with bootstrap free stopping fails before fitting", {
  dat <- simulate_gaussian_1d(n = 60, seed = 3)
  expect_error(
    fit_1d(dat, stop_control = gb_stop_ctrl("cv"),
           uq_control = gb_uq_ctrl("bootstrap",
                                   mode = "free_stopping")),
    "not available with stopping_method = \"cv\"")
})

test_that("a cv model works with summary, print, predict and the bootstrap", {
  dat <- simulate_gaussian_1d(n = 80, seed = 1)
  set.seed(2)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", iter_no_improv = 20, resampling = gb_resampling(B = 3)))
  expect_model_contract(mod, 80)
  s <- summary(mod)
  expect_equal(s$loss_name, "CV risk")
  expect_equal(s$best_loss, min(mod$cv$risk))
  expect_output(print(mod))
  expect_equal(as.vector(predict(mod, newdata = dat)),
               as.vector(mod$linear.predictors), tolerance = 1e-10)
  # Three replicates give a noisy widening factor, which may exceed 3
  b <- withCallingHandlers(
    bootstrap_gbgam(mod, n_samples = 3, mode = "free_selection"),
    warning = function(w) {
      if (grepl("widening factor", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  expect_false(is.null(b$Vp))
  expect_error(bootstrap_gbgam(mod, n_samples = 3, mode = "free_stopping"),
               "not available")
})

test_that("CV stopping checks the relative patience every chunk", {
  dat <- simulate_gaussian_1d(n = 100, seed = 1)
  chunk <- 4L
  set.seed(2)
  mod <- fit_1d(dat, stop_control = gb_stop_ctrl(
    "cv", iter_no_improv = 3, iter_no_improv_rel = 0.5,
    resampling = gb_resampling(B = 3, chunk = chunk)))
  M      <- sum(!is.na(mod$metrics$k_opts))
  m_star <- mod$mstop
  cat("cv: mstop", m_star, "run", M, "\n")
  # The rule at a check after `run` iterations (ties to the earlier one)
  stops_at <- function(run) {
    best <- which.min(mod$cv$risk[seq_len(run + 1L)]) - 1L
    run - best >= max(3, ceiling(0.5 * best))
  }
  expect_true(mod$converged)
  expect_identical(M %% chunk, 0L)
  expect_gt(ceiling(0.5 * m_star), 3)
  expect_true(stops_at(M))
  expect_false(any(vapply(seq(chunk, M - chunk, by = chunk), stops_at,
                          logical(1))))
})
