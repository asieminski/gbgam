# The starting intercept (intercept-only MLE or zero) and its
# post-iteration refit (learners_to_fit_after_each_iter = 1L).

# Fit with the intercept excluded from selection, optionally refitted
# after each iteration
fit_intercept_mode <- function(dat, family, refit_after, weights = NULL,
                               edf_method = "none") {
  ctrl <- list(step_size = 0.1, excluded_learner_indices = 1L)
  if (refit_after) ctrl$learners_to_fit_after_each_iter <- 1L
  gbgam(
    y ~ t2(x, k = 8, bs = "cr"),
    data            = dat,
    family          = family,
    weights         = weights,
    learner_control = do.call(gb_learner_ctrl, ctrl),
    stop_control    = gb_stop_ctrl("none", max_iter = 40),
    edf_control     = gb_edf_ctrl(edf_method),
    uq_control      = gb_uq_ctrl("none"),
    verbose         = FALSE
  )
}


# ---- Starting intercept ------------------------------------------------------

test_that("find_intercept_mle matches glm for links with a pole at eta = 0", {
  set.seed(201)
  n <- 120
  w <- stats::runif(n, 0.5, 2)
  cases <- list(
    list(fam = stats::Gamma(), y = function(m) stats::rgamma(n, 4, 4 / m)),
    list(fam = stats::inverse.gaussian(),
         y = function(m) stats::rgamma(n, 6, 6 / m)),
    list(fam = stats::poisson(link = "sqrt"),
         y = function(m) stats::rpois(n, 5 * m))
  )
  for (cs in cases) {
    # weighted mean of y below and above 1 (the old bracket was [-1, 1])
    for (m in c(0.4, 3)) {
      y   <- cs$y(m)
      lab <- paste(cs$fam$family, cs$fam$link, m)
      a   <- find_intercept_mle(y, new_family(cs$fam), w, rep(0, n))
      ref <- unname(stats::coef(stats::glm(y ~ 1, family = cs$fam,
                                           weights = w)))
      expect_equal(a, ref, tolerance = 1e-6, label = lab)
      # with a small offset, compare with a tightly converged glm
      off <- 0.001 * stats::rnorm(n)
      a_o <- find_intercept_mle(y, new_family(cs$fam), w, off)
      ref_o <- unname(stats::coef(suppressWarnings(stats::glm(
        y ~ 1, family = cs$fam, weights = w, offset = off,
        control = list(epsilon = 1e-14, maxit = 100)))))
      expect_equal(a_o, ref_o, tolerance = 1e-6, label = paste(lab, "offset"))
    }
  }
  # Zero-weight rows are ignored
  y <- c(stats::rgamma(n, 4, 4), NaN)
  expect_equal(find_intercept_mle(y, new_family(stats::Gamma()),
                                  c(w, 0), rep(0, n + 1)),
               find_intercept_mle(y[1:n], new_family(stats::Gamma()),
                                  w, rep(0, n)))
  # An invalid start gives a clear error
  expect_error(find_intercept_mle(rep(0, 5), new_family(stats::poisson()),
                                  rep(1, 5), rep(0, 5)),
               "valid starting intercept")
})

test_that("find_intercept_mle matches an offset-only glm from a far start", {
  # glm(y ~ 1, offset = ) with offsets over orders of magnitude and skewed
  # weights, so that the link of the weighted mean is far from the MLE. The
  # glm starts at the true intercept: its own start fails for some links.
  set.seed(301)
  n   <- 300
  w   <- stats::rexp(n) + 0.05
  off <- stats::rnorm(n, 0, 2.5)
  pos <- stats::rlnorm(n, 0, 1.2)  # for links with eta > 0
  rig <- function(mu) r_inverse_gaussian(n, mu, rep(2, n))
  cases <- list(
    list(fam = stats::poisson(), a = 1.5, off = off,
         y = function(eta) stats::rpois(n, exp(eta))),
    list(fam = stats::poisson("sqrt"), a = 0.7, off = pos,
         y = function(eta) stats::rpois(n, eta^2)),
    list(fam = stats::Gamma("log"), a = -1, off = off,
         y = function(eta) stats::rgamma(n, 1.5, 1.5 / exp(eta))),
    list(fam = stats::Gamma(), a = 0.3, off = pos,
         y = function(eta) stats::rgamma(n, 2, 2 * eta)),
    list(fam = stats::Gamma("identity"), a = 2, off = pos,
         y = function(eta) stats::rgamma(n, 2, 2 / eta)),
    list(fam = stats::inverse.gaussian("log"), a = 0.5, off = off,
         y = function(eta) rig(exp(eta))),
    list(fam = stats::inverse.gaussian(), a = 0.2, off = pos,
         y = function(eta) rig(sqrt(1 / eta))),
    list(fam = stats::binomial(), a = -1, off = off,
         y = function(eta) stats::rbinom(n, 1, stats::plogis(eta))),
    list(fam = stats::binomial("probit"), a = 0.4, off = off / 2,
         y = function(eta) stats::rbinom(n, 1, stats::pnorm(eta))),
    list(fam = stats::binomial("cloglog"), a = -0.6, off = off / 2,
         y = function(eta) stats::rbinom(n, 1, 1 - exp(-exp(eta)))),
    list(fam = stats::gaussian("log"), a = 0.5, off = off / 2,
         y = function(eta) stats::rnorm(n, exp(eta), 0.5))
  )
  for (cs in cases) {
    lab <- paste(cs$fam$family, cs$fam$link)
    y   <- cs$y(cs$a + cs$off)
    a   <- find_intercept_mle(y, new_family(cs$fam), w, cs$off)
    ref <- suppressWarnings(stats::glm(
      y ~ 1, family = cs$fam, weights = w, offset = cs$off, start = cs$a,
      control = list(epsilon = 1e-14, maxit = 500)))
    expect_true(ref$converged, label = lab)
    expect_equal(a, unname(stats::coef(ref)), tolerance = 1e-7, label = lab)
    # The solver had to move: the start is not the MLE
    start <- cs$fam$linkfun(sum(w * y) / sum(w))
    expect_gt(abs(start - a), 0.05, label = lab)
  }
})

test_that("find_intercept_mle reaches the MLE where full steps fail", {
  # Poisson with the identity link: undamped steps cycle around a = 1
  a <- find_intercept_mle(c(4, 10, 0, 0, 0),
                          new_family(stats::poisson("identity")),
                          rep(1, 5), c(0, 9, 9, 9, 9))
  expect_equal(a, 1, tolerance = 1e-8)
  # Log link with tiny exposures: the start is far below the closed form and
  # full steps overshoot
  set.seed(302)
  off <- stats::rnorm(50, -8, 1)
  y   <- stats::rpois(50, 20)
  a   <- find_intercept_mle(y, new_family(stats::poisson()), rep(1, 50), off)
  expect_equal(a, log(sum(y) / sum(exp(off))), tolerance = 1e-8)
  # gaussian(inverse): a full step jumps over the pole of a row and runs away
  y   <- c(2, -1.5)
  off <- c(0, 1)
  fam <- new_family(stats::gaussian("inverse"))
  a   <- find_intercept_mle(y, fam, c(1, 1), off)
  dev <- function(a) sum((y - 1 / (off + a))^2)
  ref <- stats::optimize(dev, c(0.2, 1), tol = 1e-9)
  expect_equal(a, ref$minimum, tolerance = 1e-6)
  # inverse.gaussian(inverse): the deviance is quadratic in a, so Newton's
  # method lands on the closed form, where Fisher steps crawl near the pole
  for (seed in c(129, 178)) {
    set.seed(seed)
    n   <- 30
    off <- stats::rlnorm(n, -0.5, 2)
    y   <- r_inverse_gaussian(n, 1 / (0.3 + off), rep(2, n))
    w   <- stats::rexp(n) + 0.02
    a   <- find_intercept_mle(y, new_family(stats::inverse.gaussian("inverse")),
                              w, off)
    expect_equal(a, (sum(w) - sum(w * y * off)) / sum(w * y),
                 tolerance = 1e-10)
  }
})

test_that("find_intercept_mle returns an optimum on the link's boundary", {
  # The row with y = 0 pushes its mean to 0, the end of the valid range
  a <- find_intercept_mle(c(0, 1), new_family(stats::poisson("identity")),
                          c(1, 1), c(0, 10))
  expect_gt(a, 0)
  expect_lt(a, 1e-6)
  # With the sqrt link only a tiny part of the Fisher step stays valid there
  a <- find_intercept_mle(c(0, 1), new_family(stats::poisson("sqrt")),
                          c(1, 1), c(0, 10))
  expect_gt(a, 0)
  expect_lt(a, 1e-6)
})

test_that("find_intercept_mle does not take a halved step for convergence", {
  # y = 0 is outside the support, so the deviance is not finite: the halved
  # steps get tiny while the Fisher step stays large
  expect_error(
    find_intercept_mle(c(0, 1, 2, 4),
                       new_family(stats::inverse.gaussian("log")), rep(1, 4),
                       c(0, 0.5, 1, 1.5)),
    "approximately zero total score")
})

test_that("find_intercept_mle converges at a perfect fit", {
  # The deviance is at its floating-point floor
  off <- seq(-2, 2, length.out = 20)
  a   <- find_intercept_mle(exp(0.3 + off), new_family(stats::Gamma("log")),
                            rep(1, 20), off)
  expect_equal(a, 0.3, tolerance = 1e-10)
  a   <- find_intercept_mle(exp(0.3 + off), new_family(stats::gaussian("log")),
                            rep(1, 20), off)
  expect_equal(a, 0.3, tolerance = 1e-10)
})

test_that("initialise_intercept = FALSE starts the intercept at zero", {
  set.seed(104)
  n <- 80
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rpois(n, exp(1 + sin(2 * pi * d$x)))
  fit_ii <- function(ii) {
    quiet_fit(y ~ t2(x, k = 8, bs = "cr"), data = d, family = stats::poisson(),
              learner_control = gb_learner_ctrl(step_size = 0.1, tau = 0,
                                                initialise_intercept = ii),
              stop_control = gb_stop_ctrl("none", max_iter = 20L),
              edf_control = gb_edf_ctrl("exact"))
  }
  m_true  <- fit_ii(TRUE)
  m_false <- fit_ii(FALSE)
  # intercept0 stays the MLE (used for the null deviance) ...
  expect_equal(m_false$intercept0, m_true$intercept0)
  expect_equal(m_false$null.deviance, m_true$null.deviance)
  # ... but the boosting start is 0, so the fits differ
  expect_equal(m_false$edfs[1], 0)
  expect_equal(m_true$edfs[1], 1)
  expect_gt(max(abs(m_false$linear.predictors - m_true$linear.predictors)),
            1e-3)
  # Running EDF (used for AIC) and the reported EDF agree
  for (m in list(m_true, m_false)) {
    expect_equal(m$edfs[m$mstop + 1L], m$edf.total, tolerance = 1e-8)
  }
  # On a path that never selects the intercept, its coefficient stays at the
  # starting value: 0 under FALSE, the MLE under TRUE
  fit_path <- function(ii) {
    quiet_fit(y ~ t2(x, k = 8, bs = "cr"), data = d, family = stats::poisson(),
              learner_control = gb_learner_ctrl(
                step_size = 0.1, tau = 0, initialise_intercept = ii,
                preselected_learner_indices = rep(3L, 5)),
              stop_control = gb_stop_ctrl("none", max_iter = 5L),
              edf_control = gb_edf_ctrl("exact"))
  }
  expect_equal(fit_path(FALSE)$learners[[1]]$beta, 0)
  expect_equal(fit_path(TRUE)$learners[[1]]$beta, m_true$intercept0)
})

test_that("initialise_intercept = FALSE errors clearly for inverse links", {
  set.seed(203)
  d <- data.frame(x = stats::runif(80))
  d$y <- stats::rgamma(80, shape = 5, rate = 5 / 3)
  for (fam in list(stats::Gamma(), stats::inverse.gaussian())) {
    expect_error(
      quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = d, family = fam,
                learner_control = gb_learner_ctrl(
                  initialise_intercept = FALSE)),
      "initialise_intercept = FALSE"
    )
  }
})

test_that("summary() reports the fitted intercept and the intercept-only MLE", {
  set.seed(209)
  d <- data.frame(x = stats::runif(80))
  d$y <- 2 + sin(2 * pi * d$x) + stats::rnorm(80, sd = 0.3)
  mod <- quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = d,
                   learner_control = gb_learner_ctrl(
                     step_size = 0.1, initialise_intercept = FALSE),
                   stop_control = sc_ff, edf_control = ec_ff)
  s <- summary(mod)
  expect_equal(s$intercept, unname(coef(mod)["(Intercept)"]))
  expect_equal(s$intercept0, mean(d$y))
  expect_false(isTRUE(all.equal(s$intercept, s$intercept0)))
  out <- cli::cli_fmt(print(s))
  expect_true(any(grepl("intercept-only MLE", out, fixed = TRUE)))
})


# ---- Post-iteration intercept refit ------------------------------------------

test_that("Gaussian post-iteration intercept refit is a no-op (as documented)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 101)
  dat$y <- dat$y + 3

  set.seed(102)
  w <- stats::runif(nrow(dat), 0.5, 2)

  for (wt in list(NULL, w)) {
    mod_after <- fit_intercept_mode(dat, stats::gaussian(), TRUE,
                                    weights = wt, edf_method = "exact")
    mod_never <- fit_intercept_mode(dat, stats::gaussian(), FALSE,
                                    weights = wt, edf_method = "exact")

    expect_equal(mod_after$linear.predictors, mod_never$linear.predictors,
                 tolerance = 1e-10)
    expect_equal(mod_after$learners[[1]]$beta, mod_never$learners[[1]]$beta,
                 tolerance = 1e-10)
    expect_equal(mod_after$edf.total, mod_never$edf.total, tolerance = 1e-8)
  }
})

test_that("Poisson post-iteration intercept refit changes the fit (not padding)", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 100, seed = 103)

  mod_after <- fit_intercept_mode(dat, stats::poisson(), TRUE)
  mod_never <- fit_intercept_mode(dat, stats::poisson(), FALSE)

  # For a non-identity link the gradient re-evaluated after the main update
  # has a non-zero weighted mean, so the intercept is actually updated.
  expect_gt(abs(mod_after$learners[[1]]$beta - mod_never$learners[[1]]$beta),
            1e-8)
  expect_gt(max(abs(mod_after$linear.predictors -
                      mod_never$linear.predictors)), 1e-8)
})
