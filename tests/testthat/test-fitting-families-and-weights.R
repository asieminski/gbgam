# Non-Gaussian families, response formats (binomial), weights, offset and
# subset.


# ---- Fits by family ----------------------------------------------------------

test_that("Poisson fits are finite and fitted means are positive", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 100, seed = 40)

  mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, stats::poisson(),
                   max_iter = 20)

  expect_model_contract(mod, nrow(dat))
  expect_true(all(mod$fitted.values > 0))
})

test_that("family can be supplied as character or family object with identical results", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 80, seed = 41)

  mod_chr <- fit_fixed(y ~ t2(x, k = 7, bs = "cr"), dat, "poisson",
                       max_iter = 18)

  mod_obj <- fit_fixed(y ~ t2(x, k = 7, bs = "cr"), dat, stats::poisson(),
                       max_iter = 18)

  expect_equal(mod_chr$family$family, "poisson")
  expect_equal(mod_obj$family$family, "poisson")
  expect_equal(mod_chr$metrics$k_opts,    mod_obj$metrics$k_opts)
  expect_equal(mod_chr$linear.predictors, mod_obj$linear.predictors, tolerance = 1e-12)
})

test_that("binomial response encodings give equivalent fitted probabilities", {
  skip_on_cran()
  dd <- simulate_binomial_formats(n = 70, seed = 42)

  mod_prop <- fit_fixed(prop ~ t2(x, k = 8, bs = "cr"), dd$agg,
                        stats::binomial(), max_iter = 20,
                        weights = dd$agg$n_trials)
  mod_mat  <- fit_fixed(cbind(successes, failures) ~ t2(x, k = 8, bs = "cr"),
                        dd$agg, stats::binomial(), max_iter = 20)
  mod_long <- fit_fixed(y_fac ~ t2(x, k = 8, bs = "cr"), dd$long,
                        stats::binomial(), max_iter = 20)

  p_prop     <- mod_prop$fitted.values
  p_mat      <- mod_mat$fitted.values
  p_long_agg <- as.vector(tapply(mod_long$fitted.values, dd$long$id, mean))

  expect_equal(p_prop, p_mat,      tolerance = 1e-6)
  expect_equal(p_prop, p_long_agg, tolerance = 1e-6)
})

test_that("zero observation weights are handled gracefully", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 43)
  w <- c(rep(0, 15), rep(1, 85))

  mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, weights = w,
                   max_iter = 15)

  expect_model_contract(mod, nrow(dat))
})

test_that("nontrivial weights alter the fit", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 100, seed = 44)
  w1 <- rep(1, nrow(dat))
  w2 <- seq_len(nrow(dat))

  mod1 <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, weights = w1,
                    max_iter = 15)
  mod2 <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, weights = w2,
                    max_iter = 15)

  expect_false(isTRUE(all.equal(
    mod1$linear.predictors,
    mod2$linear.predictors,
    tolerance = 1e-10
  )))
})

test_that("Gamma fits are finite and fitted means are positive", {
  skip_on_cran()
  dat <- simulate_gamma_1d(n = 100, seed = 45)

  mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat,
                   stats::Gamma(link = "log"), max_iter = 20)

  expect_model_contract(mod, nrow(dat))
  expect_true(all(mod$fitted.values > 0))
})

test_that("inverse.gaussian fits are finite and fitted means are positive", {
  skip_on_cran()
  dat <- simulate_invgaussian_1d(n = 100, seed = 46)

  mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat,
                   stats::inverse.gaussian(link = "log"), max_iter = 20)

  expect_model_contract(mod, nrow(dat))
  expect_true(all(mod$fitted.values > 0))
})


# ---- Families and response formats -------------------------------------------

test_that("gbgam() fits Gamma() and inverse.gaussian() with default links", {
  set.seed(202)
  n <- 150
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rgamma(n, shape = 5, rate = 5 / (0.4 + 0.3 * d$x))
  expect_lt(mean(d$y), 1)
  for (fam in list(stats::Gamma(), stats::inverse.gaussian())) {
    mod <- quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = d, family = fam,
                     learner_control = gb_learner_ctrl(step_size = 0.1),
                     stop_control = gb_stop_ctrl("aic", max_iter = 60L,
                                                 iter_no_improv = 20L))
    expect_equal(mod$intercept0, fam$linkfun(mean(d$y)), tolerance = 1e-8)
    expect_true(all(is.finite(mod$fitted.values)))
    expect_true(all(mod$fitted.values > 0))
  }
})

test_that("binomial accepts zero-trial rows and logical responses", {
  set.seed(207)
  n <- 120
  d <- data.frame(x = stats::runif(n))
  p <- stats::plogis(sin(2 * pi * d$x))
  tr <- sample(0:4, n, replace = TRUE)
  tr[1:20] <- 0
  d$s <- stats::rbinom(n, tr, p)
  d$f <- tr - d$s
  fit_b <- function(fml, dat) {
    quiet_fit(fml, data = dat, family = stats::binomial(),
              stop_control = sc_ff, edf_control = gb_edf_ctrl("exact"),
              learner_control = gb_learner_ctrl(step_size = 0.1))
  }
  full <- fit_b(cbind(s, f) ~ t2(x, k = 6, bs = "cr"), d)
  expect_true(all(is.finite(full$linear.predictors)))
  expect_equal(unname(full$y[tr == 0]), rep(0, sum(tr == 0)))
  expect_equal(unname(full$prior.weights), tr)
  # The zero-trial rows do not change the intercept-only MLE, which matches
  # glm()
  kept <- fit_b(cbind(s, f) ~ t2(x, k = 6, bs = "cr"), d[tr > 0, ])
  expect_equal(full$intercept0, kept$intercept0, tolerance = 1e-10)
  ref <- stats::glm(cbind(s, f) ~ 1, family = stats::binomial(), data = d)
  expect_equal(full$intercept0, unname(stats::coef(ref)), tolerance = 1e-6)

  # Logical response equals its 0/1 version
  d$yb <- stats::rbinom(n, 1, p) == 1
  d$yn <- as.numeric(d$yb)
  m_l <- fit_b(yb ~ t2(x, k = 6, bs = "cr"), d)
  m_n <- fit_b(yn ~ t2(x, k = 6, bs = "cr"), d)
  expect_equal(m_l$y, d$yn)
  expect_equal(m_l$linear.predictors, m_n$linear.predictors)
})

test_that("Binomial gbgam pipeline: aggregate format fits and predicts", {
  skip_on_cran()
  bindat <- simulate_binomial_formats(n = 60, seed = 81)
  agg <- bindat$agg
  mod <- fit_fixed(cbind(successes, failures) ~ t2(x, k = 6, bs = "cr"), agg,
                   stats::binomial(), max_iter = 20, step_size = 0.1)
  expect_s3_class(mod, "gbgam")
  fv <- fitted(mod)
  expect_true(all(fv >= 0 & fv <= 1))
})


# ---- Weights, offset and subset ----------------------------------------------

test_that("weights expressions mix data columns and caller variables", {
  set.seed(204)
  n <- 60
  d <- data.frame(x = stats::runif(n), w = stats::runif(n, 0.5, 2))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(n, sd = 0.3)
  fml <- y ~ t2(x, k = 6, bs = "cr")
  # data$w takes precedence over a same-named variable of the caller, and
  # the other symbols come from the caller
  w <- rep(1, n)
  k <- 2
  mixed <- gbgam(fml, data = d, weights = w * k, stop_control = sc_ff,
                 edf_control = ec_ff, verbose = FALSE)
  expect_equal(mixed$prior.weights, d$w * 2)
  # A user vector named like a base function
  scale <- stats::runif(n, 0.5, 2)
  m_sc <- gbgam(fml, data = d, weights = scale, stop_control = sc_ff,
                edf_control = ec_ff, verbose = FALSE)
  expect_equal(m_sc$prior.weights, scale)
  # A wrapper local, passed directly and forwarded via `...`
  wrap_local <- function(dat) {
    my_w <- dat$w * 3
    quiet_fit(fml, data = dat, weights = my_w, stop_control = sc_ff,
              edf_control = ec_ff)
  }
  expect_equal(wrap_local(d)$prior.weights, d$w * 3)
  wrap_dots <- function(dat, ...) {
    quiet_fit(fml, data = dat, stop_control = sc_ff, edf_control = ec_ff,
              ...)
  }
  vv <- stats::runif(n)
  expect_equal(wrap_dots(d, weights = vv)$prior.weights, vv)
  expect_equal(wrap_dots(d, weights = w)$prior.weights, d$w)
  expect_equal(wrap_dots(d, weights = d$w / 2)$prior.weights, d$w / 2)
})

test_that("weights can be a column of data, a global or a wrapper local", {
  set.seed(111)
  n <- 80
  d <- data.frame(x = stats::runif(n), ww = stats::runif(n, 0.5, 2))
  d$y <- sin(2 * pi * d$x) + stats::rnorm(n, sd = 0.3)
  sc <- gb_stop_ctrl("none", max_iter = 20L)
  ec <- gb_edf_ctrl("none")
  fml <- y ~ t2(x, k = 6, bs = "cr")
  ref <- quiet_fit(fml, data = d, weights = d$ww, stop_control = sc,
                   edf_control = ec)
  # A column of data (a same-named global must not mask it)
  ww <- rep(1, n)
  col <- quiet_fit(fml, data = d, weights = ww, stop_control = sc,
                   edf_control = ec)
  expect_equal(col$prior.weights, d$ww)
  expect_equal(col$linear.predictors, ref$linear.predictors)
  # A local variable in a wrapper, also when it is called `weights`
  wrap <- function(dat, f) {
    weights <- dat$ww
    quiet_fit(f, data = dat, weights = weights, stop_control = sc,
              edf_control = ec)
  }
  expect_equal(wrap(d, fml)$prior.weights, d$ww)
  # A global (non-data) weights vector passed by name
  gw <- stats::runif(n, 0.5, 2)
  expect_equal(quiet_fit(fml, data = d, weights = gw, stop_control = sc,
                         edf_control = ec)$prior.weights, gw)
  # subset is applied to the weights as well
  sub <- gbgam(fml, data = d, weights = ww, subset = x > 0.3,
               stop_control = sc, edf_control = ec, verbose = FALSE)
  expect_equal(sub$prior.weights, d$ww[d$x > 0.3])
})

test_that("offset and subset can be forwarded through ... like weights", {
  d  <- core_data()
  n  <- nrow(d)
  sc <- gb_stop_ctrl("none", max_iter = 5)
  f  <- function(...) gbgam(..., verbose = FALSE, stop_control = sc)
  g  <- function() {
    off <- seq(0, 0.2, length.out = n)
    ss  <- d$x > 0.3
    f(ff_core, data = d, offset = off, subset = ss)
  }
  m_wrap <- g()
  m_dir  <- gbgam(ff_core, data = d, verbose = FALSE, stop_control = sc,
                  offset = seq(0, 0.2, length.out = n), subset = d$x > 0.3)
  expect_equal(nrow(m_wrap$model), sum(d$x > 0.3))
  expect_equal(m_wrap$offset, m_dir$offset)
  expect_equal(m_wrap$coefficients, m_dir$coefficients)

  # Columns of data are used for offset and subset
  d$o   <- seq(0, 0.2, length.out = n)
  m_col <- gbgam(ff_core, data = d, verbose = FALSE, stop_control = sc,
                 offset = o, subset = x > 0.3)
  expect_equal(m_col$coefficients, m_dir$coefficients)
  expect_equal(m_col$offset, m_dir$offset)
})

test_that("weights summing to less than the EDF stop AIC stopping", {
  d <- core_data()
  expect_error(gbgam(ff_core, data = d, verbose = FALSE,
                     weights = rep(0.01, nrow(d))),
               "dispersion cannot be estimated")
})

test_that("default weights do not stop AIC stopping when the EDF passes n", {
  set.seed(2)
  n <- 25
  d <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n),
                  x4 = runif(n))
  d$y <- sin(2 * pi * d$x1) + rnorm(n, sd = 0.05)
  ff  <- y ~ t2(x1, k = 8, bs = "cr") + t2(x2, k = 8, bs = "cr") +
    t2(x3, k = 8, bs = "cr") + t2(x4, k = 8, bs = "cr")
  mod <- gbgam(ff, data = d, verbose = FALSE,
               learner_control = gb_learner_ctrl(step_size = 0.5),
               edf_control     = gb_edf_ctrl("hutchpp", z = 10),
               stop_control    = gb_stop_ctrl("aic", iter_no_improv = 300,
                                              max_iter = 1500))
  expect_gt(max(mod$edfs, na.rm = TRUE), n)
})
