# S3 methods of gbgam models and their compatibility with mgcv's predict,
# plot and summary machinery.

local_gamma_mod <- function() {
  fit_fixed(y ~ t2(x, k = 8, bs = "cr"),
            simulate_gamma_1d(n = 80, seed = 99),
            family = stats::Gamma(link = "log"),
            uq_control = gb_uq_ctrl("none"))
}

local_invgaussian_mod <- function() {
  fit_fixed(y ~ t2(x, k = 8, bs = "cr"),
            simulate_invgaussian_1d(n = 80, seed = 99),
            family = stats::inverse.gaussian(link = "log"),
            uq_control = gb_uq_ctrl("none"))
}

# One small model per family: Gaussian, Poisson, Gamma, inverse Gaussian
local_family_mods <- function() {
  list(local_gaussian_mod(), local_poisson_mod(), local_gamma_mod(),
       local_invgaussian_mod())
}

s3_sim_data <- function(n = 150L, seed = 11L) {
  set.seed(seed)
  d <- data.frame(
    x = stats::runif(n),
    z = stats::rnorm(n),
    g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  d$y <- sin(2 * pi * d$x) + 0.5 * d$z + 0.3 * as.numeric(d$g) +
    stats::rnorm(n, sd = 0.3)
  d
}

s3_fit <- function(formula, data, family = stats::gaussian(), max_iter = 30,
                   edf = "none", ...) {
  gbgam(
    formula,
    data         = data,
    family       = family,
    stop_control = gb_stop_ctrl("none", max_iter = max_iter),
    edf_control  = gb_edf_ctrl(edf),
    uq_control   = gb_uq_ctrl("none"),
    verbose      = FALSE,
    ...
  )
}


# ---- predict() ---------------------------------------------------------------

test_that("predict on training data matches fitted() (all families)", {
  skip_on_cran()
  for (mod in local_family_mods()) {
    pred <- predict(mod, type = "response")
    expect_equal(as.vector(pred), as.vector(fitted(mod)), tolerance = 1e-10,
                 label = mod$family$family)
  }
})

test_that("predict with newdata produces finite values on the link scale (all families)", {
  skip_on_cran()
  newdat <- data.frame(x = seq(0, 1, length.out = 20))
  for (mod in local_family_mods()) {
    pred <- predict(mod, newdata = newdat)
    expect_length(pred, 20L)
    expect_true(all(is.finite(pred)), label = mod$family$family)
  }
})

test_that("predict type='response' returns positive values for non-Gaussian log-link families", {
  # Gaussian is excluded: its response scale is unbounded.
  skip_on_cran()
  newdat <- data.frame(x = seq(0, 1, length.out = 15))
  for (mod in list(local_poisson_mod(), local_gamma_mod(), local_invgaussian_mod())) {
    pred <- predict(mod, newdata = newdat, type = "response")
    expect_length(pred, 15L)
    expect_true(all(is.finite(pred)), label = mod$family$family)
    expect_true(all(pred > 0),        label = mod$family$family)
  }
})

test_that("predict(newdata) works with a numeric parametric term and no factor", {
  d <- s3_sim_data()
  m <- s3_fit(y ~ z + t2(x, k = 6, bs = "cr"), d)
  expect_null(m$contrasts)
  expect_length(m$xlevels, 0L)
  expect_no_warning(p <- predict(m, newdata = d[1:5, ]))
  expect_equal(unname(as.vector(p)), unname(as.vector(predict(m)[1:5])),
               tolerance = 1e-12)
})

test_that("a factor used only inside a smooth is not stored as parametric", {
  d <- s3_sim_data()
  m <- s3_fit(y ~ t2(x, g, bs = c("cr", "re"), k = c(6, 3)), d)
  expect_null(m$contrasts)
  expect_false("g" %in% names(m$xlevels))
  expect_no_warning(p <- predict(m, newdata = d[1:5, ]))
  expect_equal(unname(as.vector(p)), unname(as.vector(predict(m)[1:5])),
               tolerance = 1e-12)
})

test_that("parametric factors keep their levels and contrasts", {
  d <- s3_sim_data()
  m <- s3_fit(y ~ z + g + t2(x, k = 6, bs = "cr"), d)
  expect_equal(m$xlevels, list(g = c("a", "b", "c")))
  expect_equal(m$contrasts, list(g = "contr.treatment"))
  nd <- d[d$g == "b", ][1:4, ]  # new data with a single observed level
  expect_equal(unname(as.vector(predict(m, newdata = nd))),
               unname(as.vector(predict(m)[as.integer(rownames(nd))])),
               tolerance = 1e-12)
})

test_that("predict(type = 'terms') gives non-zero parametric terms summing to the link", {
  d   <- s3_sim_data()
  d$y <- d$y + 2 * (d$g == "c")  # strong factor effect so g is selected
  m   <- s3_fit(y ~ z + g + t2(x, k = 6, bs = "cr"), d, max_iter = 250)
  expect_equal(m$assign, c(0L, 1L, 2L, 2L))

  pt <- predict(m, type = "terms")
  expect_true(all(c("z", "g") %in% colnames(pt)))
  expect_true(any(pt[, "z"] != 0))
  expect_true(any(pt[, "g"] != 0))
  expect_equal(unname(rowSums(pt) + attr(pt, "constant")),
               unname(as.vector(predict(m, type = "link"))),
               tolerance = 1e-10)

  # Same with new data
  nd  <- d[1:7, ]
  ptn <- predict(m, newdata = nd, type = "terms")
  expect_equal(unname(rowSums(ptn) + attr(ptn, "constant")),
               unname(as.vector(predict(m, newdata = nd, type = "link"))),
               tolerance = 1e-10)
})


# ---- fitted(), residuals() and weights() -------------------------------------

test_that("residuals() returns named types without error (Gaussian)", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  n   <- nrow(mod$model)
  for (tp in c("response", "working", "deviance", "pearson", "scaled.pearson")) {
    r <- residuals(mod, type = tp)
    expect_true(length(r) == n,    info = paste("wrong length for type =", tp))
    expect_true(all(is.finite(r)), info = paste("non-finite residuals for type =", tp))
  }
})

test_that("residuals() returns finite values for Poisson, Gamma, and inverse.gaussian", {
  skip_on_cran()
  for (mod in list(local_poisson_mod(), local_gamma_mod(), local_invgaussian_mod())) {
    n <- nrow(mod$model)
    for (tp in c("response", "working", "deviance", "pearson")) {
      r <- residuals(mod, type = tp)
      expect_true(length(r) == n,    info = paste("wrong length for type =", tp))
      expect_true(all(is.finite(r)), info = paste("non-finite residuals for type =", tp))
    }
  }
})

test_that("$residuals and $weights follow the glm working convention (Poisson)", {
  set.seed(13)
  n  <- 120
  dp <- data.frame(x = stats::runif(n), w = stats::runif(n, 0.5, 2))
  dp$y <- stats::rpois(n, exp(0.5 + sin(2 * pi * dp$x)))
  m <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), dp, family = poisson(),
              weights = dp$w)
  mu  <- m$fitted.values
  eta <- m$linear.predictors
  me  <- m$family$mu.eta(eta)
  w   <- dp$w

  expect_equal(m$user_weights, w)
  expect_equal(weights(m), m$prior.weights)
  expect_equal(weights(m, type = "working"), w * me^2 / mu)
  expect_equal(m$residuals, (m$y - mu) / me)

  expect_equal(residuals(m, "working"), (m$y - mu) / me)
  expect_equal(residuals(m, "response"), m$y - mu)
  pear <- (m$y - mu) * sqrt(w) / sqrt(mu)
  expect_equal(residuals(m, "pearson"), pear)
  expect_equal(residuals(m, "scaled.pearson"), pear)  # scale fixed at 1
  dev <- sign(m$y - mu) * sqrt(pmax(stats::poisson()$dev.resids(m$y, mu, w), 0))
  expect_equal(residuals(m, "deviance"), dev)

  # mgcv::plot.gam partial residuals use residuals * sqrt(weights), which
  # must be the Pearson residuals
  expect_equal(m$residuals * sqrt(m$weights), pear)
})

test_that("Gaussian working residuals/weights equal response residuals/prior weights", {
  d <- s3_sim_data()
  m <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d)
  expect_equal(m$residuals, m$y - m$fitted.values)
  expect_equal(m$weights, m$prior.weights)
  expect_equal(residuals(m, "scaled.pearson"),
               residuals(m, "pearson") / sqrt(m$sig2))
})

test_that("residuals() are padded under na.action = na.exclude", {
  d <- s3_sim_data()
  d$y[c(3, 17, 40)] <- NA
  m <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d, na.action = stats::na.exclude)
  expect_equal(nobs(m), nrow(d) - 3L)
  for (tp in c("response", "working", "deviance", "pearson", "scaled.pearson")) {
    r <- residuals(m, type = tp)
    expect_length(r, nrow(d))
    expect_true(all(is.na(r[c(3, 17, 40)])), info = tp)
    expect_true(all(is.finite(r[-c(3, 17, 40)])), info = tp)
  }
})

test_that("fitted() is padded under na.exclude; fitting and bootstrap work", {
  skip_on_cran()
  d <- s3_sim_data(n = 100L)
  d$y[c(3, 17, 40)] <- NA
  m <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = d,
             na.action          = stats::na.exclude,
             learner_control    = gb_learner_ctrl(step_size = 0.1),
             stop_control       = gb_stop_ctrl(iter_no_improv = 10L),
             edf_control        = gb_edf_ctrl("exact"),
             save_fitted_values = TRUE,
             verbose            = FALSE)
  n_used <- nrow(d) - 3L
  f <- fitted(m)
  expect_length(f, nrow(d))
  expect_true(all(is.na(f[c(3, 17, 40)])))
  expect_equal(unname(f[-c(3, 17, 40)]), m$fitted.values)
  expect_length(m$fitted.values, n_used)
  expect_equal(nrow(m$fitted_values), n_used)
  expect_true(is.finite(AIC(m)))

  # alpha != NULL also exercises the coverage calibration against the
  # (unpadded) fitted means; its size warnings are not the subject here
  set.seed(15)
  b <- withCallingHandlers(
    bootstrap_gbgam(m, n_samples = 5L, alpha = 0.05),
    warning = function(w) {
      if (grepl("widening factor", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_true(is.finite(b$widening_factor))
  p <- length(m$coefficients)
  expect_equal(dim(b$Vp), c(p, p))
  expect_true(all(is.finite(b$Vp)))
  expect_length(fitted(b), nrow(d))
})

test_that("weights() are padded under na.action = na.exclude", {
  set.seed(4)
  d   <- data.frame(x = stats::runif(100))
  d$y <- stats::rpois(100, exp(sin(2 * pi * d$x)))
  d$y[c(3, 40)] <- NA
  m <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d, family = stats::poisson(),
              na.action = stats::na.exclude)
  for (tp in c("prior", "working")) {
    w <- weights(m, type = tp)
    expect_length(w, nrow(d))
    expect_true(all(is.na(w[c(3, 40)])), info = tp)
    expect_true(all(is.finite(w[-c(3, 40)])), info = tp)
  }
  # The padded vectors line up: working residuals give Pearson residuals
  expect_equal(residuals(m, "working") * sqrt(weights(m, "working")),
               residuals(m, "pearson"))
})

test_that("weights() returns the prior weights vector", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  w   <- weights(mod)
  expect_length(w, nrow(mod$model))
  expect_true(all(w >= 0))
})


# ---- coef(), vcov(), logLik(), AIC() and other extractors --------------------

test_that("coef() returns named coefficients matching mod$coefficients", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  expect_equal(coef(mod), mod$coefficients)
})

test_that("X beta rebuilds the stored linear predictor (all families)", {
  skip_on_cran()
  for (mod in local_family_mods()) {
    eta_rebuilt <- rebuild_training_eta(mod)
    expect_equal(eta_rebuilt, mod$linear.predictors, tolerance = 1e-8,
                 label = mod$family$family)
  }
})

test_that("vcov() warns when Vp is absent, returns matrix when present", {
  skip_on_cran()
  mod_no_vp <- local_gaussian_mod()
  expect_warning(vcov(mod_no_vp), regexp = "No covariance")

  expect_warning(
    mod_vp <- bootstrap_gbgam(mod_no_vp, n_samples = 8L, alpha = NULL,
                              verbose = FALSE),
    "edf_compute_method"
  )
  V <- vcov(mod_vp)
  p <- length(mod_vp$coefficients)
  expect_equal(dim(V), c(p, p))
  expect_true(all(is.finite(V)))
})

test_that("logLik() returns a logLik object with correct attributes (Gaussian)", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  ll  <- logLik(mod)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_equal(attr(ll, "nobs"), nrow(mod$model))
})

test_that("logLik() returns a finite logLik for Poisson, Gamma, and inverse.gaussian", {
  skip_on_cran()
  for (mod in list(local_poisson_mod(), local_gamma_mod(), local_invgaussian_mod())) {
    ll <- logLik(mod)
    expect_s3_class(ll, "logLik")
    expect_true(is.finite(as.numeric(ll)))
    expect_equal(attr(ll, "nobs"), nrow(mod$model))
  }
})

test_that("Gaussian logLik responds to the units of y like a log-likelihood", {
  d  <- s3_sim_data()
  m1 <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d, edf = "exact")
  d$y <- 10 * d$y
  m2 <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d, edf = "exact")
  n  <- nrow(d)
  expect_equal(m2$edf.total, m1$edf.total, tolerance = 1e-8)
  expect_equal(as.numeric(logLik(m1)) - as.numeric(logLik(m2)),
               n * log(10), tolerance = 1e-8)
})

test_that("AIC() and BIC() match their definitions from logLik.gbgam", {
  skip_on_cran()
  # With EDF tracking: without it the df, and so AIC() and BIC(), are NA
  mod <- local_gaussian_mod_edf()
  ll  <- logLik(mod)
  df  <- attr(ll, "df")
  n   <- attr(ll, "nobs")

  expect_true(is.finite(AIC(mod)))
  expect_true(is.finite(BIC(mod)))

  # AIC = -2 * logLik + 2 * df (definition used by stats::AIC.default)
  expect_equal(AIC(mod), -2 * as.numeric(ll) + 2 * df, tolerance = 1e-10)

  # BIC = -2 * logLik + log(n) * df
  expect_equal(BIC(mod), -2 * as.numeric(ll) + log(n) * df, tolerance = 1e-10)

  # BIC >= AIC for n > e^2 (approx 7.4); holds for any real dataset
  expect_gte(BIC(mod), AIC(mod))
})

test_that("AIC() matches the internal stopping AIC (Gaussian: +2 for the scale)", {
  d <- s3_sim_data()
  m <- gbgam(y ~ z + t2(x, k = 6, bs = "cr"), data = d, verbose = FALSE,
             stop_control = gb_stop_ctrl(iter_no_improv = 10))
  ll <- logLik(m)
  expect_equal(as.numeric(ll), -neg2loglik(m) / 2)
  expect_equal(attr(ll, "df"), m$edf.total + 1)
  expect_equal(AIC(m) - 2, min(m$metrics$ls, na.rm = TRUE), tolerance = 1e-8)

  # Gaussian log-likelihood at the estimated scale
  expect_equal(
    as.numeric(ll),
    sum(stats::dnorm(m$y, m$fitted.values, sqrt(m$scale), log = TRUE)),
    tolerance = 1e-10
  )
})

test_that("AIC() matches the internal stopping AIC (Poisson: no scale df)", {
  set.seed(12)
  n  <- 150
  dp <- data.frame(x = stats::runif(n))
  dp$y <- stats::rpois(n, exp(0.5 + sin(2 * pi * dp$x)))
  m <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dp, family = poisson(),
             verbose = FALSE, stop_control = gb_stop_ctrl(iter_no_improv = 10))
  ll <- logLik(m)
  expect_equal(attr(ll, "df"), m$edf.total)
  expect_equal(as.numeric(ll),
               sum(stats::dpois(m$y, m$fitted.values, log = TRUE)),
               tolerance = 1e-10)
  expect_equal(AIC(m), min(m$metrics$ls, na.rm = TRUE), tolerance = 1e-8)
})

test_that("nobs() returns the number of training rows", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 77, seed = 200)
  mod <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat,
               learner_control = gb_learner_ctrl(),
               stop_control    = gb_stop_ctrl("none", max_iter = 5),
               edf_control     = gb_edf_ctrl("none"),
               uq_control      = gb_uq_ctrl("none"),
               verbose         = FALSE)
  expect_equal(nobs(mod), 77L)
})

test_that("nobs() counts the rows with non-zero weight, as glm and mgcv do", {
  dat <- simulate_poisson_1d(n = 120, seed = 201)
  w   <- rep(c(1, 2, 0, 1), length.out = nrow(dat))
  mod <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), dat, family = stats::poisson(),
                max_iter = 10, edf = "exact", weights = w)
  ref <- stats::glm(y ~ x, family = stats::poisson(), data = dat,
                    weights = w)
  expect_identical(nobs(mod), sum(w != 0))
  expect_identical(nobs(mod), nobs(ref))
  ll <- logLik(mod)
  expect_identical(attr(ll, "nobs"), sum(w != 0))
  expect_equal(BIC(mod),
               -2 * as.numeric(ll) + log(sum(w != 0)) * attr(ll, "df"))
  # binomial rows without trials do not count either
  dat$trials <- rep(c(5L, 0L, 3L), length.out = nrow(dat))
  dat$succ   <- pmin(dat$y, dat$trials)
  mb <- s3_fit(cbind(succ, trials - succ) ~ t2(x, k = 6, bs = "cr"), dat,
               family = stats::binomial(), max_iter = 10)
  expect_identical(nobs(mb), sum(dat$trials > 0))
})

test_that("deviance() and df.residual() return correct scalars (Gaussian)", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  expect_equal(deviance(mod), mod$deviance)
  expect_equal(df.residual(mod), mod$df.residual)
})

test_that("deviance() returns a finite positive scalar for Poisson, Gamma, and inverse.gaussian", {
  skip_on_cran()
  for (mod in list(local_poisson_mod(), local_gamma_mod(), local_invgaussian_mod())) {
    expect_true(is.finite(deviance(mod)) && deviance(mod) > 0)
  }
})

test_that("formula.gbgam returns the original formula", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 50)
  fml <- y ~ t2(x, k = 6, bs = "cr")
  mod <- gbgam(fml, data = dat,
               learner_control = gb_learner_ctrl(),
               stop_control    = gb_stop_ctrl("none", max_iter = 10),
               edf_control     = gb_edf_ctrl("none"),
               uq_control      = gb_uq_ctrl("none"),
               verbose         = FALSE)
  expect_equal(deparse(formula(mod)), deparse(fml))
})


# ---- plot(), print() and summary() -------------------------------------------

test_that("plot.gbgam runs without error for a basic Gaussian model", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  expect_no_error({
    grDevices::pdf(file = NULL)
    plot(mod)
    grDevices::dev.off()
  })
})

test_that("plot.gbgam labels include EDF when edf is computed", {
  skip_on_cran()
  mod <- local_gaussian_mod_edf()
  grDevices::pdf(file = NULL)
  plt <- plot(mod, pages = 1)
  grDevices::dev.off()
  expect_true(is.list(plt))
  expect_length(plt, length(mod$smooth))
  for (j in seq_along(mod$smooth)) {
    sm    <- mod$smooth[[j]]
    edf_j <- sum(mod$edf[sm$first.para:sm$last.para])
    expect_gt(edf_j, 0)
    expect_match(plt[[j]]$ylab, paste0(",", round(edf_j, 2), ")"),
                 fixed = TRUE)
  }
})

test_that("plot() works with select = and without Vp", {
  d <- s3_sim_data()
  m <- s3_fit(y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 6, bs = "cr"), d)
  expect_null(m$Vp)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(m, select = 1))
  expect_no_error(plot(m, select = 2, residuals = TRUE))
  expect_warning(plot(m, select = 1, se = TRUE), "No covariance matrix")
})

test_that("plot(residuals = TRUE) works under na.action = na.exclude", {
  d <- s3_sim_data()
  d$y[c(3, 17)] <- NA
  m_excl <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d,
                   na.action = stats::na.exclude)
  m_omit <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d, na.action = stats::na.omit)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p_excl <- plot(m_excl, residuals = TRUE)
  p_omit <- plot(m_omit, residuals = TRUE)
  expect_length(p_excl[[1L]]$p.resid, nrow(d) - 2L)
  expect_equal(p_excl[[1L]]$p.resid, p_omit[[1L]]$p.resid)
})

test_that("predict(se.fit = TRUE) without Vp gives an informative error", {
  d <- s3_sim_data()
  m <- s3_fit(y ~ t2(x, k = 6, bs = "cr"), d)
  expect_null(m$Vp)
  nd <- d[1:5, ]
  expect_error(predict(m, nd, se.fit = TRUE), "No covariance matrix")
  # positional and partially matched se.fit
  expect_error(predict(m, nd, "terms", TRUE), "No covariance matrix")
  expect_error(predict(m, nd, se = TRUE), "No covariance matrix")
  # without standard errors the result is that of mgcv::predict.gam()
  expect_identical(predict(m, nd, type = "response"),
                   mgcv::predict.gam(m, nd, type = "response"))
  expect_identical(predict(m, nd, se.fit = FALSE),
                   mgcv::predict.gam(m, nd))
})

test_that("summary.gbgam runs without error and returns a summary.gbgam object (Gaussian)", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  out <- expect_no_error(summary(mod))
  expect_s3_class(out, "summary.gbgam")
  expect_true(is.data.frame(out$learner_table))
  expect_true(nrow(out$learner_table) >= 1L)
})

test_that("summary.gbgam runs without error and returns a summary.gbgam object (non-Gaussian families)", {
  skip_on_cran()
  for (mod in list(local_poisson_mod(), local_gamma_mod(), local_invgaussian_mod())) {
    out <- expect_no_error(summary(mod))
    expect_s3_class(out, "summary.gbgam")
    expect_true(is.data.frame(out$learner_table), info = mod$family$family)
    expect_true(nrow(out$learner_table) >= 1L,    info = mod$family$family)
  }
})

test_that("print.gbgam runs without error", {
  skip_on_cran()
  mod <- local_gaussian_mod()
  expect_no_error(print(mod))
})

test_that("print(summary()) shows long formulas without cli list separators", {
  d <- s3_sim_data()
  d$x2 <- stats::runif(nrow(d))
  m <- s3_fit(
    y ~ t2(x, k = 5, bs = "cr") + t2(z, k = 5, bs = "cr") + t2(x2, k = 5, bs = "cr"),
    d, max_iter = 5
  )
  expect_gt(length(deparse(m$formula)), 1L)
  out <- paste(cli::cli_fmt(print(summary(m))), collapse = " ")
  expect_false(grepl(" and ", out, fixed = TRUE))
  expect_true(grepl("t2(x2, k = 5, bs = \"cr\")", out, fixed = TRUE))
})

test_that("print(summary()) aligns its columns for long learner labels", {
  d <- data.frame(soil_moisture_percentage = stats::runif(150),
                  mean_annual_temperature_celsius = stats::runif(150))
  d$y <- sin(2 * pi * d$soil_moisture_percentage) *
    d$mean_annual_temperature_celsius + stats::rnorm(150, sd = 0.3)
  m <- s3_fit(y ~ t2(soil_moisture_percentage, mean_annual_temperature_celsius,
                     k = c(4, 4)),
              d, max_iter = 20, edf = "exact")
  labels <- summary(m)$learner_table$learner
  token_ends <- function(x) {
    g <- gregexpr("\\S+", x)[[1L]]
    as.integer(g + attr(g, "match.length") - 1L)
  }
  table_lines <- function(width) {
    op <- options(width = width)
    on.exit(options(op))
    out <- cli::ansi_strip(cli::cli_fmt(print(summary(m))))
    out[seq(grep("^  learner ", out), length(out))]
  }
  for (width in c(80L, 200L)) {
    tl   <- table_lines(width)
    rows <- tl[-1L]
    expect_length(rows, length(labels))
    nm <- sub("^  (\\S+).*$", "\\1", rows)
    # every row keeps its own learner suffix, and the names differ
    expect_false(anyDuplicated(nm) > 0L, info = width)
    expect_true(all(endsWith(nm, sub("^.*\\)", "", labels))), info = width)
    if (width == 200L) expect_identical(nm, labels)
    if (width == 80L) expect_true(all(nchar(tl) <= 80L))
    # right-aligned columns: sel, %sel and acc_edf end where their headers end
    for (r in rows)
      expect_identical(token_ends(r)[2:4], token_ends(tl[1L])[2:4],
                       info = paste(width, r))
  }
})

test_that("print and summary separate the end of the run from mstop", {
  d <- core_data()
  m <- gbgam(ff_core, data = d, verbose = FALSE,
             stop_control = gb_stop_ctrl(iter_no_improv = 20, max_iter = 300))
  expect_output(print(m), "no AIC improvement after 20 iterations")
  e     <- extend_gbgam(m, n_iter = 30)
  m_run <- sum(!is.na(e$metrics$k_opts))
  expect_lt(e$mstop, m_run)
  expect_output(print(e),
                paste0("mstop=", e$mstop, " \\(run ended at iteration ",
                       m_run, ": max_iter = ", e$max_iter,
                       " reached\\), best AIC="))
  expect_match(summary(e)$stop_msgs, "run ended at iteration")

  m2 <- gbgam(ff_core, data = d, verbose = FALSE,
              stop_control = gb_stop_ctrl(max_iter = 10))
  expect_equal(m2$mstop, 10)
  expect_output(print(m2), "mstop=10 \\(max_iter = 10 reached\\), best AIC=")
})

test_that("print and summary show small best losses with significant digits", {
  d <- core_data()
  m <- gbgam(ff_core, data = d, verbose = FALSE,
             stop_control = gb_stop_ctrl(max_iter = 10))
  # e.g. a true_f mean squared error of a well-fitted model
  m$metrics$ls <- m$metrics$ls - min(m$metrics$ls) + 0.001246099
  expect_output(print(m), "best AIC=0.001246", fixed = TRUE)
  out <- cli::ansi_strip(cli::cli_fmt(print(summary(m))))
  expect_true(any(grepl("best AIC\\s+0\\.001246\\s", out)))
})

test_that("%d_stop_loss shares are signed and sum to 100", {
  d <- core_data()
  m <- gbgam(y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 6, bs = "cr"),
             data = transform(d, z = stats::runif(nrow(d))), verbose = FALSE,
             stop_control = gb_stop_ctrl(iter_no_improv = 20))
  sh <- summary(m)$learner_table[["%d_stop_loss"]]
  expect_equal(sum(sh), 100, tolerance = 0.5)
  expect_false(any(1 / sh[sh == 0] < 0))  # no negative zeros
})
