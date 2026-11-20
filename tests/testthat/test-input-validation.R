# Checks of the arguments of gbgam(): controls, formula, family, response.


# ---- Controls checked by gbgam() ---------------------------------------------

test_that("gbgam warns and ignores tau when edf_compute_method = 'none'", {
  dat <- simulate_gaussian_1d(n = 40, seed = 10)
  expect_warning(
    gbgam(
      y ~ t2(x, k = 6, bs = "cr"),
      data            = dat,
      learner_control = gb_learner_ctrl(tau = 1),
      stop_control    = gb_stop_ctrl("none", max_iter = 5),
      edf_control     = gb_edf_ctrl("none"),
      uq_control      = gb_uq_ctrl("none"),
      verbose         = FALSE
    ),
    "tau is set to 0"
  )
})

test_that("gbgam rejects AIC stopping without EDF computation", {
  dat <- simulate_gaussian_1d(n = 40, seed = 10)
  expect_error(
    gbgam(
      y ~ t2(x, k = 6, bs = "cr"),
      data            = dat,
      learner_control = gb_learner_ctrl(),
      stop_control    = gb_stop_ctrl("aic", max_iter = 5),
      edf_control     = gb_edf_ctrl("none"),
      uq_control      = gb_uq_ctrl("none"),
      verbose         = FALSE
    ),
    "requires EDF tracking"
  )
})

test_that("gbgam rejects out-of-range preselected_learner_indices", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 40, seed = 13)
  expect_error(
    gbgam(
      y ~ t2(x, k = 6, bs = "cr"),
      data            = dat,
      learner_control = gb_learner_ctrl(preselected_learner_indices = rep(99L, 5)),
      stop_control    = gb_stop_ctrl("none", max_iter = 5),
      edf_control     = gb_edf_ctrl("none"),
      uq_control      = gb_uq_ctrl("none"),
      verbose         = FALSE
    ),
    "preselected_learner_indices"
  )
})

test_that("gbgam rejects out-of-range excluded_learner_indices", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 40, seed = 14)
  expect_error(
    gbgam(
      y ~ t2(x, k = 6, bs = "cr"),
      data            = dat,
      learner_control = gb_learner_ctrl(excluded_learner_indices = 99L),
      stop_control    = gb_stop_ctrl("none", max_iter = 5),
      edf_control     = gb_edf_ctrl("none"),
      uq_control      = gb_uq_ctrl("none"),
      verbose         = FALSE
    ),
    "excluded_learner_indices"
  )
})

test_that("learners_to_fit_after_each_iter is bounds-checked", {
  dat <- simulate_gaussian_1d(n = 60, seed = 110)
  expect_error(
    quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = dat,
              learner_control = gb_learner_ctrl(
                learners_to_fit_after_each_iter = 99L)),
    "learners_to_fit_after_each_iter contains out-of-range"
  )
})

test_that("incompatible bootstrap settings fail before the fit", {
  d <- core_data()
  local_mocked_bindings(fit_gbgam = function(...) stop("the fit ran"))
  for (mode in c("free_smoothing", "free_selection", "free_stopping")) {
    expect_error(
      gbgam(ff_core, data = d, verbose = FALSE,
            edf_control = gb_edf_ctrl("none"),
            stop_control = gb_stop_ctrl("none", max_iter = 20),
            uq_control = gb_uq_ctrl("bootstrap", mode = mode,
                                    n_samples = 3)),
      "requires EDF tracking")
  }
  expect_error(
    gbgam(ff_core, data = d, verbose = FALSE,
          stop_control = gb_stop_ctrl("none", max_iter = 20),
          uq_control = gb_uq_ctrl("bootstrap",
                                  mode = "free_stopping",
                                  n_samples = 3)),
    "requires the main fit to stop by stopping_method = \"aic\" or")
  expect_error(
    gbgam(ff_core, data = d, verbose = FALSE,
          stop_control = gb_stop_ctrl("cv"),
          uq_control = gb_uq_ctrl(
            "bootstrap", mode = "free_stopping", n_samples = 3)),
    "not available with stopping_method = \"cv\"")
})


# ---- Formula, family and response --------------------------------------------

test_that("unsupported formula constructs are rejected", {
  dat <- simulate_gaussian_1d(n = 40, seed = 11)
  expect_error(
    new_learners(y ~ 0 + t2(x, k = 6, bs = "cr"), data = dat),
    "requires an intercept"
  )
  for (fml in list(y ~ s(x, k = 6), y ~ te(x, k = 5), y ~ ti(x, k = 5))) {
    expect_error(fit_fixed(fml, dat, max_iter = 5, step_size = 0.1),
                 "not supported")
  }
})

test_that("a formula without an intercept names the alternative", {
  d <- core_data()
  expect_error(gbgam(y ~ t2(x, k = 6) - 1, data = d, verbose = FALSE),
               "initialise_intercept = FALSE, excluded_learner_indices = 1L")
})

test_that("new_family() accepts names, functions and family objects", {
  # Names and family objects
  expect_identical(new_family("gaussian")$family, "gaussian")
  expect_identical(new_family("poisson")$family, "poisson")
  expect_identical(new_family(stats::poisson())$family, "poisson")
  expect_identical(new_family(stats::binomial())$family, "binomial")
  expect_identical(new_family(stats::binomial("probit"))$link, "probit")
  expect_identical(new_family(stats::Gamma(link = "log"))$family, "Gamma")
  expect_identical(
    new_family(stats::inverse.gaussian(link = "log"))$family,
    "inverse.gaussian"
  )
  # A generator function (not its result) is called
  expect_identical(new_family(stats::Gamma)$family, "Gamma")
  expect_identical(new_family(stats::poisson)$family, "poisson")
  expect_identical(new_family(stats::binomial)$family, "binomial")
  # Unsupported families and other types
  expect_error(new_family("quasi"), "Must be element")
  expect_error(new_family("quasipoisson"), "Must be element")
  expect_error(new_family(stats::quasipoisson()), "Must be element")
  expect_error(new_family(42L), "family must be a character string")
  expect_error(new_family(1), "family must be")
})

test_that("a negative weight is rejected", {
  dat <- simulate_gaussian_1d(n = 40, seed = 12)
  w_bad <- rep(1, nrow(dat))
  w_bad[20] <- -1
  expect_error(
    fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 5,
              step_size = 0.1, weights = w_bad),
    "Assertion on 'weights' failed: Element 20 is not >= 0"
  )
})

test_that("a binomial factor response needs exactly two levels", {
  dat <- simulate_gaussian_1d(n = 40, seed = 12)
  dat$y <- factor(sample(letters[1:3], nrow(dat), replace = TRUE))
  expect_error(
    fit_fixed(y ~ t2(x, k = 6, bs = "cr"), dat, "binomial", max_iter = 5,
              step_size = 0.1),
    "Must have exactly 2 levels"
  )
})

test_that("binomial numeric responses must lie in [0, 1]", {
  set.seed(109)
  n <- 60
  d <- data.frame(x = stats::runif(n))
  d$y3 <- stats::rbinom(n, 3, 0.5)
  d$yp <- d$y3 / 3
  expect_error(
    quiet_fit(y3 ~ t2(x, k = 6, bs = "cr"), data = d,
              family = stats::binomial()),
    "<= 1"
  )
  expect_warning(
    quiet_fit(yp ~ t2(x, k = 6, bs = "cr"), data = d,
              family = stats::binomial(),
              stop_control = gb_stop_ctrl("none", max_iter = 2L),
              edf_control = gb_edf_ctrl("none")),
    "Fractional binomial response but no weights"
  )
  expect_silent(
    quiet_fit(yp ~ t2(x, k = 6, bs = "cr"), data = d,
              family = stats::binomial(), weights = rep(3, n),
              stop_control = gb_stop_ctrl("none", max_iter = 2L),
              edf_control = gb_edf_ctrl("none"))
  )
})

test_that("response checks name the response", {
  dat <- simulate_gaussian_1d(n = 40, seed = 12)
  dat$yn <- stats::rpois(nrow(dat), 2)
  dat$yn[4] <- -1
  expect_error(
    fit_fixed(yn ~ t2(x, k = 6, bs = "cr"), dat, "poisson", max_iter = 5),
    "Assertion on 'response' failed: Element 4 is not >= 0"
  )
  dat$yf <- factor(dat$y > 0)
  expect_error(
    fit_fixed(yf ~ t2(x, k = 6, bs = "cr"), dat, max_iter = 5),
    "Assertion on 'response' failed: Must be of type 'numeric'"
  )
  dat$yb <- dat$y
  expect_error(
    fit_fixed(yb ~ t2(x, k = 6, bs = "cr"), dat, "binomial", max_iter = 5),
    "Assertion on 'response' failed"
  )
})


# ---- Covariate types ---------------------------------------------------------

test_that("a character covariate is handled like a factor", {
  dat <- simulate_gaussian_1d(n = 90, seed = 15)
  dat$ch <- rep(c("u", "v", "w"), 30)
  dat$y  <- dat$y + 0.5 * (dat$ch == "v")
  dat_f  <- dat
  dat_f$ch <- factor(dat_f$ch)
  ff <- y ~ ch + t2(x, k = 6, bs = "cr")
  m   <- fit_fixed(ff, dat, max_iter = 20, edf = "hutchpp")
  m_f <- fit_fixed(ff, dat_f, max_iter = 20, edf = "hutchpp")
  expect_identical(names(coef(m))[1:3], c("(Intercept)", "chv", "chw"))
  expect_equal(coef(m), coef(m_f))
  expect_equal(m$var.summary, m_f$var.summary)
  # predict() with the data, and with new data as character or factor
  expect_equal(predict(m, newdata = dat), m$linear.predictors,
               ignore_attr = TRUE)
  nd <- data.frame(x = c(0.2, 0.7), ch = c("w", "v"))
  p_chr <- predict(m, newdata = nd)
  nd$ch <- factor(nd$ch, levels = c("u", "v", "w"))
  expect_equal(predict(m, newdata = nd), p_chr)
  expect_equal(predict(m_f, newdata = nd), p_chr)
  expect_s3_class(summary(m), "summary.gbgam")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(m, all.terms = TRUE, pages = 1))
})

test_that("a logical covariate keeps its quantile summary", {
  dat    <- simulate_gaussian_1d(n = 40, seed = 15)
  dat$lg <- rep(c(TRUE, FALSE), 20)
  m <- fit_fixed(y ~ lg + t2(x, k = 6, bs = "cr"), dat, max_iter = 5)
  expect_equal(m$var.summary$lg,
               stats::quantile(as.numeric(dat$lg), c(0.25, 0.5, 0.75)))
  expect_type(m$model$lg, "logical")
})
