# smooth_intervals(): intervals for the intercept plus each smooth.


# ---- Intervals ---------------------------------------------------------------

# The AIC fit, its posterior and its bootstrap, made once on first use (not
# at load time, so that skip_on_cran() skips them too)
si_fits <- local({
  fits <- NULL
  function() {
    if (is.null(fits)) {
      mod <- fit_aic(f_col, sim_collinear())
      set.seed(10)
      bb <- bayesian_backfitting(mod, n_samples = 200L, n_burnin = 50L)
      set.seed(11)
      bt <- bootstrap_gbgam(mod, n_samples = 20L, correction = "smooth")
      fits <<- list(mod = mod, bb = bb, bt = bt)
    }
    fits
  }
})

test_that("smooth_intervals centres on the intercept plus each smooth", {
  skip_on_cran()
  bb <- si_fits()$bb
  ci <- smooth_intervals(bb, alpha = 0.1)
  p  <- stats::predict(bb, type = "iterms", se.fit = TRUE)
  a  <- unname(bb$coefficients[["(Intercept)"]])
  expect_equal(c(ci$fit), c(p$fit) + a)
  expect_equal(c(ci$se), c(p$se.fit))
  expect_equal(c(ci$upper - ci$fit), stats::qnorm(0.95) * c(ci$se))
})

test_that("the standard errors include the intercept's uncertainty", {
  skip_on_cran()
  bb <- si_fits()$bb
  ci <- smooth_intervals(bb, alpha = 0.1)
  # se(a + f)^2 = se(f)^2 + Var(a) + 2 Cov(a, f)
  lab <- "t2(x1)"
  pt <- stats::predict(bb, type = "terms", se.fit = TRUE)
  X  <- stats::predict(bb, type = "lpmatrix")
  sm <- bb$smooth[[which(vapply(bb$smooth, `[[`, "", "label") == lab)]]
  i  <- sm$first.para:sm$last.para
  v  <- pt$se.fit[, lab]^2 + bb$Vp[1L, 1L] +
    2 * as.vector(X[, i] %*% bb$Vp[i, 1L])
  expect_equal(unname(ci$se[, lab]^2), unname(v), tolerance = 1e-8)
})

test_that("centre = \"posterior\" centres on the posterior mean", {
  skip_on_cran()
  bb <- si_fits()$bb
  cp <- smooth_intervals(bb, centre = "posterior")
  a_b <- unname(bb$coefficients_bayes[["(Intercept)"]])
  m2 <- bb
  m2$coefficients <- bb$coefficients_bayes
  expect_equal(c(cp$fit), c(stats::predict(m2, type = "iterms")) + a_b)
})

test_that("smooth_intervals needs $Vp, and a posterior for its centre", {
  skip_on_cran()
  fits <- si_fits()
  expect_error(smooth_intervals(fits$mod), "no \\$Vp")
  expect_error(smooth_intervals(fits$bt, centre = "posterior"), "Bayesian")
})

test_that("bootstrap intervals widen with the level and carry the factors", {
  skip_on_cran()
  bt <- si_fits()$bt
  width <- function(alpha) {
    ci <- smooth_intervals(bt, alpha = alpha)
    colMeans(ci$upper - ci$lower)
  }
  expect_true(all(width(0.01) > width(0.10)))
  expect_named(smooth_intervals(bt)$factor, names(bt$widening_factor))
})

test_that("alpha is the miscoverage, and a level given as alpha is an error", {
  skip_on_cran()
  bb <- si_fits()$bb
  ci <- smooth_intervals(bb, alpha = 0.1)
  expect_identical(ci$alpha, 0.1)
  expect_equal(c(ci$upper - ci$fit), stats::qnorm(0.95) * c(ci$se))
  expect_identical(smooth_intervals(bb), smooth_intervals(bb, alpha = 0.05))
  expect_identical(smooth_intervals(bb, 0.1), ci)
  expect_error(smooth_intervals(bb, 0.95), "alpha is the miscoverage")
  expect_error(smooth_intervals(bb, alpha = 0), "positive")
  expect_error(smooth_intervals(bb, alpha = 1.5))
})

test_that("correction = chooses the widening of a bootstrap fit", {
  skip_on_cran()
  fits <- si_fits()
  bt   <- fits$bt
  # The stored correction ("smooth") is the default
  expect_identical(smooth_intervals(bt, alpha = 0.1, correction = "smooth"),
                   smooth_intervals(bt, alpha = 0.1))
  # Without widening: plus or minus z standard errors from the plain
  # covariance of the replicates
  plain <- bt
  plain$Vp <- stats::cov(bt$boot_coef)
  plain$uq_control <- gb_uq_ctrl("bayes")
  ci0 <- smooth_intervals(bt, alpha = 0.1, correction = "none")
  expect_null(ci0$factor)
  expect_equal(ci0[c("lower", "upper")],
               smooth_intervals(plain, alpha = 0.1)[c("lower", "upper")])
  # A bootstrap run without widening is widened on request, as one run with
  # that correction
  set.seed(11)
  bt0 <- bootstrap_gbgam(fits$mod, n_samples = 20L, alpha = NULL)
  expect_null(smooth_intervals(bt0, alpha = 0.1)$factor)
  for (corr in c("model", "smooth", "learner"))
    expect_identical(smooth_intervals(bt0, alpha = 0.1, correction = corr)[
                       c("lower", "upper", "factor")],
                     smooth_intervals(bt, alpha = 0.1, correction = corr)[
                       c("lower", "upper", "factor")])
  expect_error(smooth_intervals(fits$bb, correction = "model"),
               "needs a bootstrap fit")
  expect_identical(smooth_intervals(fits$bb, correction = "none"),
                   smooth_intervals(fits$bb))
  expect_error(smooth_intervals(bt, correction = "both"), "correction")
})

test_that("smooth_intervals includes the intercept for by-variable smooths", {
  skip_on_cran()
  set.seed(3)
  n <- 200
  d <- data.frame(x = stats::runif(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- sin(2 * pi * d$x) + (d$g == "b") * d$x + stats::rnorm(n, 0, 0.3)
  mod <- fit_aic(y ~ t2(x, bs = "ps", k = 8) + t2(x, by = g, bs = "ps", k = 8),
                 d)
  set.seed(4)
  bb <- bayesian_backfitting(mod, n_samples = 200L, n_burnin = 20L)
  ci <- smooth_intervals(bb)
  a  <- unname(bb$coefficients[["(Intercept)"]])
  # outside its level a by-factor smooth is zero: a + f is a, with the
  # intercept's standard error
  off <- d$g == "a"
  expect_equal(unname(ci$fit[off, "t2(x):gb"]), rep(a, sum(off)))
  expect_equal(unname(ci$se[off, "t2(x):gb"]),
               rep(sqrt(bb$Vp[1L, 1L]), sum(off)))
  expect_true(all(ci$se > 0))
})
