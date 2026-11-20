# Finite-difference checks of the tracked derivatives (dmu/dy, deta/dy,
# da/dy), the EDF built from them and their helper functions.


# ---- Trajectory FD tests: dmu/dy, deta/dy ------------------------------------

# Per family, with an offset: a reference fit (50 iterations, tau = 0, exact
# EDF), then central differences on the same selection path for 20 random j
for (fam_spec in list(
  list(label       = "gaussian / identity",
       simulate    = function() simulate_gaussian_1d(n = 30L, seed = 301L),
       family      = stats::gaussian(),
       weights_var = NULL),
  list(label       = "Gamma / log",
       simulate    = function() simulate_gamma_1d(n = 30L, seed = 302L),
       family      = stats::Gamma(link = "log"),
       weights_var = NULL),
  list(label       = "binomial / logit",
       simulate    = function() simulate_binomial_prop(n = 30L, seed = 303L),
       family      = stats::binomial(),
       weights_var = "w"),
  list(label       = "inverse.gaussian / log",
       simulate    = function() simulate_invgaussian_1d(n = 30L, seed = 304L),
       family      = stats::inverse.gaussian(link = "log"),
       weights_var = NULL)
)) {
  local({
    spec <- fam_spec
    test_that(paste("dmu/dy, deta/dy via central FD with an offset:",
                    spec$label), {
      skip_on_cran()
      # Keep the derivative matrices on the returned model (dropped by default)
      op <- options(gbgam.keep_derivs = TRUE)
      on.exit(options(op), add = TRUE)

      dat     <- spec$simulate()
      set.seed(305L)
      dat$o   <- 0.2 * stats::rnorm(nrow(dat))
      weights <- if (!is.null(spec$weights_var)) dat[[spec$weights_var]]
      fml     <- y ~ t2(x, k = 6, bs = "cr") + offset(o)
      y       <- dat$y
      n       <- length(y)
      eps     <- 1e-5

      sc <- gb_stop_ctrl("none", max_iter = 50L)
      ic <- gb_uq_ctrl("none")
      fit_one <- function(d, lc, edf) {
        quiet_fit(fml, data = d, family = spec$family, weights = weights,
                  learner_control = lc, stop_control = sc,
                  edf_control = gb_edf_ctrl(edf), uq_control = ic)
      }

      mod <- fit_one(dat, gb_learner_ctrl(step_size = 0.1, tau = 0), "exact")

      # Fixed-path controls for perturbed models (no EDF needed)
      lc_fixed <- gb_learner_ctrl(step_size = 0.1, tau = 0,
                                  preselected_learner_indices =
                                    mod$metrics$k_opts)
      expect_identical(fit_one(dat, lc_fixed, "none")$linear.predictors,
                       mod$linear.predictors)

      # Proportions of 0 or 1 cannot move both ways (binomial)
      set.seed(42L)
      ok <- which(spec$family$family != "binomial" | (eps < y & y < 1 - eps))
      js <- sample(ok, 20)
      for (j in js) {
        dp <- dat
        dm <- dat
        dp$y[j] <- y[j] + eps
        dm$y[j] <- y[j] - eps
        fp <- fit_one(dp, lc_fixed, "none")
        fm <- fit_one(dm, lc_fixed, "none")

        expect_lt(
          max_rel_err(mod$derivs$dmudy[, j],
                      (fp$fitted.values - fm$fitted.values) / (2 * eps)),
          1e-6, label = paste0("dmu/dy[,", j, "] ", spec$label))
        expect_lt(
          max_rel_err(mod$derivs$detady[, j],
                      (fp$linear.predictors - fm$linear.predictors) /
                        (2 * eps)),
          1e-6, label = paste0("deta/dy[,", j, "] ", spec$label))
      }
    })
  })
}


# ---- Full-Jacobian FD tests --------------------------------------------------

# Offsets, learner step sizes, post-iteration learners,
# initialise_intercept = FALSE and the Hutch++ sketch

# Compares columns `js` of the tracked dmu/dy with central finite differences
# of the fitted values on a fixed selection path (tau = 0, 40 iterations).
fd_jacobian_check <- function(dat, fml, family, weights = NULL,
                              lc_args = list(), js = 1:8, eps = 1e-4,
                              tol = 1e-6) {
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  sc  <- gb_stop_ctrl("none", max_iter = 40L)
  ic  <- gb_uq_ctrl("none")
  lc0 <- c(list(step_size = 0.1, tau = 0), lc_args)
  fit_one <- function(d, lc, edf) {
    suppressMessages(suppressWarnings(gbgam(
      fml, data = d, family = family, weights = weights,
      learner_control = lc, stop_control = sc,
      edf_control = gb_edf_ctrl(edf), uq_control = ic,
      verbose = FALSE)))
  }
  ref <- fit_one(dat, do.call(gb_learner_ctrl, lc0), "exact")
  lc1 <- do.call(gb_learner_ctrl, c(lc0, list(
    preselected_learner_indices = ref$metrics$k_opts)))
  for (j in js) {
    dp <- dat
    dp$y[j] <- dp$y[j] + eps
    dm <- dat
    dm$y[j] <- dm$y[j] - eps
    fd <- (fit_one(dp, lc1, "none")$fitted.values -
             fit_one(dm, lc1, "none")$fitted.values) / (2 * eps)
    expect_equal(unname(ref$derivs$dmudy[, j]), unname(fd), tolerance = tol,
                 label = paste0("dmu/dy[, ", j, "]"))
  }
  invisible(ref)
}

test_that("dmu/dy FD: poisson, Gamma(log) and probit with offset and weights", {
  skip_on_cran()
  set.seed(11)
  n <- 40
  d <- data.frame(x = sort(stats::runif(n)), o = 0.5 * stats::rnorm(n),
                  w = stats::runif(n, 0.5, 2))
  eta <- 0.3 + 0.8 * sin(2 * pi * d$x)
  fml <- y ~ t2(x, k = 6, bs = "cr") + offset(o)

  # Non-integer counts keep y - eps > 0 (also below)
  d$y <- stats::rpois(n, exp(eta + d$o)) + 0.2
  fd_jacobian_check(d, fml, stats::poisson(), weights = d$w)

  d$y <- stats::rgamma(n, shape = 3, rate = 3 / exp(eta + d$o))
  fd_jacobian_check(d, fml, stats::Gamma(link = "log"), weights = d$w)

  d$y <- stats::rbinom(n, 5, stats::pnorm(0.5 * eta + d$o)) / 5
  d$y <- pmin(pmax(d$y, 0.1), 0.9) # keep y +- eps inside [0, 1]
  fd_jacobian_check(d, fml, stats::binomial(link = "probit"),
                    weights = rep(5, n))
})

test_that("dmu/dy FD: learner step size nu_k < 1 is used in the update", {
  skip_on_cran()
  set.seed(12)
  n <- 40
  d <- data.frame(x = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rpois(n, exp(0.3 + sin(2 * pi * d$x))) + 0.2
  # sp = -3: initial learner EDF 3, not 1, so nu_k < 1
  fml <- y ~ t2(x, k = 8, bs = "cr", sp = -3) +
    t2(x2, k = 6, bs = "cr", sp = -2.5)
  ref <- fd_jacobian_check(
    d, fml, stats::poisson(),
    lc_args = list(use_learner_step_size_for_update = TRUE))
  is_dr  <- vapply(ref$learners, inherits, logical(1), "demmler_reinsch")
  nu_ks <- vapply(ref$learners[is_dr], function(lr) lr$nu_k, numeric(1))
  expect_true(all(nu_ks < 1))
})

test_that("dmu/dy FD: post-iteration learners and initialise_intercept = FALSE", {
  skip_on_cran()
  set.seed(13)
  n <- 40
  d <- data.frame(x = stats::runif(n), x2 = stats::runif(n))
  d$y <- stats::rpois(n, exp(0.3 + sin(2 * pi * d$x))) + 0.2
  fml <- y ~ t2(x, k = 6, bs = "cr") + t2(x2, k = 5, bs = "cr")
  fd_jacobian_check(d, fml, stats::poisson(),
                    lc_args = list(learners_to_fit_after_each_iter = 1L,
                                   excluded_learner_indices = 1L))
  ref <- fd_jacobian_check(d, fml, stats::poisson(),
                           lc_args = list(initialise_intercept = FALSE))
  # The tracked EDF equals the trace of the tracked Jacobian
  expect_equal(ref$edfs[ref$mstop + 1L], sum(diag(ref$derivs$dmudy)),
               tolerance = 1e-10)
})

test_that("hutchpp sketch equals the exact Jacobian times the probe matrix", {
  skip_on_cran()
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  set.seed(14)
  n <- 80
  d <- data.frame(x = stats::runif(n), x2 = stats::runif(n),
                  o = 0.3 * stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + sin(2 * pi * d$x) + d$o))
  fit_ec <- function(ec) {
    set.seed(1)
    gbgam(y ~ t2(x, k = 8, bs = "cr") + t2(x2, k = 6, bs = "cr") + offset(o),
          data = d, family = stats::poisson(),
          learner_control = gb_learner_ctrl(step_size = 0.1, tau = 0),
          stop_control = gb_stop_ctrl("none", max_iter = 30L),
          edf_control = ec, uq_control = gb_uq_ctrl("none"),
          verbose = FALSE)
  }
  ex <- fit_ec(gb_edf_ctrl("exact"))
  hp <- fit_ec(gb_edf_ctrl("hutchpp", z = 40))
  expect_identical(ex$metrics$k_opts, hp$metrics$k_opts)
  expect_equal(ex$derivs$dmudy %*% hp$derivs$GGG, hp$derivs$dmudy,
               tolerance = 1e-10)
  expect_equal(ex$derivs$detady %*% hp$derivs$GGG, hp$derivs$detady,
               tolerance = 1e-10)
})

# Family with the derivatives used by the EDF code
fix_family <- function(fam) {
  mgcv::fix.family.var(mgcv::fix.family.link(fam))
}

# Central finite difference of the intercept MLE with respect to y[j]
fd_dady <- function(y, fam, w, off, j, h = 1e-5) {
  yp <- y
  ym <- y
  yp[j] <- yp[j] + h
  ym[j] <- ym[j] - h
  (find_intercept_mle(yp, fam, w, off) -
     find_intercept_mle(ym, fam, w, off)) / (2 * h)
}


# ---- Intercept derivative, offsets and learner step sizes --------------------

test_that("da/dy of the MLE intercept matches finite differences with offsets", {
  set.seed(101)
  n   <- 30
  off <- 0.5 * stats::rnorm(n)
  w   <- stats::runif(n, 0.5, 2)
  eta <- 0.3 + off
  specs <- list(
    poisson = list(fam = stats::poisson(),
                   y = stats::rpois(n, exp(eta)) + 0.1),
    gamma_log = list(fam = stats::Gamma(link = "log"),
                     y = stats::rgamma(n, shape = 3, rate = 3 / exp(eta))),
    probit = list(fam = stats::binomial(link = "probit"),
                  y = stats::rbinom(n, 10, stats::pnorm(eta)) / 10),
    logit = list(fam = stats::binomial(),
                 y = stats::rbinom(n, 10, stats::plogis(eta)) / 10)
  )
  for (nm in names(specs)) {
    fam <- fix_family(specs[[nm]]$fam)
    y   <- specs[[nm]]$y
    a   <- find_intercept_mle(y, fam, w, off)
    d   <- compute_da_dy_for_mle_a(fam, a + off, w, y = y)
    fd  <- vapply(1:6, function(j) fd_dady(y, fam, w, off, j), numeric(1))
    expect_equal(d[1:6], fd, tolerance = 1e-6, label = nm)
  }
  # Canonical link: tr(dmu0/dy) = sum_i h'(eta_i) da/dy_i = 1 exactly
  fam <- fix_family(stats::poisson())
  y   <- specs$poisson$y
  a   <- find_intercept_mle(y, fam, w, off)
  d   <- compute_da_dy_for_mle_a(fam, a + off, w, y = y)
  expect_equal(sum(fam$mu.eta(a + off) * d), 1, tolerance = 1e-12)
})

test_that("poisson with an offset: tracked EDF equals the FD trace", {
  skip_on_cran()
  set.seed(102)
  n <- 30
  d <- data.frame(x = stats::runif(n), o = log(stats::runif(n, 1, 20)))
  d$y <- stats::rpois(n, exp(-1 + sin(2 * pi * d$x) + d$o)) + 0.5
  fit_one <- function(dd, lc, edf) {
    suppressWarnings(quiet_fit(
      y ~ t2(x, k = 6, bs = "cr") + offset(o), data = dd,
      family = stats::poisson(), learner_control = lc,
      stop_control = gb_stop_ctrl("none", max_iter = 1L),
      edf_control = gb_edf_ctrl(edf)))
  }
  ref <- fit_one(d, gb_learner_ctrl(step_size = 0.1, tau = 0), "exact")
  lc1 <- gb_learner_ctrl(step_size = 0.1, tau = 0,
                         preselected_learner_indices = ref$metrics$k_opts)
  h  <- 1e-5
  tr <- sum(vapply(seq_len(n), function(j) {
    dp <- d
    dm <- d
    dp$y[j] <- dp$y[j] + h
    dm$y[j] <- dm$y[j] - h
    (fit_one(dp, lc1, "none")$fitted.values[j] -
       fit_one(dm, lc1, "none")$fitted.values[j]) / (2 * h)
  }, numeric(1)))
  expect_equal(ref$edf.total, tr, tolerance = 1e-6)
  # One boosting step with a positive definite smoother cannot remove EDF
  expect_true(all(vapply(ref$learners, `[[`, numeric(1), "edf") >= 0))
})

test_that("use_learner_step_size_for_update: gaussian EDF equals the FD trace", {
  skip_on_cran()
  set.seed(103)
  n <- 40
  d <- data.frame(x = stats::runif(n), z = stats::runif(n))
  d$y <- sin(2 * pi * d$x) + d$z + stats::rnorm(n, sd = 0.3)
  fml <- y ~ t2(x, k = 8, bs = "cr", sp = -3) + t2(z, k = 6, bs = "cr", sp = -2)
  lc  <- gb_learner_ctrl(step_size = 0.1, tau = 0,
                         use_learner_step_size_for_update = TRUE)
  sc  <- gb_stop_ctrl("none", max_iter = 30L)
  m   <- quiet_fit(fml, data = d, learner_control = lc, stop_control = sc,
                   edf_control = gb_edf_ctrl("exact"))
  lc1 <- gb_learner_ctrl(step_size = 0.1, tau = 0,
                         use_learner_step_size_for_update = TRUE,
                         preselected_learner_indices = m$metrics$k_opts)
  # Fitted values are linear in y on a fixed path: unit perturbations are exact
  tr <- sum(vapply(seq_len(n), function(j) {
    dp <- d
    dp$y[j] <- dp$y[j] + 1
    quiet_fit(fml, data = dp, learner_control = lc1, stop_control = sc,
              edf_control = gb_edf_ctrl("none"))$fitted.values[j] -
      m$fitted.values[j]
  }, numeric(1)))
  expect_equal(m$edf.total, tr, tolerance = 1e-8)
})


# ---- Derivative helpers and stored matrices ----------------------------------

test_that("weighted du_dy and delta_dgamma_dy match the naive formulas", {
  set.seed(115)
  n   <- 12
  fam <- fix_family(stats::poisson())
  mu  <- exp(stats::rnorm(n))
  y   <- stats::rpois(n, mu)
  w   <- stats::runif(n, 0.5, 2)
  J   <- matrix(stats::rnorm(n * n, sd = 0.1), n, n)
  G   <- rrademacher(n, 6)
  om  <- omega(mu, fam)
  sc  <- domega_dmu(mu, fam) * (y - mu) - om
  naive_exact  <- w * (sc * J + diag(om))
  naive_sketch <- w * (sc * (J %*% G) + om * G)
  expect_equal(du_dy(J, mu, y, fam, GGG = NULL, weights = w), naive_exact,
               tolerance = 1e-14)
  expect_equal(du_dy(J %*% G, mu, y, fam, GGG = G, weights = w),
               naive_sketch, tolerance = 1e-14)
  # A penalised learner as built by new_learners(), at its initial rho
  lb <- new_learners(y ~ t2(x, k = 5, bs = "cr"),
                     data = data.frame(y = y, x = stats::runif(n)))
  lr <- Filter(\(l) inherits(l, "demmler_reinsch"), lb$learners)[[1L]]
  step_F <- 0.1 / (1 + lr$D * exp(lr$rho))
  expect_equal(delta_dgamma_dy(lr, naive_exact, step_size = 0.1),
               (step_F * t(lr$design)) %*% naive_exact, tolerance = 1e-14)
  # Identity link: dmu/dy is deta/dy itself
  expect_identical(dmu_dy(rep(0, n), J, fix_family(stats::gaussian())),
                   J)
})

test_that("derivative matrices are dropped unless gbgam.keep_derivs is set", {
  dat <- simulate_gaussian_1d(n = 60, seed = 113)
  fit_d <- function() {
    quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = dat,
              stop_control = gb_stop_ctrl("none", max_iter = 10L),
              edf_control = gb_edf_ctrl("exact"))
  }
  expect_null(fit_d()$derivs)
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  expect_equal(dim(fit_d()$derivs$dmudy), c(60L, 60L))
})

test_that("non-canonical link message is well formed and can be silenced", {
  dat <- simulate_gamma_1d(n = 60, seed = 114)
  expect_message(
    m <- gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat,
               family = stats::Gamma(link = "log"),
               stop_control = gb_stop_ctrl("none", max_iter = 2L),
               edf_control = gb_edf_ctrl("exact"), verbose = FALSE),
    "(Gamma/log). Evaluating EDF", fixed = TRUE
  )
  m$edf_compute_method <- "exact"
  expect_silent(init_derivs(m, nrow(dat), quiet = TRUE))
})

test_that("Hutch++ with a non-canonical link results in a warning", {
  dat <- simulate_gamma_1d(n = 60, seed = 114)
  fit <- function(family, edf = "hutchpp") {
    gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat, family = family,
          stop_control = gb_stop_ctrl("none", max_iter = 2L),
          edf_control = gb_edf_ctrl(edf, z = if (edf == "hutchpp") 12L),
          verbose = FALSE)
  }
  expect_warning(
    suppressMessages(m <- fit(stats::Gamma(link = "log"))),
    "Hutch++ EDF are unreliable with the non-canonical link Gamma/log",
    fixed = TRUE
  )
  # Neither the canonical link nor exact tracking gives it, nor a refit
  expect_no_warning(fit(stats::Gamma(link = "inverse")))
  expect_no_warning(suppressMessages(fit(stats::Gamma(link = "log"), "exact")))
  expect_no_warning(init_derivs(m, nrow(dat), quiet = TRUE))
})
