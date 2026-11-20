# Uncertainty quantification requested from gbgam() (uq_control).


# ---- Bootstrap ---------------------------------------------------------------

test_that("uq_control='bootstrap' attaches Vp and boot_coef (Gaussian)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 70)
  # edf "none" + gaussian: the bootstrap's dispersion warning is expected
  expect_warning(
    mod <- gbgam(
      y ~ t2(x, k = 6, bs = "cr"),
      data            = dat,
      learner_control = gb_learner_ctrl(),
      stop_control    = gb_stop_ctrl("none", max_iter = 15),
      edf_control     = gb_edf_ctrl("none"),
      uq_control      = gb_uq_ctrl(method = "bootstrap",
                                            n_samples = 10L, alpha = NULL),
      verbose           = FALSE
    ),
    "edf_compute_method"
  )
  expect_s3_class(mod, "gbgam")
  p <- length(mod$coefficients)
  expect_true(!is.null(mod$Vp))
  expect_true(!is.null(mod$boot_coef))
  expect_equal(dim(mod$Vp), c(p, p))
  expect_equal(dim(mod$boot_coef), c(10L, p))
})

test_that("uq_control='bootstrap' attaches Vp and boot_coef (Poisson)", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 60, seed = 70)
  mod <- suppressWarnings(gbgam(
    y ~ t2(x, k = 6, bs = "cr"),
    data            = dat,
    family          = stats::poisson(),
    learner_control = gb_learner_ctrl(),
    stop_control    = gb_stop_ctrl("none", max_iter = 15),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl(method = "bootstrap", n_samples = 10L,
                                          alpha = NULL),
    verbose           = FALSE
  ))
  p <- length(mod$coefficients)
  expect_true(!is.null(mod$Vp))
  expect_equal(dim(mod$Vp), c(p, p))
})

test_that("uq_control='bootstrap' attaches Vp and boot_coef (non-canonical log link: Gamma and inverse.gaussian)", {
  # Gamma and inverse.gaussian share the same structural contract; group them.
  skip_on_cran()
  specs <- list(
    list(family = stats::Gamma(link = "log"),            sim_fn = simulate_gamma_1d),
    list(family = stats::inverse.gaussian(link = "log"), sim_fn = simulate_invgaussian_1d)
  )
  for (sp in specs) {
    dat <- sp$sim_fn(n = 60, seed = 70)
    mod <- suppressWarnings(gbgam(
      y ~ t2(x, k = 6, bs = "cr"),
      data            = dat,
      family          = sp$family,
      learner_control = gb_learner_ctrl(),
      stop_control    = gb_stop_ctrl("none", max_iter = 15),
      edf_control     = gb_edf_ctrl("none"),
      uq_control      = gb_uq_ctrl(method = "bootstrap", n_samples = 10L,
                                            alpha = NULL),
      verbose           = FALSE
    ))
    p <- length(mod$coefficients)
    expect_true(!is.null(mod$Vp),         label = mod$family$family)
    expect_equal(dim(mod$Vp), c(p, p),    label = mod$family$family)
  }
})

test_that("predict with se.fit=TRUE after bootstrap returns finite SEs (Gaussian)", {
  skip_on_cran()
  mod  <- local_gaussian_mod()
  # local_gaussian_mod() uses edf "none": the dispersion warning is expected
  expect_warning(
    mod <- bootstrap_gbgam(mod, n_samples = 15, alpha = NULL, verbose = FALSE),
    "edf_compute_method"
  )
  pred <- predict(mod, se.fit = TRUE)
  expect_true(all(is.finite(pred$fit)))
  expect_true(all(is.finite(pred$se.fit)))
  expect_true(all(pred$se.fit >= 0))
})

test_that("predict with se.fit=TRUE after bootstrap returns finite SEs (Poisson)", {
  skip_on_cran()
  mod  <- local_poisson_mod()
  mod  <- suppressWarnings(bootstrap_gbgam(mod, n_samples = 15, alpha = NULL, verbose = FALSE))
  pred <- predict(mod, se.fit = TRUE)
  expect_true(all(is.finite(pred$fit)))
  expect_true(all(is.finite(pred$se.fit)))
  expect_true(all(pred$se.fit >= 0))
})

test_that("documented Poisson example fits a gbgam model with bootstrap Vp", {
  skip_on_cran()
  set.seed(6)
  dat_pois <- gamSim(1, n = 400, dist = "poisson", scale = .1, verbose = FALSE)

  m_pois <- gbgam(
    y ~ t2(x0, bs = "cr", k = 10) +
      t2(x1, bs = "cr", k = 10) +
      t2(x2, bs = "cr", k = 10) +
      t2(x3, bs = "cr", k = 10),
    data            = dat_pois,
    family          = poisson(),
    learner_control = gb_learner_ctrl(step_size = 0.1),
    uq_control      = gb_uq_ctrl("bootstrap", n_samples = 20),
    verbose         = FALSE
  )

  expect_model_contract(m_pois, nrow(dat_pois))
  expect_identical(m_pois$family$family, "poisson")
  expect_false(is.null(m_pois$Vp))
  expect_equal(dim(m_pois$Vp), rep(length(stats::coef(m_pois)), 2))
})


# ---- Bayesian backfitting ----------------------------------------------------

test_that("uq_control='bayes' fits with or without edf computation (Gaussian)", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 60, seed = 71)

  # Stage 1 of boost to backfit uses the boosted fit only
  mod0 <- gbgam(
    y ~ t2(x, k = 6, bs = "cr"), data = dat,
    learner_control = gb_learner_ctrl(),
    stop_control    = gb_stop_ctrl("none", max_iter = 15),
    edf_control     = gb_edf_ctrl("none"),
    uq_control      = gb_uq_ctrl(method = "bayes", n_samples = 5L,
                                          n_burnin = 2L),
    verbose           = FALSE
  )
  expect_true(all(is.finite(mod0$Vp)))

  mod <- gbgam(
    y ~ t2(x, k = 6, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(),
    stop_control    = gb_stop_ctrl("none", max_iter = 15),
    edf_control     = gb_edf_ctrl("hutchpp", z = 20L),
    uq_control      = gb_uq_ctrl(method = "bayes", n_samples = 50L,
                                          n_burnin = 20L),
    verbose           = FALSE
  )
  expect_s3_class(mod, "gbgam")
  p <- length(mod$coefficients)
  expect_true(!is.null(mod$Vp))
  expect_equal(dim(mod$Vp), c(p, p))
  expect_true(all(is.finite(mod$Vp)))
})

test_that("gb_bayesian_ctrl() in gbgam() equals bayesian_backfitting()", {
  dat  <- simulate_poisson_1d(n = 80, seed = 203)
  ctrl <- gb_bayesian_ctrl(n_samples = 20L, n_burnin = 5L, linear = "kept",
                           undersmooth = TRUE, edf_floor = 3)
  fit  <- function(ic) {
    gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat, family = stats::poisson(),
          stop_control = gb_stop_ctrl("none", max_iter = 30L),
          edf_control = gb_edf_ctrl("exact"), uq_control = ic,
          verbose = FALSE)
  }
  # the same random numbers: the posterior is drawn right after the fit
  set.seed(9)
  inside <- fit(ctrl)
  set.seed(9)
  plain  <- fit(gb_uq_ctrl("none"))
  after  <- do.call(bayesian_backfitting, c(list(plain), unclass(ctrl)))
  expect_false(is.null(after$bb_info$undersmoothed))  # the floor lifts one
  for (field in c("Vp", "coefficients_bayes", "beta_samples", "bb_info",
                  "mh_accept_rate", "uq_control"))
    expect_identical(inside[[field]], after[[field]], label = field)
})

test_that("Bayes Gibbs (Gaussian) sets mh_accept_rate to 1 for all learners", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 201)
  set.seed(1)
  mod <- gbgam(
    y ~ t2(x, k = 6, bs = "cr"),
    data            = dat,
    learner_control = gb_learner_ctrl(step_size = 0.05, tau = 0),
    stop_control    = gb_stop_ctrl("none", max_iter = 20L),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("bayes", n_samples = 50L,
                                          n_burnin = 10L),
    verbose           = FALSE
  )

  # Gibbs sampler always accepts; rate should be 1 for every learner
  expect_true(!is.null(mod$mh_accept_rate))
  expect_true(is.numeric(mod$mh_accept_rate))
  expect_equal(length(mod$mh_accept_rate), length(mod$learners))
  expect_true(all(mod$mh_accept_rate == 1))
})

test_that("Bayes MH (Poisson) stores plausible per-learner acceptance rates", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 100, seed = 202)
  set.seed(2)
  mod <- gbgam(
    y ~ t2(x, k = 5, bs = "cr"),
    data            = dat,
    family          = stats::poisson(),
    learner_control = gb_learner_ctrl(step_size = 0.05, tau = 0),
    stop_control    = gb_stop_ctrl("none", max_iter = 20L),
    edf_control     = gb_edf_ctrl("exact"),
    uq_control      = gb_uq_ctrl("bayes", n_samples = 100L,
                                          n_burnin = 20L),
    verbose           = FALSE
  )

  # MH acceptance rates are in [0, 1] for every sampled learner (NA for
  # learners never updated by boosting, which are held fixed)
  expect_true(!is.null(mod$mh_accept_rate))
  expect_true(is.numeric(mod$mh_accept_rate))
  expect_equal(length(mod$mh_accept_rate), length(mod$learners))
  ar <- mod$mh_accept_rate[!is.na(mod$mh_accept_rate)]
  expect_true(length(ar) >= 1L)
  expect_true(all(ar >= 0 & ar <= 1))
  # Every sampled learner mixes (a stuck sampler has rate ~0)
  expect_true(all(ar > 0.05),
              info = paste("accept rates:",
                           paste(round(mod$mh_accept_rate, 3), collapse = " ")))
  se <- predict(mod, type = "link", se.fit = TRUE)$se.fit
  expect_true(all(is.finite(se) & se > 0))
})


# ---- Default, time caps and failures -----------------------------------------

test_that("uq_control takes a method's own control directly", {
  dat <- simulate_gaussian_1d(n = 80, seed = 31)
  fit <- function(ic) {
    set.seed(5)
    # 4 replicates: factors below 0.9 are Monte Carlo noise
    without_low_widening_warning(
      gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat,
            stop_control = gb_stop_ctrl(max_iter = 60),
            uq_control = ic, verbose = FALSE))
  }
  pb <- gb_param_boot_ctrl(n_samples = 4L, correction = "smooth")
  b_direct <- fit(pb)
  expect_identical(b_direct$uq_control, gb_uq_ctrl(control = pb))
  expect_identical(b_direct$Vp, fit(gb_uq_ctrl(control = pb))$Vp)
  expect_identical(fit(list(pb))$Vp, b_direct$Vp)
  expect_identical(fit(list(boot = pb))$Vp, b_direct$Vp)
  bc <- gb_bayesian_ctrl(n_samples = 5L, n_burnin = 2L)
  bb_direct <- fit(bc)
  expect_identical(bb_direct$uq_control$method, "bayes")
  expect_identical(bb_direct$Vp, fit(gb_uq_ctrl(control = bc))$Vp)
  # update() reruns the method the model carries
  expect_identical(without_low_widening_warning(update(b_direct))$uq_control,
                   b_direct$uq_control)
})

test_that("uq_control with both methods' controls is an error", {
  dat <- simulate_gaussian_1d(n = 60, seed = 32)
  both <- list(gb_param_boot_ctrl(), gb_bayesian_ctrl())
  for (ic in list(both, stats::setNames(both, c("boot", "bayes"))))
    expect_error(
      quiet_fit(y ~ t2(x, k = 5), data = dat, uq_control = ic),
      "a model holds one covariance.*bootstrap_gbgam\\(\\) or bayesian")
  expect_error(quiet_fit(y ~ t2(x, k = 5), data = dat,
                         uq_control = "bootstrap"),
               "must be a list")
})

test_that("the default uncertainty quantification is none", {
  expect_identical(gb_uq_ctrl()$method, "none")
})

test_that("a failing uncertainty quantification keeps the fitted model", {
  d <- core_data()
  local_mocked_bindings(bootstrap_gbgam = function(...) stop("replicate 2"))
  e <- tryCatch(
    gbgam(ff_core, data = d, verbose = FALSE,
          stop_control = gb_stop_ctrl("none", max_iter = 10),
          uq_control = gb_uq_ctrl("bootstrap", n_samples = 3)),
    gbgam_uq_error = function(e) e)
  expect_s3_class(e, "gbgam_uq_error")
  expect_match(conditionMessage(e), "failed after the main fit: replicate 2")
  expect_s3_class(e$fit, "gbgam")
  expect_equal(e$fit$mstop, 10)
  expect_null(e$fit$Vp)
  expect_identical(e$fit$uq_control, gb_uq_ctrl("none"))
})

test_that("a used-up time limit leaves the uncertainty quantification working", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 150, seed = 112)
  fit_capped <- function(ic) {
    quiet_fit(y ~ t2(x, k = 8, bs = "cr"), data = dat,
              stop_control = gb_stop_ctrl(max_runtime_mins = 1e-6),
              uq_control = ic)
  }
  expect_warning(
    m <- fit_capped(gb_uq_ctrl("bayes", n_samples = 5, n_burnin = 1)),
    "time limit"
  )
  expect_s3_class(m, "gbgam")
  # The bootstrap warns about the time limit and about returning without Vp
  expect_warning(
    expect_warning(
      m <- fit_capped(gb_uq_ctrl("bootstrap", n_samples = 2)),
      "time limit"
    ),
    "fewer than 2 replicates"
  )
  expect_s3_class(m, "gbgam")
  expect_null(m$Vp)
  # $uq_control describes the $Vp the model carries: none
  expect_identical(m$uq_control, gb_uq_ctrl("none"))
})

test_that("apply_uq() passes every stored setting on", {
  captured <- new.env()
  testthat::local_mocked_bindings(
    bootstrap_gbgam = function(mod, ...) {
      captured$boot <- list(...)
      mod
    },
    bayesian_backfitting = function(mod, ...) {
      captured$bayes <- list(...)
      mod
    },
    .package = "gbgam"
  )
  pb <- gb_param_boot_ctrl(n_samples = 7L, mode = "free_stopping",
                           alpha = NULL, correction = "learner", phi = 2,
                           max_iter = 30L)
  apply_uq(list(), gb_uq_ctrl(control = pb), verbose = TRUE,
           max_runtime_mins = 5)
  expect_identical(captured$boot[names(pb)], unclass(pb)[names(pb)])
  expect_identical(captured$boot$max_runtime_mins, 5)
  expect_true(captured$boot$verbose)
  bc <- gb_bayesian_ctrl(n_samples = 9L, n_burnin = 3L, linear = "kept",
                         undersmooth = TRUE, edf_floor = 2)
  apply_uq(list(), gb_uq_ctrl(control = bc), verbose = FALSE,
           max_runtime_mins = NULL)
  expect_identical(captured$bayes[names(bc)], unclass(bc))
})
