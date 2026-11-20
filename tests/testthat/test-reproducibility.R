# Reproducibility: the same set.seed() gives identical fits (no dependence on
# earlier RNG use, floating-point ordering or state left from earlier fits).

# Fit twice from the same seed
fit_twice_same_seed <- function(fitter, seed = 42L) {
  set.seed(seed); mod1 <- fitter()
  set.seed(seed); mod2 <- fitter()
  list(mod1 = mod1, mod2 = mod2)
}


# ---- Gaussian: AIC stopping, exact EDF ---------------------------------------

test_that("Gaussian AIC stopping is reproducible with the same seed", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 120, seed = 7)

  fitter <- function() {
    gbgam(
      y ~ t2(x, k = 10, bs = "cr"),
      data              = dat,
      family            = stats::gaussian(),
      learner_control   = gb_learner_ctrl(step_size = 0.05),
      stop_control      = gb_stop_ctrl(stopping_method = "aic", iter_no_improv = 5),
      edf_control       = gb_edf_ctrl("exact"),
      verbose           = FALSE
    )
  }

  r <- fit_twice_same_seed(fitter, seed = 42L)
  expect_equal(r$mod1$mstop,        r$mod2$mstop)
  expect_equal(r$mod1$coefficients, r$mod2$coefficients)
  expect_equal(r$mod1$metrics$ls,   r$mod2$metrics$ls)
  expect_equal(r$mod1$metrics$k_opts, r$mod2$metrics$k_opts)
})


# ---- Poisson: AIC stopping, Hutch++ EDF --------------------------------------

# Hutch++ draws random sketch vectors, so the EDF path and the AIC stop
# depend on the RNG: a fit that did not follow the seed would differ here

test_that("Poisson AIC stopping with Hutch++ EDF is reproducible with the same seed", {
  skip_on_cran()
  dat <- simulate_poisson_1d(n = 150, seed = 8)

  fitter <- function() {
    gbgam(
      y ~ t2(x, k = 10, bs = "cr"),
      data              = dat,
      family            = stats::poisson(),
      learner_control   = gb_learner_ctrl(step_size = 0.05),
      stop_control      = gb_stop_ctrl(stopping_method = "aic", iter_no_improv = 5),
      edf_control       = gb_edf_ctrl(edf_compute_method = "hutchpp", z = 40),
      verbose           = FALSE
    )
  }

  r <- fit_twice_same_seed(fitter, seed = 43L)
  expect_equal(r$mod1$mstop,        r$mod2$mstop)
  expect_equal(r$mod1$coefficients, r$mod2$coefficients)
  expect_equal(r$mod1$metrics$ls,   r$mod2$metrics$ls)
  expect_equal(r$mod1$metrics$k_opts, r$mod2$metrics$k_opts)
})


# ---- bootstrap_gbgam(): Vp and boot_coef -------------------------------------

test_that("bootstrap_gbgam is reproducible with the same seed", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 9)

  base_mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, max_iter = 30,
                        step_size = 0.1)

  # edf "none" + gaussian: the dispersion warning is expected
  boot <- function() {
    set.seed(99L)
    bootstrap_gbgam(base_mod, n_samples = 20, alpha = NULL)
  }
  expect_warning(b1 <- boot(), "edf_compute_method")
  expect_warning(b2 <- boot(), "edf_compute_method")

  expect_equal(b1$boot_coef, b2$boot_coef)
  expect_equal(b1$Vp,        b2$Vp)
})


# ---- gbgam() with the bootstrap ----------------------------------------------

test_that("gbgam() with the bootstrap is reproducible end-to-end", {
  skip_on_cran()
  dat <- simulate_gaussian_1d(n = 80, seed = 10)

  # edf "none" + gaussian: the bootstrap's dispersion warning is expected
  fitter <- function() {
    expect_warning(
      m <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), dat, max_iter = 20,
                     step_size = 0.1,
                     uq_control = gb_uq_ctrl("bootstrap",
                                             n_samples = 10L)),
      "edf_compute_method"
    )
    m
  }

  r <- fit_twice_same_seed(fitter, seed = 77L)
  expect_equal(r$mod1$mstop,        r$mod2$mstop)
  expect_equal(r$mod1$coefficients, r$mod2$coefficients)
  expect_equal(r$mod1$boot_coef,    r$mod2$boot_coef)
  expect_equal(r$mod1$Vp,           r$mod2$Vp)
})
