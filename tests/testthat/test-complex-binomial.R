# A complex binomial model: seed-germination trials on 260 trays, counts
# cbind(germinated, failed) with 5 to 25 seeds per tray and 6 trays with no
# viable seed (zero trials). The truth, on the cloglog scale, has a site x
# season interaction, an optimum in soil moisture, a cyclic day-of-year
# effect and a seed-lot age decline that differs by storage.
#
# Covers: the non-canonical cloglog and probit links; unpenalised learners
# from factor main effects, a factor:factor interaction, poly() and the
# null spaces of factor-by smooths; a cyclic t2() with knots; factor-by
# smooths with a starting sp; zero-trial rows; deviance selection with the step
# size; tau reset to 0 without EDF tracking; the fixed-mode bootstrap with
# the per-learner correction; extend_gbgam() with the bootstrap rerun;
# predict(se.fit = TRUE); exact EDF with ASR and AIC stopping (probit).


# Seed-germination trays; `eta` is the true cloglog linear predictor
sim_germination <- function(n = 260L) {
  d <- data.frame(
    site    = factor(sample(c("north", "south"), n, TRUE)),
    season  = factor(sample(c("spring", "summer", "autumn"), n, TRUE)),
    storage = factor(sample(c("cold", "ambient"), n, TRUE)),
    moist   = stats::runif(n, 0.1, 0.6),
    doy     = stats::runif(n, 0, 365),
    age     = stats::runif(n, 0, 5),
    trials  = sample(5:25, n, TRUE)
  )
  d$trials[sample(n, 6L)] <- 0L
  d$eta <- -0.5 +
    0.4 * (d$site == "south") * (d$season == "summer") -
    0.3 * (d$season == "autumn") +
    3 * (d$moist - 0.35) - 12 * (d$moist - 0.35)^2 +
    0.5 * sin(2 * pi * d$doy / 365) -
    ifelse(d$storage == "cold", 0.1, 0.35) * d$age
  d$germ <- stats::rbinom(n, d$trials, 1 - exp(-exp(d$eta)))
  d$fail <- d$trials - d$germ
  d
}

germ_formula <- cbind(germ, fail) ~ site * season + poly(moist, 2) +
  t2(doy, bs = "cc", k = 8) +
  t2(age, by = storage, k = 6, bs = "cr", sp = 2000)

germ_learner_ctrl <- function(...) {
  gb_learner_ctrl(step_size = 0.3, selection_loss_type = "deviance",
                  use_step_size_for_selection = TRUE,
                  knots = list(doy = c(0, 365)), ...)
}


test_that("cloglog boosting with a learner-corrected bootstrap is coherent", {
  skip_on_cran()
  set.seed(401)
  d <- sim_germination()

  # tau > 0 has no effect without EDF tracking, and is reset to 0
  set.seed(1)
  expect_warning(
    mod <- gbgam(germ_formula, family = binomial(link = "cloglog"), data = d,
                 learner_control = germ_learner_ctrl(tau = 1),
                 stop_control    = gb_stop_ctrl("none", max_iter = 150),
                 edf_control     = gb_edf_ctrl("none"),
                 uq_control      = gb_uq_ctrl(
                   "bootstrap", n_samples = 20, correction = "learner"),
                 verbose = FALSE),
    "tau")
  expect_identical(mod$tau, 0)
  expect_equal(mod$mstop, 150)
  expect_complex_invariants(mod, d, truth_eta = d$eta, min_cor = 0.85)
  expect_gt(1 - mod$deviance / mod$null.deviance, 0.4)

  # One learner per unpenalised column and per penalised range space; with
  # ASR off (tau reset to 0) the factor-by smooths keep their starting sp
  labels <- vapply(mod$learners, `[[`, "", "var_label")
  expect_identical(labels, c(
    "(Intercept)", "sitesouth", "seasonspring", "seasonsummer",
    "poly(moist, 2)1", "poly(moist, 2)2", "sitesouth:seasonspring",
    "sitesouth:seasonsummer", "t2(doy)r", "t2(age):storageambientr",
    "t2(age):storageambientn1", "t2(age):storagecoldr",
    "t2(age):storagecoldn1"))
  is_pen <- vapply(mod$learners, inherits, NA, what = "demmler_reinsch")
  expect_identical(labels[is_pen], c("t2(doy)r", "t2(age):storageambientr",
                                     "t2(age):storagecoldr"))
  for (k in which(startsWith(labels, "t2(age)") & is_pen)) {
    expect_equal(exp(mod$learners[[k]]$rho), 2000, tolerance = 1e-12)
  }

  # Zero-trial rows: weight 0, y = 0 and no deviance contribution
  zero <- d$trials == 0
  expect_equal(unname(mod$prior.weights), d$trials)
  expect_true(all(mod$y[zero] == 0))
  expect_true(all(residuals(mod, type = "deviance")[zero] == 0))
  expect_true(all(mod$fitted.values[zero] > 0 & mod$fitted.values[zero] < 1))

  # Recovery of the shape of the truth on a grid
  grid <- data.frame(
    site = "north", season = "spring",
    storage = rep(c("cold", "ambient"), each = 41),
    moist = seq(0.1, 0.6, length.out = 41), doy = 100,
    age = seq(0, 5, length.out = 41))
  pt <- predict(mod, grid, type = "terms")
  cold <- grid$storage == "cold"
  # The moisture optimum is at 0.475
  expect_gte(grid$moist[cold][which.max(pt[cold, "poly(moist, 2)"])], 0.4)
  expect_lte(grid$moist[cold][which.max(pt[cold, "poly(moist, 2)"])], 0.55)
  # Seed lots decline with age, faster in ambient storage
  slope <- function(v, a) unname(stats::coef(stats::lm(v ~ a))[2])
  slope_cold <- slope(pt[cold, "t2(age):storagecold"], grid$age[cold])
  slope_amb  <- slope(pt[!cold, "t2(age):storageambient"], grid$age[!cold])
  expect_lt(slope_cold, 0)
  expect_lt(slope_amb, slope_cold)
  # The cyclic day-of-year effect joins up at the ends of the year
  ends <- predict(mod, grid[c(1, 1), ] |> transform(doy = c(0, 365)),
                  type = "terms")
  expect_equal(ends[1, "t2(doy)"], ends[2, "t2(doy)"], tolerance = 1e-10)

  # Fixed-mode bootstrap: one widening factor per learner, exactly 1 for the
  # learners never selected (they never vary over the replicates)
  expect_identical(mod$uq_control$control$mode, "fixed")
  expect_identical(mod$n_samples_done, 20L)
  expect_identical(names(mod$widening_factor), labels)
  never <- setdiff(seq_along(labels), mod$metrics$k_opts)
  expect_gt(length(never), 0)
  expect_true(all(mod$widening_factor[never] == 1))
  expect_true(all(is.finite(mod$widening_factor)))
  expect_true(all(mod$widening_factor[-never] > 0 &
                  mod$widening_factor[-never] != 1))

  # vcov() is the widened bootstrap covariance; predict() uses it
  p  <- length(mod$coefficients)
  Vp <- vcov(mod)
  expect_identical(dim(Vp), c(p, p))
  expect_true(isSymmetric(unname(Vp)))
  expect_identical(Vp, mod$Vp)
  Xp <- predict(mod, grid, type = "lpmatrix")
  pl <- predict(mod, grid, se.fit = TRUE)
  pr <- predict(mod, grid, se.fit = TRUE, type = "response")
  expect_equal(as.vector(pl$se.fit), unname(sqrt(rowSums((Xp %*% Vp) * Xp))),
               tolerance = 1e-10)
  expect_equal(as.vector(pr$se.fit),
               abs(mod$family$mu.eta(as.vector(pl$fit))) *
                 as.vector(pl$se.fit),
               tolerance = 1e-10)
  # The se of a fitted probability built from thousands of trials is below
  # the sd of a single Bernoulli trial
  expect_true(all(pr$se.fit > 0 & pr$se.fit < sqrt(pr$fit * (1 - pr$fit))))

  # Without EDF tracking the path is deterministic: extending the fit equals
  # fitting the longer path at once, and the bootstrap is rerun
  set.seed(2)
  ext <- extend_gbgam(mod, 50)
  long <- gbgam(germ_formula, family = binomial(link = "cloglog"), data = d,
                learner_control = germ_learner_ctrl(),
                stop_control    = gb_stop_ctrl("none", max_iter = 200),
                edf_control     = gb_edf_ctrl("none"),
                verbose = FALSE)
  expect_equal(ext$mstop, 200)
  expect_equal(ext$coefficients, long$coefficients, tolerance = 1e-10)
  expect_identical(ext$metrics$k_opts, long$metrics$k_opts)
  expect_identical(ext$metrics$k_opts[1:150], mod$metrics$k_opts)
  expect_identical(ext$n_samples_done, 20L)
  expect_identical(nrow(ext$boot_coef), 20L)
  expect_identical(ncol(ext$boot_coef), length(ext$coefficients))
  expect_false(identical(ext$boot_coef, mod$boot_coef))
  expect_identical(names(ext$widening_factor), labels)
  expect_false(isTRUE(all.equal(ext$Vp, mod$Vp)))
})

test_that("probit boosting with exact EDF and AIC stopping is coherent", {
  skip_on_cran()
  set.seed(401)
  d <- sim_germination()

  # A non-canonical link: exact EDF, and no comparison with Hutch++
  expect_message(
    mod <- gbgam(germ_formula, family = binomial(link = "probit"), data = d,
                 learner_control = germ_learner_ctrl(),
                 stop_control    = gb_stop_ctrl("aic", iter_no_improv = 20),
                 edf_control     = gb_edf_ctrl("exact"),
                 verbose = FALSE),
    "Non-canonical link")
  expect_true(mod$converged)
  expect_gt(mod$mstop, 0)
  # The truth on the probit scale
  expect_complex_invariants(mod, d,
                            truth_eta = stats::qnorm(1 - exp(-exp(d$eta))),
                            min_cor = 0.85)

  # edf~ exceeds the intercept-only model and stays below the full basis,
  # and mstop is the best iteration of the AIC path
  p <- length(mod$coefficients)
  expect_gt(mod$edf.total, 1)
  expect_lt(mod$edf.total, p)
  expect_equal(mod$mstop, which.min(mod$metrics$ls) - 1)
  expect_equal(unname(mod$prior.weights), d$trials)

  # With the default tau ASR moves the factor-by smooths off their starting sp
  labels <- vapply(mod$learners, `[[`, "", "var_label")
  asr    <- which(labels %in% c("t2(age):storageambientr",
                                "t2(age):storagecoldr"))
  expect_length(asr, 2)
  for (k in asr) {
    expect_false(isTRUE(all.equal(exp(mod$learners[[k]]$rho), 2000)))
  }
})
