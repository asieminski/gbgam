# A complex Gaussian (identity link) model: a field trial of crop yield
# (t/ha) on 240 plots. Yield depends on the variety (3-level factor), the
# fertiliser dose (ordered low < mid < high), irrigation (logical), the
# temperature anomaly (deg C from the 18 C site mean; a slope for variety B
# only), a seasonal effect of the sowing
# day, a gentle spatial trend in latitude and longitude, and rainfall that
# scales a second seasonal effect. The plot area sets the noise (sd
# 0.4 / sqrt(area)) and is the prior weight; a known baseline yield is the
# offset argument. Three temperatures are missing (na.exclude), and the
# plots sown on the first two days are left out (subset).
#
# Covers: factor treatment contrasts, ordered-factor polynomial contrasts, a
# logical, numeric:factor interactions, a single-column t2 penalty group, t2
# null-space columns and a numeric-by smooth as unpenalised learners; a
# starting sedf of 3 (sp < 0, nu_k < 1); step_size used for selection; rss
# and deviance selection; BIC stopping with relative patience; exact EDF; a
# free_stopping bootstrap with the smooth widening; smooth_intervals() with
# the learner widening; predict(se.fit); extend_gbgam(); update(); and the
# S3 methods.


# ---- Data and fits -----------------------------------------------------------

# Crop-yield field trial; `eta` is the true mean yield and `f_day` the true
# seasonal effect of the sowing day
sim_crop_trial <- function(n = 240L, seed = 201L) {
  set.seed(seed)
  d <- data.frame(
    variety = factor(sample(c("A", "B", "C"), n, TRUE)),
    dose    = ordered(sample(c("low", "mid", "high"), n, TRUE),
                      levels = c("low", "mid", "high")),
    irrig   = stats::runif(n) < 0.4,
    temp    = stats::rnorm(n, 0, 3),
    day     = stats::runif(n, 0, 120),
    lat     = stats::runif(n),
    lon     = stats::runif(n),
    rain    = stats::runif(n, 0.5, 1.5),
    area    = stats::runif(n, 0.5, 2),
    base    = stats::rnorm(n, 5, 0.2)
  )
  d$f_day <- sin(2 * pi * d$day / 120)
  d$eta   <- d$base +
    c(A = 0, B = 0.4, C = -0.2)[as.character(d$variety)] +
    c(low = 0, mid = 0.3, high = 0.45)[as.character(d$dose)] +
    0.3 * d$irrig +
    0.05 * d$temp * (d$variety == "B") +
    d$f_day +
    0.3 * d$rain * cos(2 * pi * d$day / 120) +
    1.6 * (d$lat - 0.5) * (d$lon - 0.5)
  d$y <- d$eta + stats::rnorm(n, sd = 0.4 / sqrt(d$area))
  d$temp[c(4L, 50L, 99L)] <- NA
  d
}

crop_formula <- y ~ variety + dose + irrig + temp:variety +
  t2(day, k = 8, bs = "cr", sp = -3) + t2(lat, lon, k = c(3, 3), bs = "cr") +
  t2(day, by = rain, k = 6, bs = "cr")

# Rows kept by subset = day > 2 and na.exclude
crop_n_fit <- function(d) sum(d$day > 2 & !is.na(d$temp))

crop_stop_ctrl <- function(d) {
  gb_stop_ctrl("aic", AIC_k = log(crop_n_fit(d)), iter_no_improv = 20,
               iter_no_improv_rel = 0.25)
}

# The BIC-stopped fit with a free_stopping bootstrap
fit_crop <- function(d, selection_loss_type = "rss",
                     uq_control = gb_uq_ctrl("none")) {
  gbgam(crop_formula, data = d, weights = area, offset = base,
        subset = day > 2, na.action = na.exclude,
        learner_control = gb_learner_ctrl(
          step_size                   = 0.2,
          use_step_size_for_selection = TRUE,
          selection_loss_type         = selection_loss_type
        ),
        stop_control = crop_stop_ctrl(d),
        edf_control  = gb_edf_ctrl("exact"),
        uq_control   = uq_control,
        verbose      = FALSE)
}

# Made once on first use (not at load time, so that skip_on_cran() skips it)
crop_fits <- local({
  fits <- NULL
  function() {
    if (is.null(fits)) {
      d <- sim_crop_trial()
      set.seed(202)
      # 10 replicates: smooth factors below 0.9 are Monte Carlo noise
      mod <- without_low_widening_warning(fit_crop(d, uq_control =
        gb_param_boot_ctrl(n_samples = 10, mode = "free_stopping",
                           correction = "smooth")))
      fits <<- list(d = d, mod = mod)
    }
    fits
  }
})


# ---- Learners and fit --------------------------------------------------------

test_that("the crop model has the expected learners and starting smoothness", {
  skip_on_cran()
  set.seed(1)
  mod <- crop_fits()$mod
  labels <- vapply(mod$learners, `[[`, "", "var_label")
  expect_identical(labels, c(
    "(Intercept)", "varietyB", "varietyC", "dose.L", "dose.Q", "irrigTRUE",
    "varietyA:temp", "varietyB:temp", "varietyC:temp",
    "t2(day)r", "t2(day)n1",
    "t2(lat,lon)rr", "t2(lat,lon)nr", "t2(lat,lon)rn",
    "t2(lat,lon)n1", "t2(lat,lon)n2", "t2(lat,lon)n3",
    "t2(day):rainr", "t2(day):rainn1", "t2(day):rainn2"
  ))
  is_dr <- vapply(mod$learners, inherits, NA, "demmler_reinsch")
  expect_identical(labels[is_dr], c("t2(day)r", "t2(lat,lon)nr",
                                    "t2(lat,lon)rn", "t2(day):rainr"))

  # Starting smoothness: sp = -3 gives sedf 3 without the step size and
  # nu_k < 1; the others start at sedf 1 with nu_k = 1
  sedf_at <- function(lr, nu) {
    s <- 1 / (1 + exp(lr$rho0) * lr$D)
    2 * nu * sum(s) - nu^2 * sum(s^2)
  }
  for (lr in mod$learners[is_dr]) {
    expect_equal(sedf_at(lr, lr$nu_k0), 1, tolerance = 1e-5,
                 label = lr$var_label)
    if (lr$var_label == "t2(day)r") {
      expect_equal(sedf_at(lr, 1), 3, tolerance = 1e-5)
      expect_lt(lr$nu_k0, 1)
    } else {
      expect_equal(lr$nu_k0, 1, tolerance = 1e-4, label = lr$var_label)
    }
    # ASR only lowers the smoothing parameter
    expect_lte(lr$rho, lr$rho0)
  }
  # tau = NULL resolves to 1 with EDF tracking
  expect_identical(mod$tau, 1)
})

test_that("the crop model satisfies the invariants and recovers the truth", {
  skip_on_cran()
  set.seed(2)
  f   <- crop_fits()
  mod <- f$mod
  d   <- f$d
  idx <- as.integer(rownames(mod$model))
  expect_identical(nrow(mod$model), crop_n_fit(d))
  expect_identical(as.integer(mod$na.action), c(4L, 50L, 95L))
  expect_equal(unname(mod$prior.weights), d$area[idx])
  expect_complex_invariants(mod, d, truth_eta = d$eta[idx], min_cor = 0.97)

  # With the temperature anomaly centred at 0, the varietyB contrast is a
  # learner of its own, not absorbed by the varietyB:temp slope
  lt <- summary(mod)$learner_table
  expect_gt(lt$sel[lt$learner == "varietyB"], 0)

  # Contrasts at temp = 0 (where variety B has no slope): B - A, C - A,
  # mid - low and high - low dose, and irrigation, against the truth (the
  # BIC stop shrinks them; mgcv::gam estimates the irrigation effect at
  # 0.19 on these data)
  ref <- data.frame(variety = factor("A", levels = c("A", "B", "C")),
                    dose = ordered("low", levels = c("low", "mid", "high")),
                    irrig = FALSE, temp = 0, day = 60, lat = 0.5,
                    lon = 0.5, rain = 1)
  nd <- ref[rep(1L, 6L), ]
  nd$variety[2:3] <- c("B", "C")
  nd$dose[4:5]    <- c("mid", "high")
  nd$irrig[6]     <- TRUE
  p <- as.vector(predict(mod, newdata = nd))
  truth <- c(0.4, -0.2, 0.3, 0.45, 0.3)
  expect_identical(sign(p[-1] - p[1]), sign(truth))
  expect_lt(max(abs(p[-1] - p[1] - truth)), 0.2)

  # The seasonal smooth follows sin(2 pi day / 120), and together with the
  # rain-by smooth (partly confounded with it, as rain has mean 1) the day
  # effect
  pt   <- predict(mod, type = "terms")
  excl <- as.integer(mod$na.action)
  expect_gt(stats::cor(pt[-excl, "t2(day)"], d$f_day[idx]), 0.95)
  expect_gt(stats::cor(
    pt[-excl, "t2(day)"] + pt[-excl, "t2(day):rain"],
    d$f_day[idx] + 0.3 * d$rain[idx] * cos(2 * pi * d$day[idx] / 120)
  ), 0.95)
})

test_that("BIC stopping uses the relative patience and the edf~ penalty", {
  skip_on_cran()
  set.seed(3)
  f   <- crop_fits()
  mod <- f$mod
  ls  <- mod$metrics$ls
  m   <- mod$mstop
  expect_equal(mod$AIC_k, log(crop_n_fit(f$d)))
  # mstop is the best iteration, and boosting ran the relative patience
  # beyond it (max(20, ceiling(0.25 m)) > 20 here)
  expect_equal(which.min(ls) - 1L, m)
  expect_gt(ceiling(0.25 * m), 20)
  expect_equal(length(ls) - 1L, m + max(20, ceiling(0.25 * m)))
  expect_length(mod$metrics$k_opts, length(ls) - 1L)
  # The loss at mstop is -2 logLik + log(n) edf~, with one more df for the
  # scale
  ll <- logLik(mod)
  expect_equal(attr(ll, "df"), mod$edf.total + 1)
  expect_equal(ls[m + 1L], -2 * as.numeric(ll) + mod$AIC_k * mod$edf.total,
               tolerance = 1e-8)
})

test_that("deviance selection equals rss selection for the Gaussian identity", {
  skip_on_cran()
  set.seed(4)
  f   <- crop_fits()
  dev <- fit_crop(f$d, selection_loss_type = "deviance")
  expect_identical(dev$selection_loss_type, "deviance")
  expect_identical(dev$metrics$k_opts, f$mod$metrics$k_opts)
  expect_equal(dev$mstop, f$mod$mstop)
  expect_equal(dev$coefficients, f$mod$coefficients, tolerance = 1e-10)
  expect_equal(dev$edf.total, f$mod$edf.total, tolerance = 1e-10)
})


# ---- Uncertainty quantification ----------------------------------------------

test_that("the free_stopping bootstrap and smooth_intervals are consistent", {
  skip_on_cran()
  set.seed(5)
  f   <- crop_fits()
  mod <- f$mod
  expect_identical(mod$n_samples_done, 10L)
  expect_length(mod$boot_mstop, 10L)
  expect_gt(length(unique(mod$boot_mstop)), 1L)
  expect_false(any(mod$boot_censored))
  expect_identical(dim(mod$boot_coef), c(10L, length(mod$coefficients)))

  sm_labels <- vapply(mod$smooth, `[[`, "", "label")
  expect_identical(sm_labels, c("t2(day)", "t2(lat,lon)", "t2(day):rain"))
  wf <- mod$widening_factor
  expect_identical(names(wf),
                   c(names(mod$coefficients)[seq_len(mod$nsdf)], sm_labels))
  expect_true(all(is.finite(wf) & wf > 0))

  si <- smooth_intervals(mod, alpha = 0.1, correction = "learner")
  n_full <- length(fitted(mod))
  expect_identical(dim(si$fit), c(n_full, 3L))
  expect_identical(colnames(si$fit), sm_labels)
  excl <- as.integer(mod$na.action)
  for (nm in c("fit", "se", "lower", "upper")) {
    expect_true(all(is.na(si[[nm]][excl, ])), label = nm)
  }
  expect_true(all(si$lower[-excl, ] < si$upper[-excl, ]))
  expect_equal(si$upper[-excl, ] - si$fit[-excl, ],
               stats::qnorm(0.95) * si$se[-excl, ])
  # The centre is the intercept plus the smooth
  pt <- predict(mod, type = "terms")
  expect_equal(unname(si$fit[-excl, ]),
               unname(pt[-excl, sm_labels] +
                        mod$coefficients[["(Intercept)"]]),
               tolerance = 1e-8)
  # One widening factor per learner, and the widening changes the se
  expect_true(all(is.finite(si$factor) & si$factor > 0))
  expect_identical(names(si$factor),
                   vapply(mod$learners, `[[`, "", "var_label"))
  si0 <- smooth_intervals(mod, alpha = 0.1, correction = "none")
  expect_null(si0$factor)
  expect_false(isTRUE(all.equal(si$se, si0$se)))

  # predict(se.fit = TRUE) at the training rows uses $Vp
  nd <- training_rows(mod, f$d)
  p  <- predict(mod, newdata = nd, se.fit = TRUE)
  expect_true(all(is.finite(p$se.fit) & p$se.fit > 0))
  X  <- predict(mod, newdata = nd, type = "lpmatrix")
  expect_equal(unname(as.vector(p$se.fit)),
               unname(sqrt(rowSums((X %*% mod$Vp) * X))), tolerance = 1e-8)
})


# ---- Extending and updating --------------------------------------------------

test_that("extend_gbgam continues the crop path and reruns the bootstrap", {
  skip_on_cran()
  set.seed(6)
  f   <- crop_fits()
  mod <- f$mod
  ext <- without_low_widening_warning(extend_gbgam(mod, 30))
  n0  <- length(mod$metrics$ls)
  expect_length(ext$metrics$ls, n0 + 30L)
  # The original path is continued, not refitted
  expect_equal(ext$metrics$ls[seq_len(n0)], mod$metrics$ls, tolerance = 1e-10)
  expect_identical(ext$metrics$k_opts[seq_len(n0 - 1L)], mod$metrics$k_opts)
  expect_equal(ext$mstop, which.min(ext$metrics$ls) - 1L)
  expect_identical(ext$n_samples_done, 10L)
  expect_length(ext$boot_mstop, 10L)
  expect_identical(ext$uq_control, mod$uq_control)
  # The main fit had no cap, so the extension gives its replicates none
  expect_null(mod$max_iter)
  expect_true(ext$uncapped)
  expect_false(any(ext$boot_censored))
  # The bootstrap is rerun (set.seed(6) above, not the original seed 202)
  expect_identical(dim(ext$boot_coef), c(10L, length(ext$coefficients)))
  expect_false(identical(ext$boot_coef, mod$boot_coef))
  idx <- as.integer(rownames(ext$model))
  expect_complex_invariants(ext, f$d, truth_eta = f$d$eta[idx],
                            min_cor = 0.97)
})

test_that("update() of the crop model keeps subset and na.exclude", {
  skip_on_cran()
  set.seed(7)
  d <- sim_crop_trial()
  # A direct gbgam() call: update() re-evaluates it in this frame
  m_dir <- gbgam(crop_formula, data = d, weights = area, offset = base,
                 subset = day > 2, na.action = na.exclude,
                 learner_control = gb_learner_ctrl(step_size = 0.2),
                 stop_control    = gb_stop_ctrl("none", max_iter = 40),
                 edf_control     = gb_edf_ctrl("exact"),
                 verbose         = FALSE)
  m_upd <- update(m_dir, stop_control = gb_stop_ctrl("none", max_iter = 20))
  expect_identical(nrow(m_upd$model), 231L)
  expect_equal(m_upd$mstop, 20)
  expect_identical(m_upd$na.action, m_dir$na.action)
  expect_s3_class(m_upd$na.action, "exclude")
  expect_length(fitted(m_upd), 234L)  # 231 kept + 3 NA-padded
  # The same path, cut at 20 iterations
  expect_identical(m_upd$metrics$k_opts, m_dir$metrics$k_opts[1:20])
  expect_complex_invariants(m_upd, d)
})


# ---- S3 methods --------------------------------------------------------------

test_that("print, summary and plot work on the crop model", {
  skip_on_cran()
  set.seed(8)
  mod <- crop_fits()$mod
  expect_output(print(mod))
  s <- summary(mod)
  expect_identical(nrow(s$learner_table), length(mod$learners))
  expect_identical(s$learner_table$learner,
                   vapply(mod$learners, `[[`, "", "var_label"))
  expect_equal(sum(s$learner_table$sel), mod$mstop)
  # print.summary.gbgam() writes through cli, as messages
  msgs <- paste(testthat::capture_messages(print(s)), collapse = "")
  for (lab in s$learner_table$learner) {
    expect_true(grepl(lab, msgs, fixed = TRUE), label = lab)
  }
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(mod, pages = 1, seWithMean = TRUE))
})
