# Complex Gamma(log) model: hospital length of stay (days; positive and
# right-skewed, shape 4). Multiplicative effects of the admission type
# (elective, emergency, transfer), a smooth age effect, a severity effect
# that differs by admission type and a ward random effect (12 wards). Three
# severity scores are missing (na.omit), 8 audited records get weight 0 and
# ward w12 is removed by the subset (drop.unused.levels = TRUE).
# Covers a non-canonical link with intercept padding, a bs = "re" learner
# refitted after every iteration and excluded from selection, nu_k in the
# update, tau = 1.5, factor-by smooths with null-space learners, AIC
# stopping with exact EDF, and Bayesian backfitting ("kept" with
# undersmoothing to edf_floor = 1.5; "free" by bayesian_backfitting()). edf~ is computed
# exactly: for the log link the Jacobian has a Frobenius norm much larger
# than its trace, so NA-Hutch++ is not asserted here.

# Length-of-stay data; `eta` is the true linear predictor and `ward_effect`
# the true random effect of each ward
simulate_length_of_stay <- function(n = 220L, n_ward = 12L, seed = 301L) {
  set.seed(seed)
  d <- data.frame(
    age  = stats::runif(n, 20, 90),
    sev  = stats::runif(n),
    ward = factor(sample(sprintf("w%02d", seq_len(n_ward)), n, TRUE)),
    adm  = factor(sample(c("elective", "emergency", "transfer"), n, TRUE,
                         prob = c(0.5, 0.4, 0.1))),
    w    = 1
  )
  ward_effect <- stats::rnorm(n_ward, sd = 0.3)
  d$eta <- 1 + 0.6 * sin(d$age / 15) +
    0.8 * d$sev^2 * (d$adm == "emergency") +
    0.3 * (d$adm == "transfer") + ward_effect[d$ward]
  d$y <- stats::rgamma(n, shape = 4, rate = 4 / exp(d$eta))
  d$w[sample(n, 8L)]   <- 0
  d$sev[c(5, 17, 40)] <- NA
  list(data = d, ward_effect = ward_effect)
}

fml_los <- y ~ adm + t2(age, k = 7, bs = "cr") +
  t2(sev, by = adm, k = 5, bs = "cr") + t2(ward, bs = "re")

# Index of the ward random-effect learner, from the learners built on the
# rows the fit uses
los_re_index <- function(d) {
  keep <- d$ward != "w12" & !is.na(d$sev)
  lb   <- new_learners(fml_los, data = droplevels(d[keep, ]))
  lab  <- vapply(lb$learners, `[[`, "", "var_label")
  k_re <- which(startsWith(lab, "t2(ward)"))
  testthat::expect_length(k_re, 1L)
  k_re
}

# The full configuration; the non-canonical link is announced
fit_los <- function(d, k_re) {
  testthat::expect_message(
    mod <- gbgam(
      fml_los, family = stats::Gamma(link = "log"), data = d, weights = w,
      subset = ward != "w12", na.action = na.omit,
      drop.unused.levels = TRUE,
      learner_control = gb_learner_ctrl(
        step_size = 0.2, intercept_padding = TRUE,
        learners_to_fit_after_each_iter  = k_re,
        excluded_learner_indices         = k_re,
        use_learner_step_size_for_update = TRUE, tau = 1.5),
      stop_control      = gb_stop_ctrl("aic", iter_no_improv = 30L,
                                       max_iter = 600L),
      edf_control = gb_edf_ctrl("exact"),
      uq_control  = gb_uq_ctrl("bayes", n_samples = 200L,
                                            n_burnin = 50L, linear = "kept",
                                            undersmooth = TRUE,
                                            edf_floor = 1.5),
      verbose = FALSE),
    "Non-canonical link")
  mod
}

test_that("Gamma(log) length-of-stay fit: data handling, path and invariants", {
  skip_on_cran()
  sim  <- simulate_length_of_stay()
  d    <- sim$data
  k_re <- los_re_index(d)
  set.seed(1)  # Bayesian draws
  mod  <- fit_los(d, k_re)
  idx  <- as.integer(rownames(mod$model))
  ms   <- mod$mstop

  expect_complex_invariants(mod, d, truth_eta = d$eta[idx], min_cor = 0.8)

  # The fit used nu_k in the update: the fixed-mode replay above reproduces
  # it with the flag on, and moves away from it with the flag off
  m_off <- mod
  m_off$use_learner_step_size_for_update <- FALSE
  r_off <- replay_on_y(m_off, mod$y, "fixed")
  expect_gt(max(abs(r_off$linear.predictors - mod$linear.predictors)), 0.05)

  # Subset, na.omit and drop.unused.levels
  keep_sub <- which(d$ward != "w12")
  expect_identical(idx, setdiff(keep_sub, c(5L, 17L, 40L)))
  expect_s3_class(mod$na.action, "omit")
  expect_identical(sort(as.integer(names(mod$na.action))),
                   intersect(c(5L, 17L, 40L), keep_sub))
  expect_identical(nlevels(mod$model$ward), 11L)
  expect_false("w12" %in% levels(mod$model$ward))

  # Zero weights at the audited rows; their deviance residuals are 0
  expect_equal(unname(weights(mod)), d$w[idx])
  zero <- which(d$w[idx] == 0)
  expect_gt(length(zero), 0L)
  expect_equal(unname(residuals(mod, type = "deviance")[zero]),
               rep(0, length(zero)))
  expect_true(all(residuals(mod, type = "deviance")[-zero] != 0))

  # Stopping path: AIC stops at the minimum of the stopping loss, after
  # iter_no_improv iterations without improvement
  ls <- mod$metrics$ls
  expect_true(mod$converged)
  expect_gt(ms, 10L)
  expect_lt(ms, 600L - 30L)
  expect_equal(which.min(ls), ms + 1)
  expect_equal(sum(!is.na(mod$metrics$k_opts)), ms + 30)

  # Neither the padded intercept nor the excluded RE learner is selected
  k_opts <- mod$metrics$k_opts[seq_len(ms + 30L)]
  expect_false(any(k_opts %in% c(1L, k_re)))

  # The RE learner is refitted after every iteration: one post_rho column,
  # finite at every iteration run; relaxed by tau away from rho0, and the
  # final learner's rho is the one recorded at mstop
  pr    <- mod$metrics$post_rhos
  re_lr <- mod$learners[[k_re]]
  expect_true(is.matrix(pr))
  expect_identical(ncol(pr), 1L)
  expect_gte(nrow(pr), ms)
  expect_true(all(is.finite(pr[seq_len(ms), 1L])))
  expect_equal(re_lr$rho, pr[ms, 1L])
  expect_gt(abs(re_lr$rho - re_lr$rho0), 1)
  expect_gt(re_lr$edf, 1)

  # ASR with tau = 1.5: the ledf of every fitted penalised learner is at
  # least 1.5 edf~_k (all its columns once lambda is zero); a penalised
  # learner never fitted keeps its starting rho
  for (lr in mod$learners) {
    if (!inherits(lr, "demmler_reinsch")) next
    s    <- if (isTRUE(lr$lambda_at_zero)) rep(1, length(lr$D)) else
      1 / (1 + lr$D * exp(lr$rho))
    ledf <- sum(s)
    if (lr$edf > 0) {
      expect_gte(ledf, min(1.5 * lr$edf, length(lr$D)) - 1e-6,
                 label = lr$var_label)
    } else {
      expect_identical(lr$rho, lr$rho0, label = lr$var_label)
    }
  }

  # AIC and logLik: free-scale family, so one extra df for the scale
  expect_equal(AIC(mod), ls[ms + 1L] + 2)
  ll <- logLik(mod)
  expect_equal(attr(ll, "df"), mod$edf.total + 1)
  expect_equal(AIC(mod), -2 * as.numeric(ll) + 2 * attr(ll, "df"))

  # Recovery of the truth (dispersion 1 / shape = 0.25)
  expect_gt(mod$scale, 0.15)
  expect_lt(mod$scale, 0.4)
  nd <- training_rows(mod, d)
  pt <- predict(mod, newdata = nd, type = "terms")
  expect_gt(cor(pt[, "t2(age)"], sin(nd$age / 15)), 0.9)
  ward_hat <- tapply(pt[, "t2(ward)"], droplevels(nd$ward), mean)
  expect_gt(cor(ward_hat, sim$ward_effect[1:11]), 0.8)
  g <- expand.grid(sev = c(0.1, 0.9), adm = levels(d$adm))
  g$age  <- 50
  g$ward <- "w01"
  pg <- predict(mod, newdata = g, type = "terms")
  rise <- function(col, lev) {
    diff(pg[g$adm == lev, col])
  }
  expect_gt(rise("t2(sev):admemergency", "emergency"), 0.3)
  expect_lt(abs(rise("t2(sev):admelective", "elective")), 0.2)
})

test_that("Gamma(log) length-of-stay fit: Bayesian backfitting and post-fit", {
  skip_on_cran()
  d    <- simulate_length_of_stay()$data
  k_re <- los_re_index(d)
  set.seed(2)  # Bayesian draws
  mod  <- fit_los(d, k_re)
  n    <- nrow(mod$model)
  p    <- length(mod$coefficients)
  ms   <- mod$mstop
  lr_edf <- vapply(mod$learners, function(lr) lr$edf, numeric(1))

  # Posterior from linear = "kept" with undersmoothing
  expect_identical(mod$uq_control$method, "bayes")
  expect_identical(mod$n_samples_done, 200L)
  expect_identical(dim(mod$Vp), c(p, p))
  expect_true(all(is.finite(mod$coefficients_bayes)))
  expect_true(all(vapply(mod$beta_samples, ncol, integer(1)) == 200L))
  info <- mod$bb_info
  expect_identical(info$linear, "kept")
  expect_true(is.na(info$linear_accept))
  expect_gt(info$phi, 0.15)
  expect_lt(info$phi, 0.4)
  # edf_floor = 1.5 lifts at least two learners, with clear margins
  us <- info$undersmoothed
  expect_s3_class(us, "data.frame")
  expect_gte(nrow(us), 2L)
  expect_true(all(us$edf_from < 1.5))
  expect_equal(us$edf_to, rep(1.5, nrow(us)), tolerance = 1e-6)
  expect_true(all(us$rho_to < us$rho_from))
  # Acceptance rates: NA exactly for the learners held at zero (never
  # fitted), in (0, 1] for the others
  acc <- mod$mh_accept_rate
  expect_length(acc, length(mod$learners))
  expect_identical(which(is.na(acc)), which(lr_edf == 0))
  expect_true(all(acc[!is.na(acc)] > 0 & acc[!is.na(acc)] <= 1))
  # The posterior mean stays close to the boosted fit
  X <- create_concat_design(mod)
  eta_bayes <- as.vector(X %*% mod$coefficients_bayes[colnames(X)])
  expect_gt(cor(eta_bayes, mod$linear.predictors), 0.98)

  # smooth_intervals(): one column per t2() smooth, both centres
  sp <- smooth_intervals(mod, centre = "posterior")
  sb <- smooth_intervals(mod)
  labs <- c("t2(age)", "t2(sev):admelective", "t2(sev):admemergency",
            "t2(sev):admtransfer", "t2(ward)")
  expect_identical(colnames(sp$fit), labs)
  expect_identical(dim(sp$fit), c(n, 5L))
  for (el in c("fit", "se", "lower", "upper")) {
    expect_true(all(is.finite(sp[[el]])), label = el)
  }
  expect_true(all(sp$se > 0))
  expect_true(all(sp$lower < sp$fit & sp$fit < sp$upper))
  pt <- predict(mod, newdata = training_rows(mod, d), type = "terms")
  expect_equal(unname(sb$fit[, "t2(age)"]),
               unname(pt[, "t2(age)"] + attr(pt, "constant")))
  expect_gt(cor(sp$fit[, "t2(age)"], sb$fit[, "t2(age)"]), 0.95)

  # predict() with standard errors
  ps <- predict(mod, newdata = training_rows(mod, d), se.fit = TRUE)
  expect_equal(as.vector(ps$fit), as.vector(mod$linear.predictors))
  expect_true(all(is.finite(ps$se.fit) & ps$se.fit > 0))

  # summary and plot
  sm <- summary(mod)
  expect_identical(sm$family, "Gamma/log")
  lt <- sm$learner_table
  expect_equal(sum(lt$sel), ms)
  expect_identical(lt$sel[c(1L, k_re)], c(0L, 0L))
  expect_equal(lt$acc_edf, round(lr_edf, 2))
  out <- paste(cli::cli_fmt(print(sm)), collapse = " ")
  expect_match(out, "best AIC", fixed = TRUE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(mod, pages = 1, seWithMean = TRUE))

  # linear = "free": the log link has no boundary, so the joint linear
  # block samples without leaving the valid region
  set.seed(303)
  bf <- bayesian_backfitting(mod, n_samples = 100L, n_burnin = 30L,
                             linear = "free")
  expect_identical(bf$coefficients, mod$coefficients)
  expect_identical(bf$n_samples_done, 100L)
  expect_identical(bf$bb_info$linear, "free")
  expect_gt(bf$bb_info$linear_accept, 0)
  expect_lte(bf$bb_info$linear_accept, 1)
  expect_true(all(is.finite(bf$Vp)))
  # Only the penalised learners never fitted stay at zero
  held <- vapply(mod$learners, inherits, logical(1), "demmler_reinsch") &
    lr_edf == 0
  expect_identical(which(is.na(bf$mh_accept_rate)), which(held))
})
