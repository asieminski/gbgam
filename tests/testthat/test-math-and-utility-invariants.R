# Small helpers and identities: EDF formulas, weight splits, matrix
# utilities, the dispersion estimate and the progress display.

test_that("learner_edf identities are internally consistent", {
  s <- c(0.1, 0.25, 0.5, 0.8, 1.0)

  expect_equal(
    learner_edf(s, "tr2A_AA"),
    2 * learner_edf(s, "trA") - learner_edf(s, "trAA"),
    tolerance = 1e-12
  )

  expect_equal(learner_edf(s, "trA"), sum(s))
  expect_equal(learner_edf(s, "trAA"), sum(s^2))
})

test_that("init_rho solves the requested EDF equation", {
  D <- c(0.05, 0.2, 0.7, 2.0)

  for (target in c(0.8, 1.0, 1.8, 2.5)) {
    rho <- init_rho(D, target_edf = target, edf_type = "trA")
    s <- 1 / (1 + exp(rho) * D)
    expect_equal(sum(s), target, tolerance = 1e-5)
  }

  for (target in c(0.5, 1.0, 1.5)) {
    rho <- init_rho(D, target_edf = target, edf_type = "trAA")
    s <- 1 / (1 + exp(rho) * D)
    expect_equal(sum(s^2), target, tolerance = 1e-4)
  }
})

test_that("Newton update moves tr(A) toward the target", {
  D <- c(0.1, 0.4, 1.5, 3.0)
  rho0 <- 0.5
  trA0 <- sum(1 / (1 + exp(rho0) * D))

  target <- trA0 + 0.5
  step <- newton_trA_update(rho_old = rho0, D = D, target_trA = target)
  rho1 <- rho0 + step

  obj0 <- abs(sum(1 / (1 + exp(rho0) * D)) - target)
  obj1 <- abs(sum(1 / (1 + exp(rho1) * D)) - target)

  expect_lt(obj1, obj0)
})

test_that("solve_rho_trA reaches the target from far above it", {
  D   <- 1 / seq_len(10)
  rho <- solve_rho_trA(200, D, target_trA = 3, tol = 1e-8, label = "learner")
  expect_equal(sum(1 / (1 + D * exp(rho))), 3, tolerance = 1e-6)
})

test_that("solve_rho_trA bisects when Newton's method cannot go on", {
  trA <- function(rho, D) sum(1 / (1 + D * exp(rho)))
  D   <- 10^seq(-6, 6, length.out = 20)
  # From here Newton's step overshoots far below the solution, onto the flat
  # part of the trace, where exp(rho) underflows and the slope is zero
  rho0  <- -log(stats::median(D)) + 30
  below <- rho0 + newton_trA_update(rho0, D, 19.9)
  expect_false(is.finite(newton_trA_update(below, D, 19.9)))
  rho <- solve_rho_trA(rho0, D, 19.9, tol = 1e-8, label = "learner")
  expect_equal(trA(rho, D), 19.9, tolerance = 1e-8)
  # No part of any step helps: bisection alone finds the solution
  local_mocked_bindings(
    halve_step = function(...) list(fraction = 0, loss = Inf))
  rho <- solve_rho_trA(0, D, target_trA = 3, tol = 1e-8, label = "learner")
  expect_equal(trA(rho, D), 3, tolerance = 1e-6)
})

test_that("bracket_rho_trA holds the solution for any spectrum", {
  trA <- function(rho, D) sum(1 / (1 + D * exp(rho)))
  set.seed(1)
  for (i in 1:200) {
    D      <- rlnorm(sample(2:30, 1L), rnorm(1L, 0, 6), runif(1L, 0.1, 6))
    target <- length(D) * runif(1L)^sample(c(1, 8), 1L)
    b      <- bracket_rho_trA(D, target)
    expect_lt(b[1L], b[2L])
    # Up to rounding: a bound is tight when one eigenvalue dominates
    expect_gte(trA(b[1L], D), target * (1 - 1e-12))
    expect_lte(trA(b[2L], D), target * (1 + 1e-12))
  }
})

test_that("solve_rho_trA solves every target from every start", {
  trA    <- function(rho, D) sum(1 / (1 + D * exp(rho)))
  # Targets within rounding of p and of 0, and spectra with a very small
  # eigenvalue, on which the trace has a long flat step
  spectra <- list(1 / seq_len(10), 10^seq(-15, 15, length.out = 40),
                  c(1e4, 1e2, 1, 1e-2, 1e-30))
  for (D in spectra) {
    p <- length(D)
    for (target in c(p - c(1e-13, 1e-9, 1e-3, 0.5), p / 2, 1, 1e-3, 1e-9)) {
      for (rho0 in -log(stats::median(D)) + c(-30, -8, 0, 8, 30)) {
        rho <- solve_rho_trA(rho0, D, target, tol = 1e-8, label = "learner")
        expect_lt(abs(trA(rho, D) - target) / min(target, p - target, 1), 1e-5)
      }
    }
  }
})

test_that("solve_rho_trA stops at the rounding error of a flat trace", {
  # The trace is nearly flat at this target, so the Newton step at the
  # solution is rounding error divided by a tiny slope, larger than tol
  D      <- 10^seq(-6, 6, length.out = 20)
  target <- 20 - 1e-9
  rho    <- solve_rho_trA(-log(stats::median(D)), D, target, tol = 1e-6,
                          label = "learner")
  expect_lt(abs(sum(1 / (1 + D * exp(rho))) - target), 1e-13)
})

test_that("solve_rho_trA judges convergence on the whole step and has a cap", {
  D <- 1 / seq_len(10)
  # Steps cut to below tol, far from the target, do not look converged
  local_mocked_bindings(
    halve_step = function(loss_at, ...) list(fraction = 1e-12,
                                             loss = loss_at(1e-12)))
  expect_error(solve_rho_trA(0, D, target_trA = 7, tol = 1e-6,
                             label = "learner", max_iter = 5L),
               "did not converge in 5 iterations")
})

test_that("update_learner_trA_smoothness sets lambda to zero at target p", {
  set.seed(1)
  dat <- data.frame(x = runif(100), y = rnorm(100))
  lr  <- new_learners(y ~ t2(x, k = 8, bs = "cr"), data = dat)$learners[[2L]]
  lr$edf <- length(lr$beta)
  lr2 <- update_learner_trA_smoothness(lr)
  expect_true(lr2$lambda_at_zero)
  expect_identical(lr2$rho, rho_lambda_zero)
  # Just below p the search stays finite, from the construction rho
  lr$edf <- length(lr$beta) - 1e-10
  lr3 <- update_learner_trA_smoothness(lr)
  expect_false(isTRUE(lr3$lambda_at_zero))
  expect_lt(abs(sum(1 / (1 + lr3$D * exp(lr3$rho))) - lr$edf), 1e-12)
})

test_that("halve_step keeps the first acceptable fraction", {
  # Whole step accepted
  expect_identical(halve_step(function(t) 1 - t, 1),
                   list(fraction = 1, loss = 0))
  # An invalid point (Inf) and NaN are halved away
  expect_equal(halve_step(function(t) if (t > 0.3) Inf else 1 - t, 1)$fraction,
               0.25)
  expect_equal(halve_step(function(t) if (t > 0.3) NaN else 1 - t, 1)$fraction,
               0.25)
  # A rise is rejected, a rise within the allowance is not
  rise <- function(t) 1 + (t - 0.1)
  expect_equal(halve_step(rise, 1)$fraction, 0.0625)
  expect_equal(halve_step(rise, 1, allowance = 2)$fraction, 1)
  # Sufficient decrease: 1 - 0.01 t decreases, but by less than a quarter of
  # the slope predicts
  slow <- function(t) 1 - 0.01 * t
  expect_equal(halve_step(slow, 1)$fraction, 1)
  expect_identical(halve_step(slow, 1, slope = -1, max_halvings = 5L)$fraction,
                   0)
  # A start without a loss cannot be compared
  expect_equal(halve_step(function(t) 5, Inf), list(fraction = 1, loss = 5))
  expect_equal(halve_step(function(t) 5, NaN), list(fraction = 1, loss = 5))
  # Nothing acceptable
  expect_identical(halve_step(function(t) 2, 1, max_halvings = 3L),
                   list(fraction = 0, loss = 1))
})

test_that("loss_allowance is the larger of its relative and absolute limits", {
  expect_equal(loss_allowance(100), 100 * sqrt(.Machine$double.eps))
  expect_equal(loss_allowance(-100), 100 * sqrt(.Machine$double.eps))
  expect_equal(loss_allowance(0, roundoff = 1e-12), 1e-12)
  expect_equal(loss_allowance(100, roundoff = 1e-12),
               100 * sqrt(.Machine$double.eps))
})

test_that("find_nu_k solves the EDF target for a penalised learner", {
  fake_learner <- list(
    D = c(0.1, 0.4, 1.2, 2.0),
    rho = 0.3
  )

  nu_k <- find_nu_k(fake_learner, target_df = 1)
  s <- 1 / (1 + exp(fake_learner$rho) * fake_learner$D)
  edf_eff <- 2 * nu_k * sum(s) - nu_k^2 * sum(s^2)

  expect_true(is.finite(nu_k))
  expect_gt(nu_k, 0)
  expect_equal(edf_eff, 1, tolerance = 1e-8)
})

test_that(paste("split_integer_by_weights returns integers summing to",
                "requested total"), {
  z <- split_integer_by_weights(17, c(1, 2, 1))

  expect_type(z, "double")
  expect_equal(sum(z), 17)
  expect_true(all(z >= 0))
  expect_true(all(abs(z - round(z)) < 1e-12))
})

test_that("split_mat_cols preserves all columns exactly", {
  X <- matrix(seq_len(30), nrow = 5, ncol = 6)
  xs <- split_mat_cols(X, c(2, 1, 3))

  expect_length(xs, 3L)
  expect_equal(do.call(cbind, xs), X)
})

test_that("find_intercept_mle returns a root of the score equation", {
  cases <- list(
    list(
      family = new_family(stats::gaussian()),
      y = c(-1.2, 0.5, 1.8, 0.1, 2.0),
      w = c(1, 2, 1, 3, 2),
      offset = c(0, 0, 0, 0, 0),
      phi = 1.7
    ),
    list(
      family = new_family(stats::poisson()),
      y = c(0, 1, 3, 2, 1),
      w = c(1, 2, 1, 1, 3),
      offset = c(-0.4, 0.1, 0.2, -0.1, 0.5),
      phi = 1
    ),
    list(
      family = new_family(stats::binomial()),
      y = c(0.2, 0.5, 0.6, 0.8, 0.1),
      w = c(5, 3, 4, 2, 6),
      offset = c(0.1, -0.2, 0.3, 0, -0.1),
      phi = 1
    ),
    list(
      family = new_family(Gamma(link = "log")),
      y = c(0.5, 1.2, 2.0, 0.8, 1.5),
      w = c(2, 1, 3, 1, 2),
      offset = c(0, 0.1, -0.1, 0, 0.2),
      phi = 0.5
    ),
    list(
      family = new_family(inverse.gaussian(link = "log")),
      y = c(0.3, 0.8, 1.5, 0.6, 1.0),
      w = c(1, 2, 1, 2, 1),
      offset = c(0, 0, 0, 0, 0),
      phi = 0.4
    )
  )

  for (cc in cases) {
    a_hat <- find_intercept_mle(
      y = cc$y,
      family = cc$family,
      prior.weights = cc$w,
      offset = cc$offset
    )

    dd    <- family_eta_mu_derivs(cc$offset + a_hat, cc$family)
    score <- sum(cc$w * (cc$y - dd$mu) / (cc$phi * dd$Vmu) * dd$mu_eta)

    expect_true(is.finite(a_hat))
    expect_lt(abs(score), 1e-7)
  }
})

test_that("split_integer_by_weights sums to the total despite rounding error", {
  for (z in 3:120) {
    for (cw in list(c(1, 1, 1), c(1, 2, 2), c(2, 1, 1), c(1, 2, 1))) {
      zs <- split_integer_by_weights(z, cw)
      expect_identical(sum(zs), as.double(z), info = paste(z, toString(cw)))
      expect_lt(max(abs(zs - z * cw / sum(cw))), 1)
    }
  }
  expect_identical(split_integer_by_weights(80, c(1, 2, 1)), c(20, 40, 20))
  expect_identical(split_integer_by_weights(85, c(5, 6, 6)), c(25, 30, 30))
  expect_identical(split_integer_by_weights(170, c(5, 6, 6)), c(50, 60, 60))
})

test_that("split_mat_cols accepts a matrix only", {
  expect_error(split_mat_cols(list(diag(2), diag(3)), c(1, 1)),
               "Must be of type 'matrix'")
})

test_that("estimate_dispersion clamps when the weights sum to at most the edf", {
  mod <- list(family = new_family(stats::gaussian()),
              linear.predictors = rep(0, 10), y = stats::rnorm(10),
              prior.weights = rep(0.1, 10), edf.total = 2)
  expect_identical(estimate_dispersion(mod), 1e-12)
  mod$prior.weights <- rep(0.2, 10)          # sum(W) == edf: was Inf or NaN
  expect_identical(estimate_dispersion(mod), 1e-12)
  mod$prior.weights <- rep(1, 10)
  expect_gt(estimate_dispersion(mod), 0)
})

test_that("short_learner_name() keeps the learner suffix of long labels", {
  t2_lab <- "t2(soil_moisture_percentage,mean_annual_temperature_celsius)"
  sfx    <- c("rr", "nr", "rn", "n1", "n2", "n3")
  long   <- paste0(t2_lab, sfx)
  short  <- c("(Intercept)", "t2(x,z)rr", "poly(x, 2)1")
  expect_identical(short_learner_name(short), short)
  for (w in c(16L, 24L, 40L)) {
    s <- short_learner_name(long, w)
    expect_true(all(nchar(s) <= w), info = w)
    expect_true(all(endsWith(s, paste0(")", sfx))), info = w)
    expect_true(all(startsWith(s, "t2(soil")), info = w)
    expect_false(anyDuplicated(s) > 0L, info = w)
  }
  # by-smooths: the factor level after the last ")" is kept with the suffix
  by_lab <- paste0("t2(soil_moisture_percentage):fhabitat_", c("wet", "dry"),
                   "r")
  s <- short_learner_name(by_lab)
  expect_identical(s, paste0("t2(soi...):fhabitat_", c("wet", "dry"), "r"))
  # an end longer than max_w - 4 characters keeps its last characters
  lab <- paste0("t2(x):f", strrep("a", 30), c("_one_r", "_two_r"))
  s   <- short_learner_name(lab)
  expect_identical(nchar(s), c(24L, 24L))
  expect_true(all(startsWith(s, "...")))
  expect_true(all(endsWith(s, c("_one_r", "_two_r"))))
})

test_that("the progress table aligns its columns and tells learners apart", {
  d <- data.frame(soil_moisture_percentage = stats::runif(150),
                  mean_annual_temperature_celsius = stats::runif(150))
  d$y <- sin(2 * pi * d$soil_moisture_percentage) *
    d$mean_annual_temperature_celsius + stats::rnorm(150, sd = 0.3)
  mod <- gbgam(y ~ t2(soil_moisture_percentage,
                      mean_annual_temperature_celsius, k = c(4, 4)),
               data = d, verbose = FALSE,
               stop_control = gb_stop_ctrl("none", max_iter = 40),
               edf_control  = gb_edf_ctrl("exact"))
  st <- list(width = 80L, loss_name = "AIC", edf_tracked = TRUE,
             loss_tracked = FALSE, fam_str = "gaussian/identity",
             has_pbar = FALSE)
  m     <- mod$mstop
  lines <- trimws(progress_lines(st, m, mod, NULL, NULL, NULL,
                                 mod$metrics$k_opts[seq_len(m)], NULL,
                                 mod$metrics$edfs[seq_len(m + 1L)], 0.1),
                  "right")
  hdr  <- grep("^  learner ", lines)
  rows <- lines[grepl("^  t2\\(", lines)]
  expect_length(hdr, 1L)
  dr   <- vapply(mod$learners, inherits, NA, "demmler_reinsch")
  expect_length(rows, sum(dr))
  # each penalised learner shows its own suffix, and the names differ
  nm <- sub("^  (\\S+).*$", "\\1", rows)
  expect_false(anyDuplicated(nm) > 0L)
  sfx <- sub("^.*\\)", "", vapply(mod$learners[dr], `[[`, "", "var_label"))
  expect_true(all(endsWith(nm, sfx)))
  # right-aligned columns: sel, %sel and edf end where their headers end
  token_ends <- function(x) {
    g <- gregexpr("\\S+", x)[[1L]]
    as.integer(g + attr(g, "match.length") - 1L)
  }
  for (r in rows)
    expect_identical(token_ends(r)[2:4], token_ends(lines[hdr])[2:4], info = r)
  # with UTF-8 glyphs and a tracked loss, every column (also the Delta-loss
  # share and an NA edf~ shown as a line) ends where its header ends
  op <- options(cli.unicode = TRUE)
  on.exit(options(op), add = TRUE)
  st$loss_tracked <- TRUE
  k_opts <- mod$metrics$k_opts[seq_len(m)]
  mod$learners[[which(dr)[1L]]]$edf <- NA_real_
  lines <- trimws(progress_lines(st, m, mod, 1, 1, m, k_opts,
                                 c(0, -seq_len(m)),
                                 mod$metrics$edfs[seq_len(m + 1L)], 0.1),
                  "right")
  hdr  <- grep("^  learner ", lines)
  rows <- lines[grepl("^  t2\\(", lines)]
  expect_true(any(grepl("\u2500", rows)))
  expect_match(lines[hdr], "\u0394AIC%", fixed = TRUE)
  for (r in rows)
    expect_identical(token_ends(r)[2:5], token_ends(lines[hdr])[2:5], info = r)
})

test_that("the progress table keeps the sign of each loss share", {
  dat    <- simulate_gaussian_1d(n = 60, seed = 7)
  dat$x2 <- stats::runif(60)
  mod <- fit_fixed(y ~ t2(x, k = 5, bs = "cr") + t2(x2, k = 5, bs = "cr"),
                   dat, max_iter = 5)
  dr  <- which(vapply(mod$learners, inherits, NA, "demmler_reinsch"))
  expect_gte(length(dr), 2L)
  # learner dr[1] lowers the loss by 3 and dr[2] raises it by 1: shares
  # 150% and -50% of the net change of -2
  k_opts <- c(dr[1L], dr[1L], dr[1L], dr[2L])
  st <- list(width = 80L, loss_name = "AIC", edf_tracked = FALSE,
             loss_tracked = TRUE, fam_str = "gaussian/identity",
             has_pbar = FALSE)
  lines <- progress_lines(st, 4L, mod, 8, 8, 4L, k_opts, c(10, 9, 8, 7, 8),
                          NULL, 0.1)
  lab  <- vapply(mod$learners[dr[1:2]], `[[`, "", "var_label")
  rows <- vapply(short_learner_name(lab), function(l)
    lines[startsWith(lines, paste0("  ", l, " "))], "")
  expect_match(rows[1L], " 150.0%", fixed = TRUE)
  expect_match(rows[2L], " -50.0%", fixed = TRUE)
})

test_that("the progress display shows the stopping loss next to edf~", {
  mod <- local_gaussian_mod_edf()
  st  <- list(width = 80L, loss_name = "AIC", edf_tracked = TRUE,
              loss_tracked = TRUE, fam_str = "gaussian/identity",
              has_pbar = FALSE, mode = "log")
  edfs <- c(1, 2.5, 3.25)
  one  <- progress_oneliner(st, 2L, mod, 86.25, 85.5, 1L, edfs, 0.1)
  expect_match(one, "EDF: 3.2", fixed = TRUE)
  expect_match(one, "AIC: 86.2500 (best 85.5000 m=1)", fixed = TRUE)
  lines <- progress_lines(st, 2L, mod, 86.25, 85.5, 1L, c(1L, 1L),
                          c(87, 86.5, 86.25), edfs, 0.1)
  expect_true(any(startsWith(lines, "  EDF ")))
  expect_true(any(grepl("^  AIC .*86\\.2500 .*best .*85\\.5000 \\(m=1\\)",
                        lines)))
})

test_that("a narrow progress line drops whole parts, and a log keeps all", {
  mod <- local_gaussian_mod_edf()
  st  <- list(width = 40L, loss_name = "AIC", edf_tracked = TRUE,
              loss_tracked = TRUE, fam_str = "gaussian/identity",
              has_pbar = FALSE, mode = "log")
  edfs <- c(1, 2.5, 12.375)
  args <- list(st, 286L, mod, 86.25, 85.5, 280L, edfs, 1)
  log_line <- do.call(progress_oneliner, args)
  for (p in c("gaussian/identity", "m=286", "EDF: 12.4", "maxdEDF: 9.88",
              "AIC: 86.2500 (best 85.5000 m=280)", "1.0m"))
    expect_match(log_line, p, fixed = TRUE)
  args[[1L]]$mode <- "cr"
  for (w in c(40L, 45L, 60L, 80L)) {
    args[[1L]]$width <- w
    one <- do.call(progress_oneliner, args)
    expect_identical(nchar(one), w, info = w)
    expect_match(one, "m=286", fixed = TRUE)
    # every number shown is whole: no part is cut
    shown <- regmatches(one, gregexpr("[0-9.]+", one))[[1L]]
    expect_true(all(shown %in% c("286", "12.4", "9.88", "86.2500",
                                 "85.5000", "280", "1.0")), info = one)
  }
  # the multi-line header keeps m= by shortening the family string
  st$width <- 30L
  hdr <- progress_lines(st, 286L, mod, 86.25, 85.5, 280L, 1L,
                        c(87, 86.25), edfs, 1)[1L]
  expect_match(hdr, " m=286 ", fixed = TRUE)
  expect_identical(nchar(hdr), 30L)
})

test_that("the progress display writes newline-terminated ASCII lines to a log", {
  dat <- data.frame(x = stats::runif(80))
  dat$y <- sin(2 * pi * dat$x) + stats::rnorm(80, 0, 0.3)
  op <- options(cli.unicode = FALSE, cli.dynamic = FALSE)
  on.exit(options(op), add = TRUE)
  out <- capture.output(gbgam(y ~ t2(x, k = 6, bs = "cr"), data = dat,
                              stop_control = gb_stop_ctrl("none",
                                                          max_iter = 5L),
                              edf_control = gb_edf_ctrl("exact"),
                              verbose = TRUE))
  expect_gt(length(out), 0L)
  expect_false(any(grepl("\r", out, fixed = TRUE)))
  expect_false(any(grepl("[^ -~]", out)))
})
