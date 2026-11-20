# Finite-difference checks of the tracked derivatives in hard settings: ASR
# with replayed rho paths, zero prior weights, offsets, non-canonical links,
# intercept padding, refits and the Hutch++ sketch. Each check differences
# fits on the reference fit's own path (selection, rhos and post-iteration
# rhos replayed), after checking that the replay reproduces the reference.


# ---- Helpers -----------------------------------------------------------------

# Two covariates, an offset and non-integer weights
adv_data <- function(n, seed, off_sd = 0.3) {
  set.seed(seed)
  data.frame(x = sort(stats::runif(n)), z = stats::runif(n),
             o = off_sd * stats::rnorm(n), w = stats::runif(n, 0.5, 2))
}

adv_fml <- y ~ t2(x, k = 6, bs = "cr") + t2(z, k = 5, bs = "cr") + offset(o)

# Fixed-length fit that keeps the derivatives and the fitted-value path
adv_fit <- function(dat, family, weights, lc, n_iter, edf = "exact",
                    z = NULL) {
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  ec <- if (edf == "hutchpp") gb_edf_ctrl("hutchpp", z = z) else
    gb_edf_ctrl(edf)
  suppressWarnings(quiet_fit(
    adv_fml, data = dat, family = family, weights = weights,
    learner_control = do.call(gb_learner_ctrl, lc),
    stop_control = gb_stop_ctrl("none", max_iter = n_iter),
    edf_control = ec, uq_control = gb_uq_ctrl("none"),
    save_fitted_values = TRUE))
}

# Refit of `ref` on the response y for its first m iterations, with the
# recorded selection, rhos and post-iteration rhos replayed (tau = 0)
adv_replay <- function(ref, y, m = ref$mstop, edf = "none") {
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  rec <- recorded_paths(ref, ref$mstop)
  b   <- reset_gbgam(ref, y, edf)
  b$preselected_learner_indices <- rec$path
  b["preselected_rhos"]      <- list(rec$rhos)
  b["preselected_post_rhos"] <- list(rec$post_rhos)
  b$stopping_method <- "none"
  b$max_iter <- m
  b$tau      <- 0
  suppressWarnings(suppressMessages(
    fit_gbgam(b, Sys.time(), save_fitted_values = TRUE, verbose = FALSE)))
}

# Central differences, over y[j] for j in js, of the replay of `ref` to
# iteration m: fitted means at every iteration (n x |js| x (m + 1)), eta
# and the intercept coefficient at m
adv_fd <- function(ref, js = seq_along(ref$y), m = ref$mstop, eps = 1e-5) {
  y  <- ref$y
  n  <- length(y)
  mu <- array(NA_real_, c(n, length(js), m + 1L))
  eta <- matrix(NA_real_, n, length(js))
  int <- numeric(length(js))
  for (a in seq_along(js)) {
    yp <- y
    ym <- y
    yp[js[a]] <- y[js[a]] + eps
    ym[js[a]] <- y[js[a]] - eps
    fp <- adv_replay(ref, yp, m)
    fm <- adv_replay(ref, ym, m)
    mu[, a, ] <- (fp$fitted_values - fm$fitted_values) / (2 * eps)
    eta[, a]  <- (fp$linear.predictors - fm$linear.predictors) / (2 * eps)
    int[a]    <- (fp$learners[[1L]]$beta - fm$learners[[1L]]$beta) /
      (2 * eps)
  }
  list(mu = mu, eta = eta, int = int)
}

# ASR moved the smoothing parameters (which the replays then fix)
expect_rhos_moved <- function(ref, label) {
  expect_gt(length(unique(stats::na.omit(ref$metrics$rhos))), 5L,
            label = paste(label, "number of distinct rhos"))
}

# The replay of `ref` on its own response reproduces it exactly
expect_replay_exact <- function(ref, label) {
  r0 <- adv_replay(ref, ref$y)
  expect_identical(r0$metrics$k_opts, ref$metrics$k_opts,
                   label = paste(label, "replay: path"))
  expect_identical(r0$linear.predictors, ref$linear.predictors,
                   label = paste(label, "replay: eta"))
  expect_identical(r0$fitted_values, ref$fitted_values,
                   label = paste(label, "replay: fitted-value path"))
  invisible(r0)
}

# Traces of the FD Jacobians at every iteration (the FD edf~ path)
fd_traces <- function(fd) {
  apply(fd$mu, 3L, function(M) sum(diag(M)))
}


# ---- ASR with replayed rho paths, zero weights and learner step sizes --------

test_that("FD: ASR (tau = 1), Gamma(log), offset, zero weights, nu_k update", {
  skip_on_cran()
  n <- 30L
  d <- adv_data(n, 401L)
  set.seed(402)
  f <- 0.6 * sin(2 * pi * d$x) + 0.3 * d$z
  d$y <- stats::rgamma(n, 4, 4 / exp(f + d$o))
  w   <- d$w * stats::rbinom(n, 1, 0.7)  # about 30% zero prior weights
  expect_gt(sum(w == 0), 0L)
  cases <- list(
    list(lab = "init", init = TRUE, post = NULL),
    list(lab = "no init", init = FALSE, post = NULL),
    # t2(z)r (learner 4) penalised and refitted after each iteration, so its
    # rho path is relaxed by ASR and replayed
    list(lab = "post learner", init = TRUE, post = 4L))
  for (cs in cases) {
    lc  <- list(step_size = 0.3, tau = 1, initialise_intercept = cs$init,
                use_learner_step_size_for_update = TRUE,
                learners_to_fit_after_each_iter = cs$post,
                excluded_learner_indices = cs$post)
    ref <- adv_fit(d, stats::Gamma(link = "log"), w, lc, n_iter = 40L)
    lab <- paste0("ASR ", cs$lab, ":")
    expect_rhos_moved(ref, lab)
    # The learner step sizes below 1 enter the update
    is_dr <- vapply(ref$learners, inherits, logical(1), "demmler_reinsch")
    expect_true(any(vapply(ref$learners[is_dr], `[[`, numeric(1),
                           "nu_k") < 1))
    # ... and some step was taken with a relaxed rho, so with the nu_k that
    # ASR recomputed: rho in effect at iteration m is the one recorded at
    # the learner's previous selection
    k <- ref$metrics$k_opts
    r <- ref$metrics$rhos
    relaxed <- vapply(seq_along(k)[-1L], function(m) {
      p <- which(k[seq_len(m - 1L)] == k[m])
      length(p) > 0L && !is.na(r[max(p)]) &&
        r[max(p)] != ref$learners[[k[m]]]$rho0
    }, logical(1))
    expect_true(any(relaxed), label = paste(lab, "step with a relaxed nu_k"))
    if (!is.null(cs$post)) {
      expect_gt(length(unique(ref$metrics$post_rhos[, 1L])), 5L)
    }
    expect_replay_exact(ref, lab)
    fd  <- adv_fd(ref)
    expect_lt(max_rel_err(ref$derivs$dmudy, fd$mu[, , 41L]), 1e-7,
              label = paste(lab, "dmu~/dy relative error"))
    # edf~ at every iteration is the trace of the FD Jacobian then
    tr_fd <- fd_traces(fd)
    expect_lt(max_rel_err(ref$edfs, tr_fd), 1e-7,
              label = paste(lab, "edf~ path relative error"))
    edf_k <- vapply(ref$learners, `[[`, numeric(1), "edf")
    expect_lt(abs(sum(edf_k) - ref$edf.total), 1e-10,
              label = paste(lab, "|sum edf~_k - edf~|"))
    if (is.null(cs$post)) {
      # Crediting rule: each iteration's change of edf~ goes to the selected
      # learner, the starting edf~ to the intercept
      k  <- ref$metrics$k_opts
      rule <- function(tr) {
        de <- diff(tr)
        r  <- vapply(seq_along(edf_k), function(j) sum(de[k == j]),
                     numeric(1))
        r[1L] <- r[1L] + tr[1L]
        r
      }
      expect_lt(max(abs(edf_k - rule(ref$edfs))), 1e-12,
                label = paste(lab, "edf~_k vs rule on tracked edf~"))
      expect_lt(max(abs(edf_k - rule(tr_fd))) / max(abs(tr_fd)), 1e-7,
                label = paste(lab, "edf~_k vs rule on FD traces"))
    }
  }
})


# ---- Link sweep with an offset and non-integer weights -----------------------

test_that("FD: dmu~/dy and deta/dy for non-canonical links", {
  skip_on_cran()
  n <- 30L
  d <- adv_data(n, 411L, off_sd = 0.05)  # small: keep mu valid
  set.seed(412)
  f   <- 0.3 * sin(2 * pi * d$x) + 0.2 * d$z
  mu0 <- exp(f)
  y_bin <- pmin(pmax(stats::rbinom(n, 5, 1 - exp(-exp(f - 0.5))) / 5,
                     0.1), 0.9)  # keeps y +- eps inside [0, 1]
  specs <- list(
    list("Gamma/identity", stats::Gamma(link = "identity"),
         stats::rgamma(n, 3, 3 / mu0), d$w),
    list("inverse.gaussian/log", stats::inverse.gaussian(link = "log"),
         stats::rgamma(n, 5, 5 / mu0), d$w),
    # 5 trials times a non-integer weight
    list("binomial/cloglog", stats::binomial(link = "cloglog"), y_bin,
         5 * d$w),
    list("gaussian/log", stats::gaussian(link = "log"),
         mu0 + 0.2 * stats::rnorm(n), d$w),
    # Non-integer counts keep y - eps > 0
    list("poisson/sqrt", stats::poisson(link = "sqrt"),
         stats::rpois(n, 3 * mu0) + 0.2, d$w))
  for (s in specs) {
    dd   <- d
    dd$y <- s[[3]]
    ref  <- adv_fit(dd, s[[2]], s[[4]], list(step_size = 0.1, tau = 0),
                    n_iter = 40L)
    expect_false(is_canonical_link(s[[2]]), label = s[[1]])
    expect_replay_exact(ref, s[[1]])
    fd <- adv_fd(ref)
    expect_lt(max_rel_err(ref$derivs$dmudy, fd$mu[, , 41L]), 1e-7,
              label = paste(s[[1]], "dmu~/dy relative error"))
    expect_lt(max_rel_err(ref$derivs$detady, fd$eta), 1e-7,
              label = paste(s[[1]], "deta/dy relative error"))
  }
})


# ---- Intercept padding with a non-canonical link and ASR ---------------------

test_that("FD: padding, Gamma(log), offset, tau = 1: dmu~/dy and dedf_int", {
  skip_on_cran()
  n <- 30L
  d <- adv_data(n, 421L)
  set.seed(422)
  f   <- 0.5 * sin(2 * pi * d$x) + 0.3 * d$z
  d$y <- stats::rgamma(n, 4, 4 / exp(f + d$o))
  n_it <- 40L
  for (init in c(TRUE, FALSE)) {
    lab <- paste0("padding init = ", init, ":")
    lc  <- list(step_size = 0.2, tau = 1, intercept_padding = TRUE,
                initialise_intercept = init)
    ref <- adv_fit(d, stats::Gamma(link = "log"), d$w, lc, n_iter = n_it)
    expect_rhos_moved(ref, lab)
    expect_replay_exact(ref, lab)
    fd  <- adv_fd(ref)
    expect_lt(max(abs(ref$derivs$dmudy - fd$mu[, , n_it + 1L])), 1e-8,
              label = paste(lab, "dmu~/dy absolute error"))
    # The intercept is credited with dedf_int = tr(diag(h'(eta)) 1
    # d Delta_int/dy) of the last iteration: the FD of its increment
    prev   <- adv_replay(ref, ref$y, n_it - 1L, edf = "exact")
    fd_prev <- adv_fd(ref, m = n_it - 1L)
    share_fd <- sum(exp(ref$linear.predictors) * (fd$int - fd_prev$int))
    share_tr <- ref$learners[[1L]]$edf - prev$learners[[1L]]$edf
    expect_gt(abs(share_fd), 1e-5)
    expect_lt(abs(share_tr - share_fd), 1e-8,
              label = paste(lab, "dedf_int absolute error"))
    # ... and the selected learner with the rest of the iteration's change
    tr_fd  <- fd_traces(fd)
    k_last <- ref$metrics$k_opts[n_it]
    expect_lt(abs(ref$learners[[k_last]]$edf - prev$learners[[k_last]]$edf -
                    (tr_fd[n_it + 1L] - tr_fd[n_it] - share_fd)), 1e-8,
              label = paste(lab, "selected learner's share, absolute error"))
  }
})


# ---- Refits: CV replicate, bootstrap refit, extend_gbgam ---------------------

test_that("FD: derivatives of a CV replicate fit (zero training weights)", {
  skip_on_cran()
  n <- 30L
  d <- adv_data(n, 431L)
  set.seed(432)
  f   <- 0.6 * sin(2 * pi * d$x) + 0.3 * d$z
  d$y <- stats::rgamma(n, 4, 4 / exp(f + d$o))
  n_it <- 40L
  lc   <- list(step_size = 0.3, tau = 1,
               use_learner_step_size_for_update = TRUE)
  main <- adv_fit(d, stats::Gamma(link = "log"), d$w, lc, n_iter = n_it)
  base <- resample_fit_template(main)
  w_tr <- stats::rbinom(n, 1, 0.75) * 1
  expect_gt(sum(w_tr == 0), 0L)
  v    <- 1 - w_tr
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  # The replicate as resampling fits run it: its own selection and lambda
  free <- resample_fit_replay("free_selection", main)
  expect_false(free$selection)
  st <- suppressMessages(new_resample_fit(base, w_tr, v, free))
  st <- suppressMessages(fit_state_run(st, n_iter = n_it,
                                       use_stopping = FALSE))
  expect_gt(length(unique(stats::na.omit(st$rhos))), 5L,
            label = "CV replicate: number of distinct rhos")
  path <- list(k_opts = st$k_opts[seq_len(n_it)],
               rhos = ifelse(is.na(st$rhos), 0, st$rhos)[seq_len(n_it)],
               post_rhos = st$post_rhos)
  rep_set <- list(selection = TRUE, smoothing = TRUE,
                  edf_compute_method = "none")
  rep_mu <- function(y) {
    b   <- base
    b$y <- y
    s <- suppressMessages(new_resample_fit(b, w_tr, v, rep_set))
    s <- suppressMessages(advance_resample_fit(s, n_it, path, rep_set))
    main$family$linkinv(s$mod$linear.predictors)
  }
  mu_st <- main$family$linkinv(st$mod$linear.predictors)
  expect_identical(rep_mu(base$y), mu_st, label = "CV replay: mu")
  eps  <- 1e-5
  J_fd <- vapply(seq_len(n), function(j) {
    yp <- base$y
    ym <- base$y
    yp[j] <- yp[j] + eps
    ym[j] <- ym[j] - eps
    (rep_mu(yp) - rep_mu(ym)) / (2 * eps)
  }, numeric(n))
  expect_lt(max_rel_err(st$mod$derivs$dmudy, J_fd), 1e-7,
            label = "CV replicate dmu~/dy relative error")
  expect_lt(abs(st$mod$edf.total - sum(diag(J_fd))) / sum(diag(J_fd)), 1e-7,
            label = "CV replicate edf~ relative error")
})

test_that("bootstrap refit by reset_gbgam() equals a fresh fit", {
  skip_on_cran()
  n <- 30L
  d <- adv_data(n, 441L)
  set.seed(442)
  f   <- 0.6 * sin(2 * pi * d$x) + 0.3 * d$z
  d$y <- stats::rgamma(n, 4, 4 / exp(f + d$o))
  n_it <- 40L
  lc   <- list(step_size = 0.3, tau = 1,
               use_learner_step_size_for_update = TRUE)
  main <- adv_fit(d, stats::Gamma(link = "log"), d$w, lc, n_iter = n_it)
  y_b  <- stats::rgamma(n, 4, 4 / main$fitted.values)
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  # A freely stopping replicate refits selection, lambda and mstop
  s <- replicate_settings("free_stopping", main$metrics$k_opts,
                          main$metrics$rhos, n_it, main)
  bfit <- suppressMessages(refit_replicate(main, y_b, s))
  db   <- d
  db$y <- y_b
  fresh <- adv_fit(db, stats::Gamma(link = "log"), d$w, lc, n_iter = n_it)
  expect_rhos_moved(fresh, "bootstrap refit:")
  expect_identical(bfit$metrics$k_opts, fresh$metrics$k_opts)
  expect_lt(max(abs(bfit$linear.predictors - fresh$linear.predictors)),
            1e-12, label = "bootstrap refit eta")
  expect_lt(max(abs(bfit$metrics$rhos - fresh$metrics$rhos), na.rm = TRUE),
            1e-12, label = "bootstrap refit rhos")
  expect_lt(max(abs(bfit$edfs - fresh$edfs)), 1e-12,
            label = "bootstrap refit edf~ path")
  expect_lt(max(abs(bfit$derivs$dmudy - fresh$derivs$dmudy)), 1e-12,
            label = "bootstrap refit dmu~/dy")
})

test_that("extend_gbgam() with Gamma(inverse), offset, tau = 1 = longer fit", {
  skip_on_cran()
  n <- 30L
  d <- adv_data(n, 451L, off_sd = 0.05)
  set.seed(452)
  f   <- 0.3 * sin(2 * pi * d$x) + 0.2 * d$z
  d$y <- stats::rgamma(n, 4, 4 / exp(f))
  lc  <- list(step_size = 0.3, tau = 1,
              use_learner_step_size_for_update = TRUE)
  op <- options(gbgam.keep_derivs = TRUE)
  on.exit(options(op), add = TRUE)
  for (edf in c("exact", "hutchpp")) {
    lab <- paste0("extend (", edf, "):")
    set.seed(5)
    short <- adv_fit(d, stats::Gamma(link = "inverse"), d$w, lc,
                     n_iter = 20L, edf = edf, z = 15L)
    ext   <- suppressWarnings(suppressMessages(
      extend_gbgam(short, n_iter = 30L)))
    set.seed(5)
    long  <- adv_fit(d, stats::Gamma(link = "inverse"), d$w, lc,
                     n_iter = 50L, edf = edf, z = 15L)
    expect_identical(ext$metrics$k_opts, long$metrics$k_opts)
    expect_rhos_moved(long, lab)
    expect_lt(max(abs(ext$metrics$rhos - long$metrics$rhos), na.rm = TRUE),
              1e-12, label = paste(lab, "rhos"))
    expect_lt(max(abs(ext$linear.predictors - long$linear.predictors)),
              1e-12, label = paste(lab, "eta"))
    expect_lt(max(abs(ext$edfs - long$edfs)), 1e-12,
              label = paste(lab, "edf~ path"))
    expect_lt(max(abs(ext$derivs$dmudy - long$derivs$dmudy)), 1e-12,
              label = paste(lab, "dmu~/dy (or its sketch)"))
  }
})


# ---- Hutch++ sketch ----------------------------------------------------------

test_that("Hutch++ sketch = exact Jacobian times probes (padding, ASR path)", {
  skip_on_cran()
  n <- 40L
  d <- adv_data(n, 461L)
  set.seed(462)
  f   <- 0.6 * sin(2 * pi * d$x) + 0.3 * d$z
  w   <- d$w * stats::rbinom(n, 1, 0.7)
  expect_gt(sum(w == 0), 0L)
  specs <- list(
    list("Gamma/log", stats::Gamma(link = "log"),
         stats::rgamma(n, 4, 4 / exp(f + d$o)), w),
    list("binomial/probit", stats::binomial(link = "probit"),
         pmin(pmax(stats::rbinom(n, 4, stats::pnorm(0.5 * f + d$o)) / 4,
                   0.1), 0.9), 4 * w))
  for (s in specs) for (pad in c(FALSE, TRUE)) {
    lab <- paste0(s[[1]], ", padding = ", pad, ":")
    dd   <- d
    dd$y <- s[[3]]
    lc   <- list(step_size = 0.3, tau = 1, intercept_padding = pad,
                 use_learner_step_size_for_update = TRUE)
    asr  <- adv_fit(dd, s[[2]], s[[4]], lc, n_iter = 50L)
    expect_rhos_moved(asr, lab)
    # Exact and sketched tracking on the same (ASR) rho path
    lc_p <- utils::modifyList(lc, list(
      tau = 0, preselected_learner_indices = asr$metrics$k_opts,
      preselected_rhos = recorded_paths(asr, 50L)$rhos))
    ex <- adv_fit(dd, s[[2]], s[[4]], lc_p, n_iter = 50L)
    set.seed(9)
    hp <- adv_fit(dd, s[[2]], s[[4]], lc_p, n_iter = 50L, edf = "hutchpp",
                  z = 30L)
    expect_identical(ex$linear.predictors, asr$linear.predictors,
                     label = paste(lab, "replayed eta"))
    expect_lt(max_rel_err(ex$derivs$dmudy %*% hp$derivs$GGG,
                          hp$derivs$dmudy), 1e-12,
              label = paste(lab, "sketch dmu~/dy relative error"))
    expect_lt(max_rel_err(ex$derivs$detady %*% hp$derivs$GGG,
                          hp$derivs$detady), 1e-12,
              label = paste(lab, "sketch deta/dy relative error"))
  }
})
