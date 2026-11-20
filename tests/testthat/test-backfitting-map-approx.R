# bayesian_backfitting(): the MAP approximation (Stage 1) and the blocks of
# the posterior (smoothers, joint linear block, aliased columns).


# Free the null-space (unpenalised) learners: the simple form of
# sample_posterior()'s
# step for linear = "free", which also skips excluded and zero-norm learners
free_null_sp <- function(sp, m) {
  Map(function(spj, sm) {
    nl <- sm$learners$kind == "null"
    spj$status[nl] <- "zero"
    spj$rho[nl]    <- NA_real_
    spj
  }, sp, m$smoothers)
}


# ---- MAP approximation -------------------------------------------------------

test_that("Stage 1 brings the mode close to the boosted fit", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  m   <- model_as_smoothers(mod)
  s1  <- map_smoothing_parameters(m)
  expect_true(s1$status %in% c("converged", "stalled"))
  f_size <- sum(vapply(m$smoothers, function(sm)
    sum(sm$beta * as.vector(sm$VtWV %*% sm$beta)), numeric(1)))
  rel    <- s1$discrepancy / f_size
  cat(sprintf("Stage 1: discrepancy / size of the fit %.2g\n", rel))
  # squared W-norms: measured 5.6e-4; 2e-3 is 4.5% in norm
  expect_lt(rel, 2e-3)
  expect_error(map_smoothing_parameters(m, max_sweeps = 1L),
               class = "map_approx_failed")
})

test_that("perturbed_chol perturbs only a matrix without a factor", {
  M <- crossprod(matrix(c(2, 1, 0, 1, 3, 1, 1, 0, 2), 3))
  expect_identical(perturbed_chol(M, "a block"),
                   list(R = chol(M), ridge = 0))
  # Rank one with very different column sizes: the ridge follows each column
  M <- tcrossprod(c(1, 1e7))
  f <- perturbed_chol(M, "a block")
  expect_equal(f$ridge / diag(M), rep(f$ridge[1], 2) / diag(M)[1])
  expect_lt(f$ridge[1], 1e-6)
  A <- M
  diag(A) <- diag(A) + f$ridge
  expect_equal(crossprod(f$R), A)
  # Nothing to perturb towards
  for (bad in list(matrix(-1), diag(0, 2), matrix(NaN), matrix(Inf))) {
    expect_error(perturbed_chol(bad, "a block"), class = "b2b_not_pd")
  }
})

test_that("a sweep is halved to stay valid and to lower the objective", {
  set.seed(31)
  n    <- 120
  d    <- data.frame(x = stats::runif(n))
  d$y  <- stats::rgamma(n, 4, 4 * (1 + 0.3 * sin(2 * pi * d$x)))
  mod  <- quiet_fit(y ~ t2(x, k = 6, bs = "cr"), data = d,
                    family = stats::Gamma(),
                    stop_control = gb_stop_ctrl("none", max_iter = 40L),
                    edf_control = gb_edf_ctrl("exact"))
  m    <- model_as_smoothers(mod)
  sp   <- starting_lambdas(m)
  bl   <- model_blocks(m, sp, FALSE)
  mode <- list(a = m$a, beta = lapply(m$smoothers, `[[`, "beta"), eta = m$eta)
  mode <- penalised_fit(m, sp, mode, 1e-8)$state
  move <- function(s) {
    s$a   <- s$a + 0.02
    s$eta <- s$eta + 0.02
    s
  }
  obj  <- function(s) penalised_deviance(m, sp, s)
  old  <- move(mode)
  # An improving sweep is kept as it is
  expect_identical(keep_sweep(m, sp, bl, old, mode), mode)
  # One that leaves the range of the link (eta < 0) is cut back
  far     <- mode
  far$a   <- far$a - 10
  far$eta <- far$eta - 10
  expect_false(in_link_range(m$family, far$eta))
  kept <- keep_sweep(m, sp, bl, old, far)
  expect_true(in_link_range(m$family, kept$eta))
  expect_lt(obj(kept), obj(old))
  expect_equal(kept$eta - old$eta, rep(kept$a - old$a, n))
  # From the mode, only a part within the rounding allowance survives
  up     <- old
  up$a   <- up$a + 0.1
  up$eta <- up$eta + 0.1
  kept   <- keep_sweep(m, sp, bl, mode, up)
  expect_lte(obj(kept), obj(mode) + loss_allowance(obj(mode)))
  # A sweep that removes columns is a projection: it is kept whole even if
  # the objective rises, and refused only if it leaves the link's range
  gone <- bl
  gone$smoothers[[1L]]$keep[which(mode$beta[[1L]] != 0)[1L]] <- FALSE
  expect_identical(keep_sweep(m, sp, gone, mode, up), up)
  expect_null(keep_sweep(m, sp, gone, mode, far))
  # A rise within the rounding allowance is kept whole
  tiny     <- mode
  tiny$a   <- tiny$a + 1e-6
  tiny$eta <- tiny$eta + 1e-6
  expect_gt(obj(tiny), obj(mode))
  expect_identical(keep_sweep(m, sp, bl, mode, tiny), tiny)
  # From outside the range of the link a sweep is kept whole
  expect_identical(keep_sweep(m, sp, bl, far, mode), mode)
  # A halved sweep moves the coefficients with eta
  V       <- m$smoothers[[1L]]$V
  v       <- rep(0.02 / max(abs(V %*% rep(1, ncol(V)))), ncol(V))
  shifted <- function(s, k) {
    s$beta[[1L]] <- s$beta[[1L]] + k * v
    s$eta        <- s$eta + k * as.vector(V %*% v)
    s
  }
  kept <- keep_sweep(m, sp, bl, shifted(mode, 1), shifted(mode, -20))
  expect_false(identical(kept, shifted(mode, -20)))
  expect_equal(kept$eta,
               m$offset + kept$a + as.vector(V %*% kept$beta[[1L]]))
  # The objective does not rise over the sweeps of a fit
  start <- list(a = m$a, beta = lapply(m$smoothers, function(sm)
    numeric(ncol(sm$V))), eta = m$offset + m$a)
  objs <- vapply(1:6, function(k)
    obj(penalised_fit(m, sp, start, 1e-12, max_sweeps = k)$state), numeric(1))
  expect_true(all(diff(objs) <= sqrt(.Machine$double.eps) * objs[-1]))
  # A fit started where a full sweep leaves the range of the link
  high     <- mode
  high$a   <- high$a + 3 * max(mode$eta)
  high$eta <- high$eta + 3 * max(mode$eta)
  fit      <- penalised_fit(m, sp, high, 1e-8)
  expect_true(fit$ok)
  expect_equal(fit$state$eta, mode$eta, tolerance = 1e-6)
  # With no part of a sweep kept, the fit ends at the state it had
  local({
    local_mocked_bindings(keep_sweep = function(...) NULL)
    fit <- penalised_fit(m, sp, high, 1e-8)
    expect_false(fit$ok)
    expect_identical(fit$state, high)
    expect_identical(fit$sweeps, 1L)
  })
  # A sweep cut back to nothing does not look converged
  local({
    local_mocked_bindings(keep_sweep = function(m, sp, blocks, old, new) old)
    expect_false(penalised_fit(m, sp, high, 1e-8, max_sweeps = 5L)$ok)
  })
})

test_that("removes_columns finds nonzero coefficients outside the blocks", {
  blocks <- function(...) list(smoothers = list(list(keep = c(...))))
  st     <- list(beta = list(c(1, 0)))
  expect_false(removes_columns(blocks(TRUE, TRUE), st))
  expect_false(removes_columns(blocks(TRUE, FALSE), st))
  expect_true(removes_columns(blocks(FALSE, TRUE), st))
  # A smoother without a block
  expect_true(removes_columns(list(smoothers = list(NULL)), st))
  expect_false(removes_columns(list(smoothers = list(NULL)),
                               list(beta = list(c(0, 0)))))
})

test_that("penalised_fit perturbs a singular block and keeps its fixed point", {
  # Two copies of a column at lambda = 0. sum(x^2) = 100 has an exact root,
  # so that rounding does not hide the singularity from chol()
  set.seed(33)
  n <- 40
  x <- rep(1:2, length.out = n)
  V <- cbind(x, x)
  w <- rep(1, n)
  y <- 1 + 0.5 * x + stats::rnorm(n)
  sm <- list(V = V, VtWV = crossprod(V), d = c(0, 0), lid = c(1L, 1L),
             beta = c(0, 0), label = "s", learners = data.frame(kind = "null"))
  m  <- list(y = y, w = w, offset = rep(0, n),
             family = new_family(stats::gaussian()), has_int = TRUE, a = 0,
             eta = rep(0, n), smoothers = list(sm))
  sp <- list(list(status = "zero", rho = NA_real_))
  bl <- model_blocks(m, sp, FALSE)
  # The sampler's factorisation stays strict
  expect_error(block_system(bl$smoothers[[1L]], w, "a block"),
               class = "b2b_not_pd")
  st <- list(a = 0, beta = list(c(0, 0)), eta = rep(0, n))
  md <- penalised_fit(m, sp, st, 1e-9)
  expect_true(md$ok)
  # A stationary point of the loss, not of a ridge-penalised one
  lm_fit <- unname(stats::fitted(stats::lm(y ~ x)))
  expect_equal(md$state$eta, lm_fit, tolerance = 1e-6)
  # Also under a large ridge, for smoother blocks and for the linear block
  local_mocked_bindings(perturbed_chol = function(M, where) {
    A       <- M
    diag(A) <- diag(A) + 0.1 * diag(M)
    list(R = chol(A), ridge = 0.1 * diag(M))
  })
  for (joint in c(FALSE, TRUE)) {
    md <- penalised_fit(m, sp, st, 1e-9, joint_linear = joint)
    expect_true(md$ok)
    expect_equal(md$state$eta, lm_fit, tolerance = 1e-6)
  }
})

test_that("a singular one-step system skips the learner for that sweep", {
  mod <- fit_aic(f_col, sim_collinear())
  m   <- model_as_smoothers(mod)
  local_mocked_bindings(one_step_proposal = function(...)
    signal_numerical("b2b_not_pd", "a block: matrix is not positive definite"))
  s1 <- map_smoothing_parameters(m)
  expect_identical(s1$status, "stalled")
  expect_identical(s1$sweeps, 1L)
})

test_that("unpenalised learners start at boosting's shrinkage (Poisson)", {
  skip_on_cran()
  set.seed(5)
  n   <- 300
  d   <- data.frame(x0 = stats::runif(n), z = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.5 + sin(2 * pi * d$x0) + 0.3 * d$z))
  m   <- model_as_smoothers(fit_aic(y ~ t2(x0, bs = "ps", k = 8) + z, d,
                                    family = stats::poisson()))
  sp0 <- starting_lambdas(m)
  wq  <- iwls_working(m, m$eta, "test")
  n_finite <- 0L
  n_zero   <- 0L
  for (j in seq_along(m$smoothers)) {
    sm <- m$smoothers[[j]]
    for (l in which(sm$learners$kind == "null" & !sm$learners$never)) {
      on_l <- sm$lid == l
      x    <- sm$V[, on_l]
      b    <- sm$beta[on_l]
      cc   <- sum(wq$wz * x^2)
      r    <- sum(wq$wz * x * wq$zr)
      if (sp0[[j]]$status[l] == "finite") {
        # one update of the learner alone returns its boosted coefficient
        expect_equal((r + cc * b) / (cc + exp(sp0[[j]]$rho[l])), b,
                     tolerance = 1e-10)
        n_finite <- n_finite + 1L
      } else {
        expect_identical(sp0[[j]]$status[l], "zero")
        expect_lte(r / b, 0)
        n_zero <- n_zero + 1L
      }
    }
  }
  expect_gt(n_finite, 0L)
})

test_that("Stage 1 starts from the penalised fit at tr(F_k) = edf~_k", {
  skip_on_cran()
  mod  <- fit_aic(f_col, sim_collinear())
  m    <- model_as_smoothers(mod)
  sp0  <- starting_lambdas(m)
  boosted <- list(a = m$a, beta = lapply(m$smoothers, `[[`, "beta"),
                  eta = m$eta)
  calls <- list()
  fit   <- penalised_fit
  local_mocked_bindings(penalised_fit = function(m, sp, state, ...) {
    calls[[length(calls) + 1L]] <<- list(sp = sp, state = state)
    fit(m, sp, state, ...)
  })
  # an unpenalised learner starts where one update of it alone, from the
  # boosted fit, returns its boosted coefficient (Gaussian: x^T W x = 1)
  n_null <- 0L
  for (j in seq_along(m$smoothers)) {
    sm <- m$smoothers[[j]]
    for (l in which(sm$learners$kind == "null" &
                    sp0[[j]]$status == "finite")) {
      on_l <- sm$lid == l
      r    <- sum(sm$V[, on_l] * m$w * (m$y - m$eta)) + sm$beta[on_l]
      expect_equal(r / (1 + exp(sp0[[j]]$rho[l])), sm$beta[on_l],
                   tolerance = 1e-10)
      n_null <- n_null + 1L
    }
  }
  expect_gt(n_null, 0L)
  # a selected penalised learner starts where its trace is its edf~_k
  n_pen <- 0L
  for (j in seq_along(m$smoothers)) {
    sm <- m$smoothers[[j]]
    for (l in which(sm$learners$kind == "pen" & !sm$learners$never)) {
      d <- sm$d[sm$lid == l]
      if (sm$learners$edf[l] >= length(d)) {
        expect_identical(sp0[[j]]$status[l], "zero")
        next
      }
      expect_identical(sp0[[j]]$status[l], "finite")
      expect_equal(sum(1 / (1 + exp(sp0[[j]]$rho[l]) * d)),
                   sm$learners$edf[l], tolerance = 1e-8)
      n_pen <- n_pen + 1L
    }
  }
  expect_gt(n_pen, 0L)
  s1 <- map_smoothing_parameters(m)
  expect_identical(calls[[1]]$sp, sp0)
  expect_identical(calls[[1]]$state, boosted)
  expect_true(s1$status %in% c("converged", "stalled"))
  local_mocked_bindings(penalised_fit = function(m, sp, state, ...)
    list(state = state, sweeps = 1000L, ok = FALSE))
  e <- tryCatch(map_smoothing_parameters(m), map_approx_failed = function(e) e)
  expect_s3_class(e, "map_approx_failed")
  expect_identical(e$status, "no_mode")
  expect_match(conditionMessage(e), "initial lambda")
})

test_that("the smoother blocks reproduce the boosted fit by construction", {
  skip_on_cran()
  set.seed(31)
  n <- 150
  d <- data.frame(x0 = stats::runif(n), x1 = stats::runif(n),
                  lin = stats::rnorm(n), expo = stats::runif(n, 0.5, 2),
                  g = factor(sample(c("a", "b", "c"), n, TRUE)),
                  nt = sample(3:8, n, TRUE), wt = stats::runif(n, 0.5, 2))
  d$y  <- sin(2 * pi * d$x0) + d$x1 + (d$g == "b") * d$x0 +
    stats::rnorm(n, 0, 0.3)
  d$yp <- stats::rpois(n, d$expo * exp(sin(2 * pi * d$x0)))
  d$s  <- stats::rbinom(n, d$nt, stats::plogis(sin(2 * pi * d$x0) + d$lin))
  fit <- function(formula, family = stats::gaussian(), lc = list(), ...) {
    suppressMessages(gbgam(formula, data = d, family = family,
      learner_control = do.call(gb_learner_ctrl,
                                c(list(step_size = 0.1), lc)),
      stop_control = gb_stop_ctrl("none", max_iter = 80L),
      edf_control  = gb_edf_ctrl("none"), verbose = FALSE, ...))
  }
  f_x  <- y ~ t2(x0, bs = "ps", k = 8) + t2(x1, bs = "ps", k = 8) + lin
  mods <- list(
    gaussian = fit(f_x),
    poisson  = fit(yp ~ t2(x0, bs = "ps", k = 8) + offset(log(expo)),
                   stats::poisson()),
    binomial = fit(cbind(s, nt - s) ~ t2(x0, bs = "ps", k = 8) + lin,
                   stats::binomial(), weights = d$wt),
    by       = fit(y ~ g + t2(x0, bs = "ps", k = 8) +
                     t2(x0, by = g, bs = "ps", k = 8)),
    tensor   = fit(y ~ t2(x0, x1, bs = "ps", k = 5)),
    padding  = fit(f_x, lc = list(intercept_padding = TRUE)),
    after    = fit(f_x, lc = list(learners_to_fit_after_each_iter = 1L,
                                  excluded_learner_indices = c(1L, 2L)))
  )
  # learner 2 (lin) is excluded and not refitted after each iteration
  lin <- model_as_smoothers(mods$after)$smoothers[["lin"]]
  expect_true(lin$learners$excluded)
  for (nm in names(mods)) {
    mod <- mods[[nm]]
    L   <- mod$learners
    m   <- model_as_smoothers(mod)
    # the invariants model_as_smoothers() relies on
    is_int <- vapply(L, function(l) identical(l$var_label, "(Intercept)"),
                     TRUE)
    expect_identical(sum(is_int), 1L, info = nm)
    expect_true(all(L[[which(is_int)]]$design == 1), info = nm)
    expect_identical(length(m$offset), length(m$y), info = nm)
    for (l in L) {
      expect_true(inherits(l, c("demmler_reinsch", "unpenalised")), info = nm)
      if (inherits(l, "demmler_reinsch")) expect_true(all(l$D > 0), info = nm)
    }
    # the blocks reproduce the boosted linear predictor
    err <- max(abs(m$eta - mod$linear.predictors))
    expect_lt(err, 1e-10 * max(1, abs(m$eta)), label = paste(nm, "eta error"))
    # each learner's columns are W-orthonormal (unpenalised ones rescaled)
    for (sm in m$smoothers) {
      for (l in unique(sm$lid)) {
        i <- sm$lid == l
        expect_lt(max(abs(sm$VtWV[i, i] - diag(sum(i)))), 1e-10,
                  label = paste(nm, sm$label, "V^T W V - I"))
      }
    }
    # the smooth blocks are the smooths of mod$smooth, in order
    sl <- vapply(mod$smooth, `[[`, "", "label")
    expect_identical(names(m$smoothers)[names(m$smoothers) %in% sl], sl,
                     info = nm)
  }
  # the pivoted QR in drop_aliased_flat_columns() never drops the
  # intercept: with all
  # linear learners flat, the by-factor model has aliased level constants
  m  <- model_as_smoothers(mods$by)
  sp <- free_null_sp(map_smoothing_parameters(m)$sp, m)
  al <- drop_aliased_flat_columns(m, sp)
  expect_gt(length(al$aliased), 0L)
  # after it, M_p is the rank of the flat columns
  flat <- cbind(1, do.call(cbind, Map(function(sm, st) {
    pr <- smoother_prior(sm, st)
    sm$V[, pr$keep, drop = FALSE][, pr$cDlam == 0, drop = FALSE]
  }, m$smoothers, al$sp)))
  expect_identical(penalty_null_dim(m, al$sp), ncol(flat))
  expect_identical(penalty_null_dim(m, al$sp), qr(sqrt(m$w) * flat)$rank)
})

test_that("M_p counts the intercept and the columns at lambda = 0", {
  mod <- fit_fixed(y ~ z + t2(x, k = 6, bs = "cr"), simulate_gaussian_2d())
  m   <- model_as_smoothers(mod)
  sp  <- starting_lambdas(m)
  n_zero <- sum(vapply(seq_along(m$smoothers), function(j) {
    sm <- m$smoothers[[j]]
    sum(sp[[j]]$status[sm$lid] == "zero")
  }, numeric(1)))
  expect_identical(penalty_null_dim(m, sp), 1L + as.integer(n_zero))
  # every unpenalised learner flat: the intercept, z and the null space
  # of t2(x)
  free  <- free_null_sp(sp, m)
  n_lin <- sum(vapply(m$smoothers, function(sm)
    sum(sm$learners$kind == "null"), integer(1)))
  expect_identical(penalty_null_dim(m, free), 1L + n_lin)
  # learners at lambda = Inf are not counted
  inf <- lapply(free, function(st) {
    st$status[] <- "inf"
    st
  })
  expect_identical(penalty_null_dim(m, inf), 1L)
})

test_that("without edf tracking, Stage 1 starts at boosting's lambda", {
  mod <- fit_fixed(y ~ t2(x, k = 8, bs = "cr"), simulate_gaussian_1d(),
                   edf = "none")
  m   <- model_as_smoothers(mod)
  sp  <- starting_lambdas(m)
  sm  <- m$smoothers[[1]]
  l   <- which(sm$learners$kind == "pen")
  expect_identical(sm$learners$edf[l], 0)
  expect_identical(sp[[1]]$status[l], "finite")
  expect_identical(sp[[1]]$rho[l], sm$learners$rho[l])
})

test_that("the joint linear block reaches the full-design mode", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  m   <- model_as_smoothers(mod)
  s1  <- map_smoothing_parameters(m)
  sp  <- free_null_sp(s1$sp, m)
  joint <- penalised_fit(m, sp, s1$state, s1$tol, 1e4L, joint_linear = TRUE)
  plain <- penalised_fit(m, sp, s1$state, s1$tol, 1e4L)
  expect_true(joint$ok)
  expect_lt(joint$sweeps, plain$sweeps)
  pri <- Map(smoother_prior, m$smoothers, sp)
  X   <- cbind(1, do.call(cbind, Map(function(sm, p)
    sm$V[, p$keep, drop = FALSE], m$smoothers, pri)))
  H   <- crossprod(X, m$w * X)
  diag(H) <- diag(H) + c(0, unlist(lapply(pri, `[[`, "cDlam")))
  b   <- solve(H, crossprod(X, m$w * (m$y - m$offset)))
  expect_lt(max(abs(joint$state$eta - (m$offset + as.vector(X %*% b)))), 1e-6)
})

# A collinear fit with the penalised learner of t2(x0) almost a straight
# line, at its posterior mode
low_edf_mode <- function() {
  mod <- fit_aic(f_col, sim_collinear())
  m   <- model_as_smoothers(mod)
  s1  <- map_smoothing_parameters(m)
  sp  <- s1$sp
  j   <- which(names(m$smoothers) == "t2(x0)")
  pen <- which(m$smoothers[[j]]$learners$kind == "pen")
  sp[[j]]$status[pen] <- "finite"
  sp[[j]]$rho[pen]    <- 15
  md  <- penalised_fit(m, sp, s1$state, s1$tol, 1e4L)
  list(mod = mod, m = m, sp = sp, state = md$state, j = j, pen = pen)
}

# edf of every finite penalised learner, one vector per term
pen_edf <- function(m, sp, state) {
  lapply(seq_along(m$smoothers), function(j) {
    sm <- m$smoothers[[j]]
    tw <- smoother_working_model(m, j, state)
    on <- sm$learners$kind == "pen" & !sm$learners$never &
      sp[[j]]$status == "finite"
    stats::setNames(smoother_shrinkage(sm, tw, sp[[j]])$edf[on],
                    sm$learners$label[on])
  })
}

test_that("the edf floor lifts low-edf learners to it and no others", {
  skip_on_cran()
  lo  <- low_edf_mode()
  m   <- lo$m
  before <- pen_edf(m, lo$sp, lo$state)
  expect_lt(before[[lo$j]][[1L]], 1)
  us  <- undersmooth_to_edf_floor(m, lo$sp, lo$state, edf_floor = 1)
  after <- pen_edf(m, us$sp, lo$state)
  cat(sprintf("edf floor: %d learners lifted, max |edf_to - 1| = %.2g\n",
              nrow(us$moved), max(abs(us$moved$edf_to - 1))))
  expect_identical(names(us$moved), c("term", "learner", "rho_from",
                                      "rho_to", "edf_from", "edf_to"))
  expect_true("t2(x0)" %in% us$moved$term)
  expect_true(all(us$moved$edf_from < 1))
  expect_true(all(us$moved$rho_to < us$moved$rho_from))
  expect_lt(max(abs(us$moved$edf_to - 1)), 1e-6)
  # a term with one lifted learner: its edf is now the floor
  expect_equal(after[[lo$j]][[1L]], 1, tolerance = 1e-6)
  # learners at or above the floor keep their lambda
  for (j in seq_along(m$smoothers)) {
    up <- names(before[[j]])[before[[j]] >= 1]
    i  <- match(up, m$smoothers[[j]]$learners$label)
    expect_identical(us$sp[[j]]$rho[i], lo$sp[[j]]$rho[i])
  }
})

test_that("an edf floor a learner cannot reach is a bb_failed", {
  skip_on_cran()
  lo  <- low_edf_mode()
  err <- tryCatch(undersmooth_to_edf_floor(lo$m, lo$sp, lo$state,
                                           edf_floor = 1e3),
                  bb_failed = function(e) e)
  expect_s3_class(err, "bb_failed")
  expect_identical(err$status, "floor_unreachable")
})

test_that("bayesian_backfitting() undersmooths by the edf floor", {
  skip_on_cran()
  lo <- low_edf_mode()
  expect_error(bayesian_backfitting(lo$mod, edf_floor = 0), "edf_floor")
  expect_error(bayesian_backfitting(lo$mod, edf_floor = -1), "edf_floor")
  expect_error(bayesian_backfitting(lo$mod, edf_floor = Inf), "edf_floor")
  expect_error(bayesian_backfitting(lo$mod, edf_floor = c(1, 2)),
               "edf_floor")
  set.seed(6)
  bb <- bayesian_backfitting(lo$mod, n_samples = 50L, n_burnin = 10L,
                             undersmooth = TRUE, edf_floor = 1.5)
  expect_true(bb$uq_control$control$undersmooth)
  expect_identical(bb$uq_control$control$edf_floor, 1.5)
  us <- bb$bb_info$undersmoothed
  expect_false(is.null(us))
  expect_lt(max(abs(us$edf_to - 1.5)), 1e-6)
  # the posterior uses the lifted lambda
  lam <- bb$bb_info$lambda
  i   <- match(paste(us$term, us$learner), paste(lam$term, lam$learner))
  expect_false(anyNA(i))
  expect_identical(lam$rho[i], us$rho_to)
  set.seed(6)
  plain <- bayesian_backfitting(lo$mod, n_samples = 50L, n_burnin = 10L)
  expect_null(plain$bb_info$undersmoothed)
  # phi is estimated at the undersmoothed mode
  m  <- model_as_smoothers(lo$mod)
  s1 <- map_smoothing_parameters(m)
  r  <- sample_posterior(m, s1, "free", TRUE, 1.5, NULL, 2L, 0L, FALSE, NULL)
  st <- penalised_fit(m, r$sp, s1$state, s1$tol, s1$max_inner,
                      joint_linear = TRUE)$state
  expect_equal(r$phi, dispersion_at_mode(m, r$sp, st), tolerance = 1e-6)
  expect_equal(bb$bb_info$phi, r$phi, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(bb$bb_info$phi, plain$bb_info$phi)))
})


# ---- Linear learners and aliased columns -------------------------------------

test_that("freed linear learners widen the intervals of collinear smooths", {
  skip_on_cran()
  mod <- fit_aic(f_col, sim_collinear())
  set.seed(1)
  free <- bayesian_backfitting(mod, n_samples = 600L, n_burnin = 100L)
  set.seed(1)
  kept <- bayesian_backfitting(mod, n_samples = 600L, n_burnin = 100L,
                               linear = "kept")
  se_free <- colMeans(smooth_intervals(free)$se)
  se_kept <- colMeans(smooth_intervals(kept)$se)
  expect_gt(se_free[["t2(x1)"]] / se_kept[["t2(x1)"]], 1.3)
  expect_gt(se_free[["t2(x2)"]] / se_kept[["t2(x2)"]], 1.3)
  expect_identical(free$bb_info$linear, "free")
  expect_true(free$bb_info$discrepancy >= 0)
  expect_identical(free$uq_control$control$linear, "free")
  lam <- free$bb_info$lambda
  expect_true(all(lam$status[grepl("nullspace|n[0-9]$", lam$learner)] ==
                    "zero"))
})

test_that("aliased flat columns are removed and recorded", {
  skip_on_cran()
  set.seed(3)
  n <- 200
  d <- data.frame(x = stats::runif(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  d$y <- sin(2 * pi * d$x) + (d$g == "b") * d$x + stats::rnorm(n, 0, 0.3)
  mod <- fit_aic(y ~ t2(x, bs = "ps", k = 8) + t2(x, by = g, bs = "ps", k = 8),
                 d)
  set.seed(4)
  bb <- bayesian_backfitting(mod, n_samples = 100L, n_burnin = 20L)
  expect_gt(length(bb$bb_info$aliased), 0L)
  expect_true(all(is.finite(bb$Vp)))
  expect_true(all(diag(bb$Vp) >= 0))

  # Metropolis-Hastings: the removed learner is held at zero, so it has no
  # acceptance rate; the joint linear block has one
  d$y <- stats::rpois(n, exp(sin(2 * pi * d$x) + (d$g == "b") * d$x))
  mp  <- fit_aic(y ~ t2(x, bs = "ps", k = 8) + t2(x, by = g, bs = "ps", k = 8),
                 d, family = stats::poisson())
  set.seed(4)
  bp  <- bayesian_backfitting(mp, n_samples = 100L, n_burnin = 20L)
  lab <- vapply(mp$learners, `[[`, "", "var_label")
  out <- lab %in% bp$bb_info$aliased
  expect_true(any(out))
  expect_true(all(is.na(bp$mh_accept_rate[out])))
  expect_false(anyNA(bp$mh_accept_rate[!out]))
  expect_true(0 < bp$bb_info$linear_accept && bp$bb_info$linear_accept <= 1)
})

test_that("a fit with no learner selected has a posterior", {
  skip_on_cran()
  set.seed(1)
  d  <- data.frame(x = stats::runif(100), y = stats::rnorm(100))
  m0 <- fit_aic(y ~ t2(x, bs = "ps", k = 8), d)
  expect_identical(m0$mstop, 0L)
  for (lin in c("free", "kept")) {
    set.seed(2)
    b0 <- bayesian_backfitting(m0, n_samples = 50L, n_burnin = 10L,
                               linear = lin)
    expect_identical(b0$bb_info$stage1, "no_learners")
    expect_gt(b0$Vp["(Intercept)", "(Intercept)"], 0)
  }
})

test_that("a fit with a large intercept has a posterior", {
  skip_on_cran()
  # the tolerance follows the size of the whole fit, intercept included
  set.seed(2)
  d2 <- data.frame(x = stats::runif(200), lin = stats::rnorm(200))
  d2$y <- 1e8 + sin(2 * pi * d2$x) + 0.5 * d2$lin + stats::rnorm(200, 0, 0.3)
  for (mi in c(1L, 5L)) {
    m1 <- gbgam(y ~ t2(x, bs = "ps", k = 8) + lin, data = d2,
                stop_control = gb_stop_ctrl("none", max_iter = mi),
                edf_control = gb_edf_ctrl("none"), verbose = FALSE)
    set.seed(3)
    b1 <- bayesian_backfitting(m1, n_samples = 50L, n_burnin = 10L)
    expect_true(b1$bb_info$stage1 %in% c("converged", "stalled"))
  }
})

test_that("a smoother whose learners are jointly aliased has a posterior", {
  # 30 distinct locations for a t2 tensor with 35 columns: each learner is
  # W-orthonormal, but the smoother's columns together have rank 28
  set.seed(2)
  n <- 150
  g <- expand.grid(x = seq(0, 1, length.out = 7), z = seq(0, 1, length.out = 7))
  d <- g[rep(sample(nrow(g), 30), length.out = n), ]
  d$y <- rpois(n, exp(1 + sin(2 * pi * d$x) * d$z))
  m  <- gbgam(y ~ t2(x, z, k = c(6, 6)), data = d, family = poisson(),
              verbose = FALSE)
  tb <- model_as_smoothers(m)
  V  <- tb$smoothers[[1]]$V
  expect_lt(qr(V * sqrt(tb$w))$rank, ncol(V))
  b <- bayesian_backfitting(m, n_samples = 50, n_burnin = 10)
  expect_true(all(is.finite(b$Vp)))
  # The root used for the discrepancy reproduces the singular V^T W V
  R0 <- matrix_root(tb$smoothers[[1]]$VtWV)
  expect_equal(crossprod(R0), tb$smoothers[[1]]$VtWV, tolerance = 1e-10)
})
