# Boost to backfit, part 2: Stage 2, Bayesian backfitting
# (sample_posterior()), and the driver boost_to_backfit(). The vocabulary,
# the data structures and the helpers of both stages are in
# uq_bayes_map_est.R.


# ---- Stage 2: Bayesian backfitting -------------------------------------------

#' Stage 2: posterior mode, undersmoothing, dispersion and draws
#'
#' At the Stage-1 smoothing parameters s1$sp, with the linear learners freed
#' for linear = "free". Every penalised fit is allowed s1$max_inner fitting
#' sweeps, as in Stage 1. phi: the dispersion, or NULL to estimate it
#' (dispersion_at_mode()) at the final mode, after any undersmoothing.
#' start: the time (proc.time() "elapsed") from which max_runtime_mins
#' counts.
#' @return list(sp, phi, aliased, undersmoothed, draws (from
#'   backfitting_sampler())).
#' @noRd
sample_posterior <- function(m, s1, linear, undersmooth, edf_floor, phi,
                             n_samples, n_burnin, verbose, max_runtime_mins,
                             start = proc.time()[["elapsed"]]) {
  force(start)
  sp    <- s1$sp
  joint <- identical(linear, "free")
  if (joint) {
    sp <- Map(function(spj, sm) {
      lr <- sm$learners
      # freed: every unpenalised learner except zero-norm ones and those the
      # user excluded (model_as_smoothers())
      freed <- lr$kind == "null" & !(lr$excluded & lr$never) & !lr$zero_norm
      spj$status[freed] <- "zero"
      spj$rho[freed]    <- NA_real_
      spj
    }, sp, m$smoothers)
  }
  al  <- drop_aliased_flat_columns(m, sp)
  sp  <- al$sp
  fit_mode <- function(sp, st, where) {
    md <- penalised_fit(m, sp, st, s1$tol, s1$max_inner, where,
                        joint_linear = joint)
    if (!md$ok)
      stop(backfitting_failed_condition(
        paste("no", where, "at the smoothing parameters"), "no_mode"))
    md$state
  }
  st  <- fit_mode(sp, s1$state, "posterior mode")
  us  <- NULL
  if (undersmooth) {
    us <- undersmooth_to_edf_floor(m, sp, st, edf_floor)
    if (!is.null(us$moved)) {
      sp <- us$sp
      st <- fit_mode(sp, st, "undersmoothed posterior mode")
    }
  }
  phi <- dispersion_at_mode(m, sp, st, phi)
  list(sp = sp, phi = phi, aliased = al$aliased, undersmoothed = us$moved,
       draws = backfitting_sampler(m, sp, st, n_samples, n_burnin, phi, joint,
                                   verbose, max_runtime_mins, start))
}

#' Remove exactly aliased columns with a flat prior
#'
#' The intercept and every column of a learner at lambda = 0 have a flat
#' prior; where they are linearly dependent (e.g. the level constants of a
#' `by = factor` smooth and the intercept) the posterior is improper. As
#' `mgcv` does for identifiability, a rank-revealing (pivoted) QR of their
#' weighted columns, intercept first, finds the aliased ones; their
#' learners are removed (lambda = Inf) and recorded. Removing them does not
#' change the model space.
#' @return list(sp, aliased) with aliased a character vector of learners.
#' @noRd
drop_aliased_flat_columns <- function(m, sp, tol = 1e-7) {
  cols <- list()
  who  <- list()
  if (m$has_int) {
    cols <- list(rep(1, length(m$y)))
    who  <- list(c(0L, 0L))
  }
  for (j in seq_along(m$smoothers)) {
    sm   <- m$smoothers[[j]]
    flat <- which(sp[[j]]$status[sm$lid] == "zero")
    for (i in flat) {
      cols[[length(cols) + 1L]] <- sm$V[, i]
      who[[length(who) + 1L]]  <- c(j, sm$lid[i])
    }
  }
  if (length(cols) < 2L) return(list(sp = sp, aliased = character()))
  X  <- sqrt(m$w) * do.call(cbind, cols)
  qx <- qr(X, tol = tol)
  if (qx$rank == ncol(X)) return(list(sp = sp, aliased = character()))
  aliased <- character()
  for (i in qx$pivot[seq(qx$rank + 1L, ncol(X))]) {
    j  <- who[[i]][1L]
    l  <- who[[i]][2L]
    sm <- m$smoothers[[j]]
    if (sm$learners$kind[l] != "null")
      stop(backfitting_failed_condition(
        paste(sm$learners$label[l], "has a flat prior and is aliased with",
              "other flat columns"), "aliased"))
    sp[[j]]$status[l] <- "inf"
    sp[[j]]$rho[l]    <- NA_real_
    aliased <- c(aliased, sm$learners$label[l])
  }
  list(sp = sp, aliased = unique(aliased))
}

#' Dispersion of the backfitting model at its mode
#'
#' 1 for Poisson and binomial; otherwise `phi` if given, else the REML
#' estimate at the fixed lambda (Wood, 2011, eq. 4, approach (a)): the root
#' of the score D_p / (2 phi^2) + l_s'(phi) + M_p / (2 phi) = 0, with D_p the
#' penalised deviance at the mode (penalised_deviance()), M_p the number of
#' columns with a flat prior (penalty_null_dim()) and l_s the saturated
#' log-likelihood (mgcv::fix.family.ls()). For the Gaussian and inverse
#' Gaussian families l_s' = -n_+ / (2 phi), so phi = D_p / (n_+ - M_p), with
#' n_+ the number of rows with non-zero weight; otherwise (Gamma) a root
#' search in log(phi) from that value (reml_dispersion()). \mathcal{O}(n + P).
#' @noRd
dispersion_at_mode <- function(m, sp, state, phi = NULL) {
  fam <- m$family
  if (fam$family %in% c("poisson", "binomial")) return(1)
  if (!is.null(phi)) return(phi)
  n   <- sum(m$w > 0)
  M_p <- penalty_null_dim(m, sp)
  if (!(n - M_p > 0))
    stop(backfitting_failed_condition(
      sprintf(paste("no residual degrees of freedom for the dispersion",
                    "(n = %d, flat columns = %d)"), n, M_p),
      "no_residual_df"))
  D_p <- penalised_deviance(m, sp, state)
  phi <- D_p / (n - M_p)
  if (!(fam$family %in% c("gaussian", "inverse.gaussian")) &&
      is.finite(phi) && phi > 0)
    phi <- reml_dispersion(fam, m$y, m$w, D_p, M_p, phi)
  if (!(is.finite(phi) && phi > 0))
    stop(backfitting_failed_condition(
      sprintf("the dispersion estimate is %g (penalised deviance %g)", phi,
              D_p),
      "zero_dispersion"))
  phi
}

#' Penalised deviance D_p at the mode: the deviance over the rows with
#' non-zero weight plus sum_k lambda_k sum(d_k beta_k^2)
#' (smoother_prior()$cDlam_all, unpenalised columns at unit W-norm)
#' @noRd
penalised_deviance <- function(m, sp, state) {
  pos <- m$w > 0
  dev <- sum(m$family$dev.resids(m$y[pos], m$family$linkinv(state$eta[pos]),
                                 m$w[pos]))
  pen <- sum(vapply(seq_along(m$smoothers), function(j) {
    cDlam <- smoother_prior(m$smoothers[[j]], sp[[j]])$cDlam_all
    sum(cDlam * state$beta[[j]]^2)
  }, numeric(1)))
  dev + pen
}

#' Dimension M_p of the penalty null space: the intercept (if any) and every
#' column kept in the posterior with a zero penalty (lambda_k = 0); removed
#' (lambda = Inf) columns are not counted. After drop_aliased_flat_columns()
#' these columns have full rank.
#' @noRd
penalty_null_dim <- function(m, sp) {
  as.integer(m$has_int) + sum(vapply(seq_along(m$smoothers), function(j)
    sum(smoother_prior(m$smoothers[[j]], sp[[j]])$cDlam == 0), integer(1)))
}

#' REML dispersion at fixed lambda: the root in log(phi) of
#' g(phi) = D_p / (2 phi^2) + l_s'(phi) + M_p / (2 phi), which is positive
#' for small and negative for large phi; searched from phi0 outwards
#' @noRd
reml_dispersion <- function(fam, y, w, D_p, M_p, phi0, tol = 1e-10) {
  ls    <- mgcv::fix.family.ls(fam)$ls
  n     <- length(y)
  score <- function(r) {
    phi <- exp(r)
    D_p / (2 * phi^2) + ls(y, w, n, phi)[2] + M_p / (2 * phi)
  }
  tryCatch(exp(stats::uniroot(score, log(phi0) + c(-1, 1),
                              extendInt = "downX", tol = tol)$root),
           error = function(e) NA_real_)
}

#' Undersmoothing by an edf floor
#'
#' Every penalised learner that boosting selected, with a finite lambda and
#' an edf below `edf_floor` at the mode (smoother_shrinkage(): its smoother
#' in the working model, the other smoothers fixed), gets its log(lambda)
#' lowered until its edf equals `edf_floor`. Learners of one smoother are
#' taken in turn, each at the lambda of those already lifted. Learners at or
#' above the floor, and those at lambda = Inf or 0, are untouched.
#' @param max_drop How far below its own value log(lambda) is searched.
#' @return list(sp, moved): moved has one row per learner lifted (term,
#'   learner, rho_from, rho_to, edf_from, edf_to), or is NULL.
#' @noRd
undersmooth_to_edf_floor <- function(m, sp, state, edf_floor, tol = 1e-8,
                                     max_drop = 60) {
  moved <- list()
  for (j in seq_along(m$smoothers)) {
    sm <- m$smoothers[[j]]
    tw <- smoother_working_model(m, j, state, paste("edf floor,", sm$label))
    for (l in which(sm$learners$kind == "pen" & !sm$learners$never &
                    sp[[j]]$status == "finite")) {
      where <- sprintf("edf floor, %s, learner %s", sm$label,
                       sm$learners$label[l])
      gap <- function(r) {                   # edf minus the floor
        st        <- sp[[j]]
        st$rho[l] <- r
        smoother_shrinkage(sm, tw, st, where)$edf[l] - edf_floor
      }
      rho0 <- sp[[j]]$rho[l]
      g0   <- gap(rho0)
      if (g0 >= 0) next                      # at or above the floor
      lo <- rho0
      repeat {                  # lower log(lambda) until edf passes the floor
        lo   <- lo - 2
        g_lo <- gap(lo)
        if (g_lo > 0 || lo < rho0 - max_drop) break
      }
      if (!(g_lo > 0))
        stop(backfitting_failed_condition(
          sprintf(paste("%s: edf %.4g at log(lambda) %.4g, still below the",
                        "floor %.4g"), where, g_lo + edf_floor, lo,
                  edf_floor),
          "floor_unreachable"))
      g_up <- gap(lo + 2)
      if (!(g_up <= 0))
        stop(backfitting_failed_condition(
          sprintf(paste("%s: no sign change of edf - floor in log(lambda)",
                        "[%.4g, %.4g]"), where, lo, lo + 2),
          "floor_no_root"))
      rho_new <- stats::uniroot(gap, c(lo, lo + 2), f.lower = g_lo,
                                f.upper = g_up, tol = tol)$root
      sp[[j]]$rho[l] <- rho_new
      moved[[length(moved) + 1L]] <- data.frame(
        term = sm$label, learner = sm$learners$label[l], rho_from = rho0,
        rho_to = rho_new, edf_from = g0 + edf_floor,
        edf_to = gap(rho_new) + edf_floor)
    }
  }
  list(sp = sp, moved = if (length(moved) > 0L) do.call(rbind, moved) else NULL)
}

#' Bayesian backfitting by smoother at fixed lambda and phi
#'
#' Gibbs for Gaussian data with the identity link (each smoother and the
#' intercept from its full conditional), Metropolis-Hastings otherwise (each
#' block from the normal at its local-scoring step, reverse proposal at the
#' candidate). With joint_linear each sweep ends with the intercept and the
#' linear learners of all smoothers drawn as one block. Starts at `state`.
#' @param start The time (proc.time() "elapsed") from which max_runtime_mins
#'   counts; a sweep starts only while the limit is not reached.
#' @return list(a, beta (per smoother p_j x T), n_done, accept); n_done <
#'   n_burnin + n_samples if max_runtime_mins stopped it early; accept: the
#'   acceptance rate after burn-in per slot (intercept, smoothers, linear
#'   block), NA for a slot never drawn.
#' @noRd
backfitting_sampler <- function(m, sp, state, n_samples, n_burnin, phi,
                                joint_linear, verbose = FALSE,
                                max_runtime_mins = NULL,
                                start = proc.time()[["elapsed"]]) {
  force(start)
  blocks <- sampler_blocks(m, sp, joint_linear)
  cur    <- sampler_start(m, state, blocks)
  a_draw <- numeric(n_samples)
  b_draw <- lapply(m$smoothers, function(sm)
    matrix(0, length(sm$beta), n_samples))
  # draws tried and accepted after burn-in, per slot
  tally   <- matrix(0L, 2L, length(m$smoothers) + 2L,
                    dimnames = list(c("try", "acc"), NULL))
  n_total <- n_burnin + n_samples
  pb      <- if (verbose)
    cli::cli_progress_bar("Bayesian backfitting", total = n_total,
                          clear = FALSE)
  n_done <- 0L
  for (t in seq_len(n_total)) {
    if (!is.null(max_runtime_mins) &&
        (proc.time()[["elapsed"]] - start) / 60 >= max_runtime_mins) break
    sw  <- backfitting_sweep(m, blocks, cur, phi, sprintf("sweep %d", t))
    cur <- sw$cur
    if (t > n_burnin) {
      drawn <- !is.na(sw$accepted)
      tally["try", drawn] <- tally["try", drawn] + 1L
      tally["acc", drawn] <- tally["acc", drawn] + sw$accepted[drawn]
      s <- t - n_burnin
      a_draw[s] <- cur$a
      for (j in seq_along(m$smoothers)) b_draw[[j]][, s] <- cur$beta[[j]]
    }
    n_done <- t
    if (verbose) cli::cli_progress_update(id = pb)
  }
  if (verbose) cli::cli_progress_done(id = pb)
  n_try <- tally["try", ]
  list(a = a_draw, beta = b_draw, n_done = n_done,
       accept = ifelse(n_try > 0, tally["acc", ] / n_try, NA_real_))
}

#' The blocks of the sampler
#'
#' model_blocks(); for Gaussian data with the identity link each block also
#' holds R, the Cholesky factor of V^T W V + diag(cDlam), which is the same
#' in every sweep (block_system()).
#' @noRd
sampler_blocks <- function(m, sp, joint_linear) {
  blocks <- model_blocks(m, sp, joint_linear)
  if (!is_gaussian_identity(m$family)) return(blocks)
  with_R <- function(bl) {
    if (is.null(bl)) return(NULL)
    bl$R <- block_system(bl, m$w, paste("sampler,", bl$label))$R
    bl
  }
  list(intercept = with_R(blocks$intercept),
       smoothers = lapply(blocks$smoothers, with_R),
       linear    = with_R(blocks$linear))
}

#' The state of the sampler at its start
#'
#' Smoothers without a block (all learners at lambda = Inf) are set to zero.
#' @param state The posterior mode, list(a, beta).
#' @return list(a, beta, eta, ll): eta the linear predictor (offset
#'   included) and ll its scaled_loglik() (NA for Gaussian data with the
#'   identity link, where the Gibbs draws do not use it).
#' @noRd
sampler_start <- function(m, state, blocks) {
  beta <- state$beta
  for (j in which(vapply(blocks$smoothers, is.null, logical(1))))
    beta[[j]][] <- 0
  eta <- m$offset + state$a +
    Reduce(`+`, Map(function(sm, b) as.vector(sm$V %*% b), m$smoothers,
                    beta), 0)
  l2  <- is_gaussian_identity(m$family)
  ll  <- if (l2) NA_real_ else scaled_loglik(m, eta)
  if (!l2 && !is.finite(ll))
    stop(backfitting_failed_condition("non-finite log-likelihood at the mode",
                                      "nonfinite"))
  list(a = state$a, beta = beta, eta = eta, ll = ll)
}

#' One sweep of the sampler: the intercept, each smoother, the linear block
#'
#' @param cur The sampler state (sampler_start()).
#' @param where Names the sweep in errors ("sweep t").
#' @return list(cur, accepted): accepted holds one logical per slot
#'   (intercept, smoothers, linear block), NA for a slot not drawn.
#' @noRd
backfitting_sweep <- function(m, blocks, cur, phi, where) {
  n_smoothers <- length(m$smoothers)
  accepted    <- rep(NA, n_smoothers + 2L)
  bl <- blocks$intercept
  if (!is.null(bl)) {
    r <- draw_block(m, bl, cur$a, cur$eta - cur$a, cur, phi, where)
    cur$a        <- r$b
    cur$eta      <- r$eta
    cur$ll       <- r$ll
    accepted[1L] <- r$accepted
  }
  for (j in seq_len(n_smoothers)) {
    bl <- blocks$smoothers[[j]]
    if (is.null(bl)) next                   # all learners at lambda = Inf
    b      <- cur$beta[[j]]
    eta_mk <- cur$eta - as.vector(m$smoothers[[j]]$V %*% b)
    r <- draw_block(m, bl, b[bl$keep], eta_mk, cur, phi, where)
    cur$beta[[j]] <- replace(numeric(length(b)), bl$keep, r$b)
    cur$eta       <- r$eta
    cur$ll        <- r$ll
    accepted[j + 1L] <- r$accepted
  }
  bl <- blocks$linear
  if (!is.null(bl)) {
    old    <- get_linear_block(m, cur$a, cur$beta, bl)
    eta_mk <- cur$eta - as.vector(bl$V %*% old)
    r      <- draw_block(m, bl, old, eta_mk, cur, phi, where)
    cur    <- set_linear_block(m, cur, bl, r$b)
    cur$eta <- r$eta
    cur$ll  <- r$ll
    accepted[n_smoothers + 2L] <- r$accepted
  }
  list(cur = cur, accepted = accepted)
}

#' One draw of block bl given the rest of the linear predictor
#'
#' gibbs_draw() for Gaussian data with the identity link, else mh_draw().
#' @param bl A block (sampler_blocks()).
#' @param old The block's current coefficients.
#' @param eta_mk The linear predictor without the block.
#' @param cur The sampler state (sampler_start()).
#' @param where Names the sweep in errors; the block's label is added.
#' @return list(b, eta, ll, accepted): the block's new coefficients, the
#'   linear predictor and scaled log-likelihood with them, and whether the
#'   draw was accepted (always TRUE for Gibbs).
#' @noRd
draw_block <- function(m, bl, old, eta_mk, cur, phi, where) {
  where <- paste0(where, ", ", bl$label)
  if (is_gaussian_identity(m$family)) return(gibbs_draw(m, bl, eta_mk, phi))
  mh_draw(m, bl, old, eta_mk, cur$eta, cur$ll, phi, where)
}

#' Gibbs draw of block bl from its full conditional, Gaussian identity
#'
#' N(R^{-1} R^{-T} V^T W (y - eta_mk), phi (R^T R)^{-1}) with R from
#' sampler_blocks().
#' @noRd
gibbs_draw <- function(m, bl, eta_mk, phi) {
  R   <- bl$R
  rhs <- as.vector(crossprod(bl$V, m$w * (m$y - eta_mk)))
  mn  <- backsolve(R, forwardsolve(t(R), rhs))
  new <- mn + sqrt(phi) * backsolve(R, stats::rnorm(ncol(bl$V)))
  list(b = new, eta = eta_mk + as.vector(bl$V %*% new), ll = NA_real_,
       accepted = TRUE)
}

#' Metropolis-Hastings draw of block bl
#'
#' The candidate is drawn from the local-scoring proposal at the current
#' linear predictor eta; the reverse proposal is taken at the candidate.
#' @param eta,ll The current linear predictor and its scaled_loglik().
#' @return As draw_block(); the current values if the candidate is rejected
#'   or outside the range of the link (posterior density 0).
#' @noRd
mh_draw <- function(m, bl, old, eta_mk, eta, ll, phi, where) {
  stay  <- list(b = old, eta = eta, ll = ll, accepted = FALSE)
  q     <- local_scoring_proposal(m, bl, eta_mk, eta, phi, where)
  cand  <- q$mean + backsolve(q$R, stats::rnorm(ncol(bl$V)))
  eta_c <- eta_mk + as.vector(bl$V %*% cand)
  if (!in_link_range(m$family, eta_c)) return(stay)
  q_c  <- local_scoring_proposal(m, bl, eta_mk, eta_c, phi,
                                 paste(where, "(reverse)"))
  ll_c <- scaled_loglik(m, eta_c)
  # log posterior ratio + log q(old | cand) - log q(cand | old)
  log_ratio <- (ll_c - ll - 0.5 * sum(bl$cDlam * cand^2) +
                  0.5 * sum(bl$cDlam * old^2)) / phi +
    proposal_log_density(old, q_c) - proposal_log_density(cand, q)
  if (is.nan(log_ratio))
    stop(backfitting_failed_condition(
      paste0(where, ": log acceptance ratio is NaN"), "nan"))
  if (log(stats::runif(1)) < log_ratio)
    return(list(b = cand, eta = eta_c, ll = ll_c, accepted = TRUE))
  stay
}

#' Local-scoring proposal of block bl at the linear predictor eta
#'
#' The normal N(mean, (R^T R)^{-1}) of one IWLS step for the block at eta,
#' with R the Cholesky factor of (V^T W_z V + diag(cDlam)) / phi.
#' @param eta_mk The linear predictor without the block.
#' @return list(R, mean, logdet), logdet = log det(R^T R).
#' @noRd
local_scoring_proposal <- function(m, bl, eta_mk, eta, phi, where) {
  V       <- bl$V
  wq      <- iwls_working(m, eta, where)
  ws      <- wq$wz / phi
  z_blk   <- eta - eta_mk + wq$zr             # working response of the block
  P       <- crossprod(V, ws * V)             # Fisher information of b
  diag(P) <- diag(P) + bl$cDlam / phi         # plus the prior precision
  R       <- checked_chol(P, paste(where, "proposal"))
  rhs     <- as.vector(crossprod(V, ws * z_blk))
  list(R = R, mean = backsolve(R, forwardsolve(t(R), rhs)),
       logdet = 2 * sum(log(diag(R))))
}

#' Log density of the proposal q at b, up to a constant
#' @noRd
proposal_log_density <- function(b, q) {
  0.5 * q$logdet - 0.5 * sum((q$R %*% (b - q$mean))^2)
}

#' Scaled log-likelihood at eta: -deviance / 2, which is the log-likelihood
#' times phi up to a constant
#' @noRd
scaled_loglik <- function(m, eta) {
  fam <- m$family
  -0.5 * sum(fam$dev.resids(m$y, fam$linkinv(eta), m$w))
}

#' Smoother j in the working model at the mode, other smoothers fixed
#'
#' VtWzV = V^T W_z V and c = V^T W_z (z - eta_-j); `where` names the caller
#' in errors.
#' @noRd
smoother_working_model <- function(m, j, state,
                                   where = m$smoothers[[j]]$label) {
  sm    <- m$smoothers[[j]]
  wq    <- iwls_working(m, state$eta, where)
  VtWzV <- crossprod(sm$V, wq$wz * sm$V)
  # zr is the full residual, so add back the smoother's own fit
  list(VtWzV = VtWzV, c = as.vector(crossprod(sm$V, wq$wz * wq$zr) +
                                      VtWzV %*% state$beta[[j]]))
}

#' Shrinkage Fmat = (VtWzV + S)^{-1} VtWzV of smoother sm and each learner's
#' edf
#' @noRd
smoother_shrinkage <- function(sm, tw, st, where = sm$label) {
  pr    <- smoother_prior(sm, st)
  keep  <- pr$keep
  VtWzV <- tw$VtWzV[keep, keep, drop = FALSE]
  Fmat  <- VtWzV
  if (any(keep)) {
    A       <- VtWzV
    diag(A) <- diag(A) + pr$cDlam
    Fmat    <- chol_solve(A, VtWzV, where)
  }
  lid   <- sm$lid[keep]
  diagF <- diag(Fmat)
  list(Fmat = Fmat, edf = vapply(seq_len(nrow(sm$learners)), function(l)
    sum(diagF[lid == l]), numeric(1)))
}

#' Condition signalled when Bayesian backfitting (Stage 2) fails
#' @noRd
backfitting_failed_condition <- function(msg, status) {
  numerical_condition("bb_failed", paste("bayesian_backfitting: Bayesian",
                                         "backfitting failed:", msg),
                      status = status)
}


# ---- Driver: both stages -----------------------------------------------------

#' Boost to backfit: Stage 1 (MAP approximation), then Stage 2 (the posterior
#' mode and the posterior draws by Bayesian backfitting)
#'
#' @param mod A fitted "gbgam" object.
#' @param linear "free": linear learners flat (lambda = 0) and drawn as one
#'   block with the intercept; "kept": at their Stage-1 lambda.
#' @param undersmooth Logical; undersmoothing of low-edf learners
#'   (undersmooth_to_edf_floor()) up to edf `edf_floor`.
#' @param phi The dispersion, or NULL to estimate it (dispersion_at_mode(),
#'   after undersmoothing).
#' @param max_runtime_mins Time limit in minutes, or NULL, counted from
#'   `start` (proc.time() "elapsed"): both stages count against it, and the
#'   sampler stops when it is reached.
#' @return list(beta_samples (per learner, p_k x T), accept (per learner),
#'   info, n_done).
#' @noRd
boost_to_backfit <- function(mod, n_samples, n_burnin, linear = "free",
                             undersmooth = FALSE, edf_floor = 1,
                             phi = NULL, verbose = FALSE,
                             max_runtime_mins = NULL,
                             start = proc.time()[["elapsed"]]) {
  force(start)
  m  <- model_as_smoothers(mod)
  s1 <- reclass_numerical_failure(map_smoothing_parameters(m),
                                  map_failed_condition)
  s2 <- reclass_numerical_failure(
    sample_posterior(m, s1, linear, undersmooth, edf_floor, phi, n_samples,
                     n_burnin, verbose, max_runtime_mins, start),
    backfitting_failed_condition)
  sp <- s2$sp
  dr <- s2$draws
  # draws per learner, on the learner's own scale
  K   <- length(mod$learners)
  n_s <- max(0L, dr$n_done - n_burnin)
  beta_samples <- lapply(mod$learners, function(lr)
    matrix(0, length(lr$beta), n_s))
  if (m$has_int) beta_samples[[m$int_k]][1L, ] <- dr$a[seq_len(n_s)]
  accept <- rep(NA_real_, K)
  if (m$has_int) accept[m$int_k] <- dr$accept[1L]
  for (j in seq_along(m$smoothers)) {
    sm <- m$smoothers[[j]]
    for (l in seq_along(sm$k)) {
      i <- sm$lid == l
      beta_samples[[sm$k[l]]] <-
        dr$beta[[j]][i, seq_len(n_s), drop = FALSE] / sm$scale[l]
      if (sp[[j]]$status[l] != "inf")             # NA if held at zero
        accept[sm$k[l]] <- dr$accept[j + 1L]
    }
  }
  lam <- do.call(rbind, lapply(seq_along(m$smoothers), function(j)
    data.frame(term = m$smoothers[[j]]$label,
               learner = m$smoothers[[j]]$learners$label,
               status = sp[[j]]$status, rho = sp[[j]]$rho,
               stringsAsFactors = FALSE)))
  list(beta_samples = beta_samples, accept = accept, n_done = dr$n_done,
       info = list(stage1 = s1$status, stage1_sweeps = s1$sweeps,
                   discrepancy = s1$discrepancy, phi = s2$phi,
                   linear = linear,
                   linear_accept = dr$accept[length(m$smoothers) + 2L],
                   aliased = s2$aliased, undersmoothed = s2$undersmoothed,
                   lambda = lam))
}
