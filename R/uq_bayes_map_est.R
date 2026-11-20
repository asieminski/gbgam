# Boost to backfit, part 1: the posterior of a boosted additive model, used
# by bayesian_backfitting(). This file holds the vocabulary, the helpers of
# both stages and Stage 1, the MAP approximation
# (map_smoothing_parameters()); Stage 2 (sample_posterior()) and the driver
# boost_to_backfit() are in uq_bayes_mcmc.R. All functions are internal.
#
# Vocabulary. Smoother j is the set of learners of one model term (all
# learners of a t2() term, or one parametric column). Learner
# k has a W-orthonormal Demmler-Reinsch basis V_k (V_k^T W V_k = I), penalty
# eigenvalues D_k and smoothing parameter lambda_k = exp(rho_k); in the code
# l indexes the learners of a smoother and sm$k maps them to mod$learners. A
# learner of kind "null" is unpenalised (a linear or parametric column,
# scaled to unit W-norm, with D_k = 1); the other kind is "pen". The
# penalty of smoother j is diagonal, S_j = diag(cDlam_j), with cDlam_j the
# concatenated lambda_k D_k of its learners (the code stores only the vector
# cDlam), so beta_k ~ N(0, phi diag(lambda_k D_k)^{-1}); lambda_k = Inf
# removes learner k and lambda_k = 0 gives it a flat prior, as for the
# intercept. beta_b are the boosted coefficients. A block is a set of
# columns fitted or drawn together (model_blocks(): the intercept, the kept
# learners of one smoother, or the linear block).
#
# Data structures.
#   sp: the smoothing parameters, per smoother list(status, rho): status per
#     learner, "finite", "zero" (lambda = 0) or "inf" (lambda = Inf); rho =
#     log(lambda), NA unless finite.
#   state: list(a, beta (per smoother), eta (linear predictor incl. offset)).
#
# Stage 1, MAP approximation: the lambda whose penalised fit
# (posterior mode) is close to the boosted smoothers (step 2).
#   1. Penalised fit at lambda by fitting sweeps over the smoothers:
#      (V_j^T W V_j + S_j) beta_j = V_j^T W (y - eta_-j), a singular system
#      perturbed and a sweep halved if the penalised deviance increases.
#   2. Discrepancy: sum_j (beta_j - beta_bj)^T V_j^T W V_j (beta_j - beta_bj).
#   3. Per learner l of smoother j, the one-step coefficients
#      (V_j^T W V_j + S_j,-l + lambda_l D_l)^{-1} V_j^T W (y - eta_-j) as a
#      function of lambda_l (generalised eigen-decomposition); the
#      log(lambda_l) minimising their discrepancy is kept if the discrepancy
#      of the new penalised fit does not increase (else the step is halved).
#   4. Start from the penalised fit at the lambda with tr(F_k) = edf~_k
#      (unpenalised learners: at boosting's shrinkage of them,
#      starting_lambdas()), then lambda sweeps over all learners until
#      the fit changes by less than tol.
#   5. Other families: W and y - eta become the IWLS working weights
#      W_z = W h'(eta)^2 / V(mu) and residuals z - eta = (y - mu) / h'(eta).
#
# Stage 2, Bayesian backfitting, at the Stage-1 lambda:
#   1. linear = "free": unpenalised learners get lambda = 0 and form one
#      block with the intercept; exactly aliased flat columns are dropped.
#   2. Posterior mode (fitting sweeps); optionally undersmoothing (learners
#      below an edf floor are lifted to it) and the mode again; then the
#      dispersion phi by REML at the fixed lambda (Wood, 2011), at that mode.
#   3. Each sweep draws the intercept, each smoother and the linear block from
#      N((V^T W_z V + S)^{-1} V^T W_z (z - eta_-j), phi (V^T W_z V + S)^{-1}):
#      exactly (Gibbs) for Gaussian identity, otherwise as the proposal of a
#      Metropolis-Hastings step.
#   4. The draws are mapped to the mgcv coefficients by
#      attach_posterior_draws().
#
# Cost (P coefficients in all, p_j in smoother j, K learners). No step forms
# the full X^T W X. For Gaussian data with the identity link each block's
# V^T W V and Cholesky factor (block_system()) are computed once per call of
# penalised_fit() (one penalised fit) and once in sampler_blocks(),
# \mathcal{O}(n sum_j p_j^2 + sum_j p_j^3), so a sampling sweep costs
# \mathcal{O}(nP + sum_j p_j^2); otherwise every sweep costs
# \mathcal{O}(n sum_j p_j^2 + sum_j p_j^3) (Metropolis-Hastings builds a
# forward and a reverse proposal per block). Stage 1 refits the penalised
# model by backfitting for every learner whose proposal it tries (moved or
# not), so each of its lambda sweeps runs up to K (1 + n_back) penalised
# fits, after the one at the initial lambda.
#
# Numerical failures signal the condition class map_approx_failed (Stage 1)
# or bb_failed (Stage 2), each with a $status; the internal classes
# b2b_nonfinite and b2b_not_pd are re-signalled as one of these
# (reclass_numerical_failure()). A singular system of a Stage-1 one-step
# proposal is not a failure: it skips that learner for that lambda sweep.


# ---- Both stages: model, prior, fitting sweeps, working quantities ----------

#' The fitted `gbgam` model as smoothers
#'
#' Groups the learners of `mod` by model term (`label` of smooth learners,
#' `var_label` of parametric ones; the smooth terms are those of
#' `mod$smooth`) with concatenated bases `V`, `VtWV` = V^T W V, penalty
#' eigenvalues `d`, learner index `lid`, boosted coefficients and the
#' learners' accumulated edf~_k (`learners$edf`).
#' Unpenalised columns are scaled to unit W-norm (penalty eigenvalue 1); a
#' column with zero W-norm (e.g. a factor level only on zero-weight rows,
#' which boosting fits as zero) is kept unscaled and flagged `zero_norm`,
#' and stays out of the posterior (lambda = Inf).
#' @param mod A fitted "gbgam" object.
#' @return A list with y, w, offset, family, has_int (FALSE without an
#'   intercept and for an intercept held at zero: excluded, not refitted
#'   after each iteration and not initialised), int_k (the intercept
#'   learner), a, eta and smoothers.
#' @noRd
model_as_smoothers <- function(mod) {
  L      <- mod$learners
  w      <- as.numeric(mod$prior.weights)
  y      <- as.numeric(mod$y)
  offset <- as.numeric(mod$offset)
  is_int <- vapply(L, function(lr) identical(lr$var_label, "(Intercept)"),
                   logical(1))
  kind   <- vapply(L, function(lr)
    if (inherits(lr, "demmler_reinsch")) "pen" else "null", character(1))
  group  <- vapply(L, function(lr) lr$label %||% lr$var_label, character(1))
  # learners the user excluded (and not refitted after each iteration) stay
  # out of the posterior even when linear learners are freed
  excl   <- setdiff(mod$excluded_learner_indices %||% integer(),
                    mod$learners_to_fit_after_each_iter %||% integer())
  labels <- unique(group[!is_int])
  smoothers <- lapply(seq_along(labels), function(j) {
    idx    <- which(group == labels[j] & !is_int)
    pieces <- lapply(idx, function(k) {
      lr <- L[[k]]
      if (kind[k] == "pen") {
        p <- list(V = lr$design, d = lr$D, beta = lr$beta, rho = lr$rho,
                  edf = lr$edf %||% 0,
                  at_zero = isTRUE(lr$lambda_at_zero), scale = 1,
                  zero_norm = sum(w * lr$design^2) == 0)
      } else {
        x         <- as.vector(lr$design)
        nx        <- sqrt(sum(w * x^2))
        zero_norm <- !(nx > 0)
        if (zero_norm) nx <- 1                     # kept unscaled, see above
        p <- list(V = matrix(x / nx), d = 1, beta = nx * lr$beta,
                  rho = NA_real_, edf = lr$edf %||% 0, at_zero = FALSE,
                  scale = nx, zero_norm = zero_norm)
      }
      p$never <- all(p$beta == 0)
      p
    })
    pick <- function(name, value) vapply(pieces, `[[`, value, name)
    V    <- do.call(cbind, lapply(pieces, `[[`, "V"))
    beta <- unlist(lapply(pieces, `[[`, "beta"))
    list(label = labels[j], V = V, VtWV = crossprod(V, w * V),
         d = unlist(lapply(pieces, `[[`, "d")),
         lid = rep(seq_along(idx), vapply(pieces, function(p) ncol(p$V), 1L)),
         k = idx, scale = pick("scale", 1), beta = beta,
         f = as.vector(V %*% beta),
         learners = data.frame(
           label     = vapply(L[idx], `[[`, "", "var_label"),
           kind      = kind[idx],
           never     = pick("never", TRUE),
           zero_norm = pick("zero_norm", TRUE),
           excluded  = idx %in% excl,
           rho       = pick("rho", 0),
           edf       = pick("edf", 0),
           at_zero   = pick("at_zero", TRUE),
           stringsAsFactors = FALSE))
  })
  names(smoothers) <- labels
  # an intercept that is excluded, not refitted and not initialised is held
  # at zero by boosting, and so in the posterior
  int_k   <- which(is_int)
  has_int <- length(int_k) > 0L &&
    !(int_k %in% excl && identical(mod$initialise_intercept, FALSE))
  a <- if (has_int) L[[int_k]]$beta else 0
  list(y = y, w = w, offset = offset, family = mod$family,
       has_int = has_int, int_k = int_k, a = a,
       eta = offset + a + Reduce(`+`, lapply(smoothers, `[[`, "f"), 0),
       smoothers = smoothers)
}

#' Columns of smoother sm kept under its smoothing parameters st and their
#' penalty diagonal cDlam
#'
#' cDlam is the prior precision times phi: exp(rho) d for "finite", 0 for
#' "zero"; "inf" drops columns. cDlam_all is cDlam over all columns (0 for
#' "inf").
#' @noRd
smoother_prior <- function(sm, st) {
  s_col <- st$status[sm$lid]
  keep  <- s_col != "inf"
  cDlam <- ifelse(s_col == "finite", exp(st$rho[sm$lid]) * sm$d, 0)
  list(keep = keep, cDlam = cDlam[keep], cDlam_all = cDlam)
}

#' The blocks of the fitting sweeps and the sampler
#'
#' A block is list(V, cDlam, label): columns V, drawn or fitted together,
#' their penalty diagonal cDlam (prior precision times phi) and a label for
#' messages.
#' @param m The model as smoothers (model_as_smoothers()).
#' @param sp The smoothing parameters.
#' @param joint_linear Whether the linear block is used (linear = "free").
#' @return list(intercept, smoothers, linear): the intercept (a column of
#'   ones with cDlam = 0; NULL without an intercept); per smoother, its kept
#'   columns with `keep` (NULL when all its learners are at lambda = Inf);
#'   the linear block (linear_block(); NULL unless joint_linear).
#' @noRd
model_blocks <- function(m, sp, joint_linear) {
  smoothers <- Map(function(sm, st) {
    pr <- smoother_prior(sm, st)
    if (any(pr$keep))
      list(V = sm$V[, pr$keep, drop = FALSE], cDlam = pr$cDlam,
           label = sm$label, keep = pr$keep)
  }, m$smoothers, sp)
  list(intercept = if (m$has_int)
         list(V = matrix(1, length(m$y)), cDlam = 0, label = "intercept"),
       smoothers = smoothers,
       linear    = if (joint_linear) linear_block(m, sp))
}

#' The penalised normal equations of block bl
#'
#' @param bl A block (model_blocks()).
#' @param wz Working weights W_z.
#' @param where Names the caller in errors.
#' @param factorise checked_chol(), or perturbed_chol() for a fitting step.
#' @return list(VtWzV, R, ridge): V^T W_z V, the Cholesky factor R of
#'   V^T W_z V + diag(cDlam) + diag(ridge), and ridge (0 unless perturbed).
#' @noRd
block_system <- function(bl, wz, where, factorise = checked_chol) {
  VtWzV   <- crossprod(bl$V, wz * bl$V)
  A       <- VtWzV
  diag(A) <- diag(A) + bl$cDlam
  f       <- factorise(A, where)
  if (is.list(f)) {
    list(VtWzV = VtWzV, R = f$R, ridge = f$ridge)
  } else {
    list(VtWzV = VtWzV, R = f, ridge = 0)
  }
}

#' block_system() of every smoother block and of the linear block, for the
#' fitting sweeps: a singular block is perturbed (perturbed_chol())
#' @return list(smoothers (NULL where the block is NULL), linear).
#' @noRd
block_systems <- function(blocks, wz, where) {
  sys <- function(bl) if (!is.null(bl))
    block_system(bl, wz, paste(where, bl$label), perturbed_chol)
  list(smoothers = lapply(blocks$smoothers, sys), linear = sys(blocks$linear))
}

#' The intercept and the linear learners of all smoothers as one block
#'
#' Columns V with penalty diagonal cDlam: the intercept and every
#' unpenalised learner not at lambda = Inf, null-space columns of smooths
#' included; Stage 2 with linear = "free". idx holds each smoother's
#' columns in the block.
#' @return A block (model_blocks()) with idx, or NULL when it has no column.
#' @noRd
linear_block <- function(m, sp) {
  idx <- lapply(seq_along(m$smoothers), function(j) {
    sm <- m$smoothers[[j]]
    which(sm$learners$kind[sm$lid] == "null" &
            sp[[j]]$status[sm$lid] != "inf")
  })
  V <- do.call(cbind, c(list(matrix(1, length(m$y), as.integer(m$has_int))),
                        Map(function(sm, i) sm$V[, i, drop = FALSE],
                            m$smoothers, idx)))
  if (ncol(V) == 0L) return(NULL)
  cDlam <- c(rep(0, as.integer(m$has_int)), unlist(Map(function(sm, st, i)
    smoother_prior(sm, st)$cDlam_all[i], m$smoothers, sp, idx)))
  list(V = V, cDlam = cDlam, label = "linear block", idx = idx)
}

#' Coefficients of the linear block lb from the intercept a and beta
#' @noRd
get_linear_block <- function(m, a, beta, lb) {
  c(if (m$has_int) a, unlist(Map(function(b, i) b[i], beta, lb$idx)))
}

#' Write the linear-block coefficients v back into state
#' @noRd
set_linear_block <- function(m, state, lb, v) {
  pos <- 0L
  if (m$has_int) {
    state$a <- v[1L]
    pos     <- 1L
  }
  for (j in seq_along(state$beta)) {
    i <- lb$idx[[j]]
    state$beta[[j]][i] <- v[pos + seq_along(i)]
    pos <- pos + length(i)
  }
  state
}

#' W-norm of V e for the coefficients e of smoother sm
#' @noRd
w_norm <- function(sm, e) sqrt(sum(e * as.vector(sm$VtWV %*% e)))

#' Penalised fit at smoothing parameters sp by fitting sweeps
#'
#' Per sweep the intercept, then each smoother jointly,
#' (V_j^T W_z V_j + S_j) b = V_j^T W_z (z - eta_-j), then (joint_linear,
#' Stage 2 only) the intercept and the unpenalised learners of all smoothers
#' (parametric columns and smooths' null-space columns, which are also in
#' their smoother's block) jointly, the joint update of modified backfitting
#' (Buja et al., 1989). Any family, through iwls_working(); used by both
#' stages.
#'
#' A sweep solves one working model, so as a whole it is a descent direction
#' of the penalised deviance (a block step alone need not be, because W_z
#' stays at the start of the sweep). Except for Gaussian data with the
#' identity link, where every step is an exact minimiser, the whole sweep is
#' therefore halved until it stays in the range of the link and the penalised
#' deviance does not increase (keep_sweep()). Convergence is judged on the
#' whole sweep, so that a halved sweep does not look converged.
#' @return list(state, sweeps, ok); ok = FALSE if max_sweeps do not reach a
#'   change below tol, the working quantities fail or no part of a sweep is
#'   kept.
#' @noRd
penalised_fit <- function(m, sp, state, tol, max_sweeps = 1000L,
                          where = "fitting sweeps", joint_linear = FALSE) {
  blocks <- model_blocks(m, sp, joint_linear)
  l2     <- is_gaussian_identity(m$family)
  for (i in seq_len(max_sweeps)) {
    wq <- tryCatch(iwls_working(m, state$eta, where),
                   b2b_nonfinite = function(e) NULL)
    if (is.null(wq)) return(list(state = state, sweeps = i, ok = FALSE))
    # For Gaussian data with the identity link W_z = W in every sweep, so the
    # systems are computed in the first sweep only
    if (i == 1L || !l2) sys <- block_systems(blocks, wq$wz, where)
    sw <- fitting_sweep(m, blocks, sys, wq, state)
    kept <- if (l2) sw$state else keep_sweep(m, sp, blocks, state, sw$state)
    if (is.null(kept)) return(list(state = state, sweeps = i, ok = FALSE))
    state <- kept
    if (sw$change < tol) return(list(state = state, sweeps = i, ok = TRUE))
  }
  list(state = state, sweeps = max_sweeps, ok = FALSE)
}

#' One fitting sweep of penalised_fit()
#'
#' A perturbed system solves for new = old + step, so its right-hand side
#' gains ridge * old, which keeps the solution of the unperturbed equations
#' as the fixed point.
#' @param sys The blocks' systems (block_systems()).
#' @param wq The working quantities at state (iwls_working()).
#' @return list(state, change): the state after the sweep and the largest
#'   change of a block in W-norm.
#' @noRd
fitting_sweep <- function(m, blocks, sys, wq, state) {
  zr     <- wq$zr
  change <- 0
  if (m$has_int) {
    da        <- sum(wq$wz * zr) / sum(wq$wz)
    state$a   <- state$a + da
    state$eta <- state$eta + da
    zr        <- zr - da
    change    <- abs(da) * sqrt(sum(m$w))            # W-norm, as below
  }
  for (j in seq_along(m$smoothers)) {
    sm  <- m$smoothers[[j]]
    bl  <- blocks$smoothers[[j]]
    old <- state$beta[[j]]
    if (is.null(bl) && all(old == 0)) next          # removed, already 0
    new <- numeric(length(old))
    if (!is.null(bl)) {
      sy  <- sys$smoothers[[j]]
      # V^T W_z (z - eta_-j): zr is the full residual, so add back the
      # smoother's own fit
      rhs <- as.vector(crossprod(bl$V, wq$wz * zr) +
                         sy$VtWzV %*% old[bl$keep]) + sy$ridge * old[bl$keep]
      new[bl$keep] <- backsolve(sy$R, forwardsolve(t(sy$R), rhs))
    }
    dV <- as.vector(sm$V %*% (new - old))
    state$eta <- state$eta + dV
    zr        <- zr - dV
    change    <- max(change, w_norm(sm, new - old))
    state$beta[[j]] <- new
  }
  lb <- blocks$linear
  if (!is.null(lb)) {
    sy  <- sys$linear
    old <- get_linear_block(m, state$a, state$beta, lb)
    rhs <- as.vector(crossprod(lb$V, wq$wz * zr) + sy$VtWzV %*% old) +
      sy$ridge * old
    new <- backsolve(sy$R, forwardsolve(t(sy$R), rhs))
    dV  <- as.vector(lb$V %*% (new - old))
    state$eta <- state$eta + dV
    change    <- max(change, sqrt(sum(m$w * dV^2)))
    state     <- set_linear_block(m, state, lb, new)
  }
  list(state = state, change = change)
}

#' Keep as much of a fitting sweep as the penalised deviance allows
#'
#' A sweep that removes columns (at lambda = Inf) is a projection and not a
#' descent step, so it is kept whole unless it leaves the link's range. A
#' sweep from a point of non-finite objective is kept whole as well.
#' @param old,new The states before and after the sweep.
#' @return The state to continue from, or NULL if none is acceptable.
#' @noRd
keep_sweep <- function(m, sp, blocks, old, new) {
  objective <- function(st) {
    if (!in_link_range(m$family, st$eta)) return(Inf)
    penalised_deviance(m, sp, st)
  }
  obj_new <- objective(new)
  if (removes_columns(blocks, old))
    return(if (is.finite(obj_new)) new)
  obj_old <- objective(old)
  between <- function(t) {
    if (t == 1) return(new)
    list(a    = old$a + t * (new$a - old$a),
         beta = Map(function(b0, b1) b0 + t * (b1 - b0), old$beta, new$beta),
         eta  = old$eta + t * (new$eta - old$eta))
  }
  srch <- halve_step(function(t) if (t == 1) obj_new else objective(between(t)),
                     obj_old,
                     loss_allowance(obj_old, deviance_roundoff(m$family,
                                                               old$eta, m$w)))
  if (srch$fraction > 0) between(srch$fraction)
}

#' TRUE if the state has a nonzero coefficient on a column outside the blocks
#' @noRd
removes_columns <- function(blocks, state) {
  any(vapply(seq_along(blocks$smoothers), function(j) {
    b  <- state$beta[[j]]
    bl <- blocks$smoothers[[j]]
    any((if (is.null(bl)) b else b[!bl$keep]) != 0)
  }, logical(1)))
}

#' Whether the family is Gaussian with the identity link (W_z = W)
#' @noRd
is_gaussian_identity <- function(family) {
  family$family == "gaussian" && family$link == "identity"
}

#' TRUE if eta and its means lie in the family's valid ranges (`valideta`,
#' `validmu`; a family without them accepts every value)
#' @noRd
in_link_range <- function(family, eta) {
  (is.null(family$valideta) || family$valideta(eta)) &&
    (is.null(family$validmu) || family$validmu(family$linkinv(eta)))
}

#' IWLS working weights W_z = W h'(eta)^2 / V(mu) and residuals z - eta
#'
#' W and y - eta for Gaussian identity.
#' @noRd
iwls_working <- function(m, eta, where) {
  if (is_gaussian_identity(m$family)) return(list(wz = m$w, zr = m$y - eta))
  fam <- m$family
  mu  <- fam$linkinv(eta)
  me  <- fam$mu.eta(eta)
  va  <- fam$variance(mu)
  if (!all(is.finite(c(mu, me, va))) || any(me == 0) || any(va <= 0))
    signal_numerical("b2b_nonfinite", paste0(where, ": non-finite mu, zero ",
                                             "h'(eta) or non-positive V(mu)"))
  list(wz = m$w * me^2 / va, zr = (m$y - mu) / me)
}

#' Cholesky factor of M; signals b2b_not_pd naming `where` if M is not
#' positive definite (an allocation failure is passed on unchanged)
#' @noRd
checked_chol <- function(M, where) {
  tryCatch(chol(M), error = function(e) {
    if (is_alloc_error(e)) stop(e)
    signal_numerical("b2b_not_pd", paste0(where, ": matrix is not positive ",
                                          "definite"))
  })
}

#' Cholesky factor of M, perturbed if M is not positive definite
#'
#' The perturbation of Wood et al. (2016, section 3.1.2) on the Jacobi-scaled
#' matrix: M + eps diag(diag(M)), with eps raised until the factor exists, so
#' that a column keeps its size relative to the others. For the steps of a
#' fit only; a matrix that defines a distribution uses checked_chol().
#' @return list(R, ridge): the factor and eps diag(M) (0 if not perturbed).
#'   Signals b2b_not_pd naming `where` if M is not finite or eps reaches its
#'   limit.
#' @noRd
perturbed_chol <- function(M, where) {
  try_chol <- function(A) tryCatch(chol(A), error = function(e) {
    if (is_alloc_error(e)) stop(e)
    NULL
  })
  if (all(is.finite(M))) {
    R <- try_chol(M)
    if (!is.null(R)) return(list(R = R, ridge = 0))
    d <- diag(M)
    for (eps in 10^seq(-14, -3)) {
      A       <- M
      diag(A) <- diag(A) + eps * d
      R       <- try_chol(A)
      if (!is.null(R)) return(list(R = R, ridge = eps * d))
    }
  }
  signal_numerical("b2b_not_pd", paste0(where, ": matrix is not positive ",
                                        "definite"))
}

#' TRUE if the error e is R's failure to allocate memory, which callers must
#' see as such (e.g. to rerun with more memory), not as a numerical failure
#' @noRd
is_alloc_error <- function(e) {
  grepl("cannot allocate|could not allocate|memory exhausted|vector memory",
        conditionMessage(e))
}

#' Solve M x = b by Cholesky, with an error naming `where`
#' @noRd
chol_solve <- function(M, b, where) {
  R <- checked_chol(M, where)
  backsolve(R, forwardsolve(t(R), b))
}

#' Error condition of class cls
#' @noRd
numerical_condition <- function(cls, msg, ...) {
  structure(class = c(cls, "error", "condition"),
            list(message = msg, call = NULL, ...))
}

#' Signal an error condition of class cls
#' @noRd
signal_numerical <- function(cls, msg) stop(numerical_condition(cls, msg))

#' Evaluate expr, re-signalling its numerical failures (b2b_nonfinite,
#' b2b_not_pd) as failed(msg, status), i.e. map_approx_failed or bb_failed
#' @noRd
reclass_numerical_failure <- function(expr, failed) {
  tryCatch(expr,
           b2b_nonfinite = function(e)
             stop(failed(conditionMessage(e), "nonfinite")),
           b2b_not_pd    = function(e)
             stop(failed(conditionMessage(e), "not_pd")))
}


# ---- Stage 1: MAP approximation (any family, IWLS) ---------------------------

#' Stage 1: smoothing parameters whose penalised fit is close to the boosted
#' smoothers
#'
#' Learners boosting never updated stay at lambda = Inf; the others start at
#' starting_lambdas(), with the penalised fit there
#' (no mode: map_approx_failed with status "no_mode"). Lambda sweeps: per
#' learner a one-step proposal, kept only if the discrepancy of the new
#' penalised fit does not increase (try_learner_proposal()); a learner whose
#' one-step system is singular gets no proposal in that sweep. Ends when a
#' sweep changes neither the intercept nor any smoother by tol: "converged", or
#' "stalled" if that sweep rejected a proposal or had none for a learner (both
#' used in the same way; no tolerance on the discrepancy left); after
#' max_sweeps lambda sweeps, signals map_approx_failed.
#' @param pad Margin of the log(lambda) range: shrinkage factors within
#'   1 / (1 + e^pad), about 1e-8, of 0 and 1 count as their limits.
#' @param n_back Halvings of a rejected step.
#' @param max_inner Fitting sweeps allowed per penalised fit (also used by
#'   Stage 2, through the returned max_inner).
#' @return list(sp, state, tol, discrepancy, sweeps, status, max_inner).
#' @noRd
map_smoothing_parameters <- function(m, eps = 1e-7, max_sweeps = 200L,
                                     step = 0.1, pad = 18.5, n_back = 4L,
                                     max_inner = 1000L) {
  tol   <- convergence_tol(m, eps)
  sp    <- starting_lambdas(m)
  state <- list(a = m$a, beta = lapply(m$smoothers, `[[`, "beta"), eta = m$eta)
  free  <- do.call(rbind, lapply(seq_along(m$smoothers), function(j) {
    lr <- m$smoothers[[j]]$learners
    if (any(!lr$never)) data.frame(j = j, l = which(!lr$never))
  }))
  if (is.null(free))
    return(list(sp = sp, state = state, tol = tol, status = "no_learners",
                sweeps = 0L, discrepancy = 0, max_inner = max_inner))
  # R0 (R0^T R0 = V^T W V) and target = R0 beta_b of the smoothers with a free
  # learner, fixed through Stage 1: one_step_proposal() measures the
  # discrepancy as ||R0 b - target||^2
  js <- unique(free$j)
  m$smoothers[js] <- lapply(m$smoothers[js], function(sm) {
    sm$R0     <- matrix_root(sm$VtWV)
    sm$target <- as.vector(sm$R0 %*% sm$beta)
    sm
  })
  md <- penalised_fit(m, sp, state, tol, max_inner, "initial lambda, mode")
  if (!md$ok)
    stop_map_failed("no penalised fit at the initial lambda", "no_mode")
  state         <- md$state
  disc_now      <- model_discrepancy(m, state$beta)
  converged     <- FALSE
  rejected_last <- character()
  sweeps        <- 0L
  for (t in seq_len(max_sweeps)) {
    start         <- state
    rejected_last <- character()
    for (i in seq_len(nrow(free))) {
      j     <- free$j[i]
      l     <- free$l[i]
      lab   <- m$smoothers[[j]]$learners$label[l]
      where <- sprintf("lambda sweep %d, %s", t, lab)
      # A singular system gives no proposal for this learner in this sweep
      p     <- tryCatch(reclass_search_failure(
        one_step_proposal(m, j, l, sp, state, step, pad, tol, where)),
        b2b_not_pd = function(e) NULL)
      if (is.null(p)) {
        rejected_last <- c(rejected_last, lab)
        next
      }
      if (w_norm(m$smoothers[[j]], p$coef - state$beta[[j]]) < tol) next
      tr <- try_learner_proposal(m, j, l, sp, state, p, p$range, disc_now,
                                 tol, n_back, max_inner, where)
      if (tr$moved) {
        state    <- tr$state
        sp       <- tr$sp
        disc_now <- tr$discrepancy
      } else {
        rejected_last <- c(rejected_last, lab)
      }
    }
    sweeps <- t
    change <- max(abs(state$a - start$a) * sqrt(sum(m$w)),
                  vapply(seq_along(m$smoothers), function(j)
                    w_norm(m$smoothers[[j]], state$beta[[j]] -
                                 start$beta[[j]]), numeric(1)))
    if (change < tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged)
    stop_map_failed(paste("no convergence in", max_sweeps, "lambda sweeps"),
                    "not_converged")
  list(sp = sp, state = state, tol = tol, discrepancy = disc_now,
       sweeps = sweeps,
       status = if (length(rejected_last) > 0L) "stalled" else "converged",
       max_inner = max_inner)
}

#' Convergence tolerance: eps times the W-norm of the boosted smoothers or of
#' the whole predictor, whichever is larger
#' @noRd
convergence_tol <- function(m, eps) {
  size <- max(sqrt(sum(vapply(m$smoothers, function(sm)
    sum(sm$beta * as.vector(sm$VtWV %*% sm$beta)), numeric(1)))),
    sqrt(sum(m$w * (m$eta - m$offset)^2)))
  if (size == 0) size <- sqrt(sum(m$w))
  eps * size
}

#' Starting smoothing parameters: never-updated and zero-norm learners at
#' lambda = Inf, penalised ones at the lambda whose trace equals their edf,
#' and unpenalised ones at boosting's shrinkage of them
#'
#' A penalised learner starts at the log(lambda) with tr(F_k) = sum_i 1 /
#' (1 + lambda d_i) = edf~_k, its accumulated edf in the boosted fit
#' (rho_matching_edf()). An unpenalised learner (one column x, D = 1)
#' updated alone from the boosted fit returns its boosted coefficient beta_b
#' at lambda = x^T W_z (z - eta) / beta_b, for any x^T W_z x; it starts
#' there, or at 0 if that is not positive.
#' @noRd
starting_lambdas <- function(m) {
  wq <- NULL
  lapply(m$smoothers, function(sm) {
    lr     <- sm$learners
    status <- ifelse(lr$never | lr$zero_norm, "inf",
                     ifelse(lr$kind == "null" | lr$at_zero, "zero", "finite"))
    rho    <- ifelse(status == "finite", lr$rho, NA_real_)
    for (l in which(lr$kind == "pen" & status != "inf")) {
      st        <- rho_matching_edf(sm, l)
      status[l] <- st$status
      rho[l]    <- st$rho
    }
    for (l in which(lr$kind == "null" & status == "zero")) {
      wq   <- wq %||% iwls_working(m, m$eta, "initial lambda")
      on_l <- sm$lid == l
      lam  <- sum(sm$V[, on_l] * wq$wz * wq$zr) / sm$beta[on_l]
      if (is.finite(lam) && lam > 0) {
        status[l] <- "finite"
        rho[l]    <- log(lam)
      }
    }
    list(status = status, rho = rho)
  })
}

#' Starting smoothing parameter of penalised learner l of smoother sm
#'
#' The log(lambda) with tr(F_k) = sum_i 1 / (1 + lambda d_i) = edf~_k
#' (solve_rho_trA(), from boosting's rho, or from 0 when boosting set lambda
#' to 0): status "zero" when edf~_k is at least the number of columns, and
#' boosting's own lambda when edf~_k is not positive (always so without EDF
#' tracking).
#' @return list(status, rho).
#' @noRd
rho_matching_edf <- function(sm, l) {
  lr  <- sm$learners
  d   <- sm$d[sm$lid == l]
  edf <- lr$edf[l]
  if (!(edf > 0))
    return(if (lr$at_zero[l]) list(status = "zero", rho = NA_real_)
           else list(status = "finite", rho = lr$rho[l]))
  if (edf >= length(d)) return(list(status = "zero", rho = NA_real_))
  rho0 <- if (lr$at_zero[l]) 0 else lr$rho[l]
  rho  <- tryCatch(
    solve_rho_trA(rho0, D = d, target_trA = edf, tol = 1e-8,
                  label = lr$label[l]),
    error = function(e) NA_real_)
  if (!is.finite(rho))
    stop_map_failed(sprintf(paste("initial lambda of %s: no log(lambda)",
                                  "with edf %.4g"), lr$label[l], edf),
                    "search_failed")
  list(status = "finite", rho = rho)
}

#' Try the one-step proposal p of learner l of smoother j
#'
#' The proposal, then up to n_back halvings of the step in log(lambda)
#' (from its current position within `range`); the first whose penalised
#' fit exists and does not increase the discrepancy is kept.
#' @return list(moved, and if moved: sp, state, discrepancy).
#' @noRd
try_learner_proposal <- function(m, j, l, sp, state, p, range, disc_now,
                                 tol, n_back, max_inner, where) {
  p0    <- log_lambda_position(sp[[j]]$status[l], sp[[j]]$rho[l], range)
  p1    <- log_lambda_position(p$status, p$rho, range)
  tries <- c(list(list(status = p$status, rho = p$rho)),
             lapply(seq_len(n_back), function(h)
               list(status = "finite", rho = p0 + (p1 - p0) / 2^h)))
  for (tr in tries) {
    sp_try <- set_lambda(sp, j, l, tr$status, tr$rho)
    md     <- penalised_fit(m, sp_try, state, tol, max_inner,
                            paste(where, "fit"))
    if (!md$ok) next
    disc_try <- model_discrepancy(m, md$state$beta)
    if (disc_try <= disc_now)
      return(list(moved = TRUE, sp = sp_try, state = md$state,
                  discrepancy = disc_try))
  }
  list(moved = FALSE)
}

#' Position of a learner's smoothing parameter on its log(lambda) range:
#' rho if finite, the range's ends for lambda = Inf ("hi") and 0 ("lo")
#' @noRd
log_lambda_position <- function(status, rho, range) {
  switch(status, finite = rho, inf = range[["hi"]], zero = range[["lo"]])
}

#' Set learner l of smoother j in sp to status (and rho, kept only if finite)
#' @noRd
set_lambda <- function(sp, j, l, status, rho) {
  sp[[j]]$status[l] <- status
  sp[[j]]$rho[l]    <- if (status == "finite") rho else NA_real_
  sp
}

#' One-step proposal for learner l of smoother j
#'
#' The one-step (backfitting / local-scoring) coefficients of smoother j as a
#' function of learner l's log(lambda), the other learners at their state:
#' (VtWzVpS + lambda D_l) b = rhs with VtWzVpS = V^T W_z V + S without
#' learner l, by the generalised eigen-decomposition of (VtWzVpS, D_l) with
#' eigenvalues delta and shrinkage factors s = 1 / (1 + lambda delta). The
#' smoothers of m carry R0 and target from map_smoothing_parameters().
#' @noRd
one_step_proposal <- function(m, j, l, sp, state, step, pad, tol, where,
                              eig_tol = 1e-10) {
  sm      <- m$smoothers[[j]]
  pr      <- smoother_prior(sm, sp[[j]])
  on_l    <- sm$lid == l
  keep    <- on_l | pr$keep
  cDlam   <- replace(pr$cDlam_all, on_l, 0)
  wq      <- iwls_working(m, state$eta, where)
  Vk      <- sm$V[, keep, drop = FALSE]
  VtWzV   <- crossprod(Vk, wq$wz * Vk)
  VtWzVpS <- VtWzV
  diag(VtWzVpS) <- diag(VtWzVpS) + cDlam[keep]
  # V^T W_z (z - eta_-j): zr is the full residual, so add back the
  # smoother's fit
  rhs     <- as.vector(crossprod(Vk, wq$wz * wq$zr) +
                         VtWzV %*% state$beta[[j]][keep])
  R       <- checked_chol(VtWzVpS, paste(where, "(V^T W_z V + S)"))
  Ri      <- backsolve(R, diag(sum(keep)))
  Dl      <- ifelse(on_l[keep], sm$d[keep], 0)
  e       <- eigen(crossprod(Ri, Dl * Ri), symmetric = TRUE)
  delta   <- e$values
  delta[delta < eig_tol * max(delta)] <- 0   # directions learner l leaves free
  U       <- Ri %*% e$vectors
  cc      <- as.vector(crossprod(U, rhs))
  coef_full <- function(s) {
    b       <- numeric(length(sm$beta))
    b[keep] <- as.vector(U %*% (cc * s))
    b
  }
  fit_at <- function(s) as.vector(sm$R0 %*% coef_full(s))
  srch <- search_log_lambda(function(r) fit_at(1 / (1 + exp(r) * delta)),
                            fit_at(as.numeric(delta == 0)),
                            fit_at(rep(1, length(delta))),
                            sm$target, delta, step, pad, tol, where)
  s <- switch(srch$status, finite = 1 / (1 + exp(srch$rho) * delta),
              inf = as.numeric(delta == 0), zero = rep(1, length(delta)))
  list(status = srch$status, rho = srch$rho, coef = coef_full(s),
       range = log_lambda_range(delta[delta > 0], pad))
}

#' One-dimensional search over log(lambda)
#'
#' The smoothing parameter minimising ||fit_at(r) - target||^2 over a grid
#' of step `step` on the log-lambda range, refined by optimize() at an
#' interior grid minimum; the limits lambda = Inf and 0 (fit_inf, fit_zero)
#' compete and win ties. The range grows by `pad` (at most max_ext times)
#' while the grid minimum is at an end not within tol of that end's limit.
#' @param opt_tol optimize()'s tolerance in log(lambda).
#' @noRd
search_log_lambda <- function(fit_at, fit_inf, fit_zero, target, delta,
                              step, pad, tol, where, max_ext = 10L,
                              opt_tol = 1e-10) {
  pen   <- delta > 0
  dis   <- function(r) sum((fit_at(r) - target)^2)
  gap   <- function(r, lim) sqrt(sum((fit_at(r) - lim)^2))
  rg    <- log_lambda_range(delta[pen], pad)
  lo    <- rg[["lo"]]
  hi    <- rg[["hi"]]
  n_ext <- 0L
  repeat {
    grid <- seq(lo, hi, by = step)
    v    <- vapply(grid, dis, numeric(1))
    if (!all(is.finite(v)))
      stop_search_failed(where, "non-finite discrepancy")
    i     <- which.min(v)
    at_hi <- i == length(grid)
    at_lo <- i == 1L
    if (at_hi && gap(grid[i], fit_inf) > tol) {
      hi <- hi + pad
    } else if (at_lo && gap(grid[i], fit_zero) > tol) {
      lo <- lo - pad
    } else {
      break
    }
    n_ext <- n_ext + 1L
    if (n_ext > max_ext)
      stop_search_failed(where, "grid end not within tol of its limit")
  }
  cand <- c(inf  = sum((fit_inf - target)^2),
            zero = sum((fit_zero - target)^2))
  rho  <- NA_real_
  if (!at_hi && !at_lo) {
    o    <- stats::optimize(dis, grid[c(i - 1L, i + 1L)], tol = opt_tol)
    rho  <- o$minimum
    cand <- c(cand, finite = o$objective)
  }
  status <- names(cand)[which.min(cand)]
  list(status = status, rho = if (status == "finite") rho else NA_real_)
}

#' Log(lambda) range keeping every 1 / (1 + lambda delta) off its limits
#'
#' Every shrinkage factor stays more than 1 / (1 + e^pad) from 0 and 1.
#' @noRd
log_lambda_range <- function(delta_pen, pad) {
  c(lo = -log(max(delta_pen)) - pad, hi = -log(min(delta_pen)) + pad)
}

#' Discrepancy of smoother sm: squared W-norm of V (b - beta_b)
#' @noRd
smoother_discrepancy <- function(sm, b) {
  e <- b - sm$beta
  sum(e * as.vector(sm$VtWV %*% e))
}

#' Discrepancy of the model: sum of the smoother discrepancies
#' @noRd
model_discrepancy <- function(m, beta) {
  sum(vapply(seq_along(m$smoothers), function(j)
    smoother_discrepancy(m$smoothers[[j]], beta[[j]]), numeric(1)))
}

#' A root R with R^T R = M: Cholesky, or from the eigen-decomposition for a
#' singular M
#' @noRd
matrix_root <- function(M) {
  R <- tryCatch(chol(M), error = function(e)
    if (is_alloc_error(e)) stop(e) else NULL)
  if (!is.null(R)) return(R)
  # The learners of a smoother can be jointly aliased on the data (e.g. a
  # tensor product with fewer distinct locations than columns)
  e <- eigen(M, symmetric = TRUE)
  sqrt(pmax(e$values, 0)) * t(e$vectors)
}

#' Evaluate expr, re-signalling a failed log(lambda) search as
#' map_approx_failed with status "search_failed"
#' @noRd
reclass_search_failure <- function(expr) {
  tryCatch(expr, map_approx_search_failed = function(e)
    stop_map_failed(conditionMessage(e), "search_failed"))
}

#' Signal map_approx_failed
#' @noRd
stop_map_failed <- function(msg, status) stop(map_failed_condition(msg, status))

#' Condition signalled when the MAP approximation (Stage 1) fails
#' @noRd
map_failed_condition <- function(msg, status) {
  numerical_condition("map_approx_failed", paste("bayesian_backfitting: MAP",
                                                 "approximation failed:", msg),
                      status = status)
}

#' Signal a failed one-dimensional log(lambda) search
#' @noRd
stop_search_failed <- function(where, msg) {
  signal_numerical("map_approx_search_failed", paste0(where, ": ", msg))
}
