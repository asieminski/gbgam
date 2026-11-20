# Fitting and updating a base learner, and its smoothness: adaptive
# smoothness relaxation (ASR), the smoothing parameter solver, the learner's
# edf and its step size nu_k.


# LEARNER COMPUTATION -----------------------------------------------------


#' Fit learner obj to the weighted pseudo-residuals Wg
#'
#' The (shrunken) least-squares coefficients b of Wg on the learner's design
#' and the fitted values ghat.
#' @noRd
fit_learner <- function(obj, Wg) {
  resp <- Wg
  Xp   <- obj$design
  if (is.matrix(Xp) && ncol(Xp) == 1) Xp <- as.vector(Xp)

  if (inherits(obj, "demmler_reinsch")) {
    b_up <- as.vector(crossprod(Xp, resp))
    if (obj$lambda_at_zero) {
      b <- b_up
    } else {
      lam <- exp(obj$rho)
      b   <- as.vector(b_up / (1 + lam * obj$D))
    }
  } else if (inherits(obj, "unpenalised")) {
    # A column with zero weighted norm (e.g. a factor level without weight
    # in a resampling fit) carries no information: its fit is zero
    b <- if (obj$xtWx > 0) as.vector(crossprod(Xp, resp) / obj$xtWx) else 0
  }
  ghat <- if (is.matrix(Xp)) as.vector(Xp %*% b) else Xp * b

  return(list(b = b, ghat = ghat))
}


# LEARNER UPDATES ---------------------------------------------------------

#' Add one boosting step to a learner and relax its smoothness (ASR)
#' @noRd
update_learner <- function(object, mod, dbeta, deta, dedf=NULL){

  object$beta <- object$beta + as.vector(dbeta)
  object$eta  <- object$eta + as.vector(deta)
  if (!is.null(dedf)) {
    object$edf <- object$edf + dedf

    # Relax the smoothness (ASR)
    if (inherits(object, "demmler_reinsch")) {
      # Per-update progress stays silent; the one-off lambda-to-zero message
      # follows gbgam(verbose = )
      object <- update_learner_trA_smoothness(
        learner = object, tol = 1e-6,
        target_tr2AminAA = mod$target_tr2AminAA, tau = mod$tau,
        verbose = FALSE, notify_lambda_zero = isTRUE(mod$verbose)
      )
    }
  }
  return(object)
}


# log(lambda) of a learner whose lambda is set to zero
rho_lambda_zero <- -1e16

#' Relax the smoothness of a penalised learner
#'
#' @description
#' One step of adaptive smoothness relaxation (ASR) for one penalised
#' (Demmler-Reinsch) learner: lowers its smoothing parameter so that its
#' ledf rises to `tau` times its accumulated edf, then recomputes its step
#' size nu_k so that its sedf stays at `target_tr2AminAA`. [gbgam()] calls
#' it after every iteration that selects the learner.
#'
#' @param learner A penalised learner (class `"demmler_reinsch"`), from
#'   `mod$learners` or `new_learners()$learners`, whose field `edf` holds
#'   its accumulated edf \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k}.
#' @param tol Positive number. Convergence tolerance of Newton's method on
#'   `rho`. Default `1e-6`, the value [gbgam()] uses.
#' @param target_tr2AminAA Positive number. The sedf the learner keeps
#'   after nu_k is recomputed. Default `1`, the value [gbgam()] uses.
#' @param tau Non-negative number. Multiplier of the accumulated edf that
#'   gives the target ledf; values below one slow the relaxation down.
#'   Default `1`.
#' @param verbose Logical. Report every smoothness update with [message()].
#'   Default `FALSE`.
#' @param notify_lambda_zero Logical. Give a message when the smoothing
#'   parameter is set to zero. Default `verbose`.
#'
#' @return The learner with updated `rho` and `nu_k`, and with
#'   `lambda_at_zero = TRUE` once the smoothing parameter is set to zero. Its
#'   other fields are unchanged.
#'
#' @details
#' The hat matrix \eqn{\mathbf{A}_k}{A_k} of the learner has the eigenvalues
#' 1 / (1 + lambda_k d), where lambda_k = exp(`rho`) is its smoothing
#' parameter and d runs over the entries of its field `D`. Its ledf is
#' \eqn{\operatorname{tr}(\mathbf{A}_k)}{tr(A_k)}, and the sedf of the
#' learner scaled by its step size nu_k is
#' \eqn{2 \nu_k \operatorname{tr}(\mathbf{A}_k) - \nu_k^2
#' \operatorname{tr}(\mathbf{A}_k^2)}{2 nu_k tr(A_k) - nu_k^2 tr(A_k^2)}.
#'
#' ## What happens
#'
#' With the target ledf `tau * edf` and \eqn{p} the number of columns of the
#' learner, the first of these that applies decides the result:
#' 1. `lambda_at_zero` is `TRUE`: the learner is returned as it is.
#' 1. The target is at least \eqn{p}: the smoothing parameter is set to zero
#'    (`rho` is then a large negative number, and `lambda_at_zero` is `TRUE`)
#'    and `nu_k` is recomputed.
#' 1. The ledf is more than `1e-6` below the target: `rho` is lowered until
#'    the ledf equals the target, and `nu_k` is recomputed.
#' 1. Otherwise the ledf already reaches the target and the learner is
#'    returned as it is, so `target_tr2AminAA` has no effect and the smoothing
#'    parameter is never raised.
#'
#' ## Finding `rho`
#'
#' The ledf falls from \eqn{p} to 0 as `rho` rises, so every target between
#' them has one solution, which lies between two bounds known in closed form.
#' `rho` is found by Newton's method, each step halved until the distance of
#' the ledf from the target does not grow; a step that cannot be used is
#' replaced by bisection between the bounds, which close in on the solution
#' with every evaluation. The search has converged when a whole Newton step
#' or the distance between the bounds is smaller than `tol`, or when the
#' distance of the ledf from the target is at the rounding error of the ledf
#' and no step reduces it. More than 1000 iterations result in an error, as
#' does a `target_tr2AminAA` that no step size reaches.
#'
#' ## Fields of the learner
#'
#' The function reads `edf`, `D`, `rho`, `lambda_at_zero`, the length of
#' `beta` (\eqn{p}) and, for its messages, `var_label`. It writes `rho`,
#' `nu_k` and `lambda_at_zero`. It does not refit the learner: `beta`, `eta`,
#' `design`, `edf`, `rho0` and `nu_k0` keep their values, and the new `rho`
#' and `nu_k` act from the next boosting update of the learner. See
#' [new_learners()] for what each field holds.
#'
#' @seealso [new_learners()] for the fields of a learner.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = runif(100), y = rnorm(100))
#' lb  <- new_learners(y ~ t2(x, k = 8, bs = "cr"), data = dat)
#' lr  <- lb$learners[[2]]                  # the penalised learner of t2(x)
#' lr$edf <- 3                              # accumulated edf~_k
#' lr2 <- update_learner_trA_smoothness(lr)
#' c(rho = lr$rho, rho_new = lr2$rho)
#'
#' # The eigenvalues of the hat matrix give the ledf and the sedf
#' s <- 1 / (1 + lr2$D * exp(lr2$rho))
#' sum(s)                                   # ledf tr(A_k) = 3
#' 2 * lr2$nu_k * sum(s) - lr2$nu_k^2 * sum(s^2)  # sedf 1
#'
#' # tau slows the relaxation down: here the target ledf is 1.5
#' lr3 <- update_learner_trA_smoothness(lr, tau = 0.5)
#' sum(1 / (1 + lr3$D * exp(lr3$rho)))
#'
#' # A target above the number of columns sets lambda to zero
#' lr$edf <- length(lr$beta) + 1
#' update_learner_trA_smoothness(lr)$lambda_at_zero
#' @noRd
update_learner_trA_smoothness <- function(learner, tol = 1e-6,
                                          target_tr2AminAA = 1,
                                          tau = 1, verbose = FALSE,
                                          notify_lambda_zero = verbose) {
  check_trA_smoothness_args(learner, tol, target_tr2AminAA, tau, verbose,
                            notify_lambda_zero)

  # Once lambda is already at zero there is nothing to update: skip silently.
  if (isTRUE(learner$lambda_at_zero)) return(learner)

  target_trA <- learner$edf * tau  # target ledf: tau times edf~_k

  # Only lambda = 0 reaches a target at the basis dimension or above it
  p <- length(learner$beta)
  if (target_trA >= p) {
    if (notify_lambda_zero) {
      message(sprintf(
        paste("%s has %.2f edf, but its basis has only %d columns.",
              "Setting lambda to zero instead."),
        learner$var_label, target_trA, p
      ))
    }

    learner$rho            <- rho_lambda_zero
    learner$lambda_at_zero <- TRUE
    # Keep the selection edf at its target after setting lambda to zero:
    # with lambda = 0 this gives nu_k = 1 - sqrt(1 - target / p)
    learner$nu_k <- find_nu_k(learner = learner,
                              target_df = target_tr2AminAA)
  } else {
    gap <- sum(1 / (1 + learner$D * exp(learner$rho))) - target_trA
    # Relax only when the learner trace is below its target
    if (gap < -1e-6) {
      learner$rho <- solve_rho_trA(learner$rho, D = learner$D,
                                   target_trA = target_trA, tol = tol,
                                   label = learner$var_label)
      if (verbose) {
        message(sprintf("%s updated to tr(A) = %.2f", learner$var_label,
                        target_trA))
      }
      # Recompute nu_k for unbiased selection
      learner$nu_k <- find_nu_k(learner = learner,
                                target_df = target_tr2AminAA)
    }
  }

  return(learner)
}

#' Check the arguments of update_learner_trA_smoothness()
#'
#' @inheritParams update_learner_trA_smoothness
#' @noRd
check_trA_smoothness_args <- function(learner, tol, target_tr2AminAA, tau,
                                      verbose, notify_lambda_zero) {
  checkmate::assert_class(learner, "demmler_reinsch", .var.name = "learner")
  checkmate::assert_number(tol, lower = 0, finite = TRUE)
  if (tol <= 0) stop("tol must be positive.", call. = FALSE)
  checkmate::assert_number(target_tr2AminAA, lower = 0, finite = TRUE)
  if (target_tr2AminAA <= 0) {
    stop("target_tr2AminAA must be positive.", call. = FALSE)
  }
  checkmate::assert_number(tau, lower = 0, finite = TRUE)
  checkmate::assert_flag(verbose)
  checkmate::assert_flag(notify_lambda_zero)
  invisible(NULL)
}


#' Find the log smoothing parameter giving the desired edf
#'
#' The edf decreases in rho; the search interval c(-32, 32) is extended when
#' it does not bracket the root (e.g. for very large weights, which shift
#' the root by log of their scale).
#' @param D Diagonal smoothing matrix.
#' @param edf_type One of "tr2A_AA" (the sedf), "trA", or "trAA". Selects
#'   the desired edf formulation.
#' @param target_edf what edf do we want the learner to have
#' @param label Learner name for the error message.
#'
#' @returns rho
#' @noRd
init_rho <- function(D, target_edf, edf_type = "tr2A_AA", label = "learner") {
  edf_gap <- function(rho) {
    learner_edf(1 / (1 + exp(rho) * D), type = edf_type) - target_edf
  }
  tryCatch(
    stats::uniroot(edf_gap, interval = c(-32, 32), extendInt = "downX")$root,
    error = function(e) {
      stop("Could not find the smoothing parameter giving edf ", target_edf,
           " for ", label, ": ", conditionMessage(e), call. = FALSE)
    })
}

#' rho giving the learner trace tr(A) = target_trA
#'
#' tr(A) falls from p to 0 as rho rises, so the root is unique and lies in
#' the bracket of bracket_rho_trA(), which every evaluated rho narrows. A
#' Newton step inside the bracket is halved until |tr(A) - target_trA| does
#' not grow. A step that is not finite, leaves the bracket or cannot be
#' halved to an improvement is replaced by the midpoint of the bracket: past
#' the root, tr(A) can be flat to rounding, with no slope for Newton.
#'
#' Stops once a whole Newton step or the bracket is smaller than tol, or once
#' the gap is at the rounding error of the trace and no step reduces it,
#' where the Newton step is rounding error divided by a tiny slope. Results
#' in an error after max_iter iterations.
#' @param target_trA Target trace, between 0 and length(D), both excluded.
#' @param label Learner name for the error message.
#' @noRd
solve_rho_trA <- function(rho, D, target_trA, tol, label, max_iter = 1000L) {
  trA_gap  <- function(r) sum(1 / (1 + D * exp(r))) - target_trA
  bracket  <- bracket_rho_trA(D, target_trA)
  lower    <- bracket[1L]
  upper    <- bracket[2L]
  gap      <- trA_gap(rho)
  rounding <- 4 * length(D) * .Machine$double.eps
  for (i in seq_len(max_iter)) {
    if (gap > 0) lower <- max(lower, rho) else upper <- min(upper, rho)
    drho   <- newton_trA_update(rho, D = D, target_trA = target_trA)
    inside <- is.finite(drho) && rho + drho >= lower && rho + drho <= upper
    # Halve until the step is no worse than drho = 0
    srch <- if (inside) {
      halve_step(function(t) abs(trA_gap(rho + t * drho)), abs(gap),
                 max_halvings = 1000L)
    } else {
      list(fraction = 0)
    }
    if (srch$fraction == 0) {
      rho <- (lower + upper) / 2
      if (upper - lower < tol) return(rho)
      gap <- trA_gap(rho)
      next
    }
    rho     <- rho + srch$fraction * drho
    stalled <- srch$loss >= abs(gap) && srch$loss <= rounding
    gap     <- trA_gap(rho)
    if (abs(drho) < tol || stalled) return(rho)
  }
  stop("update_learner_trA_smoothness: the search for rho of ", label,
       " did not converge in ", max_iter, " iterations (rho = ", rho,
       ", target tr(A) = ", target_trA, ").", call. = FALSE)
}

#' Interval of rho that holds the root of tr(A) = target_trA
#'
#' With a = lambda D_i, tr(A) = sum_i 1 / (1 + a_i). Jensen's inequality gives
#' tr(A) >= p / (1 + lambda mean(D)), and 1 / (1 + a_i) <= 1 / a_i gives
#' tr(A) <= sum_i 1 / (lambda D_i); setting each bound to the target gives
#' the two ends.
#' @param D Positive eigenvalues of the penalty.
#' @param target_trA Target trace, between 0 and length(D), both excluded.
#' @return c(lower, upper).
#' @noRd
bracket_rho_trA <- function(D, target_trA) {
  p <- length(D)
  c(log((p - target_trA) / target_trA / mean(D)),
    log(sum(1 / D) / target_trA))
}

#' Newton step in rho towards tr(A) = target_trA (\eqn{\mathcal{O}(p)})
#' @noRd
newton_trA_update <- function(rho_old, D, target_trA) {
  lam    <- exp(rho_old)
  lamD   <- lam * D
  s      <- 1 / (1 + lamD)
  frho   <- sum(s) - target_trA
  dfdrho <- sum(-(s^2) * lamD)
  drho   <- -frho / dfdrho
  return(drho)
}

#' Edf of a learner with shrinkage factors s: tr(2A - A^2) (the sedf),
#' tr(A) or tr(A^2)
#' @noRd
learner_edf <- function(shrinkage_factors, type = c("tr2A_AA", "trA", "trAA")) {
  type <- match.arg(type)
  s    <- shrinkage_factors
  switch(type,
         tr2A_AA = 2 * sum(s) - sum(s^2),  # sedf = tr(2A - A^2)
         trA      = sum(s),
         trAA    = sum(s^2)
  )
}


#' Learner step size nu_k giving selection edf target_df
#'
#' Solves 2 a tr(A) - a^2 tr(A^2) = target_df for a; the smaller root (the
#' larger gives the same sedf with a step past the peak).
#' @noRd
find_nu_k <- function(learner, target_df = 1) {
  D     <- learner$D
  rho   <- learner$rho
  s     <- 1 / (1 + exp(rho) * D)
  trA   <- sum(s)
  trAA <- sum(s^2)
  Delta <- 4 * (trA^2 - trAA * target_df)
  # Delta < 0 (or 0 / 0 below): the target is above the largest sedf,
  # tr(A)^2 / tr(A^2), that any step size gives
  unreachable <- function() {
    stop("find_nu_k: no step size gives ", learner$var_label, " sedf ",
         target_df, "; the largest reachable is ", trA^2 / trAA, ".",
         call. = FALSE)
  }
  if (isTRUE(Delta < 0)) unreachable()
  nu_k <- (2 * trA - sqrt(Delta)) / (2 * trAA)
  if (is.nan(nu_k)) unreachable()
  return(nu_k)
}

#' The learner's step size nu_k at construction (selection edf 1, capped at 1)
#' Computed afresh for a learner without `nu_k0`.
#' @noRd
initial_nu_k <- function(lr) {
  lr$nu_k0 %||% min(1, find_nu_k(lr, target_df = 1))
}
