#' Derivative of the MLE intercept w.r.t. each observation
#'
#' Implicit differentiation of the intercept score equation
#' \eqn{\sum_i w_i \omega_i (y_i - \mu_i) = 0}, with
#' \eqn{\eta_i = o_i + a} and \eqn{\omega = h'(\eta) / V(\mu)}:
#' \deqn{\frac{da}{dy_j} = \frac{w_j \omega_j}{\sum_i w_i h'(\eta_i)
#'   (\omega_i - (y_i - \mu_i) \omega'_i)}.}
#' Exact for any link and any offset. For canonical links \eqn{\omega = 1}
#' and \eqn{\omega' = 0}; without an offset it reduces to
#' \eqn{w_j / (\sum_i w_i \, h'(a))}.
#'
#' @param family A family object (augmented with `dvar` and `d2link`).
#' @param eta0 Linear predictor vector at the MLE intercept (intercept + offset).
#' @param prior.weights A numeric vector of prior weights.
#' @param y The response vector.
#' @noRd
compute_da_dy_for_mle_a <- function(family, eta0, prior.weights, y) {
  den <- observed_intercept_information(family, eta0, prior.weights, y)
  prior.weights * omega(family$linkinv(eta0), family) / den
}

#' Observed information of the intercept
#'
#' Half the second derivative of the deviance in \eqn{a},
#' \eqn{\sum_i w_i h'(\eta_i) (\omega_i - (y_i - \mu_i) \omega'_i)}. It equals
#' the Fisher information for canonical links; for other links it can be
#' zero or negative away from a minimum of the deviance.
#'
#' @param family A family object (augmented with `dvar` and `d2link`).
#' @param eta Linear predictor vector (intercept + offset).
#' @param prior.weights A numeric vector of prior weights.
#' @param y The response vector.
#' @noRd
observed_intercept_information <- function(family, eta, prior.weights, y) {
  mu     <- family$linkinv(eta)
  mu_eta <- family$mu.eta(eta)
  om     <- omega(mu, family)
  dom    <- domega_dmu(mu, family)
  sum(prior.weights * mu_eta * (om - (y - mu) * dom))
}


#' Find intercept via maximum likelihood estimation
#'
#' Newton's method for the intercept \eqn{a} of the model
#' \eqn{\eta_i = o_i + a}, started at \eqn{g(\bar{y}_w)}, the link of the
#' weighted mean response. Without an offset this start is the exact MLE, on
#' the valid side of any pole of the link (e.g. \eqn{\eta = 0} for the
#' `inverse` and `1/mu^2` links).
#'
#' A step divides the score by the observed information, or by the Fisher
#' information where the observed one is not positive (the Fisher information
#' always is, so that step descends the deviance). Observed information is
#' used because Fisher steps become tiny near a pole, where the Fisher
#' information far exceeds the curvature of the deviance. Every step is
#' halved until \eqn{\eta} and \eqn{\mu} stay valid (`valideta()` and
#' `validmu()` of the family) and the deviance decreases sufficiently
#' (`halve_step()`); if no part of a step passes, the Fisher step is tried.
#' Where the predicted decrease is below the rounding error of the deviance,
#' or the deviance is not finite (a response outside the family's support),
#' the size of the score decides instead.
#'
#' Convergence is judged on the whole Fisher step, whichever step is taken.
#' An optimum on the boundary of the valid region, where only a tiny part of
#' the Fisher step stays valid, is returned as well.
#'
#' @param y A numeric vector of response values.
#' @param family A family object (augmented with `dvar` and `d2link`).
#' @param prior.weights A numeric vector of prior weights.
#' @param offset A numeric vector of offsets.
#' @noRd
find_intercept_mle <- function(y, family, prior.weights, offset) {
  offset    <- rep_len(offset, length(y))
  max_steps <- 200L
  tol       <- 1e-12  # converged: relative Fisher step
  tol_slow  <- 1e-8   # accepted once the iteration cannot move further
  min_step  <- 1e-15  # smallest step tried at a boundary
  # Rows with zero prior weight do not enter the score
  keep  <- prior.weights > 0
  w     <- prior.weights[keep]
  yk    <- y[keep]
  off_k <- offset[keep]
  canonical <- is_canonical_link(family)
  # Validity is checked on every row, zero-weight rows included
  valid <- function(a) {
    eta <- offset + a
    all(is.finite(eta)) &&
      (is.null(family$valideta) || isTRUE(family$valideta(eta))) &&
      (is.null(family$validmu) ||
         isTRUE(family$validmu(family$linkinv(eta))))
  }
  score_at <- function(a) {
    dd <- family_eta_mu_derivs(off_k + a, family)
    sum(w * (yk - dd$mu) * dd$mu_eta / dd$Vmu)
  }
  # Losses of a step: Inf at an invalid point, where they are undefined
  deviance_at <- function(a) {
    if (!valid(a)) return(Inf)
    sum(family$dev.resids(yk, family$linkinv(off_k + a), w))
  }
  score_size_at <- function(a) if (valid(a)) abs(score_at(a)) else Inf

  # Without an offset: the exact MLE, on the valid side of any pole
  a <- family$linkfun(sum(w * yk) / sum(w))
  if (length(a) != 1L || !is.finite(a) || !valid(a)) {
    stop("Could not find a valid starting intercept: the link of the ",
         "weighted mean response (", signif(a, 4), ") is outside the valid ",
         "range of the \"", family$link, "\" link (with the offset).",
         call. = FALSE)
  }
  # Newton's method, each step halved to stay valid and to reduce the deviance
  fisher <- Inf
  for (i in seq_len(max_steps)) {
    dd     <- family_eta_mu_derivs(off_k + a, family)
    score  <- sum(w * (yk - dd$mu) * dd$mu_eta / dd$Vmu)
    info   <- sum(w * dd$mu_eta^2 / dd$Vmu)
    fisher <- score / info
    if (!is.finite(fisher)) break
    if (abs(fisher) <= tol * (1 + abs(a))) return(a)
    steps <- fisher
    if (!canonical) {
      obs <- observed_intercept_information(family, off_k + a, w, yk)
      if (is.finite(obs) && obs > 0) steps <- c(score / obs, fisher)
    }
    dev   <- deviance_at(a)
    slack <- loss_allowance(dev, deviance_roundoff(family, off_k + a, w))
    for (step in steps) {
      srch <- if (!is.finite(dev) || abs(2 * score * step) <= slack) {
        halve_step(function(t) score_size_at(a + t * step), abs(score))
      } else {
        halve_step(function(t) deviance_at(a + t * step), dev, slack,
                   slope = -2 * score * step)
      }
      if (srch$fraction > 0) break
    }
    if (srch$fraction == 0 || a + srch$fraction * step == a) break
    a <- a + srch$fraction * step
  }
  # Accept a slowly converging solution, or one on the boundary of the valid
  # region, where the Fisher step stays valid only once it is tiny
  if (is.finite(fisher)) {
    while (!valid(a + fisher) && abs(fisher) > min_step) fisher <- fisher / 2
    if (abs(fisher) <= tol_slow * (1 + abs(a))) return(a)
  }
  stop("Could not find an intercept with approximately zero total score.",
       call. = FALSE)
}


# EDF derivatives ---------------------------------------------------------

#' Initialise EDF-tracking derivative objects on a fully assembled model list
#'
#' Must be called AFTER the learner/edf/stop controls have been merged into
#' `mod` (so that `mod$edf_compute_method`, `mod$z`, `mod$cweights`, and
#' `mod$initialise_intercept` are all present) and after `mod$y` and
#' `mod$linear.predictors` have been set to the starting values.
#'
#' Also sets \eqn{\widetilde{\operatorname{edf}}^{[0]}}{edf~^[0]},
#' `mod$edf.total` and the intercept learner's `edf`, to
#' \eqn{\mathrm{tr}(d\tilde{\mu}^{[0]}/dy)}: 0 when
#' `initialise_intercept = FALSE`, 1 for canonical links, and
#' \eqn{\sum_i h'(\eta_i) \, da/dy_i} otherwise (which differs from 1 only
#' for non-canonical links with an offset).
#'
#' @param mod A model list.
#' @param n Number of observations.
#' @param quiet If `TRUE`, suppress the non-canonical link message and the
#'   warning about Hutch++ with such a link (e.g. for bootstrap refits).
#' @return `mod` with `mod$derivs` set, or `mod` unchanged when
#'   `edf_compute_method == "none"`.
#' @noRd
init_derivs <- function(mod, n, quiet = FALSE) {
  if (mod$edf_compute_method == "none") return(mod)

  canonical <- is_canonical_link(mod$family)
  if (!canonical && !quiet) {
    message(
      "Non-canonical link function used (",
      mod$family$family, "/",
      mod$family$link,
      "). ",
      "Evaluating EDF at observed y, ",
      "which no longer corresponds to AIC (see Efron 2004)."
    )
    if (mod$edf_compute_method == "hutchpp") warn_hutchpp_non_canonical(mod)
  }

  derivs   <- list()
  dmu_deta <- mod$family$mu.eta(mod$linear.predictors)

  if (isTRUE(mod$initialise_intercept)) {
    dady <- compute_da_dy_for_mle_a(
      mod$family,
      eta0          = mod$linear.predictors,
      prior.weights = mod$prior.weights,
      y             = mod$y
    )
    # tr(dmu~0/dy) is exactly 1 for canonical links
    edf0 <- if (canonical) 1 else sum(dmu_deta * dady)
  } else {
    # The starting intercept (0) does not depend on y
    dady <- rep(0, n)
    edf0 <- 0
  }
  mod$edf.total <- edf0
  if (identical(mod$learners[[1L]]$var_label, "(Intercept)")) {
    mod$learners[[1L]]$edf <- edf0
  }

  if (mod$edf_compute_method == "hutchpp") {
    GGG             <- rrademacher(n, mod$z)
    dady            <- as.vector(matrix(dady, nrow = 1L) %*% GGG)
    derivs$GGG      <- GGG
    derivs$z        <- mod$z
    derivs$cweights <- mod$cweights
    derivs$zs       <- split_integer_by_weights(mod$z, mod$cweights)
    # The probe blocks are constant: split them once here
    derivs$Gs       <- split_mat_cols(GGG, derivs$zs)
  }

  derivs$detady   <- matrix(dady, n, length(dady), byrow = TRUE)
  derivs$dmudy    <- outer(dmu_deta, dady)

  mod$derivs <- derivs
  mod
}

#' Warn that Hutch++ is unreliable with a non-canonical link
#'
#' With such a link the Jacobian is far from symmetric, and the estimate of
#' edf~ can be dominated by noise (see the Details of gb_edf_ctrl()).
#' @param mod A model list.
#' @return `NULL`, invisibly.
#' @noRd
warn_hutchpp_non_canonical <- function(mod) {
  warning(
    "Hutch++ EDF are unreliable with the non-canonical link ",
    mod$family$family, "/", mod$family$link, ": the estimate of the model ",
    "edf can be far from the true value, even negative, so that stopping ",
    "by AIC and adaptive smoothness relaxation follow noise and the fit ",
    "can be poor. Use gb_edf_ctrl(\"exact\") (practical up to about ",
    "n = 2000). For larger data, use a canonical link, or ",
    "gb_edf_ctrl(\"none\") with stopping by cross-validation.",
    call. = FALSE
  )
  invisible(NULL)
}
