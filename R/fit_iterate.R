#' Name of the stopping loss, for printing
#'
#' @param mod Object with class `gbgam`.
#' @noRd
stop_loss_name <- function(mod) {
  switch(mod$stopping_method,
         aic    = "AIC",
         true_f = "MSE",
         cv     = if (is.null(mod$resampling)) "CV risk"
                  else resampling_risk_label(mod$resampling),
         "Loss")
}

#' Stopping loss of the current model
#'
#' @param mod Object with class `gbgam`.
#' @param stopping_method Stopping method from [gb_stop_ctrl()].
#' @param AIC_k Penalty multiplier of edf~ in the AIC.
#' @return Scalar loss, or `NULL` when no loss is tracked by the main fit.
#' @noRd
compute_loss <- function(mod, stopping_method, AIC_k = 2) {
  res <- switch(stopping_method,
    aic    = model_aic(mod, k = AIC_k),
    true_f = true_f_mse(mod),
    # cross-validated stopping: the resampling estimate is computed by
    # fit_gbgam_cv() across fits, not by the main fit
    cv     = NULL,
    none   = NULL)
  return(as.vector(res))
}

#' Stop when the stopping loss is not finite
#'
#' A non-finite loss would make every iteration look like an improvement
#' (`Inf > Inf` is `FALSE`), so loss-based stopping would never trigger.
#'
#' @param loss Scalar loss from `compute_loss()` or `NULL` (no loss tracked).
#' @param m Iteration number (0 for the initial model).
#' @noRd
check_loss_finite <- function(loss, m) {
  if (!is.null(loss) && !is.finite(loss)) {
    stop("Stopping loss became non-finite (", loss, ") at iteration ", m,
         "; check the response against the family (e.g. non-integer ",
         "Poisson counts or responses outside the family's support).",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Set the smoothing parameter of a learner during a replayed path
#'
#' Used for `preselected_rhos`. For Demmler-Reinsch learners, `nu_k` and
#' `lambda_at_zero` are updated together with `rho`, mirroring
#' `update_learner_trA_smoothness()`, so that a replay reproduces the
#' original fit also when `use_learner_step_size_for_update = TRUE`. `nu_k`
#' is left alone when `rho` does not change. As in the ASR update, `nu_k`
#' keeps the selection EDF at `target_tr2AminAA` also when `lambda` is set to
#' zero, giving
#' \eqn{\nu_k = 1 - \sqrt{1 - \mathrm{target} / p}}.
#'
#' @param lr A learner.
#' @param new_rho Scalar; the new log smoothing parameter.
#' @param target_tr2AminAA Target selection EDF used to recompute `nu_k`.
#' @return The updated learner.
#' @noRd
set_learner_rho <- function(lr, new_rho, target_tr2AminAA) {
  if (!inherits(lr, "demmler_reinsch") || identical(lr$rho, new_rho)) {
    lr$rho <- new_rho
    return(lr)
  }
  lr$rho            <- new_rho
  # check if new rho implies lambda=0
  lr$lambda_at_zero <- new_rho <= rho_lambda_zero / 10
  lr$nu_k <- if (identical(new_rho, lr$rho0)) {
    initial_nu_k(lr)
  } else {
    find_nu_k(lr, target_df = target_tr2AminAA)
  }
  lr
}



#' Compute -2 log L(y, mu, phi, wt) using Fletcher's dispersion estimate
#'
#' For binomial, uses Bernoulli likelihood and assumes y is a fraction and
#' prior.weights indicate the number of trials*weights.
#'
#' For Gaussian, Gamma, and inverse Gaussian, computes the log-likelihood
#' with phi estimated via `estimate_dispersion()`.
#'
#'
#' @param mod A `gbgam` model object.
#' @return Scalar: -2 * log-likelihood.
#' @noRd
neg2loglik <- function(mod) {
  mu  <- fitted_mu(mod)
  y   <- mod$y
  wt  <- mod$prior.weights
  fam <- mod$family$family
  edf <- mod$edf.total %||% 0
  # Only weights rescaled below the row count are rejected: with default or
  # frequency weights a noisy (Hutch++) or long-tailed edf~ may pass sum(wt),
  # and phi, clamped to 1e-12, then gives a huge AIC that is never the minimum
  if (!fam %in% c("poisson", "binomial") && sum(wt) < length(wt) &&
      sum(wt) - edf <= 0) {
    stop("The dispersion cannot be estimated: the prior weights sum to ",
         signif(sum(wt), 4), ", which does not exceed the model EDF (",
         round(edf, 2), "). The weights are frequency weights, so their ",
         "sum acts as the number of observations; rescale them, e.g. to ",
         "sum to the number of observations.", call. = FALSE)
  }
  phi <- estimate_dispersion(mod)

  ll <- switch(fam,
    gaussian = stats::dnorm(y, mean = mu, sd = sqrt(phi), log = TRUE),
    poisson  = {
      # dpois() for integer y; its continuous extension
      # y log(mu) - mu - log Gamma(y + 1) for non-integer y (dpois() would
      # return -Inf there and AIC stopping would never trigger)
      ll_p <- stats::dpois(round(y), lambda = mu, log = TRUE)
      frac <- y != round(y)
      ll_p[frac] <- y[frac] * log(mu[frac]) - mu[frac] - lgamma(y[frac] + 1)
      ll_p
    },
    binomial = y*log(mu) + (1 - y)*log(1 - mu),
    Gamma    = stats::dgamma(y, shape = 1/phi, rate = 1/(mu*phi), log = TRUE),
    # Closed-form inverse Gaussian log density with shape lambda = 1/phi
    inverse.gaussian = -0.5 * log(2 * pi * phi * y^3) -
      (y - mu)^2 / (2 * phi * mu^2 * y),
    stop("neg2loglik: unsupported family \"", fam, "\".", call. = FALSE)
  )
  return(-2 * sum(wt * ll))
}

#' AIC of the current model
#'
#' @param mod Object with class `gbgam`.
#' @param k Penalty multiplier of edf~.
#' @noRd
model_aic <- function(mod, k = 2) {
  neg2loglik(mod) + k * mod$edf.total
}

#' Weighted MSE of the linear predictor against `true_f`
#'
#' @param mod Object with class `gbgam`.
#' @noRd
true_f_mse <- function(mod) {
  r2 <- (mod$true_f - mod$linear.predictors)^2
  mse <- stats::weighted.mean(r2, w = mod$prior.weights)
  return(mse)
}


#' Construct GAM smooth terms
#'
#' @param mf model.frame with all variables appearing in `mgcv_terms`
#' @param mgcv_terms list of [`mgcv`][mgcv::mgcv-package] terms created with its
#'   constructor ([mgcv::t2()] in our case)
#' @param knots a named list (or NULL) for knot locations to be used in basis
#'   construction (see [mgcv::smoothCon()]).
#'
#' @return list of `smooth.spec` objects
#' @noRd
construct_smooths <- function(mf, mgcv_terms, knots=NULL){
  smooths <- list()
  for(j in seq_along(mgcv_terms)){
    term <- mgcv_terms[[j]]
    sms <- mgcv::smoothCon(
      term,
      data = mf,
      # We set our own constraints based on weighted averages
      absorb.cons = FALSE,
      # This argument doesn't change the predicted matrix in any way
      diagonal.penalty = FALSE,
      scale.penalty = FALSE,  # keep sp= on the scale of the user's penalty
      knots=knots
    )
    smooths <- append(smooths, sms)
  }
  return(smooths)
}

#' Initialise a family object
#'
#' @param family A family object, family function or family name.
#' @return The family object with the derivatives added by
#'   `mgcv::fix.family.var()` and `mgcv::fix.family.link()`.
#' @noRd
new_family <- function(family) {
  supported <- c("gaussian", "poisson", "binomial", "Gamma", "inverse.gaussian")
  if (checkmate::test_character(family, len = 1)) {
    checkmate::assert_choice(family, supported)
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) family <- family()
  if (!checkmate::test_class(family, "family")) {
    stop("family must be a character string, a family function, or a ",
         "\"family\" object.", call. = FALSE)
  }
  checkmate::assert_choice(family$family, supported)
  # Add derivatives
  family <- family |>
    mgcv::fix.family.var() |>
    mgcv::fix.family.link()
  return(family)
}



#' u = omega (y - mu): the gradient of the log-likelihood (the negative
#' gradient of the loss) w.r.t. eta
#'
#' Computes \eqn{(y - \mu) / (V(\mu) g'(\mu))}.
#'
#' @param y Response vector.
#' @param eta Linear predictor.
#' @param family A family object.
#' @noRd
pseudo_residuals <- function(y, eta, family){
  fam <- family
  mu <- fam$linkinv(eta)
  omega <- omega(mu, family)
  return(omega*(y-mu))
}

#' Training deviance at a linear predictor
#'
#' Sum of the family's deviance residuals at the prior weights over the
#' observations with positive weight (zero-weight observations, e.g. those
#' held out by cross-validation, are not fitted). A mean outside the
#' family's support gives `NaN` (without a warning).
#'
#' @param mod Object with class `gbgam`.
#' @param eta Linear predictor.
#' @noRd
train_deviance <- function(mod, eta) {
  pos <- mod$prior.weights > 0
  mu  <- mod$family$linkinv(eta[pos])
  suppressWarnings(
    sum(mod$family$dev.resids(mod$y[pos], mu, wt = mod$prior.weights[pos])))
}

#' Stop when the negative gradient is not finite
#'
#' @param u Negative gradient (working residuals).
#' @param w Prior weights; observations with zero weight are not checked.
#' @param iter Iteration number, or `NULL`.
#' @noRd
check_finite_residuals <- function(u, w, iter = NULL) {
  if (!all(is.finite(u[w > 0]))) stop_infinite_residuals(iter)
  invisible(NULL)
}

#' Stop with `mboost`'s "Infinite residuals" error
#'
#' @param iter Iteration number, or `NULL`.
#' @noRd
stop_infinite_residuals <- function(iter = NULL) {
  stop_diverged(paste0(
    "Infinite residuals at ", iteration_label(iter), ": the step size is too ",
    "large for this data; decrease step_size in gb_learner_ctrl()."), iter)
}

#' Stop with an error of class `gbgam_diverged`
#'
#' Scripts can catch the class, e.g. to refit with a smaller step size.
#'
#' @param msg Error message.
#' @param iter Iteration number, or `NULL`.
#' @noRd
stop_diverged <- function(msg, iter = NULL) {
  stop(structure(
    class = c("gbgam_diverged", "error", "condition"),
    list(message = msg, call = NULL, iter = iter)
  ))
}

#' Label of an iteration for error messages
#'
#' @param iter Iteration number, or `NULL`.
#' @noRd
iteration_label <- function(iter) {
  if (is.null(iter)) "this iteration" else paste("iteration", iter)
}

#' Absolute limit of the rounding error in the deviance
#'
#' A change in the deviance is normally judged relative to the deviance. That
#' fails for a model that reproduces the data almost exactly: the deviance is
#' then close to 0, but the rounding error in it is not, so any tiny increase
#' would look large. This function gives an absolute limit for that case; see
#' `loss_allowance()`.
#'
#' The limit is the deviance caused by a relative error of
#' `sqrt(.Machine$double.eps)` in every mean. An error \eqn{\delta \mu_i} in
#' the mean adds about \eqn{w_i (\delta \mu_i)^2 / V(\mu_i)} to the deviance,
#' so with \eqn{\delta^2} = `.Machine$double.eps` the limit is
#' `.Machine$double.eps` times the sum of \eqn{w_i \mu_i^2 / V(\mu_i)}, over
#' the observations with positive weight. Terms that are not finite (a
#' binomial mean of exactly 1, where \eqn{V(\mu_i) = 0}) are left out.
#'
#' @param family A family object.
#' @param eta Linear predictor.
#' @param w Prior weights.
#' @return The limit, a non-negative number.
#' @noRd
deviance_roundoff <- function(family, eta, w) {
  pos <- w > 0
  mu  <- family$linkinv(eta[pos])
  r   <- w[pos] * mu^2 / family$variance(mu)
  .Machine$double.eps * sum(r[is.finite(r)])
}

#' Stop when a boosting update diverges
#'
#' With a stable step size every update moves the linear predictor along a
#' descent direction of the training deviance (the learners fit the negative
#' gradient by weighted least squares), so the deviance cannot increase. The
#' function stops with an informative error when the updated linear
#' predictor, its mean or its negative gradient is not finite, or when the
#' training deviance increases beyond round-off (`loss_allowance()` with
#' `deviance_roundoff()`). It is called before the
#' derivative update, so that a diverging step never reaches the EDF
#' computation.
#'
#' @param mod Object with class `gbgam` before the update.
#' @param eta_new Updated linear predictor.
#' @param iter Iteration number, or `NULL`.
#' @noRd
check_boost_step <- function(mod, eta_new, iter = NULL) {
  fam    <- mod$family
  w      <- mod$prior.weights
  mu_new <- fam$linkinv(eta_new)
  if (!all(is.finite(eta_new)) || !all(is.finite(mu_new))) {
    stop_infinite_residuals(iter)
  }
  check_finite_residuals(pseudo_residuals(mod$y, eta_new, fam), w, iter)

  dev_old <- train_deviance(mod, mod$linear.predictors)
  dev_new <- train_deviance(mod, eta_new)
  if (!is.finite(dev_new)) {
    stop_diverged(paste0(
      "The training deviance became non-finite at ", iteration_label(iter),
      ": the linear predictor has left the valid range of the '", fam$link,
      "' link. Decrease step_size in gb_learner_ctrl() or use a link such ",
      "as 'log'."), iter)
  }
  incr <- dev_new - dev_old
  if (incr > 0 && incr > loss_allowance(
    dev_old, deviance_roundoff(fam, mod$linear.predictors, w))) {
    # The increase can be tiny relative to the deviance, so it is printed
    # itself instead of the two (then equal-looking) deviances
    rel <- if (dev_old > 0) {
      paste0(", ", format(signif(incr / dev_old, 2)), " of its value")
    } else {
      ""
    }
    stop_diverged(paste0(
      "Boosting diverged at ", iteration_label(iter),
      ": the training deviance increased by ", format(signif(incr, 3)), rel,
      " (from ",
      format(signif(dev_old, 6)), "), so the step size is too large for ",
      "this data. Decrease step_size in gb_learner_ctrl()."), iter)
  }
  invisible(NULL)
}

#' Loss for learner selection
#'
#' @param learner A base learner.
#' @param learner_fit Its fit to the negative gradient (with `$ghat`).
#' @param mod Object with class `gbgam`.
#' @param u Negative gradient of the loss (`pseudo_residuals()`).
#' @param ghat_int Fitted values of the intercept-only learner for the same
#'   gradient (intercept padding) or `NULL`. It shifts the provisional
#'   deviance update; the `"rss"` criterion ignores it, because the padded
#'   weighted RSS differs from the unpadded one by a learner-independent
#'   constant.
#' @noRd
learner_selection_loss <- function(learner, learner_fit, mod, u,
                                   ghat_int = NULL) {
  nu_k <- learner$nu_k %||% 1
  w    <- mod$prior.weights
  nu   <- if (mod$use_step_size_for_selection) mod$step_size else 1

  if (mod$selection_loss_type == "deviance") {
    eta_hat <- (nu * nu_k * learner_fit$ghat) + mod$linear.predictors
    if (!is.null(ghat_int)) eta_hat <- eta_hat + nu * ghat_int
    fam <- mod$family
    if (fam$family == "poisson" && fam$link == "log") {
      # Poisson deviance up to a learner-independent constant and a factor
      # of 2 (avoids log/which in dev.resids); only the argmin over learners
      # is used, so this is NOT an actual deviance value.
      total_loss <- sum(w * (exp(eta_hat) - mod$y * eta_hat))
    } else {
      mu_hat     <- fam$linkinv(eta_hat)
      total_loss <- sum(fam$dev.resids(mod$y, mu = mu_hat, wt = w))
    }
  } else {
    total_loss <- sum(w * (u - (nu * nu_k) * learner_fit$ghat)^2)
  }
  total_loss
}

#' Perform a single boosting iteration
#' @param mod Object with class `gbgam`
#' @param k_opt integer index of base learner to select.
#' If NULL, selects base learner based on the data.
#' @param exclude integer vector of base learners to exclude from
#'   consideration, or `NULL` (default). A test hook: fits use
#'   `mod$excluded_learner_indices`.
#' @param iter Iteration number for error messages, or `NULL`.
#' @noRd
iterate_mod <- function(mod, k_opt = NULL, exclude=NULL, iter = NULL){
  learners <- mod$learners
  step_size <- mod$step_size
  # Compute gradient --------------------------------------------------------
  u <- pseudo_residuals(mod$y, mod$linear.predictors, family = mod$family)
  check_finite_residuals(u, mod$prior.weights, iter)
  Wu <- mod$prior.weights*u

  # Intercept padding: the intercept-only learner (learner 1) is fitted to
  # the same gradient as every candidate and never competes on its own
  pad     <- isTRUE(mod$intercept_padding)
  fit_int <- if (pad) fit_learner(learners[[1L]], Wu) else NULL

  # Select learner ----------------------------------------------------------
  if(!is.null(k_opt)){
    fit_opt <- fit_learner(learners[[k_opt]], Wu)
  }else{
     # Initialise
    K <- length(learners)
    learner_fits <- vector("list", K)
    learner_loss <- rep(NA_real_, K)
    skip <- c(exclude, mod$excluded_learner_indices, if (pad) 1L)
    # Fit
    for(k in 1:K){
      if (k %in% skip) next
      # Fit the learner
      lr <- learners[[k]]
      learner_fits[[k]] <- fit_learner(lr, Wu)
      learner_loss[k] <- learner_selection_loss(
        learner=lr, learner_fit=learner_fits[[k]],
        mod=mod, u=u, ghat_int = fit_int$ghat)

    }

    if (!any(is.finite(learner_loss))) {
      if (all(vapply(learner_fits, is.null, logical(1)))) {
        stop("No base learner is left to select at this iteration.",
             call. = FALSE)
      }
      stop_diverged(paste0(
        "All candidate base learners give a non-finite selection loss ",
        "at this iteration; the linear predictor has left the valid ",
        "range of the '", mod$family$link, "' link. Reduce step_size ",
        "or use a link such as 'log'."), iter)
    }
    k_opt <- which.min(learner_loss) # Select
    fit_opt <- learner_fits[[k_opt]]
  }
  # Update ------------------------------------------------------------------
  learner_opt <- learners[[k_opt]]
  apply_nu_k <- mod$use_learner_step_size_for_update &&
    !is.null(learner_opt$nu_k)
  total_step_size <- if(apply_nu_k) step_size*learner_opt$nu_k else step_size
  dbeta <- total_step_size * fit_opt$b
  deta <- total_step_size * fit_opt$ghat

  # Pass the step actually applied so the derivative update matches it
  new_mod <- update_model(mod, k_opt, dbeta, deta, step_size = total_step_size,
                          fit_int = fit_int, iter = iter)

  return(new_mod)
}

#' Apply the update of the selected learner to the model
#'
#' @param mod Object with class `gbgam`.
#' @param k_opt Index of the selected learner.
#' @param dbeta,deta Coefficient and linear predictor updates.
#' @param step_size The step actually applied to the selected learner
#'   (`step_size * nu_k` when `use_learner_step_size_for_update = TRUE`).
#' @param fit_int Fit of the intercept-only learner (learner 1) to the same
#'   gradient (intercept padding), or `NULL`. It is added with the global
#'   `mod$step_size` in the same iteration. The iteration's change of edf~ is
#'   split: the trace share of the intercept update goes to learner 1, the
#'   remainder to `k_opt` (and drives its ASR update).
#' @param iter Iteration number for error messages, or `NULL`.
#' @noRd
update_model <- function(mod, k_opt, dbeta, deta, step_size = mod$step_size,
                         fit_int = NULL, iter = NULL){

  pad <- !is.null(fit_int)
  if (pad) {
    dbeta_int <- mod$step_size * fit_int$b
    deta_int  <- mod$step_size * fit_int$ghat
    deta_tot  <- deta + deta_int
  } else {
    deta_tot  <- deta
  }

  # Stop on a diverging step before it reaches the derivatives and edf~
  eta_new  <- mod$linear.predictors+deta_tot
  check_boost_step(mod, eta_new, iter)

  # Update derivatives and edf~
  dedf_int <- NULL
  if(mod$edf_compute_method == "none"){
    derivs <- dedf <- NULL
  }else{
    # Update dmudy and detady (sketched, i.e. %*% GGG, under NA-Hutch++)
    derivs <- update_derivatives(
      mod, eta_new, k_opt, step_size = step_size,
      int_step_size = if (pad) mod$step_size else NULL)
    # Trace share of the padded intercept update (NULL when not padded)
    dedf_int        <- derivs$dedf_int
    derivs$dedf_int <- NULL
    # Compute edf~
    edf <- model_edf(mod, derivs)
    dedf <- edf-mod$edf.total

  }
  # edf~ credited to k_opt: the iteration's change of edf~ minus the intercept's
  # share, so that ASR of k_opt only sees its own part
  dedf_k <- if (is.null(dedf_int)) dedf else dedf - dedf_int

  # Update learner: edf~_k, coefs, fitted values, smoothing parameter
  mod$learners[[k_opt]] <- update_learner(
    mod$learners[[k_opt]],
    mod=mod,
    dbeta=dbeta,
    deta=deta,
    dedf=dedf_k)
  # Padded intercept: coefficient, fitted values and its edf~ share (no ASR,
  # the intercept is unpenalised)
  if (pad) {
    mod$learners[[1L]] <- update_learner(
      mod$learners[[1L]],
      mod   = mod,
      dbeta = dbeta_int,
      deta  = deta_int,
      dedf  = dedf_int)
  }
  # Update model
  update_gbgam(mod, dedf = dedf, deta = deta_tot, derivs = derivs,
               k_opt = k_opt)
}


#' Store the updated model-level state
#'
#' @param object Object with class `gbgam`.
#' @param dedf Change of edf~, or `NULL` (not tracked).
#' @param deta Change of the linear predictor.
#' @param derivs Updated EDF-tracking derivatives.
#' @param k_opt Index of the selected learner.
#' @noRd
update_gbgam <- function(object, dedf, deta, derivs, k_opt) {
  if (!is.null(dedf)) object$edf.total <- object$edf.total + dedf
  object$linear.predictors <- object$linear.predictors + deta
  object$derivs            <- derivs
  object$k_opt             <- k_opt
  object
}
