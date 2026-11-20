# Finalise ----------------------------------------------------------------

#' Attach the `mgcv`-compatible fields to a fitted model
#'
#' Run after the boosting loop. Separated from gbgam() to keep the main
#' function readable and to make the finalisation step testable on its own.
#'
#' @param final_model List returned by fit_gbgam().
#' @param formula The original user formula (stored verbatim).
#' @param gp Output of mgcv::interpret.gam(formula).
#'
#' @return final_model with class c("gbgam","gam","glm","lm") and all fields
#'   required by `mgcv` predict / plot / summary methods.
#' @noRd
finalise_gbgam <- function(final_model, formula, gp) {

  final_model$formula      <- formula
  final_model$pred.formula <- gp$pred.formula  # used by predict.gam for new data
  final_model$fitted.values <- as.vector(
    final_model$family$linkinv(final_model$linear.predictors))

  # glm/mgcv convention: $residuals are the working residuals
  # (y - mu) / mu'(eta) and $weights are the working weights
  # w * mu'(eta)^2 / V(mu). mgcv::plot.gam() builds partial residuals from
  # residuals * sqrt(weights). The user weights are kept in $user_weights.
  mu_eta <- final_model$family$mu.eta(final_model$linear.predictors)
  final_model$residuals <-
    (final_model$y - final_model$fitted.values) / mu_eta
  final_model$weights <- as.vector(
    final_model$prior.weights * mu_eta^2 /
      final_model$family$variance(final_model$fitted.values))

  # Per-coefficient EDF: each learner's edf~_k is spread evenly over its
  # original basis columns (nrow(A) for a Demmler-Reinsch learner, which may
  # have fewer coefficients), and 0 is put back at the dropped ind0var
  # columns, as in learners_to_mgcv_coef() (all NA without EDF tracking).
  edf_compact <- unlist(lapply(final_model$learners, function(lr) {
    p <- if (is.null(lr$A)) length(lr$beta) else nrow(lr$A)
    rep(lr$edf / p, p)
  }))
  ind0var <- final_model$ind0var
  if (length(ind0var) > 0L) {
    edf_coef <- rep(0, length(final_model$coefficients))
    edf_coef[-ind0var] <- edf_compact
  } else {
    edf_coef <- edf_compact
  }
  # Without EDF tracking every EDF field is NA rather than 0
  edf_tracked <- final_model$edf_compute_method != "none"
  if (!edf_tracked) edf_coef[] <- NA_real_
  final_model$edf1 <- final_model$edf <- edf_coef
  final_model$edf.total <- sum(edf_coef)

  # Per-smooth EDF: sum coefficient-level EDFs over each smooth's parameter
  # range so that mgcv's plot.gam can label axes as "t2(x,3.43)"
  for (j in seq_along(final_model$smooth)) {
    sm <- final_model$smooth[[j]]
    final_model$smooth[[j]]$edf <-
      sum(final_model$edf[sm$first.para:sm$last.para])
  }
  final_model$cmX <- model_col_means(final_model)

  final_model$rank <- if (edf_tracked) {
    sum(final_model$edf > sqrt(.Machine$double.eps))
  } else {
    NA_integer_
  }
  final_model$df.residual <- nrow(final_model$model) - final_model$edf.total

  final_model$deviance <- sum(
    final_model$family$dev.resids(
      final_model$y, final_model$fitted.values, final_model$prior.weights))

  eta0 <- final_model$intercept0 + final_model$offset
  mu0  <- final_model$family$linkinv(eta0)
  final_model$null.deviance <- sum(
    final_model$family$dev.resids(final_model$y, mu0, final_model$prior.weights))

  # Dispersion: Fletcher's (2012) estimator with n = sum(prior weights) and
  # edf~ as the model degrees of freedom (0 when EDF is not tracked)
  # for gaussian, Gamma and inverse.gaussian; fixed at 1 for poisson and
  # binomial. edf~ ignores the search degrees of freedom, so
  # the scale and the quantities built on it are approximate.
  disp_mod <- final_model
  if (!edf_tracked) disp_mod$edf.total <- 0
  final_model$scale <- estimate_dispersion(disp_mod)
  # sig2 is the mgcv alias for the estimated scale parameter; set it for all
  # families whose dispersion is free (gaussian, Gamma, inverse.gaussian).
  # Poisson and binomial have phi fixed at 1 and do not need sig2.
  if (final_model$family$family %in%
        c("gaussian", "Gamma", "inverse.gaussian")) {
    final_model$sig2 <- final_model$scale
  }

  final_model$iter      <- final_model$mstop %||% NA_integer_
  final_model$optimizer <- "boosting"
  final_model$method    <- "boosting"
  # Converged: loss-based stopping found no improvement within its patience
  # (not a fit ended by max_iter or the time limit)
  final_model$converged <-
    isTRUE(final_model$metrics$stopping_conditions["no_improvements"])

  # var.summary: 3-quantile (Q1, median, Q3) summary of every predictor
  # variable in the model frame, and the first value of every factor
  # (character predictors included, as they are stored as factors below).
  # Matrix predictors (functional terms) are flattened before summarising.
  resp_nm <- deparse(formula(final_model)[[2L]])
  pred_nms <- setdiff(names(final_model$model), resp_nm)
  # Character predictors are stored as factors (with the levels in
  # $xlevels), as model.matrix() reads them, so that methods reading the
  # model frame (e.g. termplot() in plot(all.terms = TRUE)) treat them alike
  for (nm in pred_nms) {
    v <- final_model$model[[nm]]
    if (is.character(v)) {
      final_model$model[[nm]] <-
        factor(v, levels = final_model$xlevels[[nm]] %||% sort(unique(v)))
    }
  }
  final_model$var.summary <- setNames(
    lapply(pred_nms, function(nm) {
      v <- final_model$model[[nm]]
      if (is.factor(v)) v[1L]
      else stats::quantile(as.vector(v), c(0.25, 0.5, 0.75), na.rm = TRUE)
    }),
    pred_nms
  )

  class(final_model) <- c("gbgam", "gam", "glm", "lm")
  final_model
}


#' Weighted column means of the model matrix
#'
#' Computes `$cmX`, which `mgcv::plot.gam(seWithMean = TRUE)` and
#' `predict.gam(type = "iterms")` use to add the uncertainty of the overall
#' mean to the term-wise standard errors of centred smooths. The means are
#' weighted by the prior weights, matching the centring of the learners, so
#' they are exactly zero for the columns of `"gbgam_centred"` smooths.
#'
#' @param mod A model list with `$pterms`, `$model`, `$prior.weights` and
#'   `$smooth`.
#' @return A numeric vector with one entry per coefficient.
#' @noRd
model_col_means <- function(mod) {
  w      <- mod$prior.weights
  w_mean <- function(X) colSums(X * w) / sum(w)
  pX     <- stats::model.matrix(mod$pterms, mod$model)
  sm_means <- lapply(mod$smooth, function(sm) {
    if (inherits(sm, "gbgam_centred")) rep(0, ncol(sm$X)) else w_mean(sm$X)
  })
  unname(c(w_mean(pX), unlist(sm_means)))
}
