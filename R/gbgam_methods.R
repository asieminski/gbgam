# Methods -----------------------------------------------------------------

#' Extract fitted values from a `gbgam` model
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#' @return A numeric vector of fitted values on the response scale
#'   (i.e. after applying the inverse link function). Under
#'   `na.action = na.exclude` it is padded with `NA` for the excluded rows.
#' @examples
#' set.seed(1)
#' n <- 100
#' df <- data.frame(x = runif(n))
#' eta <- 0.5 + sin(2 * pi * df$x)
#' df$y <- rpois(n, lambda = exp(eta))
#' mod <- gbgam(
#'   y ~ t2(x, k = 6, bs = "cr"),
#'   data = df,
#'   family = poisson(),
#'   stop_control = gb_stop_ctrl(stopping_method = "none", max_iter = 20),
#'   edf_control  = gb_edf_ctrl(edf_compute_method = "none")
#' )
#' head(fitted(mod))
#' @export
fitted.gbgam <- function(object, ...) {
  stats::napredict(object$na.action, fitted_mu(object))
}

#' Fitted means of the rows used in the fit (never NA-padded)
#'
#' Internal counterpart of [fitted.gbgam()] for use during and after fitting,
#' where vectors must have one entry per row of the model frame.
#'
#' @param mod A model list or fitted `"gbgam"` object.
#' @return Numeric vector `linkinv(mod$linear.predictors)`.
#' @noRd
fitted_mu <- function(mod) {
  as.vector(mod$family$linkinv(mod$linear.predictors))
}


#' Summarise a fitted `gbgam` model
#'
#' Returns a `"summary.gbgam"` object containing model-level stopping
#' information and learner-level metrics. The object is printed via
#' [print.summary.gbgam()], so assigning the result to a variable does not
#' produce any output. Unlike [mgcv::summary.gam()], it reports no p-values
#' or approximate significance tests for the terms.
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#'
#' @returns An object of class `"summary.gbgam"` (returned visibly, so
#'   `summary(object)` prints it). It is a list with:
#'   \describe{
#'     \item{`family`}{family name and link}
#'     \item{`formula`}{the model formula}
#'     \item{`mstop`}{selected iteration}
#'     \item{`stop_msgs`}{character vector of stopping reason(s)}
#'     \item{`step_size`}{boosting step size}
#'     \item{`intercept0`}{MLE of the intercept-only model (the boosting
#'       start when `initialise_intercept = TRUE`)}
#'     \item{`intercept`}{fitted intercept coefficient}
#'     \item{`edf`}{the model's \eqn{\widetilde{\operatorname{edf}}}{edf~},
#'       or `NULL` when EDF is not tracked
#'       (`edf_compute_method = "none"`)}
#'     \item{`edf_label`}{EDF method annotation string}
#'     \item{`loss_name`}{name of the tracked loss (`"AIC"`, `"MSE"`, etc.)}
#'     \item{`best_loss`}{minimum tracked loss value, or `NULL`}
#'     \item{`learner_table`}{data frame with one row per learner: its
#'       label, whether it is a penalised (Demmler-Reinsch) learner, the
#'       number (`sel`) and percentage (`%sel`) of selections up to `mstop`,
#'       with loss-based stopping the signed share of the change of the
#'       stopping loss up to `mstop` over the iterations that selected it
#'       (`%d_stop_loss`, in percent of the total change; negative when its
#'       selections raised the loss; forced post-iteration and padded
#'       intercept updates count towards the selected learner), and with
#'       EDF tracking its accumulated
#'       \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} (`acc_edf`)}
#'   }
#' @export
summary.gbgam <- function(object, ...) {
  loss_name <- stop_loss_name(object)

  # Model-level metrics -----------------------------------------------------
  stop_msgs <- stop_reasons(object)

  # EDF
  edf <- if (object$edf_compute_method != "none") {
    round(object$edf.total, 2)
  }
  edf_label <- edf_method_label(object)
  # Loss
  best_loss <- if (object$stopping_method != "none")
    min(object$metrics$ls, na.rm = TRUE) else NULL

  # Learner table -----------------------------------------------------------
  learners   <- object$learners
  n_learners <- length(learners)
  m_stop     <- object$mstop
  k_opts     <- object$metrics$k_opts[seq_len(m_stop)]
  sel_count  <- tabulate(k_opts, nbins = n_learners)
  tot_sel    <- max(1L, sum(sel_count))  # mstop = 0: avoid 0/0
  is_dr      <- vapply(learners, inherits, logical(1L), "demmler_reinsch")

  df <- data.frame(
    learner          = vapply(learners, function(lr) lr$var_label,
                              character(1L)),
    is_dr            = is_dr,
    sel              = sel_count,
    `%sel`           = round(100 * sel_count / tot_sel, 1),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  if (object$stopping_method != "none") {
    # metrics$ls starts at iteration 0
    dloss       <- diff(object$metrics$ls)[seq_len(m_stop)]
    loss_change <- vapply(seq_len(n_learners), function(k) {
      idx <- which(k_opts == k)
      if (length(idx) > 0L) sum(dloss[idx], na.rm = TRUE) else 0
    }, numeric(1L))
    tot_change <- sum(loss_change)
    # Signed shares (+ 0 turns -0 into 0 for printing)
    df[["%d_stop_loss"]] <- if (abs(tot_change) > 1e-10)
      round(100 * loss_change / tot_change, 1) + 0 else 0
  }

  if (object$edf_compute_method != "none") {
    df[["acc_edf"]] <- round(
      vapply(learners, function(lr) lr$edf %||% NA_real_, numeric(1L)), 2)
  }

  structure(
    list(
      family      = paste0(object$family$family, "/", object$family$link),
      formula     = object$formula,
      mstop       = object$mstop,
      stop_msgs   = stop_msgs,
      step_size   = object$step_size,
      intercept0  = object$intercept0,
      intercept   = unname(object$coefficients["(Intercept)"]),
      edf         = edf,
      edf_label   = edf_label,
      loss_name   = loss_name,
      best_loss   = best_loss,
      learner_table = df
    ),
    class = "summary.gbgam"
  )
}


#' Print a `summary.gbgam` object
#'
#' @param x A `"summary.gbgam"` object returned by [summary.gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#'
#' @returns `x`, invisibly.
#' @export
print.summary.gbgam <- function(x, ...) {
  df <- x$learner_table

  # Header ------------------------------------------------------------------
  cli::cli_rule(left = paste0("gbgam  ", x$family))
  # One string: cli would join a multi-line deparse() with ", and "
  fml_str <- paste(deparse(x$formula, width.cutoff = 500L), collapse = " ")
  cli::cli_text("{.field Formula:} {fml_str}")
  # Fit summary -------------------------------------------------------------
  stop_msg <- paste(x$stop_msgs, collapse = "; ")
  cli::cli_text("{.field mstop}      {x$mstop}  ({stop_msg})")
  cli::cli_text("{.field step size}  {x$step_size}")
  # Fitted intercept, and the intercept-only MLE (the boosting start when
  # initialise_intercept = TRUE)
  int_fit <- round(x[["intercept"]] %||% NA_real_, 4)
  int_mle <- round(x$intercept0, 4)
  cli::cli_text(
    "{.field intercept}  {int_fit}  (intercept-only MLE: {int_mle})")
  if (!is.null(x$edf))
    cli::cli_text("{.field edf}        {x$edf}{x$edf_label}")
  if (!is.null(x$best_loss)) {
    best <- format(x$best_loss, digits = 4)
    cli::cli_text(
      "{.field best {x$loss_name}}   {best}  (m = {x$mstop})")
  }

  # Learner table -----------------------------------------------------------
  cli::cli_rule(left = "Learner metrics")

  is_dr    <- df$is_dr
  loss_col <- if ("%d_stop_loss" %in% names(df)) "%d_stop_loss" else NULL
  edf_col  <- if ("acc_edf" %in% names(df)) "acc_edf" else NULL

  # The learner column fits the longest label, but leaves the other columns
  # within the console width (at least 24 characters); longer labels are
  # shortened, keeping their learner suffix (short_learner_name())
  other_w <- 2L + 7L + 8L + (if (!is.null(loss_col)) 14L else 0L) +
    (if (!is.null(edf_col)) 9L else 0L)
  col_w   <- min(max(24L, nchar(df$learner)),
                 max(24L, getOption("width", 80L) - other_w))
  labels  <- short_learner_name(df$learner, col_w)

  hdr_parts <- sprintf("  %-*s  %5s  %6s", col_w, "learner", "sel", "%sel")
  if (!is.null(loss_col))
    hdr_parts <- paste0(hdr_parts, sprintf("  %12s", loss_col))
  if (!is.null(edf_col))
    hdr_parts <- paste0(hdr_parts, sprintf("  %7s", edf_col))
  cli::cli_verbatim(hdr_parts)

  for (k in seq_len(nrow(df))) {
    nm_plain <- labels[k]
    pad      <- max(0L, col_w - nchar(nm_plain))
    nm_str   <- if (is_dr[k])
                  paste0(cli::style_bold(nm_plain), strrep(" ", pad))
                else
                  sprintf("%-*s", col_w, nm_plain)
    row <- sprintf("  %s  %5d  %5.1f%%", nm_str, df$sel[k], df$`%sel`[k])
    if (!is.null(loss_col))
      row <- paste0(row, sprintf("  %11.1f%%", df[[loss_col]][k]))
    if (!is.null(edf_col))
      row <- paste0(row, sprintf("  %7.2f", df[[edf_col]][k]))
    cli::cli_verbatim(row)
  }

  invisible(x)
}


#' Print a `gbgam` object
#'
#' @param x A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#'
#' @returns
#' `x`, invisibly, after printing model information to the console.
#' @export
print.gbgam <- function(x, ...) {
  print(x$family)
  cat("Formula:\n")
  print(x$formula)

  edf <- if (x$edf_compute_method == "none") {
    "not tracked (edf_compute_method = \"none\")"
  } else {
    paste0(round(x$edf.total, 2), edf_method_label(x))
  }
  cat("Estimated model degrees of freedom = ", edf, "\n", sep = "")

  reasons <- stop_reasons(x)
  metric  <- stop_loss_name(x)
  cat(
    "mstop=", x$mstop,
    if (length(reasons) > 0L)
      paste0(" (", paste(reasons, collapse = "; "), ")"),
    if (x$stopping_method != "none")
      paste0(", best ", metric, "=",
             format(min(x$metrics$ls, na.rm = TRUE), digits = 4)),
    "\n",
    sep = ""
  )

  invisible(x)
}

#' Label of the EDF method, for printing
#'
#' @param x A fitted `gbgam` object.
#' @return `" (Hutch++, z=<z>)"` for Hutch++ EDF tracking, else `""`.
#' @noRd
edf_method_label <- function(x) {
  if (identical(x$edf_compute_method, "hutchpp"))
    paste0(" (Hutch++, z=", x$z, ")") else ""
}

#' Why a boosting run ended
#'
#' The end of the run and the selected iteration are reported separately:
#' with loss-based stopping, `mstop` can be earlier than the iteration at
#' which `max_iter` or the time limit ended the run (e.g. after
#' [extend_gbgam()]). The time limit reported is that of the run that ended
#' the path: `x$metrics$max_runtime_mins` after an extension stopped by its
#' own time limit, otherwise `x$max_runtime_mins`.
#'
#' @param x A fitted `gbgam` object.
#' @return A character vector with one message per stopping condition met.
#' @noRd
stop_reasons <- function(x) {
  flags <- x$metrics$stopping_conditions
  m_run <- sum(!is.na(x$metrics$k_opts))
  ended <- if (x$mstop < m_run) paste0("run ended at iteration ", m_run, ": ")
           else ""
  limit <- x$metrics[["max_runtime_mins"]] %||% x$max_runtime_mins
  c(
    if (isTRUE(flags["no_improvements"]))
      paste0("no ", stop_loss_name(x), " improvement after ",
             m_run - x$mstop, " iterations"),
    if (isTRUE(flags["mstop_reached"]))
      paste0(ended, "max_iter = ", x$max_iter, " reached"),
    if (isTRUE(flags["time_limit_exceeded"]))
      paste0(ended, "time limit of ", limit, " mins reached")
  )
}

#' Plot a `gbgam` object
#'
#' Calls [mgcv::plot.gam()]. A model with the covariance `x$Vp`, from
#' [bootstrap_gbgam()] or [bayesian_backfitting()], is plotted with intervals
#' without any further step.
#'
#' @param x A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Additional arguments passed on to [mgcv::plot.gam()].
#'
#' @returns
#' The result of [mgcv::plot.gam()]. If `se` is not supplied, it is set to
#' `TRUE` when `x$Vp` is present and otherwise `FALSE`. Requesting `se = TRUE`
#' without `x$Vp` gives a warning and plots without standard errors. With
#' `residuals = TRUE`, the partial residuals are those of the rows used in
#' the fit, also under `na.action = na.exclude`.
#'
#' @details
#' The defaults are those of [mgcv::plot.gam()]: plus or minus 2 standard
#' errors that leave out the uncertainty of the intercept. For intervals of
#' level 1 - `alpha`, with coverage across the function (see
#' [smooth_intervals()]), set `se = qnorm(1 - alpha / 2)` and
#' `seWithMean = TRUE`. `shift = coef(x)[1]` adds the intercept to the curves.
#' The covariance of a bootstrap fit is widened for the bootstrap's own
#' `alpha`, so use that `alpha` here; [smooth_intervals()] gives its intervals
#' at other levels.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x   <- runif(150)
#' dat <- data.frame(x = x, y = sin(2 * pi * x) + rnorm(150, sd = 0.3))
#' mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = dat)
#' plot(mod)
#'
#' # The covariance is all that the intervals need
#' mod <- bootstrap_gbgam(mod, n_samples = 50)
#' plot(mod)
#'
#' # 95% intervals, widened for this level by the bootstrap's default
#' # alpha = 0.05, including the uncertainty of the intercept
#' plot(mod, seWithMean = TRUE, se = qnorm(0.975))
#' }
#' @export
plot.gbgam <- function(x, ...) {
  args <- list(...)
  # exact matching: args$se would partially match e.g. `select`
  se <- args[["se", exact = TRUE]]
  if (is.null(se)) {
    args$se <- !is.null(x$Vp)
  } else if (isTRUE(as.logical(se[1L])) && is.null(x$Vp)) {
    warning(
      "No covariance matrix available (x$Vp is NULL); plotting without ",
      "standard errors. Run bootstrap_gbgam() or bayesian_backfitting() ",
      "first.",
      call. = FALSE
    )
    args$se <- FALSE
  }
  # Without na.action, predict() is not NA-padded, so that the partial
  # residuals line up with the unpadded x$residuals and x$weights
  x$na.action <- NULL
  do.call(mgcv::plot.gam, c(list(x = x), args))
}


# Additional S3 methods ---------------------------------------------------

#' Extract residuals from a `gbgam` model
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param type Character. Type of residuals, defined as in [stats::glm()]
#'   and [mgcv::gam()] (\eqn{w} are the prior weights). One of:
#'   \describe{
#'     \item{`"deviance"`}{Signed square-root deviance residuals.}
#'     \item{`"response"`}{Raw response residuals \eqn{y - \hat{\mu}}.}
#'     \item{`"working"`}{Working residuals \eqn{(y-\hat\mu)/\mu'(\hat\eta)},
#'       as used by IRLS (stored in `object$residuals`).}
#'     \item{`"pearson"`}{Pearson residuals
#'       \eqn{(y-\hat\mu)\sqrt{w}/\sqrt{V(\hat\mu)}}.}
#'     \item{`"scaled.pearson"`}{Pearson residuals divided by
#'       \eqn{\sqrt{\hat\phi}}, the estimated scale (`object$sig2`; 1 for
#'       Poisson and binomial).}
#'   }
#'   Default `"deviance"`.
#' @param ... Ignored; present for S3 compatibility.
#' @return A numeric vector of residuals of the requested type. Under
#'   `na.action = na.exclude` it is padded with `NA` for the excluded rows.
#' @export
residuals.gbgam <- function(object, type = "deviance", ...) {
  checkmate::assert_choice(type, c("response", "working", "deviance",
                                   "pearson", "scaled.pearson"))

  mu     <- object$fitted.values
  r_resp <- object$y - mu
  res <- switch(type,
    response = r_resp,
    working  = object$residuals,
    deviance = {
      dr <- object$family$dev.resids(object$y, mu, object$prior.weights)
      sign(r_resp) * sqrt(pmax(dr, 0))
    },
    pearson = {
      r_resp * sqrt(object$prior.weights) / sqrt(object$family$variance(mu))
    },
    scaled.pearson = {
      r_resp * sqrt(object$prior.weights) /
        sqrt((object$sig2 %||% 1) * object$family$variance(mu))
    }
  )
  stats::naresid(object$na.action, res)
}


#' Variance-covariance matrix of a `gbgam` model
#'
#' Returns the covariance matrix `$Vp` of the coefficients: the (widened)
#' bootstrap covariance from [bootstrap_gbgam()] or
#' `gbgam(uq_control = gb_uq_ctrl("bootstrap"))`, or the
#' posterior covariance from [bayesian_backfitting()] or
#' `gb_uq_ctrl("bayes")`. Unlike `mgcv`'s `Vp`, it is not the Bayesian
#' posterior covariance at the smoothing parameters of the fit: it is either
#' a bootstrap covariance, widened for smoothing bias, or a backfitting
#' posterior conditional on the smoothing parameters chosen in Stage 1 of
#' [bayesian_backfitting()]. [predict(se.fit = TRUE)][predict.gbgam()],
#' [plot()][plot.gbgam()] and [smooth_intervals()] use it.
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#' @return A \eqn{p \times p} numeric matrix, or `NULL` with a warning when
#'   the model has no `$Vp`.
#' @export
vcov.gbgam <- function(object, ...) {
  if (!is.null(object$Vp)) return(object$Vp)
  warning(no_vp_message(), call. = FALSE)
  NULL
}

#' Message for a model without a covariance matrix
#'
#' @return A character string naming the ways to attach `$Vp`.
#' @noRd
no_vp_message <- function() {
  paste0(
    "No covariance matrix available. Fit with uq_control = ",
    "gb_uq_ctrl(\"bootstrap\") or gb_uq_ctrl(\"bayes\"), ",
    "or call bootstrap_gbgam() or bayesian_backfitting() on the fitted ",
    "model."
  )
}


#' Predictions from a `gbgam` model
#'
#' Calls [mgcv::predict.gam()]. Standard errors (`se.fit = TRUE`) use the
#' covariance matrix `object$Vp` (see [vcov.gbgam()]); requesting them from
#' a model without `$Vp` results in an error.
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Arguments passed on to [mgcv::predict.gam()], such as
#'   `newdata`, `type` and `se.fit`.
#' @return The result of [mgcv::predict.gam()].
#' @export
predict.gbgam <- function(object, ...) {
  if (is.null(object$Vp) && se_fit_requested(...)) {
    stop(no_vp_message(), call. = FALSE)
  }
  mgcv::predict.gam(object, ...)
}

#' Whether a call of `mgcv::predict.gam()` requests standard errors
#'
#' Matches `...` to the arguments of [mgcv::predict.gam()] (by name, partial
#' name or position) and returns its `se.fit`.
#'
#' @param ... Arguments of [mgcv::predict.gam()] after `object`.
#' @return `TRUE` or `FALSE`.
#' @noRd
se_fit_requested <- function(...) {
  # A call predict.gam(object, ...) with the evaluated arguments, matched
  # to the formals of predict.gam() as R would match them
  cl <- match.call(mgcv::predict.gam,
                   as.call(c(quote(predict.gam), quote(object), list(...))))
  isTRUE(as.logical(cl$se.fit)[1L])
}


#' Log-likelihood of a fitted `gbgam` model
#'
#' Computes the full log-likelihood at the fitted values. It is the same
#' likelihood as in the AIC stopping criterion: for the gaussian, Gamma and
#' inverse.gaussian families the scale is estimated by Fletcher's (2012)
#' method, and for binomial the Bernoulli kernel weighted by the prior
#' weights is used. The result is used by [AIC()], [BIC()], and
#' likelihood-ratio tests.
#'
#' The `df` attribute is the model's
#' \eqn{\widetilde{\operatorname{edf}}}{edf~} (`object$edf.total`),
#' plus one for the estimated scale parameter of the gaussian, Gamma and
#' inverse.gaussian families (as in [stats::glm()] and [mgcv::gam()]). For
#' these families `AIC(object)` therefore exceeds the stopping AIC (with
#' `AIC_k = 2`) by 2. If the model was fitted with
#' `edf_compute_method = "none"`, the EDF is not tracked: `df` is `NA` (so
#' [AIC()] is `NA`) and the scale is estimated with 0 model degrees of
#' freedom.
#'
#' @references
#' Fletcher, D. J. (2012).
#' Estimating overdispersion when fitting a generalized linear model to sparse
#' data.
#' *Biometrika, 99*(1), 230--237.
#' \href{https://doi.org/10.1093/biomet/asr083}{doi:10.1093/biomet/asr083}
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#' @return A `"logLik"` object with attributes `df` (estimated degrees of
#'   freedom consumed) and `nobs`.
#' @export
logLik.gbgam <- function(object, ...) {
  # Same -2 log L as the AIC stopping rule (neg2loglik()); the dispersion
  # comes from estimate_dispersion() for the free-scale families.
  scale_est   <- object$family$family %in%
    c("gaussian", "Gamma", "inverse.gaussian")
  edf_tracked <- object$edf_compute_method != "none"
  # Without EDF tracking the dispersion uses 0 model degrees of freedom, as
  # in `$scale`
  if (!edf_tracked) object$edf.total <- 0
  val         <- -neg2loglik(object) / 2
  attr(val, "df")   <- if (edf_tracked) {
    object$edf.total + as.numeric(scale_est)
  } else {
    NA_real_
  }
  attr(val, "nobs") <- nobs(object)
  class(val) <- "logLik"
  val
}


#' Number of observations used to fit a `gbgam` model
#'
#' The number of rows with a non-zero prior weight, as for [stats::glm()]
#' and [mgcv::gam()] models: rows with weight zero (e.g. binomial rows with
#' no trials) do not enter the fit. It is the `nobs` attribute of
#' [logLik()][logLik.gbgam()], and so the sample size of [BIC()].
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#' @return Integer scalar.
#' @export
nobs.gbgam <- function(object, ...) {
  w <- object$prior.weights
  if (is.null(w)) nrow(object$model) else sum(w != 0)
}


#' Prior or working weights from a fitted `gbgam` model
#'
#' With `type = "prior"`, returns the prior weights used during boosting.
#' Following [stats::glm], for a `cbind(successes,
#' failures)` binomial response these are
#' `n_trials * user_supplied_weights` (where
#' `n_trials = successes + failures`); for all other response types they
#' equal the user-supplied `weights` argument (defaulting to a vector of
#' ones). The user-supplied weights are stored in `object$user_weights`.
#'
#' With `type = "working"`, returns the working weights
#' \eqn{w \mu'(\hat\eta)^2 / V(\hat\mu)} at the fitted values (stored in
#' `object$weights`, as in [stats::glm()] and [mgcv::gam()]).
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param type Character. `"prior"` or `"working"`. Default `"prior"`.
#' @param ... Ignored; present for S3 compatibility.
#' @return A numeric vector of weights. Under `na.action = na.exclude` it is
#'   padded with `NA` for the excluded rows.
#' @export
weights.gbgam <- function(object, type = c("prior", "working"), ...) {
  type <- match.arg(type)
  res  <- if (type == "prior") object$prior.weights else object$weights
  stats::naresid(object$na.action, res)
}


# Formula and update methods ----------------------------------------------

#' Extract the formula from a fitted `gbgam` model
#'
#' @param x A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param ... Ignored; present for S3 compatibility.
#' @return The formula used to fit the model.
#' @export
formula.gbgam <- function(x, ...) x$formula


#' Update a fitted `gbgam` model
#'
#' Refits the model, with a new formula or other arguments of [gbgam()] if
#' given, as [stats::update.default()] does. Unless `uq_control` is
#' given, the refit runs the uncertainty quantification recorded in the
#' model (`$uq_control`, the one that attached its `$Vp`).
#'
#' @param object A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param formula. Formula. A new formula or a one-sided update expression
#'   such as `~ . + t2(z, k = 6)`. If omitted, the original formula is used.
#' @param ... Named arguments of [gbgam()] to change or add in the refitting
#'   call (e.g. `stop_control`, `learner_control`, `family`).
#' @param evaluate Logical. Evaluate the updated call and return the
#'   refitted model; `FALSE` returns the unevaluated call. Default `TRUE`.
#'
#' @return A refitted `"gbgam"` object (when `evaluate = TRUE`) or an
#'   unevaluated call (when `evaluate = FALSE`).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x   <- runif(80)
#' dat <- data.frame(x = x, y = sin(2 * pi * x) + rnorm(80, sd = 0.3))
#' mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = dat,
#'              stop_control = gb_stop_ctrl(stopping_method = "none",
#'                                          max_iter = 30),
#'              edf_control  = gb_edf_ctrl(edf_compute_method = "none"))
#'
#' # More iterations
#' mod2 <- update(mod, stop_control = gb_stop_ctrl(stopping_method = "none",
#'                                                 max_iter = 60))
#' }
#' @export
update.gbgam <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  extras <- as.list(match.call(expand.dots = FALSE))[["..."]]

  call[["formula"]] <- if (!missing(formula.)) {
    stats::update.formula(object$formula, formula.)
  } else {
    object$formula
  }

  # Without a data argument, the model frame can stand in for the data only
  # when it has the same rows; after NA removal or a subset it has fewer, so
  # the call is re-evaluated in the calling environment instead (as
  # stats::update.default() does)
  if (is.null(call[["data"]]) && is.null(object$na.action) &&
      is.null(call[["subset"]])) {
    call[["data"]] <- object$model
  }

  if (is.null(extras$uq_control) &&
      !is.null(object$uq_control)) {
    call[["uq_control"]] <- object$uq_control
  }

  if (length(extras) > 0L) {
    existing <- names(extras) %in% names(call)
    for (nm in names(extras)[existing]) call[[nm]] <- extras[[nm]]
    if (any(!existing)) call <- as.call(c(as.list(call), extras[!existing]))
  }

  if (evaluate) eval(call, parent.frame()) else call
}
