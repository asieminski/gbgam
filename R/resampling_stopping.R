# Stopping by resampling estimates of the prediction error
# (gb_stop_ctrl(stopping_method = "cv")).
#
# Vocabulary used by every name in this file:
#   main fit           the boosting fit to all rows; it is the fit returned.
#   resample fits      B boosting fits, one per resample, run alongside the
#                      main fit to the same iteration.
#   in-bag weights     n x B matrix: how often row i is in resample b. They
#                      multiply the prior weights of resample fit b; 0 means
#                      row i is held out of fit b (for k-fold: it is in fold b).
#   evaluation weights n x B matrix: the out-of-bag error is the sum over
#                      rows i and fits b of evaluation weight times the
#                      deviance of row i under fit b. By default they average
#                      each row's deviance over the fits that hold it out, and
#                      then average over rows by prior weight (Efron &
#                      Tibshirani, 1997).
#   apparent error     the main fit's training deviance per unit prior weight.
#   out-of-bag error   Err(1) of Efron and Tibshirani (1997); the
#                      cross-validation error for k-fold.
#   no-information     the deviance expected if responses and fitted means
#   error              were paired at random (gamma; .632+ only).
#   resampling plan    everything fixed before boosting starts (see
#                      plan_resampling()).
#
# The algorithm (fit_gbgam_cv()):
#   1. plan_resampling():       draw the in-bag weights, derive the evaluation
#                               weights, record the RNG state.
#   2. build_resample_fits():   set up the B resample fits at iteration 0.
#   3. fit_state_init():        set up the main fit, recording its history.
#   4. boost_with_resampling(): repeat {advance the main fit by `chunk`
#                               iterations; advance the resample fits to the
#                               same iteration; estimate the prediction error
#                               at every iteration so far; stop when the
#                               patience, max_iter or the time limit is
#                               reached}.
#   5. resampling_result():     rewind the main fit to the iteration with the
#                               lowest estimate and return it.
# Sections 1-5 below hold the functions of steps 1-5, in this order.
# extend_gbgam_cv() continues steps 4-5 for extend_gbgam(). The estimators
# themselves are in resampling_estimators.R.


#' Control parameters for stopping by resampling
#'
#' @description
#' Settings for stopping by resampling (`gb_stop_ctrl("cv")`): the scheme,
#' the number of resampling fits, the steps each fit redoes on its own, and
#' the estimator. The main fit and the resampling fits run side by side, and
#' the stopping rule of [gb_stop_ctrl()] is checked every `chunk`
#' iterations. The returned model is the main fit at the iteration with the
#' lowest estimate of the prediction error, measured by the family's
#' deviance. Supplying an argument that the chosen scheme does not use
#' results in an error.
#'
#' @param type Character. The resampling scheme: `"kfold"`, k-fold
#'   cross-validation; `"bootstrap"`, the nonparametric bootstrap; or
#'   `"subsample"`, samples of a `fraction` of the observations drawn
#'   without replacement. Supplying it together with `weights` results in an
#'   error. Default `"kfold"`.
#' @param B Integer, at least 2, or NULL. The number of resampling fits
#'   (folds for `"kfold"`). Supplying it together with `folds` or `weights`
#'   results in an error. With `"kfold"`, `B` must not exceed the number of
#'   observations. Default `NULL`: `10` for `"kfold"` and `25` otherwise.
#' @param mode Character. Which steps of the fit each resampling fit redoes
#'   on its own; the steps it does not redo are copied from the main fit.
#'   One of:
#'   \describe{
#'     \item{`"free_selection"`}{each fit selects its own learners and
#'       adapts its own smoothing parameters by adaptive smoothness
#'       relaxation (ASR).}
#'     \item{`"free_smoothing"`}{replays the main fit's learner selection;
#'       each fit adapts its own smoothing parameters by ASR.}
#'     \item{`"fixed"`}{replays the main fit's learner selection and
#'       smoothing parameters. It needs no tracking of the effective
#'       degrees of freedom (EDF; see [gb_edf_ctrl()]), so it is the
#'       fastest.}
#'   }
#'   When the main fit does not use ASR (no EDF tracking, or `tau = 0`),
#'   the smoothing parameters stay at their starting values. Default
#'   `"free_selection"`.
#' @param estimator Character, or NULL. How the deviances of the resampling
#'   fits are combined into an estimate (see Details): `"oob"`, the
#'   out-of-sample estimate (the cross-validation error for `"kfold"`);
#'   `"632"`, the .632 estimator; or `"632plus"`, the .632+ estimator. The
#'   .632 estimators are only available for `"bootstrap"` and for custom
#'   `weights`. Default `NULL`: `"632plus"` for `"bootstrap"`, `"oob"`
#'   otherwise.
#' @param folds Vector, or NULL. Only with `type = "kfold"`. Fold labels,
#'   one per observation used in the fit (after `subset` and `na.action`),
#'   e.g. blocks of consecutive days for time series. Default `NULL`: random
#'   folds of (nearly) equal size.
#' @param weights Matrix, or NULL. Custom training weights: a non-negative
#'   matrix with one row per observation used in the fit and one column per
#'   resampling fit, which multiplies the prior weights (see Details).
#'   Default `NULL`: the weights of `type`.
#' @param eval_weights Matrix, or NULL. Only with `weights`. Evaluation
#'   weights of the same size as `weights` (see Details). Default `NULL`:
#'   the out-of-sample averaging found from the zeros of `weights`.
#' @param fraction Number in (0, 1), or NULL. Only with
#'   `type = "subsample"`. Fraction of the observations in each subsample.
#'   Default `NULL`: `0.5`.
#' @param chunk Positive integer. Number of iterations the main fit and the
#'   resampling fits run between two checks of the stopping rule. Default
#'   `100`.
#' @param keep_fits Logical. If `TRUE` and the fit is resumable
#'   (`gb_stop_ctrl(resumable = TRUE)`), the resampling fits are stored in
#'   the model (`$resume$cv$fits`), so that [extend_gbgam()] continues them
#'   instead of refitting them from the start. They take about `B` times
#'   the memory of the main fit's learners. Default `FALSE`.
#'
#' @return A list of class `"gb_resampling"`.
#'
#' @details
#' ## The resampling fits
#'
#' Each resampling fit is a boosting fit to all observations, with the prior
#' weights multiplied by its training weights: 1 for the observations in
#' its training part and 0 for the others with `"kfold"` and
#' `"subsample"`, and the number of times each observation is drawn with
#' `"bootstrap"`. This is the weight-matrix design of `mboost`'s `cv()`.
#' A resampling fit uses the bases of the main fit's smooths, centred and
#' reparametrised again for its own weights. An observation with training
#' weight 0 is left out of the fit.
#'
#' ## The estimates
#'
#' All estimates are deviances per unit of prior weight. The out-of-sample
#' estimate (`"oob"`; Efron & Tibshirani, 1997) averages each observation's deviance
#' over the fits that leave it out, and then averages these over the
#' observations, weighted by their prior weights; the observations that no
#' fit leaves out do not enter it. For k-fold cross-validation, which leaves
#' each observation out once, this is the usual cross-validation error. It
#' is the estimate that `mboost::cvrisk()` uses for every scheme.
#'
#' The .632 estimator puts weight 0.368 on the apparent error (the main
#' fit's training deviance) and 0.632 on the out-of-sample estimate. The
#' .632+ estimator (Efron & Tibshirani, 1997) moves weight towards the
#' out-of-sample estimate when the fit overfits, that is, when the
#' out-of-sample estimate is well above the apparent error relative to the
#' no-information error. The no-information error is the mean deviance over
#' all pairs of a response and a fitted mean (Efron & Tibshirani, 1997, eq.
#' 26). With an offset (e.g. an exposure), each
#' fitted mean in a pair keeps its own observation's offset.
#'
#' ## Custom weights
#'
#' Other schemes, such as blocked bootstrap samples for dependent data, can
#' be given as `weights`. The estimate is then the out-of-sample estimate
#' found from the zeros of `weights`, or, with `eval_weights`, the sum over
#' all fits and observations of each observation's deviance times its
#' evaluation weight; the .632 estimators can be added to either. A column
#' of `weights` without zeros leaves no observation out, so without
#' `eval_weights` its fit does not enter the estimate (with a warning).
#'
#' @section The fitted model's $cv:
#' A model fitted with `stopping_method = "cv"` holds the resampling
#' results in `$cv`, a list with the fields below. Its paths have one
#' element (or row) per iteration, starting at iteration 0.
#' \describe{
#'   \item{`type`, `estimator`, `B`}{The resampling scheme (`"custom"` for
#'     `weights`), the estimator and the number of resampling fits.}
#'   \item{`mstop`}{The iteration with the lowest `risk`; equal to
#'     `mod$mstop`.}
#'   \item{`risk`}{The estimate of the prediction error, which is the
#'     stopping loss (also stored as `mod$metrics$ls`).}
#'   \item{`apparent`}{The apparent (training) error of the main fit.}
#'   \item{`oob`}{The out-of-sample estimate (the cross-validation error
#'     for `"kfold"`): the row sums of `rep_risk`.}
#'   \item{`noinfo`}{The no-information error (`"632plus"` only, else
#'     `NULL`).}
#'   \item{`weight`}{The weight 0.632 / (1 - 0.368 R) of Efron and
#'     Tibshirani (1997, eq. 29), with R the relative overfitting rate:
#'     `oob` truncated at `noinfo`, minus `apparent`, over `noinfo` minus
#'     `apparent`, and 0 when the truncated `oob` is at most `apparent`
#'     (`"632plus"` only, else `NULL`). Where `oob` is at most `noinfo`,
#'     `risk` is (1 - `weight`) `apparent` + `weight` `oob`; where `oob` is
#'     above `noinfo`, the correction uses `noinfo` in place of `oob`, so
#'     `risk` does not have this form.}
#'   \item{`rep_risk`}{A matrix with one column per resampling fit: its
#'     share of the out-of-sample estimate, that is, the sum of its
#'     deviances times their evaluation weights.}
#'   \item{`weights`}{The training weights, one row per observation and one
#'     column per resampling fit.}
#'   \item{`n_never_out`}{The number of observations (with positive prior
#'     weight) that no fit leaves out, which do not enter the estimate;
#'     `NA` with `eval_weights`.}
#' }
#'
#' @seealso [gb_stop_ctrl()], [extend_gbgam()]
#'
#' @references
#'
#' Efron, B., & Tibshirani, R. (1997).
#' Improvements on cross-validation: The .632+ bootstrap method.
#' *Journal of the American Statistical Association, 92*(438), 548--560.
#' \href{https://doi.org/10.1080/01621459.1997.10474007}{doi:10.1080/01621459.1997.10474007}
#'
#' @examples
#' # 10-fold cross-validation of the deviance (the default)
#' gb_resampling()
#'
#' # .632+ bootstrap with 50 samples
#' gb_resampling("bootstrap", B = 50)
#'
#' \donttest{
#' set.seed(1)
#' dat <- mgcv::gamSim(1, n = 200, verbose = FALSE)
#' mod <- gbgam(y ~ t2(x0) + t2(x1) + t2(x2), data = dat,
#'              stop_control = gb_stop_ctrl("cv",
#'                resampling = gb_resampling(B = 5)),
#'              verbose = FALSE)
#' mod$mstop
#' plot(mod$cv$risk, type = "l", xlab = "iteration", ylab = "CV risk")
#' }
#' @export
gb_resampling <- function(type = c("kfold", "bootstrap", "subsample"),
                          B = NULL,
                          mode = c("free_selection", "free_smoothing",
                                   "fixed"),
                          estimator = NULL,
                          folds = NULL,
                          weights = NULL,
                          eval_weights = NULL,
                          fraction = NULL,
                          chunk = 100L,
                          keep_fits = FALSE) {
  type_given <- !missing(type)
  type <- match.arg(type)
  check_choice(mode, eval(formals()$mode), "mode")
  mode <- match.arg(mode)
  if (!is.null(B) && (!is.null(folds) || !is.null(weights)))
    stop("B is not used with ", if (!is.null(folds)) "folds" else "weights",
         ": the number of fits is taken from it.", call. = FALSE)
  if (!is.null(weights)) {
    checkmate::assert_matrix(weights, mode = "numeric", any.missing = FALSE,
                             min.cols = 2L, .var.name = "weights")
    if (any(weights < 0) || any(!is.finite(weights)))
      stop("weights must be finite and non-negative.", call. = FALSE)
    if (!is.null(folds))
      stop("supply either folds or weights, not both.", call. = FALSE)
    if (type_given)
      stop("type is not used with weights: the weights define the ",
           "resampling scheme.", call. = FALSE)
    type <- "custom"
    B    <- ncol(weights)
  }
  if (!is.null(eval_weights)) {
    if (is.null(weights))
      stop("eval_weights can only be supplied together with weights.",
           call. = FALSE)
    checkmate::assert_matrix(eval_weights, mode = "numeric",
                             any.missing = FALSE, nrows = nrow(weights),
                             ncols = ncol(weights),
                             .var.name = "eval_weights")
    if (any(!is.finite(eval_weights)))
      stop("eval_weights must be finite.", call. = FALSE)
  }
  if (!is.null(folds)) {
    if (type != "kfold")
      stop("folds can only be supplied for type = \"kfold\".", call. = FALSE)
    checkmate::assert_atomic_vector(folds, any.missing = FALSE,
                                    .var.name = "folds")
    B <- length(unique(folds))
    if (B < 2L) stop("folds must contain at least 2 folds.", call. = FALSE)
  }
  B <- B %||% if (type == "kfold") 10L else 25L
  checkmate::assert_count(B, .var.name = "B")
  if (B < 2L) stop("B must be at least 2.", call. = FALSE)

  estimator <- estimator %||% if (type == "bootstrap") "632plus" else "oob"
  checkmate::assert_choice(estimator, c("oob", "632", "632plus"),
                           .var.name = "estimator")
  if (estimator != "oob" && !type %in% c("bootstrap", "custom"))
    stop("estimator = \"", estimator, "\" is only available for ",
         "type = \"bootstrap\" or custom weights.", call. = FALSE)

  if (type == "subsample") {
    fraction <- fraction %||% 0.5
    checkmate::assert_number(fraction, lower = 0, upper = 1,
                             .var.name = "fraction")
    if (fraction <= 0 || fraction >= 1)
      stop("fraction must lie strictly between 0 and 1.", call. = FALSE)
  } else if (!is.null(fraction)) {
    stop("fraction is only used with type = \"subsample\".", call. = FALSE)
  }

  checkmate::assert_count(chunk, positive = TRUE, .var.name = "chunk")
  checkmate::assert_flag(keep_fits, .var.name = "keep_fits")

  structure(
    list(type         = type,
         B            = as.integer(B),
         estimator    = estimator,
         folds        = folds,
         weights      = weights,
         eval_weights = eval_weights,
         fraction     = fraction,
         mode         = mode,
         chunk        = as.integer(chunk),
         keep_fits    = keep_fits),
    class = "gb_resampling")
}


# ---- 1. The resampling plan (plan_resampling()) ------------------------------

#' In-bag weights of the resample fits
#'
#' @param resampling A `"gb_resampling"` object.
#' @param n Number of rows used in the fit.
#' @return An n x B matrix: how often each row is in each resample.
#' @noRd
make_inbag_weights <- function(resampling, n) {
  if (!is.null(resampling$weights)) {
    if (nrow(resampling$weights) != n)
      stop("gb_resampling(weights = ) must have one row per observation ",
           "used in the fit (", n, "), got ", nrow(resampling$weights), ".",
           call. = FALSE)
    return(unname(resampling$weights * 1))  # logical/integer -> double
  }
  if (resampling$type == "kfold") {
    kfold_inbag_weights(n, resampling$B, resampling$folds)
  } else if (resampling$type == "bootstrap") {
    bootstrap_inbag_weights(n, resampling$B)
  } else {
    subsample_inbag_weights(n, resampling$B, resampling$fraction)
  }
}

#' In-bag weights of k-fold cross-validation: 0 for the rows of fold b
#'
#' @param n Number of rows.
#' @param B Number of folds.
#' @param folds User fold labels, or NULL for random folds of (nearly) equal
#'   size.
#' @noRd
kfold_inbag_weights <- function(n, B, folds) {
  if (is.null(folds)) {
    if (B > n)
      stop("gb_resampling(B = ", B, "): the number of folds must not ",
           "exceed the number of observations used in the fit (", n, ").",
           call. = FALSE)
    fold <- sample(rep_len(seq_len(B), n))
  } else {
    if (length(folds) != n)
      stop("gb_resampling(folds = ) must have one label per observation ",
           "used in the fit (", n, "), got ", length(folds), ".",
           call. = FALSE)
    fold <- as.integer(factor(folds))
  }
  outer(fold, seq_len(B), `!=`) * 1
}

#' In-bag weights of the nonparametric bootstrap: draw counts of n draws
#' with replacement
#' @noRd
bootstrap_inbag_weights <- function(n, B) {
  inbag_weights <- matrix(0, n, B)
  for (b in seq_len(B)) {
    inbag_weights[, b] <- tabulate(sample.int(n, n, replace = TRUE), n)
  }
  inbag_weights
}

#' In-bag weights of subsampling: 1 for a `fraction` of the rows drawn
#' without replacement
#' @noRd
subsample_inbag_weights <- function(n, B, fraction) {
  size <- max(1L, floor(fraction * n))
  inbag_weights <- matrix(0, n, B)
  for (b in seq_len(B)) {
    inbag_weights[sample.int(n, size), b] <- 1
  }
  inbag_weights
}

#' Evaluation weights of the out-of-bag error
#'
#' Each row's deviance is averaged over the fits that hold it out, and the
#' row averages are averaged by prior weight (the leave-one-out bootstrap of
#' Efron & Tibshirani, 1997; k-fold cross-validation when each row is held out once).
#'
#' @param inbag_weights In-bag weights (n x B).
#' @param prior_weights Prior weights of the main fit.
#' @return A list with `evaluation_weights` (n x B) and `n_never_out`, the
#'   number of rows with positive prior weight that no fit holds out.
#' @noRd
oob_evaluation_weights <- function(inbag_weights, prior_weights) {
  held_out       <- inbag_weights == 0
  times_held_out <- rowSums(held_out)
  counted        <- times_held_out > 0 & prior_weights > 0
  if (!any(counted))
    stop("No observation is left out by any resampling fit, so the ",
         "prediction error cannot be estimated.", call. = FALSE)
  row_scale <- ifelse(counted, 1 / pmax(times_held_out, 1), 0) /
    sum(prior_weights[counted])
  list(evaluation_weights = held_out * row_scale,
       n_never_out        = sum(!counted & prior_weights > 0))
}

#' Warn when resample fits or rows do not enter the out-of-bag error
#'
#' @param plan A resampling plan (see `plan_resampling()`).
#' @noRd
warn_out_of_bag_coverage <- function(plan) {
  n <- nrow(plan$inbag_weights)
  B <- ncol(plan$inbag_weights)
  holds_none_out <- colSums(plan$inbag_weights == 0) == 0
  if (is.null(plan$resampling$eval_weights) && any(holds_none_out))
    warning(sprintf(paste0(
      "%d of %d resampling fits leave no observation out and do not enter ",
      "the estimate."), sum(holds_none_out), B), call. = FALSE)
  max_never_out <- 0.05  # share of rows
  if (isTRUE(plan$n_never_out > max_never_out * n))
    warning(sprintf(paste0(
      "%d of %d observations are not left out by any resampling fit and ",
      "do not enter the estimate; increase B."), plan$n_never_out, n),
      call. = FALSE)
  invisible(NULL)
}

#' The resampling plan: everything fixed before boosting starts
#'
#' @param mod The main model (assembled by `gbgam()`, or a fitted model when
#'   extending it).
#' @param inbag_weights In-bag weights, or NULL to draw them from
#'   `mod$resampling`.
#' @param rng_state A `.Random.seed` value, or NULL to record the current one.
#' @return A list of class `"gb_resampling_plan"` with fields
#'   \describe{
#'     \item{template}{the main model without its fitted state
#'       (`resample_fit_template()`); every resample fit starts from it.}
#'     \item{resampling}{the `"gb_resampling"` settings.}
#'     \item{inbag_weights}{n x B in-bag weights; B is its number of
#'       columns.}
#'     \item{evaluation_weights}{n x B evaluation weights.}
#'     \item{n_never_out}{rows with positive prior weight that no fit holds
#'       out (`NA` with custom evaluation weights).}
#'     \item{rng_state}{the RNG state from which the B resample fits draw
#'       their Hutch++ probes, one fit after the other, so that
#'       `extend_gbgam()` can rebuild them.}
#'     \item{replay}{what the resample fits copy from the main fit
#'       (`resample_fit_replay()`).}
#'   }
#' @noRd
plan_resampling <- function(mod, inbag_weights = NULL, rng_state = NULL) {
  resampling <- mod$resampling
  n <- length(mod$y)
  if (is.null(inbag_weights)) {
    inbag_weights <- make_inbag_weights(resampling, n)
  }
  if (is.null(resampling$eval_weights)) {
    evaluation <- oob_evaluation_weights(inbag_weights, mod$prior.weights)
  } else {
    evaluation <- list(
      evaluation_weights = unname(resampling$eval_weights * 1),  # -> double
      n_never_out        = NA_integer_)
  }
  plan <- structure(list(
    template           = resample_fit_template(mod),
    resampling         = resampling,
    inbag_weights      = inbag_weights,
    evaluation_weights = evaluation$evaluation_weights,
    n_never_out        = evaluation$n_never_out,
    rng_state          = rng_state %||% current_rng_state(),
    replay             = resample_fit_replay(resampling$mode, mod)
  ), class = "gb_resampling_plan")
  stopifnot(nrow(inbag_weights) == n,
            identical(dim(plan$evaluation_weights), dim(inbag_weights)))
  plan
}

#' The main model without the fields a resample fit does not start from
#'
#' @param mod The main model (before or after fitting).
#' @return The template of the resample fits.
#' @noRd
resample_fit_template <- function(mod) {
  drop <- c("learners", "derivs", "resume", "fitted_values", "metrics", "cv",
            "edfs", "Vp", "Ve", "boot_coef", "boot_learner_coef",
            "boot_mstop", "boot_censored", "widening_factor",
            "beta_samples", "coefficients_bayes", "mh_accept_rate",
            "bb_info", "residuals", "weights", "fitted.values",
            "linear.predictors", "resampling")
  mod[setdiff(names(mod), drop)]
}

#' What the resample fits copy from the main fit
#'
#' @param mode The resampling mode (see [gb_resampling()]).
#' @param mod The main model.
#' @return A list: `selection` (replay the main fit's learner selection),
#'   `smoothing` (replay its smoothing parameters) and `edf_compute_method`
#'   of the resample fits (`"none"` when the smoothing is replayed, since
#'   ASR then needs no EDF).
#' @noRd
resample_fit_replay <- function(mode, mod) {
  smoothing <- mode == "fixed"
  list(
    selection          = mode != "free_selection",
    smoothing          = smoothing,
    edf_compute_method = if (smoothing) "none" else mod$edf_compute_method
  )
}


# ---- 2. The resample fits at iteration 0 (build_resample_fits()) -------------

#' Monitor of a resample fit: its share of the out-of-bag error
#'
#' @param evaluation_weights The fit's column of the evaluation weights.
#' @param family,y,prior_weights Of the main fit.
#' @return A function of the running resample model returning the sum over
#'   the held-out rows of evaluation weight times deviance (0 for a fit
#'   without evaluation weight).
#' @noRd
heldout_deviance_monitor <- function(evaluation_weights, family, y,
                                     prior_weights) {
  force(family)
  held_out <- which(evaluation_weights != 0)
  if (length(held_out) == 0L) return(function(mod) 0)
  v <- evaluation_weights[held_out]
  y <- y[held_out]
  w <- prior_weights[held_out]
  function(mod) {
    mu <- family$linkinv(mod$linear.predictors[held_out])
    sum(v * family$dev.resids(y, mu, w))
  }
}

#' Set up one resample fit at iteration 0
#'
#' Draws its Hutch++ probes (if any) from the current RNG stream.
#'
#' @param template From `resample_fit_template()`.
#' @param inbag_weights,evaluation_weights The fit's columns of the plan's
#'   weights.
#' @param replay From `resample_fit_replay()`.
#' @return A fitting state (see `fit_state_init()`).
#' @noRd
new_resample_fit <- function(template, inbag_weights, evaluation_weights,
                             replay) {
  fit <- reweight_model(template, inbag_weights, replay$edf_compute_method)
  fit <- without_own_stopping(fit, replay)
  fit <- init_derivs(fit, length(fit$y), quiet = TRUE)
  fit_state_init(fit, start_time = Sys.time(), save_fitted_values = FALSE,
                 verbose = FALSE,
                 monitor = heldout_deviance_monitor(
                   evaluation_weights, template$family, template$y,
                   template$prior.weights))
}

#' The template refitted for its in-bag weights: starting intercept and
#' learners as in `gbgam()`, on the bases of the main fit's smooths
#'
#' @param template From `resample_fit_template()`.
#' @param inbag_weights The fit's in-bag weights.
#' @param edf_compute_method EDF method of the resample fit.
#' @return The model of the resample fit at iteration 0 (no derivatives).
#' @noRd
reweight_model <- function(template, inbag_weights, edf_compute_method) {
  fit <- template
  n   <- length(fit$y)
  pw  <- template$prior.weights * inbag_weights
  fit$prior.weights <- pw
  fit$user_weights  <- (template$user_weights %||% rep(1, n)) * inbag_weights

  init <- isTRUE(fit$initialise_intercept)
  fit$intercept0 <- find_intercept_mle(
    y = fit$y, family = fit$family, prior.weights = pw, offset = fit$offset)
  a0 <- if (init) fit$intercept0 else 0
  fit$edf_compute_method <- edf_compute_method
  compute_edf <- edf_compute_method != "none"
  lb <- tryCatch(
    new_learners_internal(smooths = template$smooth, data = template$model,
                          pterms = template$pterms, weights = pw, ahat = a0,
                          edf0 = as.numeric(compute_edf && init),
                          knots = template$knots, clamp_sp = TRUE),
    error = function(e) stop(
      "A resampling fit could not be set up: ", conditionMessage(e),
      " Too few observations of some term may be left in the fit; use ",
      "fewer folds, larger subsamples or other folds.", call. = FALSE))
  fit$learners          <- lb$learners
  fit$cen_const         <- lb$cen_const
  fit$linear.predictors <- a0 + fit$offset
  fit$edf.total         <- if (init) 1 else 0
  fit
}

#' Remove a resample fit's own stopping rule and replay settings
#'
#' A resample fit runs until the iteration the main fit has reached; it
#' copies the main fit's selection (and smoothing parameters) when advanced
#' (`advance_resample_fit()`), so no path is preset here.
#'
#' @param fit From `reweight_model()`.
#' @param replay From `resample_fit_replay()`.
#' @return `fit` with stopping method `"none"`, the stopping and replay
#'   fields `NULL` (kept as names, which `fit_state_init()` requires), no
#'   progress output, and ASR off (`tau = 0`) when the smoothing is replayed.
#' @noRd
without_own_stopping <- function(fit, replay) {
  fit$stopping_method <- "none"
  unset <- c("max_iter", "iter_no_improv", "iter_no_improv_rel",
             "max_runtime_mins", "preselected_learner_indices",
             "preselected_rhos", "preselected_post_rhos")
  fit[unset] <- list(NULL)
  fit$verbose <- FALSE
  if (replay$smoothing) fit$tau <- 0
  fit
}

#' Set up the B resample fits at iteration 0, in order
#'
#' The fits draw their Hutch++ probes one after the other from the current
#' RNG stream; `plan$rng_state` is that stream's state before the first
#' draw.
#'
#' @param plan From `plan_resampling()`.
#' @return A list of B fitting states.
#' @noRd
build_resample_fits <- function(plan) {
  B    <- ncol(plan$inbag_weights)
  fits <- vector("list", B)
  for (b in seq_len(B)) {
    fits[[b]] <- new_resample_fit(plan$template, plan$inbag_weights[, b],
                                  plan$evaluation_weights[, b], plan$replay)
  }
  fits
}

#' Give stored resample fits back their monitors (see
#' `strip_resample_fits()`)
#'
#' @param fits Stored fitting states.
#' @param plan From `plan_resampling()`.
#' @return The fitting states, ready to be advanced.
#' @noRd
reattach_heldout_monitors <- function(fits, plan) {
  template <- plan$template
  for (b in seq_along(fits)) {
    fits[[b]]$monitor <- heldout_deviance_monitor(
      plan$evaluation_weights[, b], template$family, template$y,
      template$prior.weights)
  }
  fits
}

#' The resample fits without their monitors and best models, to store in
#' the model (`$resume$cv$fits`)
#' @noRd
strip_resample_fits <- function(fits) {
  for (b in seq_along(fits)) {
    fits[[b]]$monitor <- NULL
    fits[[b]]$opt_mod <- NULL
  }
  fits
}


# ---- 3. The main fit (fit_state_init()) --------------------------------------

#' Monitor of the main fit: its apparent and no-information errors
#'
#' @param family,y,prior_weights Of the main fit.
#' @param with_no_information Logical; compute the no-information error
#'   (needed by `"632plus"` only), else `NA`.
#' @return A function of the running main model returning
#'   `c(apparent = , no_information = )`; `fit_state_init()` keeps its values
#'   in `$trace`.
#' @noRd
apparent_error_monitor <- function(family, y, prior_weights,
                                   with_no_information) {
  force(family)
  force(y)
  force(with_no_information)
  total_weight <- sum(prior_weights)
  function(mod) {
    mu <- family$linkinv(mod$linear.predictors)
    c(apparent = sum(family$dev.resids(y, mu, prior_weights)) / total_weight,
      no_information = if (with_no_information) {
        no_information_deviance(y, mu, prior_weights, family)
      } else {
        NA_real_
      })
  }
}


# ---- 4. Boosting with resampling (boost_with_resampling()) -------------------

#' The main fit's selection path up to its current iteration
#'
#' @param main Fitting state of the main fit.
#' @return A list with the selected learners `k_opts`, their smoothing
#'   parameters `rhos` and those of the post-iteration learners `post_rhos`.
#' @noRd
main_selection_path <- function(main) {
  iterations <- seq_len(main$m)
  list(k_opts    = main$k_opts[iterations],
       rhos      = main$rhos[iterations],
       post_rhos = main$post_rhos[iterations, , drop = FALSE])
}

#' Advance a resample fit to an iteration
#'
#' @param fit Its fitting state.
#' @param to_iteration Iteration to reach.
#' @param main_path From `main_selection_path()`, up to at least
#'   `to_iteration`.
#' @param replay From `resample_fit_replay()`.
#' @return The updated fitting state.
#' @noRd
advance_resample_fit <- function(fit, to_iteration, main_path, replay) {
  if (to_iteration <= fit$m) return(fit)
  if (replay$selection) {
    fit$mod$preselected_learner_indices <- main_path$k_opts
    if (replay$smoothing) {
      fit$mod$preselected_rhos <- main_path$rhos
      if (ncol(main_path$post_rhos) > 0L)
        fit$mod$preselected_post_rhos <- main_path$post_rhos
    }
  }
  fit_state_run(fit, n_iter = to_iteration - fit$m, use_stopping = FALSE)
}

#' Advance every resample fit to the main fit's iteration
#'
#' @param fits Fitting states of the resample fits.
#' @param main Fitting state of the main fit.
#' @param replay From `resample_fit_replay()`.
#' @return The updated fitting states.
#' @noRd
advance_resample_fits <- function(fits, main, replay) {
  main_path <- main_selection_path(main)
  for (b in seq_along(fits)) {
    fits[[b]] <- tryCatch(
      advance_resample_fit(fits[[b]], main$m, main_path, replay),
      error = function(e) stop(
        "Resampling fit ", b, " failed (the main fit is not affected): ",
        conditionMessage(e), call. = FALSE))
  }
  fits
}

#' Out-of-bag deviance paths of the resample fits
#'
#' @param fits Fitting states of the resample fits, all at the same
#'   iteration m.
#' @return An (m + 1) x B matrix: column b is fit b's share of the
#'   out-of-bag error at iterations 0 to m.
#' @noRd
resample_risk_paths <- function(fits) {
  risk_paths <- matrix(NA_real_, length(fits[[1L]]$trace), length(fits))
  for (b in seq_along(fits)) {
    risk_paths[, b] <- unlist(fits[[b]]$trace)
  }
  risk_paths
}

#' Estimate of the prediction error at every iteration
#'
#' @param estimator `"oob"`, `"632"` or `"632plus"`.
#' @param apparent,oob,no_information Paths of the apparent, out-of-bag and
#'   no-information errors (the last for `"632plus"` only).
#' @return A list with the estimate `risk` and, for `"632plus"`, the
#'   relative overfitting rate and the weight on the out-of-bag error (see
#'   `err_632plus()`; `NULL` otherwise).
#' @noRd
resampling_estimate <- function(estimator, apparent, oob, no_information) {
  if (estimator == "oob") {
    list(risk = oob, relative_overfitting = NULL, weight = NULL)
  } else if (estimator == "632") {
    list(risk = err_632(apparent, oob), relative_overfitting = NULL,
         weight = NULL)
  } else {
    err_632plus(apparent, oob, no_information)
  }
}

#' Name of the estimate, for the progress display and summary()
#'
#' @param resampling A `"gb_resampling"` object.
#' @noRd
resampling_risk_label <- function(resampling) {
  if (resampling$estimator == "632plus") {
    ".632+ risk"
  } else if (resampling$estimator == "632") {
    ".632 risk"
  } else if (resampling$type == "kfold") {
    "CV risk"
  } else {
    "OOB risk"
  }
}

#' Number of main-fit iterations until the next check of the stopping rule
#'
#' @param main Fitting state of the main fit.
#' @param chunk Iterations between two checks.
#' @return `chunk`, or fewer to stop at `max_iter` (0 once it is reached).
#' @noRd
iterations_to_next_check <- function(main, chunk) {
  max_iter <- main$mod$max_iter
  if (is.null(max_iter)) return(chunk)
  max(0, min(chunk, max_iter - main$m))
}

#' Stopping conditions of resampling stopping, checked after a chunk
#'
#' @param mod The main model (its `iter_no_improv`, `iter_no_improv_rel`,
#'   `max_iter` and `max_runtime_mins`; `NULL` turns a condition off).
#' @param iteration Iteration the fits have reached.
#' @param best_iteration Iteration with the lowest estimate so far.
#' @param elapsed_mins Minutes since the start of the fit.
#' @return A named logical vector, as recorded by `fit_state_run()`.
#' @noRd
resampling_stopping_conditions <- function(mod, iteration, best_iteration,
                                           elapsed_mins) {
  c(no_improvements = !is.null(mod$iter_no_improv) &&
      (iteration - best_iteration >= patience(mod, best_iteration)),
    mstop_reached = !is.null(mod$max_iter) && iteration >= mod$max_iter,
    time_limit_exceeded = !is.null(mod$max_runtime_mins) &&
      elapsed_mins > mod$max_runtime_mins)
}

#' Advance the main fit and the resample fits until the stopping rule
#'
#' @param main Fitting state of the main fit (with `record_history = TRUE`
#'   and `apparent_error_monitor()`).
#' @param fits Fitting states of the resample fits.
#' @param plan From `plan_resampling()`.
#' @param verbose Logical; report the estimate at every check.
#' @return A list with the final states `main` and `fits`, the paths
#'   `risk_paths` (from `resample_risk_paths()`) and `main_errors` (an
#'   (m + 1) x 2 matrix with columns `apparent` and `no_information`), the
#'   `estimate` (from `resampling_estimate()`), `best_iteration` and the
#'   `stopping_conditions`.
#' @noRd
boost_with_resampling <- function(main, fits, plan, verbose) {
  resampling <- plan$resampling
  repeat {
    n_iter <- iterations_to_next_check(main, resampling$chunk)
    if (n_iter > 0)
      main <- fit_state_run(main, n_iter = n_iter, use_stopping = FALSE)
    fits <- advance_resample_fits(fits, main, plan$replay)

    risk_paths  <- resample_risk_paths(fits)
    main_errors <- do.call(rbind, main$trace)
    estimate    <- resampling_estimate(
      resampling$estimator,
      apparent       = main_errors[, "apparent"],
      oob            = rowSums(risk_paths),
      no_information = main_errors[, "no_information"])
    best_iteration <- which.min(estimate$risk) - 1L

    elapsed_mins <- as.numeric(difftime(Sys.time(), main$start_time,
                                        units = "mins"))
    stopping_conditions <- resampling_stopping_conditions(
      main$mod, main$m, best_iteration, elapsed_mins)
    if (verbose) {
      cli::cli_text(
        "iteration {main$m}: best m = {best_iteration}, ",
        "{resampling_risk_label(resampling)} = ",
        "{signif(estimate$risk[best_iteration + 1L], 6)} ",
        "({format_minutes(elapsed_mins)})")
    }
    if (any(stopping_conditions)) break
  }
  list(main = main, fits = fits, risk_paths = risk_paths,
       main_errors = main_errors, estimate = estimate,
       best_iteration = best_iteration,
       stopping_conditions = stopping_conditions)
}


# ---- 5. The returned model (resampling_result()) -----------------------------

#' The fitted model's `$cv` (documented in [gb_resampling()])
#'
#' @param run From `boost_with_resampling()`.
#' @param plan From `plan_resampling()`.
#' @noRd
resampling_summary <- function(run, plan) {
  resampling <- plan$resampling
  list(
    type        = resampling$type,
    estimator   = resampling$estimator,
    B           = ncol(plan$inbag_weights),
    mstop       = run$best_iteration,
    risk        = run$estimate$risk,
    apparent    = run$main_errors[, "apparent"],
    oob         = rowSums(run$risk_paths),
    noinfo      = if (resampling$estimator == "632plus")
                    run$main_errors[, "no_information"] else NULL,
    weight      = run$estimate$weight,
    rep_risk    = run$risk_paths,
    weights     = plan$inbag_weights,
    n_never_out = plan$n_never_out
  )
}

#' The main fit at the best iteration, with the resampling results
#'
#' @param run From `boost_with_resampling()`.
#' @param plan From `plan_resampling()`.
#' @param resumable Logical; attach the running states as `$resume`.
#' @return The fitted model.
#' @noRd
resampling_result <- function(run, plan, resumable) {
  main <- run$main
  best <- run$best_iteration
  risk <- run$estimate$risk

  # The main fit at the best iteration; its loss path is the estimate
  main$opt_mod             <- fit_state_rewind(main, best)
  main$stopping_conditions <- run$stopping_conditions
  main$ls0                 <- risk[1L]
  main$ls                  <- risk[-1L]
  main$min_loss            <- risk[best + 1L]
  model <- fit_state_result(main, resumable = resumable)

  model$cv <- resampling_summary(run, plan)
  if (resumable) {
    model$resume$cv <- list(
      rng_state = plan$rng_state,
      fits      = if (isTRUE(plan$resampling$keep_fits))
                    strip_resample_fits(run$fits) else NULL)
  }
  model
}


# ---- Entry points ------------------------------------------------------------

#' Fit with stopping by a resampling estimate
#'
#' Called by gbgam() for `stopping_method = "cv"`.
#'
#' @inheritParams fit_gbgam
#' @noRd
fit_gbgam_cv <- function(mod, start_time, save_fitted_values, verbose,
                         resumable = TRUE) {
  if (!is.null(mod$preselected_learner_indices))
    stop("stopping_method = \"cv\" cannot be combined with ",
         "preselected_learner_indices.", call. = FALSE)
  plan <- plan_resampling(mod)
  warn_out_of_bag_coverage(plan)
  fits <- build_resample_fits(plan)
  main <- fit_state_init(
    mod, start_time, save_fitted_values, verbose = FALSE,
    record_history = TRUE,
    monitor = apparent_error_monitor(
      mod$family, mod$y, mod$prior.weights,
      with_no_information = plan$resampling$estimator == "632plus"))
  run <- boost_with_resampling(main, fits, plan, verbose)
  resampling_result(run, plan, resumable)
}

#' Extend a fit stopped by a resampling estimate (for extend_gbgam())
#'
#' Runs exactly `n_iter` further iterations, or until the time limit: the
#' patience is not applied to an extension (`extend_gbgam()` restores the
#' model's patience afterwards). The resample fits are continued if stored,
#' and otherwise rebuilt from the recorded RNG state and rerun.
#'
#' @noRd
extend_gbgam_cv <- function(mod, n_iter, start_time, verbose,
                            max_runtime_mins) {
  stored <- mod$resume$cv
  if (is.null(stored$rng_state))
    stop("extend_gbgam: the resampling fits of this model cannot be ",
         "rebuilt: it was fitted with an older gbgam. Refit it.",
         call. = FALSE)
  plan <- plan_resampling(mod, inbag_weights = mod$cv$weights,
                          rng_state = stored$rng_state)
  fits <- if (is.null(stored$fits)) {
    with_rng_state(plan$rng_state, build_resample_fits(plan))
  } else {
    reattach_heldout_monitors(stored$fits, plan)
  }
  main <- fit_state_resume(
    mod, start_time, verbose = FALSE,
    monitor = apparent_error_monitor(
      mod$family, mod$y, mod$prior.weights,
      with_no_information = plan$resampling$estimator == "632plus"))
  main$mod$max_iter <- main$m + n_iter
  main$mod[c("iter_no_improv", "iter_no_improv_rel")] <- list(NULL)
  main$mod["max_runtime_mins"] <- list(max_runtime_mins)

  run <- boost_with_resampling(main, fits, plan, verbose)
  resampling_result(run, plan, resumable = TRUE)
}
