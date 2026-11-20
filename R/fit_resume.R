# Resumable fits: the running state of a boosting fit, its restoration and
# extend_gbgam().
#
# A learner's running state is small ("light"): its coefficients, EDF and
# smoothing state. Its fitted values are recomputed from the coefficients,
# and its structural fields (design, Demmler-Reinsch basis, penalty) are
# taken from the learners of the returned model, so that a resumable model
# does not store the design matrices twice.

# Minutes of the time limit kept back from the uncertainty
# quantification after boosting
finalise_margin_mins <- 0.01


#' Light running state of a learner
#'
#' @param lr A learner.
#' @return A list with the fields that boosting updates.
#' @noRd
learner_light <- function(lr) {
  list(beta           = lr$beta,
       edf            = lr$edf,
       rho            = lr$rho,
       nu_k           = lr$nu_k,
       lambda_at_zero = lr$lambda_at_zero)
}


#' Light states of the learners updated in one iteration
#'
#' @param mod The running model after the iteration.
#' @param k Indices of the learners updated in the iteration.
#' @return A list with the indices `k` and their light states `s`.
#' @noRd
learner_history_entry <- function(mod, k) {
  k <- k[!is.na(k)]
  list(k = k, s = lapply(mod$learners[k], learner_light))
}


#' Restore a learner from its light state
#'
#' @param lr A learner with the structural fields (design, A, D, ...).
#' @param light Its light state from `learner_light()`.
#' @param beta0 Its coefficients at the start of boosting; its fitted values
#'   are `design %*% (beta - beta0)`, as accumulated by `update_learner()`.
#' @return The learner in the given state.
#' @noRd
restore_learner <- function(lr, light, beta0) {
  for (nm in names(light)) lr[[nm]] <- light[[nm]]
  lr$eta <- as.vector(lr$design %*% (light$beta - beta0))
  lr
}


#' Linear predictor implied by the learners
#'
#' @param mod A model whose learners carry their fitted values `eta`.
#' @param beta0 Coefficients of the learners at the start of boosting.
#' @return The offset plus the starting fit plus the sum of the learners'
#'   fitted values.
#' @noRd
learners_linear_predictor <- function(mod, beta0) {
  start <- Reduce(`+`, Map(function(lr, b0) as.vector(lr$design %*% b0),
                           mod$learners, beta0))
  as.vector(mod$offset + start +
              Reduce(`+`, lapply(mod$learners, `[[`, "eta")))
}


#' Pack the running part of a fitting state
#'
#' The best model and the structural fields of the learners are not stored:
#' `fit_state_resume()` takes them from the returned model.
#'
#' @param st A fitting state.
#' @return A list, stored as `$resume` on the returned model.
#' @noRd
fit_state_pack <- function(st) {
  mod <- st$mod
  list(
    learners          = lapply(mod$learners, learner_light),
    linear.predictors = mod$linear.predictors,
    edf.total         = mod$edf.total,
    k_opt             = mod$k_opt,
    derivs            = mod$derivs,
    m                 = st$m,
    counter           = st$counter,
    ls0               = st$ls0,
    min_loss          = st$min_loss,
    edf0              = st$edf0,
    k_opts            = st$k_opts,
    rhos              = st$rhos,
    ls                = st$ls,
    edfs              = st$edfs,
    post_rhos         = st$post_rhos,
    record_history    = st$record_history,
    history           = st$history,
    init_light        = st$init_light,
    beta0             = st$beta0,
    trace             = st$trace
  )
}


#' Rebuild a fitting state from a resumable model
#'
#' @param object A model with `$resume` (from `gbgam(..., stop_control =
#'   gb_stop_ctrl(resumable = TRUE))`).
#' @param start_time Start of the extension, for time limits and progress.
#' @param verbose Logical; show progress.
#' @param monitor As in `fit_state_init()`.
#' @return A fitting state positioned at the last iteration of the fit.
#' @noRd
fit_state_resume <- function(object, start_time, verbose, monitor = NULL) {
  rs             <- object$resume
  opt_mod        <- object
  opt_mod$resume <- NULL

  mod                   <- opt_mod
  mod$learners          <- Map(restore_learner, mod$learners, rs$learners,
                               rs$beta0)
  mod$linear.predictors <- rs$linear.predictors
  mod$edf.total         <- rs$edf.total
  mod$k_opt             <- rs$k_opt
  mod$derivs            <- rs$derivs
  mod$verbose           <- verbose

  # The best model as it was before finalisation: edf.total is the running
  # total at mstop (finalise_gbgam() replaces it by the sum over learners)
  if (mod$edf_compute_method != "none") {
    opt_mod$edf.total <- c(rs$edf0, rs$edfs)[opt_mod$mstop + 1]
  }

  save_fitted   <- !is.null(object$fitted_values)
  fitted_values <- if (save_fitted) {
    lapply(seq_len(ncol(object$fitted_values)),
           function(j) object$fitted_values[, j])
  } else NULL

  list(
    mod                 = mod,
    opt_mod             = opt_mod,
    m                   = rs$m,
    counter             = rs$counter,
    ls0                 = rs$ls0,
    min_loss            = rs$min_loss,
    edf0                = rs$edf0,
    keep_derivs         = isTRUE(getOption("gbgam.keep_derivs", FALSE)),
    save_fitted_values  = save_fitted,
    start_time          = start_time,
    progress            = gb_progress_init(verbose, mod, start_time),
    post_idx            = mod$learners_to_fit_after_each_iter %||% integer(0),
    k_opts              = rs$k_opts,
    rhos                = rs$rhos,
    ls                  = rs$ls,
    edfs                = rs$edfs,
    post_rhos           = rs$post_rhos,
    fitted_values       = fitted_values,
    stopping_conditions = object$metrics$stopping_conditions,
    record_history      = rs$record_history,
    history             = rs$history,
    init_light          = rs$init_light,
    beta0               = rs$beta0,
    monitor             = monitor,
    trace               = rs$trace
  )
}


#' The model of a fitting state at an earlier iteration
#'
#' Requires a state run with `record_history = TRUE`.
#'
#' @param st A fitting state.
#' @param m_star Iteration (0 to `st$m`).
#' @return The model at iteration `m_star` (without EDF-tracking
#'   derivatives), with `$mstop = m_star`.
#' @noRd
fit_state_rewind <- function(st, m_star) {
  stopifnot(isTRUE(st$record_history), m_star >= 0, m_star <= st$m)
  mod   <- st$mod
  light <- st$init_light
  for (j in seq_len(m_star)) {
    h <- st$history[[j]]
    light[h$k] <- h$s
  }
  mod$learners <- Map(restore_learner, mod$learners, light, st$beta0)
  mod$linear.predictors <- learners_linear_predictor(mod, st$beta0)
  if (mod$edf_compute_method != "none") {
    mod$edf.total <- c(st$edf0, st$edfs)[m_star + 1]
  }
  mod$derivs <- NULL
  mod["k_opt"] <- list(if (m_star > 0) st$k_opts[m_star] else NULL)
  mod$mstop <- m_star
  mod
}


#' Continue a boosting fit
#'
#' @description
#' Runs further boosting iterations on a fitted model, from its last
#' iteration (not from `mstop`), so that a fit that ended too early, e.g. at
#' `max_iter`, does not have to be refitted from the start. The model must
#' have been fitted with `gb_stop_ctrl(resumable = TRUE)` (the default).
#'
#' @param mod A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param n_iter Integer, at least 1. Number of further iterations.
#' @param verbose Logical. Show the progress display of the boosting loop,
#'   as [gbgam()] does, and the progress bar of the recomputed uncertainty
#'   quantification. Default `FALSE`.
#' @param max_runtime_mins Non-negative number, or NULL. Time limit in
#'   minutes for the extension and the recomputed uncertainty
#'   quantification; when it is reached, the extension ends at the last
#'   iteration run. Default `NULL` (no limit).
#'
#' @return The extended model. Its `$metrics` cover the whole path, and
#'   `$mstop` is the best iteration over that path. `iter_no_improv`,
#'   `iter_no_improv_rel` and `max_runtime_mins` stay those of the original
#'   fit; `max_iter` is the last iteration run. A fit without a `max_iter`
#'   stays without a cap for its freely stopping bootstrap replicates: its
#'   `$uncapped` is `TRUE`. When the extension ended at its own time limit,
#'   `$metrics$max_runtime_mins` holds that limit, which
#'   [print()][print.gbgam()] and [summary()][summary.gbgam()] report. If the
#'   recomputed uncertainty quantification fails, the error has class
#'   `"gbgam_uq_error"` and keeps the extended model, without
#'   uncertainty quantification, in its `fit` field, as in [gbgam()].
#'
#' @details
#' The extension runs exactly `n_iter` iterations (fewer at the time
#' limit); the patience (`iter_no_improv`, `iter_no_improv_rel`) is not
#' applied to them. The best iteration is
#' then chosen over the whole path by the model's stopping method: the AIC
#' (or the `true_f` mean squared error) as in [gbgam()], or the resampling
#' estimate for `stopping_method = "cv"`, whose resampling fits are extended
#' as well.
#'
#' The uncertainty quantification results of `mod` (`$Vp` and the bootstrap
#' or Bayesian backfitting fields) are removed, and the method recorded in
#' `mod$uq_control`, the one that attached its `$Vp`, is run again
#' on the extended model with the same settings. If the time limit stops it
#' early, it keeps the replicates or draws completed so far; with fewer than
#' two, the extended model has no `$Vp` and its `$uq_control` is
#' `gb_uq_ctrl("none")`.
#'
#' @seealso [gbgam()] and [gb_stop_ctrl()] (`resumable`).
#'
#' @examples
#' set.seed(1)
#' dat <- mgcv::gamSim(1, n = 200, verbose = FALSE)
#' mod <- gbgam(y ~ t2(x0) + t2(x1) + t2(x2), data = dat,
#'              stop_control = gb_stop_ctrl(max_iter = 50), verbose = FALSE)
#' mod$mstop
#' mod2 <- extend_gbgam(mod, n_iter = 150)
#' mod2$mstop
#' @export
extend_gbgam <- function(mod, n_iter, verbose = FALSE,
                         max_runtime_mins = NULL) {
  checkmate::assert_class(mod, "gbgam")
  checkmate::assert_count(n_iter, positive = TRUE)
  check_verbose(verbose)
  check_runtime_mins(max_runtime_mins)
  if (is.null(mod$resume)) {
    stop("extend_gbgam: the model cannot be extended: it was fitted with ",
         "gb_stop_ctrl(resumable = FALSE). Refit it with ",
         "gb_stop_ctrl(resumable = TRUE).", call. = FALSE)
  }
  check_patience_fields(mod, "extend_gbgam")
  if (!is.null(mod$preselected_learner_indices) &&
      mod$resume$m + n_iter > length(mod$preselected_learner_indices)) {
    stop("extend_gbgam: the model replays preselected_learner_indices of ",
         "length ", length(mod$preselected_learner_indices), "; it cannot ",
         "be extended beyond that.", call. = FALSE)
  }
  start_time <- Sys.time()
  orig <- mod[c("iter_no_improv", "iter_no_improv_rel", "max_runtime_mins")]
  ic   <- mod$uq_control
  mod  <- remove_uq_fields(mod)

  final_model <- if (identical(mod$stopping_method, "cv")) {
    extend_gbgam_cv(mod, n_iter, start_time, verbose, max_runtime_mins)
  } else {
    # Run n_iter iterations (or until the time limit): patience is not
    # applied to the extension
    st <- fit_state_resume(mod, start_time, verbose)
    st$mod$max_iter            <- st$m + n_iter
    st$mod["iter_no_improv"]   <- list(NULL)
    st$mod["iter_no_improv_rel"] <- list(NULL)
    st$mod["max_runtime_mins"] <- list(max_runtime_mins)
    st$progress                <- gb_progress_init(verbose, st$mod, start_time)
    st                         <- fit_state_run(st, n_iter = n_iter)
    fit_state_result(st, resumable = TRUE)
  }
  # The extended path ends at the new max_iter; the user's patience and time
  # limit stay those of the original fit. A time limit that ended the
  # extension is kept in the metrics, for print() and summary()
  final_model$uncapped <- is.null(mod$max_iter) || isTRUE(mod$uncapped)
  final_model$max_iter <- final_model$resume$m %||% final_model$max_iter
  final_model[names(orig)] <- orig
  if (isTRUE(final_model$metrics$stopping_conditions["time_limit_exceeded"]))
    final_model$metrics$max_runtime_mins <- max_runtime_mins
  final_model <- finalise_gbgam(final_model, mod$formula,
                                list(pred.formula = mod$pred.formula))

  elapsed_mins   <- as.numeric(Sys.time() - start_time, units = "mins")
  remaining_mins <- if (is.null(max_runtime_mins)) NULL else
    max(0, max_runtime_mins - elapsed_mins - finalise_margin_mins)
  # An error here keeps the extended fit: it is attached to the condition
  # as `$fit`, as in gbgam()
  tryCatch(
    apply_uq(final_model, ic,
             verbose = verbose,
             max_runtime_mins = remaining_mins),
    error = function(e) stop_uq_failed(e, final_model, ic$method)
  )
}
