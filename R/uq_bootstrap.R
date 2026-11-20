#' Parametric bootstrap for a fitted `gbgam` model
#'
#' @description
#' Attaches a bootstrap covariance matrix `$Vp` to a fitted [gbgam()] model.
#' Each replicate simulates a new response at the fitted means and refits
#' the model, redoing on its own the steps that `mode` names; the covariance
#' of the replicates' coefficients, widened for smoothing bias unless
#' `alpha = NULL`, becomes `$Vp`. [plot()][plot.gbgam()] and
#' [predict(se.fit = TRUE)][predict.gbgam()] then use it, with no further
#' step; [smooth_intervals()] gives intervals at another level or with
#' another `correction`, and the limits as numbers.
#'
#' @inheritParams gb_param_boot_ctrl
#' @param mod A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param max_runtime_mins Non-negative number, or NULL. Time limit in
#'   minutes; when it is reached, the replicates completed so far are used
#'   (with a warning; their number is in `$n_samples_done`), and with fewer
#'   than 2 of them `mod` is returned unchanged (with a warning). Default
#'   `NULL` (no limit).
#'
#' @return
#' The original `mod` object with additional fields:
#' \describe{
#'   \item{`$Vp`}{Bootstrap covariance matrix of `mod$coefficients`,
#'     widened by the factor(s) in `$widening_factor`. Dimension `p x p`
#'     where `p = length(mod$coefficients)`.}
#'   \item{`$boot_coef`}{`replicates x p` matrix of the replicates'
#'     coefficients, one row per completed replicate.}
#'   \item{`$boot_learner_coef`}{List with one `replicates x p_k` matrix per
#'     learner: the replicates' coefficients of each learner on its own
#'     scale, one row per replicate. This is the transpose of the
#'     orientation of `$beta_samples` from [bayesian_backfitting()].}
#'   \item{`$boot_mstop`}{The stopping iteration of each replicate.}
#'   \item{`$boot_censored`}{Logical, per replicate: with mode
#'     `"free_stopping"`, whether the replicate reached `max_iter` before
#'     its stopping method stopped it (censored; with a warning). Always
#'     `FALSE` in the other modes.}
#'   \item{`$widening_factor`}{One number for `correction = "model"`; for
#'     `"smooth"` and `"learner"`, one per smooth and parametric coefficient,
#'     or one per learner, named by them. `NULL` if `alpha = NULL`.}
#'   \item{`$n_samples_done`}{The number of replicates completed:
#'     `n_samples`, or fewer when `max_runtime_mins` stopped the bootstrap.
#'     [bayesian_backfitting()] records its number of draws in the same
#'     field.}
#'   \item{`$uq_control`}{`gb_uq_ctrl("bootstrap")` with the
#'     settings that were run (`n_samples` as requested).}
#' }
#' Uncertainty quantification already attached to `mod` (by an earlier
#' bootstrap or by [bayesian_backfitting()]) is replaced. If no boosting
#' iteration improved on the intercept-only model (`mod$mstop = 0`), or the
#' time limit leaves fewer than 2 replicates, `mod` is returned unchanged,
#' with its earlier results if it had any, and a warning. Either way,
#' `$uq_control` describes the `$Vp` that the returned model carries.
#'
#' @details
#' ## The modes
#'
#' `mode` sets which steps of the fit each replicate redoes on its own (see
#' [gb_param_boot_ctrl()]). The steps it does not redo are copied from the
#' main fit: the sequence of selected learners, the smoothing parameter of
#' every learner at every iteration, and the stopping iteration
#' `mod$mstop`. A replicate that adapts its own smoothing parameters does
#' so by adaptive smoothness relaxation (ASR), with the main fit's
#' settings. With `"free_stopping"`, each replicate stops by the main fit's
#' stopping method, so the intervals include the uncertainty of the
#' stopping iteration (on choosing it, see Mayr et al., 2012); with
#' `stopping_method = "true_f"`, this is the mean squared error against the
#' true linear predictor `true_f` given in [gb_stop_ctrl()], the same vector
#' as for the main fit. Every replicate reuses the main fit's
#' learners (their design matrices and Demmler-Reinsch decompositions), so
#' a replicate costs only its boosting iterations.
#'
#' ## Widening
#'
#' The replicates are simulated around the fitted means, which carry the
#' fit's smoothing bias, so plain bootstrap intervals tend to be too narrow.
#' They are widened by factors found by Wood's (2020) rule:
#' `correction = "model"` finds one factor for the whole model, as Wood
#' (2020) does, and `"smooth"` and `"learner"` apply the same rule to each
#' smooth (and parametric coefficient) or to each learner. At each
#' observation \eqn{i}, the
#' normal-approximation interval is the main fit's value
#' \eqn{\tilde t_i}{t~_i} plus or minus
#' \eqn{\Phi^{-1}(1 - \frac{\alpha}{2}) \tilde s_i}{Phi^-1(1 - alpha/2) s~_i},
#' where \eqn{\Phi^{-1}}{Phi^-1} is the standard
#' normal quantile function and \eqn{\tilde s_i}{s~_i} the bootstrap
#' standard deviation of the value; the widening factor \eqn{c} stretches it
#' to \eqn{\tilde t_i \pm c \, \Phi^{-1}(1 - \frac{\alpha}{2}) \tilde s_i}{
#' t~_i +/- c Phi^-1(1 - alpha/2) s~_i}. The factor is chosen so that a
#' share 1 - `alpha` of all replicate values, over all observations, falls
#' inside the stretched intervals around the main fit, so the intervals
#' target coverage across the function, not at each observation (see
#' [smooth_intervals()]): \eqn{c} is the
#' 1 - `alpha` quantile of the replicates' absolute deviations from the main
#' fit, each divided by its \eqn{\tilde s_i}{s~_i}, divided by
#' \eqn{\Phi^{-1}(1 - \frac{\alpha}{2})}{Phi^-1(1 - alpha/2)}. The
#' values are the fitted means for `"model"` and the partial effects on the
#' link scale of each smooth or learner otherwise; the bootstrap standard
#' errors are multiplied by the factors. A factor is 1 when the values do
#' not vary over the replicates. A factor above 3 or below 0.9 gives a
#' warning that the intervals may not reach their nominal coverage; with
#' `"smooth"` and `"learner"`, one warning per side names the smooths or
#' learners concerned (`$widening_factor` holds every factor).
#' `$Vp` holds the widening for `alpha` only: [smooth_intervals()] recomputes
#' the factors at its own `alpha`, in the same way.
#'
#' ## Dispersion
#'
#' For families with a free dispersion (`gaussian`, `Gamma`,
#' `inverse.gaussian`), the replicate responses are simulated with
#' dispersion `phi`, by default `mod$scale`: Fletcher's (2012) estimate with
#' the sum of the prior weights minus
#' \eqn{\widetilde{\operatorname{edf}}}{edf~} as its residual degrees of
#' freedom (see [gbgam()]). When the model was fitted with
#' `edf_compute_method = "none"`, `mod$scale` is estimated with zero model
#' degrees of freedom and is biased low; when no residual degrees of
#' freedom are left, it is set to a floor of 1e-12. In both cases a warning
#' says so unless `phi` is supplied.
#'
#' ## Families
#'
#' The supported families are `gaussian`, `poisson`, `binomial`, `Gamma` and
#' `inverse.gaussian`. For `binomial`, the prior weights are taken as the
#' number of trials when all of them are integer-valued (the `glm()`
#' convention), so `cbind(successes, failures)` and a proportion response
#' with `weights = trials` give the same bootstrap distribution; otherwise
#' the number of trials is the prior weight divided by the user weight,
#' rounded (at least 1). For `poisson`, a row with prior weight \eqn{w} is
#' simulated as \eqn{Y / w}, with \eqn{Y} a Poisson draw of mean \eqn{w}
#' times the fitted mean, so its variance is the fitted mean divided by
#' \eqn{w}; for an integer \eqn{w} this is the mean of \eqn{w} Poisson
#' draws, and a weighted fit matches a fit on \eqn{w} copies of the row.
#' Rows with zero prior weight do not affect the fit and are kept at their
#' fitted mean.
#'
#' @seealso [gb_param_boot_ctrl()], [gb_uq_ctrl()],
#'   [smooth_intervals()], [gbgam()]
#'
#' @references
#' Fletcher, D. J. (2012).
#' Estimating overdispersion when fitting a generalized linear model to sparse
#' data.
#' *Biometrika, 99*(1), 230--237.
#' \href{https://doi.org/10.1093/biomet/asr083}{doi:10.1093/biomet/asr083}
#'
#' Mayr, A., Hofner, B., & Schmid, M. (2012).
#' The importance of knowing when to stop.
#' *Methods of Information in Medicine, 51*(2), 178--186.
#' \href{https://doi.org/10.3414/ME11-02-0030}{doi:10.3414/ME11-02-0030}
#'
#' Wood, S. N. (2020).
#' Rejoinder on: Inference and computation with generalized additive models and
#' their extensions.
#' *TEST, 29*(2), 354--358.
#' \href{https://doi.org/10.1007/s11749-020-00716-0}{doi:10.1007/s11749-020-00716-0}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' x <- runif(n)
#' df <- data.frame(x = x, y = sin(2 * pi * x) + rnorm(n, sd = 0.4))
#' mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = df)
#'
#' # Default mode "fixed": replicates replay the main fit's selection and
#' # smoothing parameters (use more replicates, e.g. the default 200, in
#' # practice)
#' mod <- bootstrap_gbgam(mod, n_samples = 50)
#' mod$widening_factor
#'
#' # plot() now draws intervals. With these settings they are 95% intervals,
#' # at the level of the widening, including the uncertainty of the intercept
#' plot(mod, seWithMean = TRUE, se = qnorm(0.975))
#'
#' # Each replicate adapts its own smoothing parameters
#' mod_fs <- bootstrap_gbgam(mod, n_samples = 10, mode = "free_smoothing")
#' }
#' @export
bootstrap_gbgam <- function(mod, n_samples = 200L,
                            mode = c("fixed", "free_smoothing",
                                     "free_selection", "free_stopping"),
                            alpha = 0.05,
                            correction = c("model", "smooth", "learner"),
                            phi = NULL,
                            max_iter = NULL,
                            verbose = FALSE,
                            max_runtime_mins = NULL) {

  checkmate::assert_class(mod, "gbgam")
  ctrl <- gb_param_boot_ctrl(n_samples = n_samples, mode = mode,
                             alpha = alpha, correction = correction,
                             phi = phi, max_iter = max_iter)
  n_samples  <- ctrl$n_samples
  mode       <- ctrl$mode
  correction <- ctrl$correction
  check_verbose(verbose)
  check_runtime_mins(max_runtime_mins)
  check_boot_requirements(mode, mod$stopping_method, mod$edf_compute_method)
  if (mode == "free_stopping")
    check_patience_fields(mod, "bootstrap_gbgam")

  m <- mod$mstop
  if (is.null(m))
    stop("bootstrap_gbgam: mod$mstop is missing; was the model fitted?",
         call. = FALSE)
  if (m < 1L) {
    warning("bootstrap_gbgam: mstop = 0 (no boosting iteration improved on ",
            "the intercept-only model); returning the model unchanged, ",
            "without bootstrap uncertainty quantification.",
            call. = FALSE)
    return(mod)
  }
  mod_in <- mod  # returned if the time limit leaves too few replicates
  mod    <- remove_uq_fields(mod)

  # Recorded selection path and smoothing-parameter paths, replayed in the
  # modes that do not redo these steps
  rec <- recorded_paths(mod, m)

  # Unpadded fitted means mu-tilde (fitted() is NA-padded under na.exclude)
  tmu     <- fitted_mu(mod)
  fam     <- mod$family
  phi_eff <- phi %||% mod$scale
  user_w  <- mod$user_weights %||% rep(1, length(tmu))

  free_phi <- is.null(phi) && !(fam$family %in% c("poisson", "binomial"))
  if (free_phi && identical(mod$edf_compute_method, "none"))
    warning("bootstrap_gbgam: the model was fitted with edf_compute_method = ",
            "\"none\", so mod$scale uses 0 model degrees of freedom ",
            "and underestimates the dispersion by a factor (n - edf)/n. ",
            "Refit with gb_edf_ctrl(\"hutchpp\") or \"exact\" for a ",
            "df-corrected dispersion.", call. = FALSE)
  else if (free_phi &&
           !(sum(mod$prior.weights) - (mod$edf.total %||% 0) > 0))
    warning("bootstrap_gbgam: no residual degrees of freedom are left for ",
            "the dispersion (sum of the prior weights <= edf), so mod$scale ",
            "is the floor 1e-12; supply phi.", call. = FALSE)

  s <- replicate_settings(mode, rec$path, rec$rhos, m, mod,
                          post_rhos = rec$post_rhos,
                          max_iter = max_iter)

  one_rep <- replicate_function(mod, s, tmu, user_w, phi_eff)

  if (verbose) {
    # format_done reads n_done (the replicates run, set below) when the bar
    # is done: cli sets pb_current to the total then, also after a time stop
    pb_boot <- cli::cli_progress_bar(
      "Bootstrap",
      total       = n_samples,
      format      = paste0("  {cli::pb_name}  {cli::pb_bar} ",
                           "{cli::pb_current}/{cli::pb_total} ",
                           "({cli::pb_percent})  ETA: {cli::pb_eta}"),
      format_done = paste0("  {cli::pb_name}  {n_done} replicates ",
                           "in {cli::pb_elapsed}"),
      clear       = FALSE
    )
  }

  boot_start <- proc.time()[["elapsed"]]
  n_done     <- 0L
  results    <- vector("list", n_samples)
  for (b in seq_len(n_samples)) {
    if (!is.null(max_runtime_mins) &&
        (proc.time()[["elapsed"]] - boot_start) / 60 >= max_runtime_mins)
      break
    results[[b]] <- one_rep(b)
    n_done <- n_done + 1L
    if (verbose) cli::cli_progress_update(id = pb_boot)
  }
  if (verbose) cli::cli_progress_done(id = pb_boot)

  if (n_done < n_samples) {
    warning(sprintf(
      "bootstrap_gbgam: time limit reached after %d / %d replicates.",
      n_done, n_samples), call. = FALSE)
    if (n_done < 2L) {
      warning("bootstrap_gbgam: fewer than 2 replicates completed; ",
              "returning the model unchanged, without bootstrap results.",
              call. = FALSE)
      return(mod_in)
    }
    results <- results[seq_len(n_done)]
  }

  coef_mat <- do.call(rbind, lapply(results, `[[`, "coefs"))
  dimnames(coef_mat) <- list(NULL, names(mod$coefficients))

  mod$boot_coef <- coef_mat
  # Replicate learner coefficients, one n_done x p_k matrix per learner
  mod$boot_learner_coef <- lapply(seq_along(mod$learners), function(k)
    do.call(rbind, lapply(results, function(r) r$beta[[k]])))
  # Stopping iteration of each replicate, and whether a freely stopping
  # replicate ran into max_iter before its stopping method stopped it
  mod$boot_mstop    <- vapply(results, function(r) as.integer(r$mstop), 1L)
  mod$boot_censored <- vapply(results, function(r) isTRUE(r$capped), TRUE)
  if (any(mod$boot_censored))
    warning(sprintf(paste0(
      "bootstrap_gbgam: %d of %d replicates reached max_iter before their ",
      "stopping method stopped them (see $boot_censored)."),
      sum(mod$boot_censored), length(results)), call. = FALSE)
  mod$Vp <- stats::cov(coef_mat)

  if (!is.null(alpha)) {
    inf <- widen_covariance(mod, alpha, correction)
    mod$Vp              <- inf$Vp
    mod$widening_factor <- inf$factor
    warn_widening_factors(inf$factor, n_done, correction)
  } else {
    mod$widening_factor <- NULL
  }

  mod$n_samples_done    <- n_done
  mod$uq_control <- gb_uq_ctrl("bootstrap", control = ctrl)

  mod
}


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

#' Compute the bootstrap settings of each replicate
#'
#' Returns a named list of control overrides that `refit_replicate()` applies to
#' each bootstrap replicate, from the bootstrap mode.
#'
#' | field | fixed | free_smoothing | free_selection | free_stopping |
#' |-------|-------|----------------|----------------|---------------|
#' | `edf_compute_method` | `"none"` | inherited | inherited | inherited |
#' | `preselected_learner_indices` | `path` | `path` | `NULL` | `NULL` |
#' | `preselected_rhos` | `rhos` | `NULL` | `NULL` | `NULL` |
#' | `preselected_post_rhos` | `post_rhos` | `NULL` | `NULL` | `NULL` |
#' | `stopping_method` | `"none"` | `"none"` | `"none"` | inherited |
#' | `max_iter` | `m` | `m` | `m` | `max_iter`, or inherited |
#' | `iter_no_improv` | `NULL` | `NULL` | `NULL` | inherited |
#' | `iter_no_improv_rel` | `NULL` | `NULL` | `NULL` | inherited |
#'
#' @param mode The bootstrap mode, a name of `boot_modes`.
#' @param path Integer vector; the recorded learner selection path.
#' @param rhos Numeric vector or `NULL`; the recorded `rho` of the selected
#'   learner at each iteration.
#' @param m Integer; `mod$mstop`.
#' @param mod The fitted `"gbgam"` model (used to inherit settings).
#' @param post_rhos Numeric matrix or `NULL`; the recorded `rho` of the
#'   learners in `learners_to_fit_after_each_iter` (one row per iteration,
#'   one column per learner), as returned by `recorded_paths()`.
#' @param max_iter The cap of freely stopping replicates, or `NULL` for the
#'   main fit's `max_iter` (none when `mod$uncapped`).
#'
#' @return A named list of control overrides.
#' @noRd
replicate_settings <- function(mode, path, rhos, m, mod,
                               post_rhos = NULL, max_iter = NULL) {
  steps       <- boot_modes[[mode]]
  free_lambda <- "lambda" %in% steps
  free_sel    <- "selection" %in% steps
  free_stop   <- "stopping" %in% steps
  replay_rho  <- !free_sel && !free_lambda
  # extend_gbgam() sets max_iter of a fit that had no cap (`mod$uncapped`)
  main_cap    <- if (isTRUE(mod$uncapped)) NULL else mod$max_iter
  list(
    edf_compute_method          = if (free_lambda || free_stop)
                                    mod$edf_compute_method else "none",
    preselected_learner_indices = if (free_sel) NULL else path,
    preselected_rhos            = if (replay_rho) rhos else NULL,
    preselected_post_rhos       = if (replay_rho) post_rhos else NULL,
    stopping_method             = if (free_stop) mod$stopping_method
                                  else "none",
    max_iter                    = if (free_stop) max_iter %||% main_cap else m,
    iter_no_improv              = if (free_stop) mod$iter_no_improv else NULL,
    iter_no_improv_rel          = if (free_stop) mod$iter_no_improv_rel
                                  else NULL
  )
}


#' Recorded selection and smoothing-parameter paths of a fitted model
#'
#' Extracts the first `m` iterations of the learner selection path
#' (`metrics$k_opts`), of the `rho` of the selected learners
#' (`metrics$rhos`) and of the `rho` of the post-iteration learners
#' (`metrics$post_rhos`, one column per entry of
#' `learners_to_fit_after_each_iter`). Missing rhos (learners without a
#' smoothing parameter) are replaced by the learner's final rho, or 0.
#'
#' @param mod A fitted `"gbgam"` object.
#' @param m Integer; the number of iterations to replay (`mod$mstop`).
#' @return A list with `path`, `rhos` (`NULL` if no rho was recorded) and
#'   `post_rhos` (`NULL` if there are no post-iteration learners or their
#'   `rho` was not recorded).
#' @noRd
recorded_paths <- function(mod, m) {
  final_rho <- function(k) mod$learners[[k]]$rho %||% 0

  path     <- mod$metrics$k_opts[seq_len(m)]
  rhos_raw <- mod$metrics$rhos[seq_len(m)]
  rhos     <- if (all(is.na(rhos_raw))) NULL else rhos_raw
  for (i in which(is.na(rhos))) rhos[i] <- final_rho(path[i])

  post_rhos <- mod$metrics$post_rhos
  post_idx  <- mod$learners_to_fit_after_each_iter
  if (is.null(post_rhos) || ncol(post_rhos) == 0L ||
      ncol(post_rhos) != length(post_idx)) {
    post_rhos <- NULL
  } else {
    post_rhos <- post_rhos[seq_len(m), , drop = FALSE]
    for (j in seq_along(post_idx)) {
      na_j <- is.na(post_rhos[, j])
      post_rhos[na_j, j] <- final_rho(post_idx[j])
    }
  }
  list(path = path, rhos = rhos, post_rhos = post_rhos)
}


#' One bootstrap replicate as a function of its index
#'
#' Built from the model without the fields the replicates do not need: the
#' running state, the resampling settings and results, the fitted-value
#' path and the uncertainty quantification settings.
#'
#' @param mod The fitted model, without uncertainty quantification.
#' @param s Control overrides from `replicate_settings()`.
#' @param tmu,user_w,phi_eff Fitted means, user weights and dispersion of
#'   the simulation.
#' @return A function of the replicate index `b`.
#' @noRd
replicate_function <- function(mod, s, tmu, user_w, phi_eff) {
  force(s)
  force(tmu)
  force(user_w)
  force(phi_eff)
  mod[c("resume", "cv", "resampling", "fitted_values",
        "uq_control")] <- list(NULL)
  fam <- mod$family
  function(b) {
    y_b  <- simulate_response(fam, tmu, mod$prior.weights, user_w, phi_eff)
    bmod <- tryCatch(refit_replicate(mod, y_b, s), error = function(e) stop(
      "Bootstrap replicate ", b, " failed (the main fit is not affected): ",
      conditionMessage(e), call. = FALSE))
    list(
      coefs  = bmod$coefficients,
      beta   = lapply(bmod$learners, `[[`, "beta"),
      mstop  = bmod$mstop,
      capped = isTRUE(bmod$metrics$stopping_conditions[["mstop_reached"]]) &&
                 !identical(s$stopping_method, "none")
    )
  }
}


#' Remove the uncertainty quantification attached to a model
#'
#' Used before the uncertainty quantification is recomputed, so that the
#' model carries no covariance that belongs to other coefficients or to
#' another method.
#'
#' @param mod A `"gbgam"` model.
#' @return `mod` without `$Vp` and the bootstrap and Bayesian backfitting
#'   results, and with `$uq_control` set to
#'   `gb_uq_ctrl("none")`, which describes it.
#' @noRd
remove_uq_fields <- function(mod) {
  fields <- c("Vp", "Ve", "boot_coef", "boot_learner_coef", "boot_mstop",
              "boot_censored", "widening_factor", "coefficients_bayes",
              "beta_samples", "mh_accept_rate", "bb_info", "n_samples_done")
  mod[intersect(fields, names(mod))] <- NULL
  mod$uq_control <- gb_uq_ctrl("none")
  mod
}


#' Warn about widening factors far from 1
#'
#' With `"model"`, one factor: a warning when it is above 3 or below 0.9.
#' With `"smooth"` or `"learner"`, one factor per group: at most one warning
#' per side, naming the groups beyond the threshold with their factors (the
#' 5 farthest from it, then how many more).
#'
#' @param fac The widening factor(s), named by group unless `"model"`.
#' @param n_done The number of replicates they are computed from.
#' @param correction `"model"`, `"smooth"` or `"learner"`.
#' @noRd
warn_widening_factors <- function(fac, n_done, correction) {
  c_warn_hi <- 3    # far above 1
  c_warn_lo <- 0.9  # values slightly below 1 are Monte Carlo noise
  why_hi <- paste0("this typically indicates insufficient undersmoothing ",
                   "during fitting, and the intervals may not reach their ",
                   "nominal coverage.")
  why_lo <- sprintf(paste0("this may indicate a poor model fit (or too few ",
                           "replicates, n_samples = %d), and the intervals ",
                           "may not reach their nominal coverage."), n_done)
  if (correction == "model") {
    if (fac > c_warn_hi)
      warning(sprintf("bootstrap_gbgam: widening factor = %.2f > %g: %s",
                      fac, c_warn_hi, why_hi), call. = FALSE)
    if (fac < c_warn_lo)
      warning(sprintf("bootstrap_gbgam: widening factor = %.2f < %g: %s",
                      fac, c_warn_lo, why_lo), call. = FALSE)
    return(invisible(NULL))
  }
  groups <- if (correction == "smooth") {
    "smooths and parametric coefficients"
  } else {
    "learners"
  }
  # "name (factor)" for the groups beyond a threshold, farthest first
  list_groups <- function(f) {
    max_shown <- 5L
    f     <- f[order(-abs(log(f)))]
    shown <- sprintf("%s (%.2f)", names(f), f)[seq_len(min(length(f),
                                                           max_shown))]
    more  <- length(f) - length(shown)
    paste0(paste(shown, collapse = ", "),
           if (more > 0L) sprintf(" and %d more", more))
  }
  hi <- fac[fac > c_warn_hi]
  lo <- fac[fac < c_warn_lo]
  if (length(hi) > 0L)
    warning(sprintf(
      "bootstrap_gbgam: widening factors > %g for %d of %d %s: %s; %s",
      c_warn_hi, length(hi), length(fac), groups, list_groups(hi), why_hi),
      call. = FALSE)
  if (length(lo) > 0L)
    warning(sprintf(
      "bootstrap_gbgam: widening factors < %g for %d of %d %s: %s; %s",
      c_warn_lo, length(lo), length(fac), groups, list_groups(lo), why_lo),
      call. = FALSE)
  invisible(NULL)
}


#' Model matrix at the observations used in the fit
#'
#' `predict(type = "lpmatrix")` is padded with `NA` rows for the observations
#' excluded by `na.action = na.exclude`; this drops them, so that the rows
#' match `mod$offset` and the fitted values.
#'
#' @param mod A fitted `"gbgam"` model.
#' @return The `n x p` model matrix.
#' @noRd
train_lpmatrix <- function(mod) {
  X <- stats::predict(mod, type = "lpmatrix")
  if (inherits(mod$na.action, "exclude")) {
    X <- X[-mod$na.action, , drop = FALSE]
  }
  X
}


#' Refit a model to a new response with mode-specific overrides
#'
#' Resets `mod` to its pre-fit state for `y_new` (see `reset_gbgam()`),
#' applies the control overrides from `replicate_settings()` and refits
#' quietly. This is one bootstrap replicate of [bootstrap_gbgam()].
#'
#' @param mod A fitted `"gbgam"` object.
#' @param y_new New response vector.
#' @param s A list of control overrides from `replicate_settings()`.
#' @return The refitted model.
#' @noRd
refit_replicate <- function(mod, y_new, s) {
  bmod <- reset_gbgam(mod, y_new, edf_compute_method = s$edf_compute_method)
  bmod["preselected_learner_indices"] <- list(s$preselected_learner_indices)
  bmod["preselected_rhos"]            <- list(s$preselected_rhos)
  bmod["preselected_post_rhos"]       <- list(s$preselected_post_rhos)
  bmod$stopping_method                <- s$stopping_method
  bmod["max_iter"]                    <- list(s$max_iter)
  bmod["iter_no_improv"]              <- list(s$iter_no_improv)
  bmod["iter_no_improv_rel"]          <- list(s$iter_no_improv_rel)
  bmod["max_runtime_mins"]            <- list(NULL)
  # Replicates are silent (e.g. the ASR lambda-to-zero message)
  bmod$verbose                        <- FALSE
  fit_gbgam(bmod, start_time = Sys.time(), save_fitted_values = FALSE,
            verbose = FALSE)
}


#' Reset a fitted `gbgam` model to its pre-fit state for a new response vector
#'
#' Reuses the existing learner structure (design matrices, DR decompositions,
#' etc.) without rebuilding any basis or decompositions, keeping bootstrap
#' refitting cheap. The returned object is ready to pass directly to
#' `fit_gbgam()`; mode-specific overrides (selection path, `rho` values,
#' stopping settings) are applied by the caller via
#' `replicate_settings()`.
#'
#' @param mod A fitted `"gbgam"` object.
#' @param y_new New response vector (e.g. a bootstrap-simulated response).
#' @param edf_compute_method EDF computation method to use for this replicate.
#'   Pass `"none"` to skip derivative tracking when replaying a fixed path with
#'   fixed rhos (the `"fixed"` bootstrap mode); otherwise pass
#'   `mod$edf_compute_method` to inherit the original setting.
#'
#' @return A model list ready to pass to `fit_gbgam()`.
#' @noRd
reset_gbgam <- function(mod, y_new, edf_compute_method) {
  n <- length(y_new)

  # Reset learner running state (beta, edf, eta, and the smoothing state rho,
  # nu_k, lambda_at_zero that ASR updates); preserve structural fields
  # (design, A, D, xtWx, rho0, nu_k0, ...).
  mod$learners <- lapply(mod$learners, function(lr) {
    lr$beta <- rep(0, length(lr$beta))
    lr$edf  <- 0
    lr$eta  <- rep(0L, nrow(lr$design))
    if (inherits(lr, "demmler_reinsch")) {
      lr$rho            <- lr$rho0
      lr$nu_k           <- initial_nu_k(lr)
      lr$lambda_at_zero <- FALSE
    }
    lr
  })

  # New response and intercept initialisation, mirroring gbgam():
  # initialise_intercept = TRUE starts the intercept learner at the MLE ahat
  # of y_new; FALSE starts it at zero.
  mod$y    <- y_new
  init     <- isTRUE(mod$initialise_intercept)
  ahat     <- find_intercept_mle(
    y             = y_new,
    family        = mod$family,
    prior.weights = mod$prior.weights,
    offset        = mod$offset
  )
  a0       <- if (init) ahat else 0
  mod$intercept0          <- ahat
  mod$learners[[1L]]$beta <- a0
  mod$linear.predictors   <- a0 + mod$offset

  # Starting edf~ when EDF is not tracked, as in gbgam(): the intercept
  # learner gets 0 and the (unused) running total 1 or 0. When EDF is tracked,
  # init_derivs() overwrites both with tr(dmu~0/dy).
  mod$learners[[1L]]$edf <- 0
  mod$edf.total          <- if (init) 1 else 0

  # Apply EDF method override and re-initialise derivatives (quietly: the
  # non-canonical link message was already shown for the main fit).
  mod$edf_compute_method <- edf_compute_method
  mod <- init_derivs(mod, n, quiet = TRUE)

  mod
}


#' Simulate one response vector from a family at given fitted means.
#'
#' @param family The model family (one of the supported families).
#' @param mu Fitted means at which to simulate (length n).
#' @param prior_weights mod$prior.weights  (= n_trials * user_weights for
#'  cbind binomial, = user_weights otherwise). For binomial, integer-valued
#'  prior weights are used as the number of trials (glm/`mgcv` convention).
#' @param user_weights `mod$user_weights`, the raw user-supplied weights
#'  vector (`mod$weights` holds the working weights); only used for
#'  binomial with non-integer prior weights, where the number of trials is
#'  prior_weights / user_weights.
#' @param phi Dispersion of the families with a free dispersion.
#' @return The simulated response (length n), on the scale of the fitted
#'   response (proportions for binomial).
#' @noRd
simulate_response <- function(family, mu, prior_weights, user_weights, phi) {
  n     <- length(mu)
  tiny  <- 1e-12  # floor of weights, phi and mu where they divide
  p_lim <- 1e-10  # keeps binomial probabilities inside (0, 1)
  y <- switch(family$family,
    gaussian = {
      sd_vec <- sqrt(phi / pmax(prior_weights, tiny))
      stats::rnorm(n, mean = mu, sd = sd_vec)
    },
    poisson = {
      # Poisson(w_i mu) / w_i: E(y) = mu and Var(y) = mu / w_i (the model's
      # working variance) for any w_i > 0; for an integer w_i, the mean of
      # w_i draws. Zero-weight rows are reset below.
      w_sim <- ifelse(prior_weights > 0, prior_weights, 1)
      stats::rpois(n, lambda = w_sim * mu) / w_sim
    },
    binomial = {
      # glm/mgcv convention: integer-valued prior weights are the number of
      # trials (cbind(s, f): n_i * w_i; proportion y with weights = n_i;
      # 0/1 y with frequency weights). Otherwise fall back to the
      # prior/user ratio (e.g. binary data with fractional sampling weights).
      n_trials <- if (checkmate::test_integerish(prior_weights, tol = 1e-8)) {
        pmax(round(prior_weights), 1)
      } else {
        round(pmax(prior_weights / pmax(user_weights, tiny), 1))
      }
      prob <- pmin(pmax(mu, p_lim), 1 - p_lim)
      stats::rbinom(n, size = n_trials, prob = prob) / n_trials
    },
    Gamma = {
      # Var(Y) = mu^2 * phi / prior_weights, so shape = prior_weights / phi
      shape <- pmax(prior_weights, tiny) / pmax(phi, tiny)
      stats::rgamma(n, shape = shape, rate = shape / pmax(mu, tiny))
    },
    inverse.gaussian = {
      # Var(Y) = mu^3 * phi / prior_weights, so lambda = prior_weights / phi
      # (as in stats::inverse.gaussian()$simulate)
      lambda <- pmax(prior_weights, tiny) / pmax(phi, tiny)
      r_inverse_gaussian(n, mu, lambda)
    },
    stop("bootstrap_gbgam: family \"", family$family,
         "\" is not supported.", call. = FALSE)
  )
  # Zero-weight rows carry no information; keep them at the fitted mean so
  # that the simulated response stays inside the family's support (a floored
  # shape of 1e-12 otherwise gives y = 0 for Gamma / inverse Gaussian, and a
  # NaN stopping loss in the "free_stopping" mode).
  zero_w    <- prior_weights <= 0
  y[zero_w] <- mu[zero_w]
  y
}

#' Draw from the inverse Gaussian distribution
#'
#' The algorithm of Michael, Schucany and Haas (1976, The American
#' Statistician, 30(2), 88--90): with \eqn{y = v^2}, \eqn{v} standard
#' normal, the smaller root \eqn{x_1} of
#' \eqn{\lambda (x - \mu)^2 / (\mu^2 x) = y} is returned with probability
#' \eqn{\mu / (\mu + x_1)}, and the larger root \eqn{\mu^2 / x_1} otherwise.
#' The smaller root is computed as \eqn{\mu / (1 + r + \sqrt{r (2 + r)})},
#' \eqn{r = \mu y / (2 \lambda)}, which equals the textbook
#' \eqn{\mu (1 + r - \sqrt{r (2 + r)})} without its cancellation when
#' \eqn{\mu / \lambda} is large.
#'
#' @param n Number of draws.
#' @param mu Means (length n).
#' @param lambda Shape parameters (length n): the variance is
#'   \eqn{\mu^3 / \lambda}.
#' @return The draws (length n).
#' @noRd
r_inverse_gaussian <- function(n, mu, lambda) {
  y  <- stats::rnorm(n)^2
  r  <- mu * y / (2 * lambda)
  x1 <- mu / (1 + r + sqrt(r * (2 + r)))
  u  <- stats::runif(n)
  ifelse(u <= mu / (mu + x1), x1, mu^2 / x1)
}


#' Widening factor of one fitted quantity, by Wood's (2020) rule
#'
#' The quantile over observations and replicates of
#' \eqn{|t^{(r)}_i - \tilde t_i| / \tilde s_i}, divided by
#' \eqn{\Phi^{-1}(1 - \frac{\alpha}{2})}; 1 when no observation varies over
#' the replicates. The intervals
#' \eqn{\tilde t_i \pm c \, \Phi^{-1}(1 - \frac{\alpha}{2}) \tilde s_i} then
#' hold a share 1 - alpha of the replicates' values, over all observations.
#'
#' @param t0 The main fit's values at the data (length n): fitted means or
#'   a partial effect.
#' @param tr The replicates' values (replicates x n).
#' @param alpha The miscoverage, in (0, 0.5].
#' @return The factor, a scalar.
#' @noRd
wood_widening_factor <- function(t0, tr, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  s <- apply(tr, 2L, stats::sd)
  if (all(s == 0)) return(1)
  keep <- s > 0
  d    <- abs(sweep(tr[, keep, drop = FALSE], 2L, t0[keep], "-"))
  d    <- sweep(d, 2L, s[keep], "/")
  stats::quantile(d, 1 - alpha, names = FALSE) / z
}

#' Widen the bootstrap covariance
#'
#' For `"model"`, Wood's (2020) factor from the replicates' fitted means
#' (computed from `$boot_coef`), and `Vp` is the replicates' covariance
#' times its square. Otherwise the same rule gives one factor per group of
#' coefficients: each smooth (and each parametric coefficient, the intercept
#' included) for `"smooth"`, each learner for `"learner"`. The factor of
#' group \eqn{g} at level \eqn{1 - \alpha} is
#' \deqn{c_g = q_{1-\alpha}\left(|t_{g,i}^{(r)} - \tilde t_{g,i}| /
#'   \tilde s_{g,i}\right) / \Phi^{-1}(1 - \frac{\alpha}{2}),}
#' the quantile over observations \eqn{i} and replicates \eqn{r} of the
#' distance of each replicate's partial effect of the group (on the link
#' scale, at the data) from the main fit's, in units of its bootstrap
#' standard deviation, relative to the normal quantile; 1 for a group
#' without bootstrap variation. The replicate deviations from the main fit
#' are scaled by these factors and `Vp` is their covariance on the `mgcv`
#' scale, so that `predict(se.fit = TRUE)` uses the widened standard errors.
#'
#' @param mod A model returned by [bootstrap_gbgam()] (with `$boot_coef` and
#'   `$boot_learner_coef`).
#' @param alpha The miscoverage, in (0, 0.5].
#' @param correction `"model"`, `"smooth"` or `"learner"`.
#' @return list(Vp, factor); `factor` is a number for `"model"`, otherwise
#'   named by group.
#' @noRd
widen_covariance <- function(mod, alpha, correction) {
  coef0 <- mod$coefficients
  B     <- nrow(mod$boot_coef)
  if (correction == "model") {
    X     <- train_lpmatrix(mod)
    off   <- mod$offset %||% 0
    tmu   <- mod$family$linkinv(as.vector(X %*% coef0) + off)
    tmu_r <- mod$family$linkinv(sweep(mod$boot_coef %*% t(X), 2L, off, "+"))
    fac   <- wood_widening_factor(tmu, tmu_r, alpha)
    Vp    <- fac^2 * stats::cov(mod$boot_coef)
  } else if (correction == "smooth") {
    X     <- train_lpmatrix(mod)
    n_par <- if (length(mod$smooth) > 0L) mod$smooth[[1L]]$first.para - 1L
             else length(coef0)
    groups <- c(as.list(seq_len(n_par)),
                lapply(mod$smooth, function(sm) sm$first.para:sm$last.para))
    names(groups) <- c(names(coef0)[seq_len(n_par)],
                       vapply(mod$smooth, `[[`, "", "label"))
    fac <- vapply(groups, function(i) {
      Xi <- X[, i, drop = FALSE]
      wood_widening_factor(as.vector(Xi %*% coef0[i]),
                           mod$boot_coef[, i, drop = FALSE] %*% t(Xi),
                           alpha)
    }, numeric(1))
    fac_coef <- rep(fac, lengths(groups))[order(unlist(groups))]
    scaled   <- sweep(sweep(mod$boot_coef, 2L, coef0, "-"), 2L, fac_coef, "*")
    Vp       <- stats::cov(sweep(scaled, 2L, coef0, "+"))
  } else {
    L   <- mod$learners
    fac <- vapply(seq_along(L), function(k) {
      Xk <- L[[k]]$design
      wood_widening_factor(as.vector(Xk %*% L[[k]]$beta),
                           mod$boot_learner_coef[[k]] %*% t(Xk), alpha)
    }, numeric(1))
    names(fac) <- vapply(L, `[[`, "", "var_label")
    coefs <- t(vapply(seq_len(B), function(r) {
      learners_r <- L
      for (k in seq_along(L))
        learners_r[[k]]$beta <- L[[k]]$beta +
          fac[k] * (mod$boot_learner_coef[[k]][r, ] - L[[k]]$beta)
      learners_to_mgcv_coef(learners_r, mod$coef_names, mod$ind0var,
                            mod$cen_const, cen_absorb = mod$cen_absorb)
    }, numeric(length(coef0))))
    Vp <- stats::cov(coefs)
  }
  dimnames(Vp) <- list(names(coef0), names(coef0))
  list(Vp = Vp, factor = fac)
}
