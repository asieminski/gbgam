#' Fit a `gbgam`
#'
#' Sets up the fitting state, runs it until a stopping condition is met, and
#' returns the best model.
#'
#' @param mod A fully assembled model list (see `gbgam()`).
#' @param start_time Start of the fit, for time limits and progress.
#' @param save_fitted_values Logical; keep the fitted means of every
#'   iteration.
#' @param verbose Logical; show progress.
#' @param resumable Logical; attach the running state as `$resume`, so that
#'   `extend_gbgam()` can continue the fit.
#' @noRd
fit_gbgam <- function(mod, start_time, save_fitted_values, verbose,
                      resumable = FALSE) {
  st <- fit_state_init(mod, start_time, save_fitted_values, verbose)
  st <- fit_state_run(st)
  fit_state_result(st, resumable = resumable)
}


#' Initialise the state of a boosting fit
#'
#' The state holds the running model `mod`, the best model so far `opt_mod`
#' and the per-iteration metrics. `fit_state_run()` advances it and
#' `fit_state_result()` turns it into the returned model.
#'
#' @param mod A fully assembled model list (see `gbgam()`).
#' @param start_time Start of the fit, for time limits and progress.
#' @param save_fitted_values Logical; keep the fitted means of every
#'   iteration.
#' @param verbose Logical; show progress.
#' @param record_history Logical; record the state of the learners updated
#'   in each iteration, so that `fit_state_rewind()` can return the model at
#'   any earlier iteration (used by cross-validated stopping).
#' @param monitor `NULL` or a function of the running model returning a
#'   numeric vector, evaluated at the start and after every iteration; the
#'   values are kept in `$trace` (one element per iteration, the first for
#'   iteration 0).
#' @return A list (the fitting state).
#' @noRd
fit_state_init <- function(mod, start_time, save_fitted_values, verbose,
                           record_history = FALSE, monitor = NULL) {

  required_objects <-
    c("stopping_method",
      "learners")
  potentially_null_objects <-
    c("preselected_learner_indices",
      "preselected_rhos",
      "learners_to_fit_after_each_iter",
      "AIC_k",
      "max_iter",
      "iter_no_improv",
      "iter_no_improv_rel",
      "max_runtime_mins"
     )
  stopifnot(all(c(required_objects, potentially_null_objects) %in% names(mod)))
  psel <- mod$preselected_learner_indices
  if (!is.null(psel) &&
      (is.null(mod$max_iter) || mod$max_iter > length(psel))) {
    stop("preselected_learner_indices has length ", length(psel),
         " but max_iter is ",
         if (is.null(mod$max_iter)) "not set" else mod$max_iter,
         ". Set gb_stop_ctrl(max_iter = ...) to at most ",
         "length(preselected_learner_indices).", call. = FALSE)
  }
  edf0 <- mod$edf.total # Save this to return later

  # The n x n (exact) or n x z (hutchpp) derivative matrices are only needed
  # by the running `mod`; the cached best model drops them unless the
  # internal option gbgam.keep_derivs = TRUE (used in tests) is set.
  keep_derivs <- isTRUE(getOption("gbgam.keep_derivs", FALSE))

  # Initialise return model
  opt_mod        <- mod
  if (!keep_derivs) opt_mod$derivs <- NULL
  opt_mod$mstop  <- 0L  # will be overwritten once any iteration improves

  # Additional results
  ls0 <- compute_loss(mod, mod$stopping_method, AIC_k=mod$AIC_k)
  check_loss_finite(ls0, m = 0L)
  u0 <- pseudo_residuals(mod$y, mod$linear.predictors, mod$family)
  if (!all(is.finite(u0[mod$prior.weights > 0]))) {
    stop("The negative gradient is not finite at the starting linear ",
         "predictor: the starting mean lies on the boundary of the \"",
         mod$family$link, "\" link (e.g. mu = 0 when initialise_intercept = ",
         "FALSE). Use initialise_intercept = TRUE or a valid offset.",
         call. = FALSE)
  }

  # Progress display state (NULL when verbose = FALSE)
  progress <- gb_progress_init(verbose, mod, start_time)

  # Post-iteration learners and (optional, internal) replay of their rho
  # path: mod$preselected_post_rhos has one row per iteration and one column
  # per entry of learners_to_fit_after_each_iter (NA leaves rho unchanged)
  post_idx <- mod$learners_to_fit_after_each_iter %||% integer(0)
  pp_rhos  <- mod$preselected_post_rhos
  if (!is.null(pp_rhos) &&
      (!is.matrix(pp_rhos) || ncol(pp_rhos) != length(post_idx) ||
         is.null(mod$max_iter) || nrow(pp_rhos) < mod$max_iter)) {
    stop("preselected_post_rhos must be a matrix with one column per entry ",
         "of learners_to_fit_after_each_iter and at least max_iter rows.",
         call. = FALSE)
  }

  # Initialise iteration metrics (preallocated when max_iter is known; they
  # grow as needed otherwise, e.g. when a fit is extended)
  n_alloc <- mod$max_iter %||% 0L
  fitted_values <- NULL
  if (save_fitted_values) {
    fitted_values <- vector("list", n_alloc + 1L)
    fitted_values[[1L]] <- fitted_mu(opt_mod)
  }

  list(
    mod                 = mod,
    opt_mod             = opt_mod,
    m                   = 0,
    counter             = 0,
    ls0                 = ls0,
    min_loss            = ls0,
    edf0                = edf0,
    keep_derivs         = keep_derivs,
    save_fitted_values  = save_fitted_values,
    start_time          = start_time,
    progress            = progress,
    post_idx            = post_idx,
    k_opts              = rep(NA_integer_, n_alloc),
    rhos                = rep(NA_real_,    n_alloc),
    ls                  = rep(NA_real_,    n_alloc),
    edfs                = rep(NA_real_,    n_alloc),
    post_rhos           = matrix(NA_real_, nrow = n_alloc,
                                 ncol = length(post_idx)),
    fitted_values       = fitted_values,
    stopping_conditions = NULL,
    record_history      = record_history,
    history             = list(),
    init_light          = lapply(mod$learners, learner_light),
    beta0               = lapply(mod$learners, `[[`, "beta"),
    monitor             = monitor,
    trace               = if (is.null(monitor)) NULL else list(monitor(mod))
  )
}


#' Run boosting iterations on a fitting state
#'
#' @param st A fitting state from `fit_state_init()`.
#' @param n_iter Maximum number of iterations to run in this call.
#' @param use_stopping Logical; stop as soon as a stopping condition of the
#'   model (the patience, `max_iter`, `max_runtime_mins`) is met. With
#'   `FALSE`, exactly `n_iter` iterations are run and the conditions are only
#'   recorded.
#' @return The updated state.
#' @noRd
fit_state_run <- function(st, n_iter = Inf, use_stopping = TRUE) {
  mod           <- st$mod
  opt_mod       <- st$opt_mod
  min_loss      <- st$min_loss
  counter       <- st$counter
  ls0           <- st$ls0
  edf0          <- st$edf0
  k_opts        <- st$k_opts
  rhos          <- st$rhos
  ls            <- st$ls
  edfs          <- st$edfs
  post_rhos     <- st$post_rhos
  fitted_values <- st$fitted_values
  post_idx      <- st$post_idx
  pp_rhos       <- mod$preselected_post_rhos
  stopping_conditions <- st$stopping_conditions

  target <- st$m + n_iter
  m      <- st$m
  while (m < target) {
    m <- m + 1

    # Iterate model -----------------------------------------------------------
    # Fit main base learner. mod$excluded_learner_indices (if any) are
    # respected inside iterate_mod during data-driven selection.
    mod <- iterate_mod(mod, mod$preselected_learner_indices[m], iter = m)

    # (Optional) Set rho (and the dependent nu_k / lambda_at_zero) to a
    # prespecified value
    if(!is.null(mod$preselected_rhos)){
      mod$learners[[mod$k_opt]] <- set_learner_rho(
        mod$learners[[mod$k_opt]],
        new_rho           = mod$preselected_rhos[m],
        target_tr2AminAA = mod$target_tr2AminAA
      )
    }

    # Record metrics for the main selection before any forced post-iter updates.
    k_opts[m] <- mod$k_opt
    rhos[m]   <- mod$learners[[mod$k_opt]]$rho %||% NA_real_

    post      <- fit_post_learners(mod, m, post_idx, pp_rhos, post_rhos)
    mod       <- post$mod
    post_rhos <- post$post_rhos

    # Compute loss after all updates
    new_loss <- compute_loss(mod, mod$stopping_method, AIC_k=mod$AIC_k)
    check_loss_finite(new_loss, m = m)

    # Save model-level metrics --------------------------------------------
    if(st$save_fitted_values){
      fitted_values[[m+1]] <- fitted_mu(mod)
    }
    if(!is.null(new_loss)) ls[m] <- new_loss
    if(mod$edf_compute_method!="none") edfs[m] <- mod$edf.total
    if (st$record_history) {
      st$history[[m]] <- learner_history_entry(
        mod, unique(c(k_opts[m], post_idx,
                      if (isTRUE(mod$intercept_padding)) 1L)))
    }
    if (!is.null(st$monitor)) st$trace[[m + 1]] <- st$monitor(mod)

    # Progress display -----------------------------------------------------
    elapsed_time <- as.numeric(difftime(Sys.time(), st$start_time,
                                        units = "mins"))
    gb_progress_update(
      state        = st$progress,
      m            = m,
      mod          = mod,
      new_loss     = new_loss,
      min_loss     = min_loss,
      opt_m        = opt_mod$mstop,
      k_opts_m     = k_opts[seq_len(m)],
      ls_all       = c(ls0, ls[seq_len(m)]),
      edfs_all     = if (mod$edf_compute_method != "none")
                       c(edf0, edfs[seq_len(m)]) else NULL,
      elapsed_mins = elapsed_time
    )

    # Cache the model if it's the new best (or the latest if loss not tracked)
    best     <- update_best(mod, m, new_loss, min_loss, opt_mod, counter,
                            st$keep_derivs)
    opt_mod  <- best$opt_mod
    min_loss <- best$min_loss
    counter  <- best$counter

    # Check stopping conditions -------------------------------------------
    stopping_conditions <- c(
      no_improvements = !is.null(mod$iter_no_improv) &&
        (counter >= patience(mod, opt_mod$mstop)),
      # TRUE when max_iter was hit
      mstop_reached = !is.null(mod$max_iter) && (m==mod$max_iter),
      time_limit_exceeded = !is.null(mod$max_runtime_mins) &&
        (elapsed_time > mod$max_runtime_mins)
    )
    if (use_stopping && any(stopping_conditions)) break
  }

  st$mod                 <- mod
  st$opt_mod             <- opt_mod
  st$min_loss            <- min_loss
  st$counter             <- counter
  st$k_opts              <- k_opts
  st$rhos                <- rhos
  st$ls                  <- ls
  st$edfs                <- edfs
  st$post_rhos           <- post_rhos
  st$fitted_values       <- fitted_values
  st$stopping_conditions <- stopping_conditions
  st$m                   <- m
  st
}


#' Force-fit the post-iteration learners of iteration m
#'
#' Fits each learner in `post_idx` (e.g. intercept, random effects) after
#' the main selection, optionally replays its rho path (`pp_rhos`), and
#' records its rho in row m of `post_rhos` (NA without a smoothing
#' parameter).
#'
#' @param mod Object with class `gbgam`.
#' @param m Iteration number.
#' @param post_idx Indices of the post-iteration learners.
#' @param pp_rhos `mod$preselected_post_rhos` or `NULL`.
#' @param post_rhos Matrix of recorded rhos, one row per iteration; grown by
#'   a row when m exceeds its rows.
#' @return A list with the updated `mod` and `post_rhos`.
#' @noRd
fit_post_learners <- function(mod, m, post_idx, pp_rhos, post_rhos) {
  if (m > nrow(post_rhos)) {
    post_rhos <- rbind(post_rhos,
                       matrix(NA_real_, nrow = 1L, ncol = length(post_idx)))
  }
  for (j in seq_along(post_idx)) {
    k_post <- post_idx[j]
    mod    <- iterate_mod(mod, k_opt = k_post, iter = m)
    if (!is.null(pp_rhos) && !is.na(pp_rhos[m, j])) {
      mod$learners[[k_post]] <- set_learner_rho(
        mod$learners[[k_post]],
        new_rho           = pp_rhos[m, j],
        target_tr2AminAA = mod$target_tr2AminAA
      )
    }
    post_rhos[m, j] <- mod$learners[[k_post]]$rho %||% NA_real_
  }
  list(mod = mod, post_rhos = post_rhos)
}


#' Patience of loss-based stopping
#'
#' The number of iterations without a new best loss after which boosting
#' stops: `max(iter_no_improv, ceiling(iter_no_improv_rel * m_best))`.
#'
#' @param mod Object with class `gbgam` (its `iter_no_improv` and
#'   `iter_no_improv_rel`, both set).
#' @param m_best The best iteration so far.
#' @return A whole number.
#' @noRd
patience <- function(mod, m_best) {
  stopifnot(!is.null(mod$iter_no_improv), !is.null(mod$iter_no_improv_rel))
  # Rounding before ceiling() keeps e.g. 0.07 * 100 at 7, not 8
  rel <- ceiling(round(mod$iter_no_improv_rel * m_best, 8))
  max(mod$iter_no_improv, rel)
}


#' Stop on a model fitted before `iter_no_improv_rel` existed
#'
#' @param mod Object with class `gbgam`.
#' @param caller Name of the calling function, for the message.
#' @return `NULL`, invisibly; stops if `mod` has no `iter_no_improv_rel`.
#' @noRd
check_patience_fields <- function(mod, caller) {
  if (!"iter_no_improv_rel" %in% names(mod)) {
    stop(caller, ": the model was fitted with an older gbgam that has no ",
         "iter_no_improv_rel; refit it.", call. = FALSE)
  }
  invisible(NULL)
}


#' Keep the model of iteration m if it is the best so far
#'
#' @param mod Object with class `gbgam` after iteration m.
#' @param m Iteration number.
#' @param new_loss,min_loss Its stopping loss and the best so far (`NULL`
#'   when no loss is tracked: every iteration is then the best).
#' @param opt_mod The best model so far.
#' @param counter Iterations since the last best.
#' @param keep_derivs Logical; keep the derivatives in the cached model.
#' @return A list with the updated `opt_mod`, `min_loss` and `counter`.
#' @noRd
update_best <- function(mod, m, new_loss, min_loss, opt_mod, counter,
                        keep_derivs) {
  # New best unless the loss rose: a tie moves mstop to the later iteration
  if (!is.null(new_loss) && (new_loss > min_loss)) {
    return(list(opt_mod = opt_mod, min_loss = min_loss,
                counter = counter + 1))
  }
  opt_mod <- mod
  if (!keep_derivs) opt_mod$derivs <- NULL
  opt_mod$mstop <- m
  list(opt_mod = opt_mod, min_loss = new_loss, counter = 0)
}


#' Return the best model of a fitting state
#'
#' @param st A fitting state.
#' @param resumable Logical; attach the running state as `$resume`.
#' @return The best model (`st$opt_mod`) with its metrics and coefficients.
#' @noRd
fit_state_result <- function(st, resumable = FALSE) {
  opt_mod  <- st$opt_mod
  mod      <- st$mod
  ls0      <- st$ls0
  ls       <- st$ls
  edfs     <- st$edfs
  k_opts   <- st$k_opts
  min_loss <- st$min_loss
  edf0     <- st$edf0
  # A preallocated loss path (max_iter known, no loss tracked) keeps one
  # entry per iteration when the fit is extended
  if (length(ls) > 0L && length(ls) < length(k_opts)) {
    length(ls) <- length(k_opts)
  }

  # Save all metrics
  opt_mod$metrics <- list(
    stopping_conditions = st$stopping_conditions,
    k_opts = k_opts,
    rhos = st$rhos,
    # rho of each post-iteration learner after each iteration (rows
    # aligned with rhos; NA for learners without a smoothing parameter)
    post_rhos = st$post_rhos,
    ls = c(ls0, ls)
  )
  if(mod$edf_compute_method!="none"){
    opt_mod$edfs <- c(edf0, edfs) # Insert 0th iter edf
  }
  if(st$save_fitted_values){
    opt_mod$fitted_values <- do.call(what = cbind, st$fitted_values)
  }
  opt_mod$elapsed_secs <- as.numeric(difftime(Sys.time(), st$start_time,
                                              units = "secs"))
  opt_mod$coefficients <- learners_to_mgcv_coef(
    learners=opt_mod$learners,
    coef_names=opt_mod$coef_names,
    ind0var=opt_mod$ind0var,
    cen_const=opt_mod$cen_const,
    cen_absorb=opt_mod$cen_absorb)

  # Force-display the final state at the optimal iteration so the user
  # always sees the true stopping point, not the last throttled snapshot.
  gb_progress_update(
    state        = st$progress,
    m            = opt_mod$mstop,
    mod          = opt_mod,
    new_loss     = min_loss,
    min_loss     = min_loss,
    opt_m        = opt_mod$mstop,
    k_opts_m     = k_opts[seq_len(opt_mod$mstop)],
    ls_all       = c(ls0, ls[seq_len(opt_mod$mstop)]),
    edfs_all     = if (mod$edf_compute_method != "none")
                     c(edf0, edfs[seq_len(opt_mod$mstop)]) else NULL,
    elapsed_mins = as.numeric(difftime(Sys.time(), st$start_time,
                                       units = "mins")),
    force        = TRUE
  )
  gb_progress_close(st$progress)
  opt_mod$resume <- if (resumable) fit_state_pack(st) else NULL
  opt_mod
}
