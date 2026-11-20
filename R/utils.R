#' Simulate example data for GAMs
#'
#' Re-exported from [`mgcv`][mgcv::mgcv-package]. See [gamSim()][mgcv::gamSim]
#' for full documentation.
#'
#' @inheritParams mgcv::gamSim
#' @importFrom mgcv gamSim
#' @export
gamSim <- mgcv::gamSim

#' Tensor product smooths for `gbgam` formulae
#'
#' Re-exported from [`mgcv`][mgcv::mgcv-package]. See [t2()][mgcv::t2] for full
#' documentation.
#'
#' @inheritParams mgcv::t2
#' @importFrom mgcv t2
#' @export
t2 <- mgcv::t2


#' IWLS weight factor omega(mu) = 1 / (V(mu) g'(mu)) = h'(eta) / V(mu)
#' @noRd
omega <- function(mu, family) {
  eta <- family$linkfun(mu)
  family$mu.eta(eta) / family$variance(mu)
}


#' Compute Fletcher's (2012) bias-corrected dispersion estimator.
#'
#' phi_Fletcher = phi_Pearson / (1 + s_bar)
#'
#' where s_bar = weighted.mean(V'(mu_i) * (y_i - mu_i) / V(mu_i)),
#' V'(mu) = dV/dmu (available via mgcv::fix.family.var() as family$dvar),
#' and V(mu) is the variance function.
#'
#' Works during the boosting loop (before finalise_gbgam, so
#' $fitted.values and $df.residual are not yet available).
#'
#' Prior weights are treated as frequency weights: n = sum(W), so an
#' observation with weight w counts as w repeated observations. This matches
#' the weighting convention of the boosting fit, where a weight is
#' equivalent to repeating the observation. Consequently the estimate
#' depends on the scale of the weights; rescale them to sum to the
#' number of observations if they are not frequency weights.
#'
#' Weights summing to no more than the model edf~ (e.g. weights normalised to
#' sum to 1) leave no residual degrees of freedom; the estimate is then
#' 1e-12 (`phi_floor`).
#'
#' Reference: Fletcher, D. J. (2012). Estimating overdispersion when fitting
#' a generalized linear model to sparse data. *Biometrika, 99*(1), 230--237.
#'
#' @param mod A `gbgam` model object with `family`, `linear.predictors`,
#'   `prior.weights`, `y`, and `edf.total` fields.
#' @return A positive scalar dispersion estimate (1 for Poisson/Binomial).
#' @noRd
estimate_dispersion <- function(mod) {
  fam_name <- mod$family$family
  if (fam_name %in% c("poisson", "binomial")) return(1)
  phi_floor <- 1e-12
  W         <- mod$prior.weights
  # Residual df under the frequency-weight convention (see above)
  df_resid  <- sum(W) - (mod$edf.total %||% 0)
  if (isTRUE(df_resid <= 0)) return(phi_floor)

  mu          <- mod$family$linkinv(mod$linear.predictors)
  V           <- mod$family$variance(mu)
  V_floor     <- 1e-12  # guards the division by V(mu)
  X2          <- sum(W * (mod$y - mu)^2 / pmax(V, V_floor))  # Pearson
  phi_pearson <- X2 / df_resid

  # Fletcher s_bar
  dV    <- mod$family$dvar(mu)     # V'(mu), added by mgcv::fix.family.var()
  s_bar <- stats::weighted.mean(dV * (mod$y - mu) / pmax(V, V_floor), W)
  pmax(phi_pearson / (1 + s_bar), phi_floor)
}


#' TRUE when the link is the canonical link of the family
#' @noRd
is_canonical_link <- function(family) {
  switch(family$family,
         gaussian          = family$link == "identity",
         binomial          = family$link == "logit",
         poisson           = family$link == "log",
         Gamma             = family$link == "inverse",
         inverse.gaussian  = family$link == "1/mu^2",
         FALSE
  )
}


#' Fallback operator: x %||% y returns x unless it is NULL
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' The session's RNG state
#'
#' @return The value of `.Random.seed`, which is first created (by one
#'   uniform draw) if the session has not used its RNG yet.
#' @noRd
current_rng_state <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1)
  }
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}

#' Run code from a recorded RNG state, then restore the caller's
#'
#' @param state A `.Random.seed` value (it also records the RNG kind).
#' @param code Code to evaluate (lazily, after the state is set).
#' @return The value of `code`. The caller's `.Random.seed` (and with it
#'   the RNG kind) is restored, or removed if it did not exist.
#' @noRd
with_rng_state <- function(state, code) {
  env <- globalenv()
  if (exists(".Random.seed", envir = env, inherits = FALSE)) {
    old_state <- get(".Random.seed", envir = env, inherits = FALSE)
    on.exit(assign(".Random.seed", old_state, envir = env))
  } else {
    on.exit(rm(".Random.seed", envir = env))
  }
  assign(".Random.seed", state, envir = env)
  code
}

#' Split the integer total_int in proportion to weights (largest remainder)
#'
#' The parts are integers summing to total_int. The floor tolerates rounding
#' error in total_int * weights / sum(weights).
#' @noRd
split_integer_by_weights <- function(total_int, weights) {
  zs_frac <- weights / sum(weights) * total_int
  zs      <- floor(zs_frac + 1e-9)
  n_extra <- round(total_int - sum(zs))
  ind     <- order(zs_frac - zs, decreasing = TRUE)[seq_len(n_extra)]
  zs[ind] <- zs[ind] + 1
  zs
}

#' tr(x^T y) = sum(x * y), without forming the product
#' @noRd
tr_xprod <- function(x, y = NULL) {
  if (is.null(y)) y <- x
  sum(x * y)
}

#' tr(x y) = sum(t(x) * y), without forming the product
#' @noRd
tr_prod <- function(x, y = NULL) {
  if (is.null(y)) y <- x
  sum(t(x) * y)
}

#' Largest increase of a loss that counts as round-off
#'
#' @param loss The loss before the step.
#' @param roundoff Absolute floor of the rounding error, for a loss at its
#'   floating-point floor (see `deviance_roundoff()`).
#' @return `max(sqrt(.Machine$double.eps) * abs(loss), roundoff)`.
#' @noRd
loss_allowance <- function(loss, roundoff = 0) {
  max(sqrt(.Machine$double.eps) * abs(loss), roundoff)
}

#' Halve a step until the loss does not increase
#'
#' Tries the fractions 1, 1/2, 1/4, ... of a step, in that order, and keeps
#' the first with a finite loss of at most
#' `loss0 + allowance + 0.25 * fraction * slope`. With a negative `slope` this
#' is Armijo's sufficient decrease; with `slope = 0`, a loss that does not
#' increase. A step from a point of non-finite loss cannot be compared and is
#' taken whole.
#'
#' @param loss_at Function of the fraction, returning the loss there and
#'   `Inf` at an invalid point.
#' @param loss0 The loss before the step.
#' @param allowance Increase allowed for rounding (see `loss_allowance()`).
#' @param slope Derivative of the loss along the whole step, at most 0.
#' @param max_halvings Number of halvings to try.
#' @return list(fraction, loss); fraction 0 if no fraction was accepted.
#' @noRd
halve_step <- function(loss_at, loss0, allowance = 0, slope = 0,
                       max_halvings = 60L) {
  if (!is.finite(loss0)) return(list(fraction = 1, loss = loss_at(1)))
  for (h in 0:max_halvings) {
    fraction <- 2^-h
    loss     <- loss_at(fraction)
    if (is.finite(loss) &&
        loss <= loss0 + allowance + 0.25 * fraction * slope) {
      return(list(fraction = fraction, loss = loss))
    }
  }
  list(fraction = 0, loss = loss0)
}


# Progress display of the boosting loop ----------------------------------
#
# Live status of fit_gbgam(); the bootstrap and MCMC progress bars use cli
# directly. gb_progress_init() returns NULL when verbose = FALSE, and every
# other function then returns at once. Otherwise it picks one of three modes:
# "ansi" (a multi-line block redrawn in place by cursor movement), "cr" (the
# RStudio console, which ignores cursor movement: one line overwritten via
# \r) or "log" (no dynamic console: one newline-terminated line per update,
# at most every 5 s). gb_progress_update() redraws at most every 0.25 s,
# except with force = TRUE (the final update after the loop);
# gb_progress_close() ends the last line.

#' The box-drawing line or the Delta sign, as ASCII without UTF-8 output
#' @noRd
progress_glyph <- function(name) {
  utf8 <- cli::is_utf8_output()
  switch(name,
         line  = if (utf8) "\u2500" else "-",
         delta = if (utf8) "\u0394" else "d")
}

#' Text progress bar string: `[=====>   ]`
#' @noRd
progress_bar <- function(frac, width = 30L) {
  frac <- max(0.0, min(1.0, frac))
  n_eq <- floor(frac * width)
  tip  <- if (n_eq < width) 1L else 0L
  paste0("[",
         strrep("=", n_eq),
         if (tip) ">" else "",
         strrep(" ", max(0L, width - n_eq - tip)),
         "]")
}

#' Format a duration (minutes) as a compact string
#' @noRd
format_minutes <- function(mins) {
  if (is.null(mins) || is.na(mins)) return("?")
  if (mins < 1 / 60)                return("< 1s")
  if (mins < 1)                     return(paste0(round(mins * 60), "s"))
  if (mins < 60)                    return(paste0(sprintf("%.1f", mins), "m"))
  paste0(sprintf("%.1f", mins / 60), "h")
}

#' Learner labels shortened to at most max_w characters
#'
#' A longer label keeps its end from its last ")" (the learner suffix, e.g.
#' `rr` or `n1`, after the factor level of a by-smooth) and as much of its
#' start as fits, with "..." in between: the learners of one long `t2()`
#' term stay distinguishable. An end longer than max_w - 4 characters keeps
#' only its last characters.
#' @param labels Character vector of learner labels.
#' @param max_w Maximum number of characters, at least 4.
#' @noRd
short_learner_name <- function(labels, max_w = 24L) {
  vapply(labels, function(nm) {
    if (nchar(nm) <= max_w) return(nm)
    close <- regexpr(")[^)]*$", nm)            # the last ")", if any
    end   <- if (close > 0L) substring(nm, close) else ""
    keep  <- max_w - 3L - nchar(end)           # characters of the start
    if (keep >= 1L) {
      paste0(substr(nm, 1L, keep), "...", end)
    } else {
      paste0("...", substring(nm, nchar(nm) - max_w + 4L))
    }
  }, character(1L), USE.NAMES = FALSE)
}

#' State of the boosting-loop progress display (an environment)
#'
#' NULL when verbose = FALSE, which makes every other display function a
#' no-op.
#' @noRd
gb_progress_init <- function(verbose, mod, start_time) {
  if (!verbose) return(NULL)
  e <- new.env(parent = emptyenv())
  e$last_proc    <- proc.time()[["elapsed"]] - 1.0  # print at once
  e$n_lines      <- 0L
  e$width        <- min(getOption("width", 80L), 120L)
  # RStudio console supports ANSI colours but ignores cursor-movement codes
  # (ESC[nA).  Detect it explicitly and fall back to the \r single-line path;
  # without a dynamic console (a log file) print one line per update.
  e$mode         <- if (!cli::is_dynamic_tty()) {
    "log"
  } else if (Sys.getenv("RSTUDIO") == "") {
    "ansi"
  } else {
    "cr"
  }
  e$interval     <- if (e$mode == "log") 5 else 0.25
  e$loss_name    <- stop_loss_name(mod)
  e$edf_tracked  <- mod$edf_compute_method != "none"
  e$loss_tracked <- mod$stopping_method    != "none"
  e$fam_str      <- paste0(mod$family$family, "/", mod$family$link)
  e$has_pbar     <- !is.null(mod$max_iter) || !is.null(mod$max_runtime_mins)
  e$max_iter     <- mod$max_iter
  e$max_mins     <- mod$max_runtime_mins
  e
}

#' The multi-line progress block, one element per terminal line
#'
#' Each line is padded to exactly state$width characters.
#' @noRd
progress_lines <- function(state, m, mod, new_loss, min_loss, opt_m,
                           k_opts_m, ls_all, edfs_all, elapsed_mins) {
  W   <- state$width
  lnm <- state$loss_name
  lines <- character(0L)

  # 1. Header: a narrow console shortens the family string to keep m= ----
  rhs  <- paste0(" m=", m, " ")
  fam  <- state$fam_str
  room <- W - nchar(rhs) - 11L        # 11 = 2 + nchar(" gbgam  ") + 1
  if (nchar(fam) > room)
    fam <- if (room >= 4L) paste0(substr(fam, 1L, room - 3L), "...") else ""
  hdr <- if (nzchar(fam)) paste0(" gbgam  ", fam, " ") else " gbgam "
  if (nchar(hdr) + nchar(rhs) + 2L > W) hdr <- ""
  pad <- max(0L, W - nchar(hdr) - nchar(rhs) - 2L)
  lines <- c(lines, paste0(strrep(progress_glyph("line"), 2L), hdr,
                           strrep(progress_glyph("line"), pad), rhs))

  # 2. Progress bar (only when a fixed iteration or time budget is set) ----
  if (state$has_pbar) {
    pb_w <- max(20L, W - 26L)
    if (!is.null(state$max_iter)) {
      frac  <- m / state$max_iter
      lines <- c(lines,
                 sprintf("  %s  %d / %d  (%d%%)",
                         progress_bar(frac, pb_w), m, state$max_iter,
                         round(100 * frac)))
    } else {
      frac  <- min(1.0, elapsed_mins / state$max_mins)
      lines <- c(lines,
                 sprintf("  %s  %s / %s",
                         progress_bar(frac, pb_w),
                         format_minutes(elapsed_mins),
                         format_minutes(state$max_mins)))
    }
  }

  # 3. EDF line and loss line (when tracked; the first one gives the
  # elapsed time), or an iteration/elapsed line ---------------------------
  edf       <- edf_summary(edfs_all)
  edf_cur   <- edf$cur
  max_d_edf <- edf$max_d
  show_edf  <- state$edf_tracked && !is.na(edf_cur)
  show_loss <- state$loss_tracked && !is.null(new_loss) && !is.na(new_loss)

  if (show_edf) {
    delta_str <- if (!is.na(max_d_edf))
                   sprintf("   max %sEDF %6.2f", progress_glyph("delta"),
                           max_d_edf) else ""
    lines <- c(lines,
               sprintf("  EDF  %7.2f%s   %s",
                       edf_cur, delta_str, format_minutes(elapsed_mins)))
  }
  if (show_loss) {
    lines <- c(lines,
               sprintf("  %-5s  cur %12.4f (m=%d)   best %12.4f (m=%d)%s",
                       lnm, new_loss, m,
                       min_loss %||% NA_real_, opt_m %||% m,
                       if (show_edf) "" else
                         paste0("   ", format_minutes(elapsed_mins))))
  }
  if (!show_edf && !show_loss) {
    lines <- c(lines,
               sprintf("  m = %-6d   elapsed %s", m,
                       format_minutes(elapsed_mins)))
  }

  # 4. Blank separator ----------------------------
  lines <- c(lines, "")

  # 5. Learner table -- all demmler_reinsch learners in their natural order
  learners  <- mod$learners
  K         <- length(learners)
  sel_count <- tabulate(k_opts_m, nbins = K)

  # Per-learner cumulative loss change (negative = improvement in AIC/MSE)
  dloss    <- if (state$loss_tracked && length(ls_all) >= 2L) diff(ls_all)
  loss_red <- if (!is.null(dloss)) {
    vapply(seq_len(K), function(k) {
      idx <- which(k_opts_m == k)
      if (length(idx) > 0L) sum(dloss[idx], na.rm = TRUE) else 0.0
    }, numeric(1L))
  } else NULL

  # DR learners only, in the order they appear in $learners
  dr_idx   <- which(vapply(learners, inherits, logical(1L), "demmler_reinsch"))
  tot_sel  <- max(1L, sum(sel_count))
  # Signed shares of the loss change: improvements positive, increases
  # negative (+ 0 turns -0 into 0 for printing)
  loss_pct <- if (!is.null(loss_red)) {
    tot_loss <- abs(sum(loss_red))
    if (tot_loss > 1e-10) round(-100 * loss_red / tot_loss, 1) + 0 else
      numeric(K)
  }

  # Column header; cells with UTF-8 glyphs are padded by display width
  edf_hdr  <- if (state$edf_tracked) sprintf("%6s", "edf") else ""
  loss_hdr <- if (!is.null(loss_red))
    cli::ansi_align(paste0(progress_glyph("delta"), lnm, "%"), 8L,
                    "right") else ""
  lines <- c(lines,
             sprintf("  %-24s  %5s  %6s  %6s  %8s",
                     "learner", "sel", "%sel", edf_hdr, loss_hdr))

  lr_lines <- vapply(dr_idx, function(k) {
    lr       <- learners[[k]]
    ns       <- sel_count[k]
    edf_str  <- if (state$edf_tracked) {
      ev <- lr$edf %||% NA_real_
      if (is.na(ev)) cli::ansi_align(progress_glyph("line"), 6L, "right") else
        sprintf("%6.2f", ev)
    } else ""
    loss_str <- if (!is.null(loss_pct))
      sprintf("%7.1f%%", loss_pct[k]) else ""
    sprintf("  %-24s  %5d  %5.1f%%  %6s  %8s",
            short_learner_name(lr$var_label %||% "(unknown)"), ns,
            100.0 * ns / tot_sel, edf_str, loss_str)
  }, character(1L))
  lines <- c(lines, lr_lines)

  # 6. Footer ---------------------------------
  lines <- c(lines, strrep(progress_glyph("line"), W))

  # Pad / truncate every line to exactly W chars for clean in-place overwriting
  vapply(lines, function(l) {
    nc <- cli::ansi_nchar(l, type = "chars")
    if (nc < W) paste0(l, strrep(" ", W - nc)) else cli::ansi_substr(l, 1L, W)
  }, character(1L), USE.NAMES = FALSE)
}

#' Current edf~ and its largest increment along the path edfs_all
#'
#' NA where the path is too short (or NULL).
#' @noRd
edf_summary <- function(edfs_all) {
  n_edf <- length(edfs_all)
  list(cur   = if (n_edf >= 1L) edfs_all[n_edf] else NA_real_,
       max_d = if (n_edf >= 2L) max(diff(edfs_all), na.rm = TRUE) else
         NA_real_)
}

#' Single-line summary for the "cr" mode (overwritten in place via \r, as in
#' txtProgressBar()) and the "log" mode
#'
#' The "log" mode writes the whole line. The "cr" mode pads it to
#' state$width characters; a longer line drops whole parts, in this order,
#' until it fits: the family string, the largest edf~ increment, the best
#' loss, the "gbgam" label, the elapsed time, edf~ and the current loss.
#' The iteration is always shown.
#' @noRd
progress_oneliner <- function(state, m, mod, new_loss, min_loss, opt_m,
                              edfs_all, elapsed_mins) {
  lnm      <- state$loss_name
  edf      <- edf_summary(edfs_all)
  edf_cur  <- edf$cur
  max_d    <- edf$max_d

  # Each part carries its leading separator
  parts <- c(gbgam = "  gbgam", fam = paste0("  ", state$fam_str),
             m = sprintf("   m=%-5d", m))
  if (state$edf_tracked && !is.na(edf_cur)) {
    parts["edf"] <- sprintf("   EDF: %.1f", edf_cur)
    if (!is.na(max_d))
      parts["max_d"] <- sprintf("  max%sEDF: %.2f", progress_glyph("delta"),
                                max_d)
  }
  if (state$loss_tracked && !is.null(new_loss) && !is.na(new_loss)) {
    parts["loss"] <- sprintf("   %s: %.4f", lnm, new_loss)
    parts["best"] <- sprintf(" (best %.4f m=%d)", min_loss %||% new_loss,
                             opt_m %||% m)
  }
  parts["time"] <- paste0("   ", format_minutes(elapsed_mins))

  line_width <- function(p) nchar(paste(p, collapse = ""), type = "width")
  if (identical(state$mode, "log")) return(paste(parts, collapse = ""))
  W <- state$width
  for (drop in c("fam", "max_d", "best", "gbgam", "time", "edf", "loss")) {
    if (line_width(parts) <= W) break
    parts <- parts[names(parts) != drop]
  }
  line <- paste(parts, collapse = "")
  paste0(line, strrep(" ", max(0L, W - line_width(parts))))
}

#' Draw or update the boosting-loop progress display
#'
#' At most every state$interval seconds; force = TRUE bypasses the throttle
#' (the final display).
#' @noRd
gb_progress_update <- function(state, m, mod, new_loss, min_loss, opt_m,
                                k_opts_m, ls_all, edfs_all = NULL,
                                elapsed_mins, force = FALSE) {
  if (is.null(state)) return(invisible(NULL))
  now <- proc.time()[["elapsed"]]
  if (!force && now - state$last_proc < state$interval)
    return(invisible(NULL))
  state$last_proc <- now

  if (state$mode == "ansi") {
    lines <- progress_lines(state, m, mod, new_loss, min_loss, opt_m,
                            k_opts_m, ls_all, edfs_all, elapsed_mins)
    n <- length(lines)
    if (state$n_lines > 0L) {
      cat(sprintf("\033[%dA", state$n_lines))
      cat(paste0("\033[2K", lines, "\n", collapse = ""), sep = "")
    } else {
      cat(paste(lines, collapse = "\n"), "\n", sep = "")
    }
    state$n_lines <- n
  } else {
    one <- progress_oneliner(state, m, mod, new_loss, min_loss, opt_m,
                             edfs_all, elapsed_mins)
    if (state$mode == "log") {
      cat(trimws(one, "right"), "\n", sep = "")
    } else {
      if (state$n_lines > 0L) cat("\r")
      cat(one)
    }
    flush.console()
    state$n_lines <- 1L
  }
  invisible(NULL)
}

#' Close the progress display: finish the current line (log lines already
#' end with a newline)
#' @noRd
gb_progress_close <- function(state) {
  if (is.null(state) || state$mode == "log") return(invisible(NULL))
  cat("\n")
  invisible(NULL)
}


