# Check functions ---------------------------------------------------------



#' Check the response against the family
#'
#' @param family A family object.
#' @param model.resp Output of [stats::model.response()].
#' @param weights User-supplied prior weights, or `NULL`.
#' @noRd
check_family_inputs <- function(family, model.resp, weights) {
  fam <- family$family
  Y   <- model.resp
  if (fam == "binomial") return(check_binomial_response(Y, weights))

  # Every other family: a finite numeric vector, non-negative for Poisson
  checkmate::assert_vector(Y, .var.name = "response")
  if (fam == "poisson") {
    checkmate::assert_numeric(Y, lower = 0, finite = TRUE,
                              .var.name = "response")
    if (!checkmate::test_integerish(Y)) {
      warning("Response is not an integer despite using Poisson family.")
    }
  } else {
    checkmate::assert_numeric(Y, finite = TRUE, .var.name = "response")
  }
  if (fam %in% c("Gamma", "inverse.gaussian") && any(Y <= 0)) {
    stop(fam, " family requires a strictly positive response.",
         call. = FALSE)
  }
  invisible(NULL)
}

#' Check a binomial response
#'
#' One of: a two-column matrix of counts, a two-level factor, or a numeric
#' (or logical, read as 0/1 as in `glm()`) vector between 0 and 1.
#'
#' @param Y Output of [stats::model.response()].
#' @param weights User-supplied prior weights, or `NULL`; without them a
#'   fractional response gets a warning.
#' @noRd
check_binomial_response <- function(Y, weights) {
  if (is.logical(Y)) Y <- as.numeric(Y)
  if (is.matrix(Y)) {
    checkmate::assert_matrix(Y, ncols = 2, .var.name = "response")
    checkmate::assert_integerish(as.vector(Y), lower = 0,
                                 .var.name = "response")
  } else if (is.factor(Y)) {
    checkmate::assert_factor(Y, n.levels = 2, .var.name = "response")
  } else if (is.numeric(Y)) {
    checkmate::assert_numeric(Y, lower = 0, upper = 1, finite = TRUE,
                              .var.name = "response")
    if (is.null(weights) && !checkmate::test_integerish(Y)) {
      warning("Fractional binomial response but no weights supplied.")
    }
  } else {
    stop("Response type and/or class seems to be incompatible with ",
         "binomial family.", call. = FALSE)
  }
  invisible(NULL)
}

#' Check that the control arguments are compatible
#'
#' @param learner_control,stop_control,edf_control,uq_control Control
#'   lists from the `gb_*_ctrl()` functions.
#' @noRd
check_control_args_compatible <- function(learner_control,
                                          stop_control,
                                          edf_control,
                                          uq_control = NULL) {
  s_method <- stop_control$stopping_method
  e_method <- edf_control$edf_compute_method
  if (e_method == "none" && s_method == "aic") {
    stop("stopping_method = \"", s_method, "\" requires EDF tracking: use ",
         "gb_edf_ctrl(\"hutchpp\") or gb_edf_ctrl(\"exact\").",
         call. = FALSE)
  }
  if (identical(uq_control$method, "bootstrap")) {
    check_boot_requirements(
      uq_control$control$mode,
      stopping_method    = s_method,
      edf_compute_method = e_method)
  }
  invisible(NULL)
}



#' Check that learner indices are within range
#'
#' @param preselected,excluded,after `preselected_learner_indices`,
#'   `excluded_learner_indices` and `learners_to_fit_after_each_iter`.
#' @param n_learners Number of base learners.
#' @noRd
check_learner_index_bounds <- function(preselected, excluded, after = NULL,
                                       n_learners) {
  valid <- seq_len(n_learners)
  args  <- list(preselected_learner_indices     = preselected,
                excluded_learner_indices        = excluded,
                learners_to_fit_after_each_iter = after)
  for (nm in names(args)) {
    # setdiff(NULL, valid) is empty, so NULL arguments are skipped
    bad <- setdiff(args[[nm]], valid)
    if (length(bad) > 0)
      stop(nm, " contains out-of-range index/indices: ",
           paste(bad, collapse = ", "),
           ". Valid range is 1:", n_learners, ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Check that a learner list supports intercept padding
#'
#' @param learners List of base learners; the first must be the
#'   `"(Intercept)"` learner.
#' @noRd
check_intercept_padding <- function(learners) {
  if (!identical(learners[[1L]]$var_label, "(Intercept)")) {
    stop("intercept_padding = TRUE requires learner 1 to be the ",
         "\"(Intercept)\" learner.", call. = FALSE)
  }
  if (length(learners) < 2L) {
    stop("intercept_padding = TRUE requires at least one learner ",
         "besides the intercept.", call. = FALSE)
  }
  invisible(NULL)
}

#' Check that data-driven selection has at least one candidate learner
#'
#' @param n_learners Number of base learners.
#' @param excluded `excluded_learner_indices` (or `NULL`).
#' @param intercept_padding Flag; learner 1 is then not a candidate.
#' @noRd
check_selection_candidates <- function(n_learners, excluded,
                                       intercept_padding) {
  candidates <- setdiff(seq_len(n_learners), excluded)
  if (intercept_padding) candidates <- setdiff(candidates, 1L)
  if (length(candidates) == 0L) {
    stop("No base learner is left to select: excluded_learner_indices ",
         "removes every candidate",
         if (intercept_padding)
           " (learner 1 is not a candidate under intercept_padding = TRUE)",
         ".", call. = FALSE)
  }
  invisible(NULL)
}

#' Check that a formula is supported
#'
#' @param formula Model formula.
#' @param require_intercept Logical; reject a formula without an intercept.
#' @noRd
check_formula <- function(formula, require_intercept = TRUE) {
  # Only t2() smooth terms allowed (applies everywhere)
  trm_labs <- attr(stats::terms(formula), "term.labels")
  is_not_supported <- grepl("^(s|te|ti)\\(.*\\)$", trm_labs)
  if (any(is_not_supported)) {
    stop("s(.), ti(.), and te(.) terms are not supported. Use t2(.) instead.",
         call. = FALSE)
  }
  # gbgam() and new_learners() both build the intercept learner
  if (require_intercept) {
    intercept_included <- as.logical(attr(stats::terms(formula), "intercept"))
    if (isFALSE(intercept_included)) {
      stop("the formula requires an intercept: remove \"+ 0\" or \"- 1\" ",
           "from the ",
           "formula. To keep the intercept at zero, use gb_learner_ctrl(",
           "initialise_intercept = FALSE, excluded_learner_indices = 1L).",
           call. = FALSE)
    }
  }
  invisible(NULL)
}

#' Check that `true_f` matches the number of observations used in the fit
#'
#' @param true_f Numeric vector or `NULL` (from [gb_stop_ctrl()]).
#' @param n Number of rows in the model frame (after `subset`/`na.action`).
#' @noRd
check_true_f_length <- function(true_f, n) {
  if (!is.null(true_f) && length(true_f) != n) {
    stop("true_f must have one value per observation used in the fit (",
         n, "), got ", length(true_f), ".", call. = FALSE)
  }
  invisible(NULL)
}


#' Check a time limit in minutes
#'
#' @param x `NULL`, or a finite number of at least 0 (0 is valid: a budget
#'   used up before the uncertainty quantification).
#' @noRd
check_runtime_mins <- function(x) {
  checkmate::assert_number(x, lower = 0, finite = TRUE, null.ok = TRUE,
                           .var.name = "max_runtime_mins")
  invisible(NULL)
}


#' Check a `verbose` flag
#'
#' @param x A single logical value.
#' @noRd
check_verbose <- function(x) {
  checkmate::assert_flag(x, .var.name = "verbose")
  invisible(NULL)
}


#' Check an argument before match.arg()
#'
#' The argument must be a single string, unless it is the default vector of
#' choices, which [match.arg()] reads as its first element. The check is
#' explicit because the message of [match.arg()] for other lengths differs
#' between R versions.
#'
#' @param arg The argument.
#' @param choices Character vector. Its default vector of choices.
#' @param .var.name String. The argument name used in the message.
#' @noRd
check_choice <- function(arg, choices, .var.name) {
  if (!identical(arg, choices)) {
    checkmate::assert_string(arg, .var.name = .var.name)
  }
  invisible(NULL)
}
