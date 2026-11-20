#' Control parameters for EDF tracking
#'
#' @description
#' Settings for tracking the effective degrees of freedom (EDF) of the fit
#' during boosting. The EDF drive adaptive smoothness relaxation (ASR) and
#' stopping by AIC, and the bootstrap modes other than `"fixed"` need them
#' (see [gb_param_boot_ctrl()]). Supplying an argument that the chosen
#' method does not use results in an error.
#'
#' @param edf_compute_method Character. One of:
#'   \describe{
#'     \item{`"hutchpp"`}{estimates the EDF with the stochastic trace
#'       estimator NA-Hutch++ (see Details). It scales well to large data
#'       sets; its accuracy is set by `z`. It is unreliable with non-canonical
#'       links (see Details).}
#'     \item{`"exact"`}{tracks the full \eqn{n \times n} Jacobian of the
#'       fitted values with respect to the response and takes its trace, at
#'       \eqn{\mathcal{O}(n^2 p_{\hat{k}})}{O(n^2 p_k)} time per iteration
#'       (\eqn{p_{\hat{k}}}{p_k} is the number of columns of the selected
#'       learner) and \eqn{\mathcal{O}(n^2)}{O(n^2)} memory; practical up to
#'       about \eqn{n = 2000}.}
#'     \item{`"none"`}{tracks no EDF, so iterations are faster. Stopping by
#'       AIC and the bootstrap modes `"free_smoothing"`, `"free_selection"`
#'       and `"free_stopping"` are then not available, and ASR is turned off
#'       (`tau = 0`). Bayesian backfitting then starts its search for the
#'       smoothing parameters from their starting values.}
#'   }
#'   Default `"hutchpp"`.
#' @param z Integer, at least 3, or NULL. Only with
#'   `edf_compute_method = "hutchpp"`. Number of probe vectors. More probe
#'   vectors give more accurate EDF at a higher cost: the cost of each EDF
#'   update grows linearly in `z`. The error of the estimate depends on the
#'   Jacobian whose trace is estimated, not directly on \eqn{n}; increase
#'   `z` for more precise EDF (e.g. for reliable stopping by AIC on noisy
#'   data), and consider `"exact"` if `z` would need to approach \eqn{n}.
#'   Default `NULL`: `80`.
#' @param cweights Positive numeric vector of length 3, or NULL. Only with
#'   `edf_compute_method = "hutchpp"`. For advanced users: the proportions
#'   in which the `z` probe vectors are split among the three stages of
#'   NA-Hutch++ (see Details). Each stage needs at least one probe vector.
#'   Default `NULL`: `c(1, 2, 1)`, which should rarely need changing.
#'
#' @return A list of class `"gb_edf_ctrl"`.
#'
#' @details
#' Non-adaptive Hutch++ (NA-Hutch++; Meyer et al., 2021) estimates the trace
#' of a matrix in two parts: the trace of a low-rank approximation of the
#' matrix, computed exactly, and a Hutchinson estimate of the trace of the
#' rest. The `z` probe vectors are split in the proportions `cweights`
#' among the three stages: the two sketches that build the low-rank
#' approximation, and the Hutchinson estimate.
#'
#' The matrix here is the Jacobian of the fitted means with respect to the
#' response, \eqn{d\tilde{\mu}/dy}{d mu~/dy}, whose trace is the model edf
#' \eqn{\widetilde{\operatorname{edf}}}{edf~}. `gbgam` tracks only its
#' products with the probe vectors (right sketches), so the left sketch of
#' NA-Hutch++ is replaced by the transpose of a right sketch. The estimate
#' stays unbiased.
#'
#' ## Non-canonical links
#'
#' Do not rely on `"hutchpp"` with a non-canonical link, such as
#' `Gamma("log")`; [gbgam()] gives a warning for this combination. Unbiased
#' does not mean accurate: the Jacobian of such a model is far from symmetric,
#' its Frobenius norm can be much larger than its trace, and a single estimate
#' of \eqn{\widetilde{\operatorname{edf}}}{edf~} can then be far from the true
#' value, even negative. In one Gamma(log) fit with an exposure offset, its
#' true value was 1.08 after one iteration and the estimates had a standard
#' deviation of about 3, which grew to about 32 by the end of the fit. Stopping
#' by AIC and ASR then follow noise. Use `"exact"` where the data allow it.
#' With larger data, use a canonical link, or `"none"` with stopping by
#' cross-validation.
#'
#' @references
#' Meyer, R. A., Musco, C., Musco, C., & Woodruff, D. P. (2021).
#' Hutch++: Optimal stochastic trace estimation.
#' In *Symposium on Simplicity in Algorithms (SOSA)* (pp. 142--155). Society
#' for Industrial and Applied Mathematics.
#' \href{https://doi.org/10.1137/1.9781611976496.16}{doi:10.1137/1.9781611976496.16}
#'
#' @seealso [gb_stop_ctrl()], [gb_learner_ctrl()]
#'
#' @examples
#' # Default: NA-Hutch++ with z = 80 probe vectors
#' gb_edf_ctrl()
#'
#' # Exact trace for small data sets (n up to about 2000)
#' gb_edf_ctrl(edf_compute_method = "exact")
#'
#' # No EDF tracking, when neither AIC stopping nor ASR is needed
#' gb_edf_ctrl(edf_compute_method = "none")
#'
#' # More probe vectors for more accurate EDF
#' gb_edf_ctrl(edf_compute_method = "hutchpp", z = 200)
#'
#' @export
gb_edf_ctrl <- function(
    edf_compute_method = c("hutchpp", "exact", "none"),
    z = NULL,
    cweights = NULL) {

  # Validate inputs
  edf_compute_method <- match.arg(edf_compute_method)
  if (edf_compute_method != "hutchpp") {
    unused <- c(z = !is.null(z), cweights = !is.null(cweights))
    if (any(unused))
      stop(names(unused)[unused][1L], " is only used with ",
           "edf_compute_method = \"hutchpp\".", call. = FALSE)
  }

  if (!is.null(z)) {
    checkmate::assert_integerish(z, lower = 3, len = 1, any.missing = FALSE,
                                 .var.name = "z (Hutch++ samples)")
  }

  if (!is.null(cweights)) {
    checkmate::assert_numeric(cweights, lower = 0, finite = TRUE, len = 3,
                              any.missing = FALSE, .var.name = "cweights")
    if (any(cweights <= 0)) {
      stop("cweights must be strictly positive.", call. = FALSE)
    }
  }

  # Set defaults logic
  if (edf_compute_method == "hutchpp") {
    z        <- z %||% 80
    cweights <- cweights %||% c(1, 2, 1)
    zs <- split_integer_by_weights(z, cweights)
    if (any(zs < 1)) {
      stop("Each Hutch++ block needs at least one probe vector: z = ", z,
           " with cweights = c(", toString(cweights), ") gives blocks of ",
           "size ", toString(zs), ". Increase z or balance cweights.",
           call. = FALSE)
    }
  }
  structure(list(
    edf_compute_method = edf_compute_method,
    z = z,
    cweights = cweights
  ), class = "gb_edf_ctrl")
}


#' Control parameters for stopping
#'
#' @description
#' Settings for when boosting stops. Three rules can be combined: a
#' stopping method, which tracks a loss and stops once it has not improved
#' for a number of iterations (the patience, set by `iter_no_improv` and
#' `iter_no_improv_rel`), returning the model at its best iteration; an
#' iteration cap `max_iter`; and a time limit `max_runtime_mins`. With
#' `stopping_method = "none"`, supply `max_iter` or `max_runtime_mins`.
#' Supplying an argument that the chosen method does not use results in an
#' error.
#'
#' @param stopping_method Character. The loss that is tracked. One of:
#'   \describe{
#'     \item{`"aic"`}{the conditional AIC,
#'       \eqn{-2 \log L + c \, \widetilde{\operatorname{edf}}}{
#'       -2 log L + c edf~}, with \eqn{c} = `AIC_k`. The log-likelihood
#'       \eqn{L} is the full one, with the dispersion phi estimated again at
#'       each iteration by Fletcher's (2012) estimator for the Gaussian,
#'       Gamma and inverse Gaussian families, and fixed at 1 for the Poisson
#'       and binomial families. Needs EDF tracking (see [gb_edf_ctrl()]).}
#'     \item{`"cv"`}{a resampling estimate of the prediction error, such as
#'       k-fold cross-validation or the .632+ bootstrap, with the settings
#'       of `resampling` (see [gb_resampling()]). It needs no EDF tracking,
#'       and no `max_iter` has to be guessed.}
#'     \item{`"true_f"`}{the mean squared error against the true linear
#'       predictor `true_f`; for simulation studies only.}
#'     \item{`"none"`}{no loss; boosting runs until `max_iter` or
#'       `max_runtime_mins` is reached.}
#'   }
#'   Default `"aic"`.
#' @param iter_no_improv Integer, at least 1, or NULL. Not used with
#'   `stopping_method = "none"`. The smallest number of iterations
#'   without a new best loss after which boosting stops. With `"aic"` and
#'   `"true_f"` the rule is checked after every iteration, and a loss equal
#'   to the best one counts as a new best (the later iteration is kept).
#'   With `"cv"` it is checked every `chunk` iterations (see
#'   [gb_resampling()]), so up to `chunk - 1` more iterations may run, and
#'   ties go to the earlier iteration. Default `NULL`: `50`.
#' @param iter_no_improv_rel Non-negative number, or NULL. Not used with
#'   `stopping_method = "none"`. The patience relative to the best
#'   iteration so far, m*: boosting stops once
#'   `max(iter_no_improv, ceiling(iter_no_improv_rel * m*))` iterations
#'   have run without a new best loss, so a fit whose best iteration comes
#'   late is given more iterations to improve. With `0`, boosting stops
#'   after `iter_no_improv` iterations whatever m*. Default `NULL`: `0.1`.
#' @param true_f Numeric vector, or NULL. Only with
#'   `stopping_method = "true_f"`, where it is required. The true linear
#'   predictor, including any offset, one value per observation used in the
#'   fit (after `subset` and `na.action`). Default `NULL`.
#' @param AIC_k Non-negative number, or NULL. Only with
#'   `stopping_method = "aic"`. The penalty \eqn{c} per effective degree of
#'   freedom in the AIC; `log(n)` gives the BIC. Default `NULL`: `2`.
#' @param max_iter Integer, at least 1, or NULL. Iteration cap: boosting
#'   stops after this many iterations, whatever the stopping method. Default
#'   `NULL` (no cap).
#' @param max_runtime_mins Non-negative number, or NULL. Time limit in
#'   minutes; when it is reached, boosting stops as at the iteration cap.
#'   Default `NULL` (no limit).
#' @param resampling A list from [gb_resampling()], or NULL. Only with
#'   `stopping_method = "cv"`. The resampling settings. Default `NULL`:
#'   `gb_resampling()`, 10-fold cross-validation of the deviance.
#' @param resumable Logical. If `TRUE`, the fitted model keeps the state of
#'   its last iteration (`$resume`), so that [extend_gbgam()] can continue
#'   the fit. The state is small apart from the EDF-tracking matrices
#'   (\eqn{n \times z}{n x z} for `"hutchpp"`, \eqn{n \times n}{n x n} for
#'   `"exact"`). Default `TRUE`.
#'
#' @return A list of class `"gb_stop_ctrl"`.
#'
#' @references
#' Fletcher, D. J. (2012).
#' Estimating overdispersion when fitting a generalized linear model to sparse
#' data.
#' *Biometrika, 99*(1), 230--237.
#' \href{https://doi.org/10.1093/biomet/asr083}{doi:10.1093/biomet/asr083}
#'
#' @seealso [gb_resampling()], [gb_edf_ctrl()], [extend_gbgam()]
#'
#' @examples
#' # Stopping by AIC (default): stops after max(50, ceiling(0.1 * m*))
#' # iterations without a new best AIC, m* the best iteration so far
#' gb_stop_ctrl()
#'
#' # Stopping by AIC after 100 iterations without a new best, whatever m*
#' gb_stop_ctrl(iter_no_improv = 100, iter_no_improv_rel = 0)
#'
#' # Patience of at least 20 iterations and a quarter of m*
#' gb_stop_ctrl(iter_no_improv = 20, iter_no_improv_rel = 0.25)
#'
#' # Stopping by BIC (stops earlier than AIC)
#' n <- 200 # number of observations
#' gb_stop_ctrl(AIC_k = log(n))
#'
#' # Stopping by BIC with an iteration cap of 1000
#' gb_stop_ctrl(AIC_k = log(n), iter_no_improv = 100, max_iter = 1000)
#'
#' # A fixed number of iterations, with no loss tracked
#' gb_stop_ctrl(stopping_method = "none", max_iter = 200)
#'
#' # A time limit only
#' gb_stop_ctrl(stopping_method = "none", max_runtime_mins = 2)
#'
#' # 10-fold cross-validation of the deviance
#' gb_stop_ctrl("cv")
#'
#' # .632+ bootstrap with 25 samples
#' gb_stop_ctrl("cv", resampling = gb_resampling("bootstrap", B = 25))
#'
#' @export
gb_stop_ctrl <- function(
    stopping_method = c("aic", "none", "true_f", "cv"),
    iter_no_improv = NULL,
    iter_no_improv_rel = NULL,
    true_f = NULL,
    AIC_k = NULL,
    max_iter = NULL,
    max_runtime_mins = NULL,
    resampling = NULL,
    resumable = TRUE) {

  # Validate inputs
  stopping_method <- match.arg(stopping_method)

  checkmate::assert_integerish(iter_no_improv, lower = 1, len = 1,
                               any.missing = FALSE, null.ok = TRUE,
                               .var.name = "iter_no_improv")

  checkmate::assert_number(iter_no_improv_rel, lower = 0, finite = TRUE,
                           null.ok = TRUE, .var.name = "iter_no_improv_rel")

  checkmate::assert_numeric(true_f, finite = TRUE, any.missing = FALSE,
                            null.ok = TRUE, .var.name = "true_f")

  checkmate::assert_number(AIC_k, lower = 0, finite = TRUE, null.ok = TRUE,
                           .var.name = "AIC_k")

  checkmate::assert_integerish(max_iter, lower = 1, len = 1,
                               any.missing = FALSE, null.ok = TRUE,
                               .var.name = "max_iter")

  check_runtime_mins(max_runtime_mins)

  checkmate::assert_flag(resumable, .var.name = "resumable")

  # Logical consistency checks
  if (!is.null(resampling)) {
    if (stopping_method != "cv") {
      stop("resampling is only used with stopping_method = \"cv\".",
           call. = FALSE)
    }
    if (!inherits(resampling, "gb_resampling")) {
      stop("resampling must be created with gb_resampling().", call. = FALSE)
    }
  }
  if (stopping_method == "cv" && is.null(resampling)) {
    resampling <- gb_resampling()
  }
  if (stopping_method == "true_f" && is.null(true_f)) {
    stop("true_f is required with stopping_method = \"true_f\".",
         call. = FALSE)
  }
  unused <- c(true_f         = stopping_method != "true_f" && !is.null(true_f),
              AIC_k          = stopping_method != "aic" && !is.null(AIC_k),
              iter_no_improv = stopping_method == "none" &&
                               !is.null(iter_no_improv),
              iter_no_improv_rel = stopping_method == "none" &&
                                   !is.null(iter_no_improv_rel))
  if (any(unused)) {
    nm <- names(unused)[unused][1L]
    stop(nm, " is not used with stopping_method = \"", stopping_method,
         "\".", call. = FALSE)
  }

  # Ensure valid stopping condition exists
  no_hard_stop <- is.null(max_iter) && is.null(max_runtime_mins)
  if (stopping_method == "none" && no_hard_stop) {
    stop("No stopping condition provided. Specify max_iter, ",
         "max_runtime_mins, or a stopping_method other than \"none\".",
         call. = FALSE)
  }

  # Set defaults logic
  if (is.null(iter_no_improv) && stopping_method != "none") {
    iter_no_improv <- 50
  }
  if (is.null(iter_no_improv_rel) && stopping_method != "none") {
    iter_no_improv_rel <- 0.1
  }
  if (is.null(AIC_k) && stopping_method == "aic") {
    AIC_k <- 2
  }

  structure(list(
    stopping_method = stopping_method,
    iter_no_improv = iter_no_improv,
    iter_no_improv_rel = iter_no_improv_rel,
    true_f = true_f,
    AIC_k = AIC_k,
    max_iter = max_iter,
    max_runtime_mins = max_runtime_mins,
    resampling = resampling,
    resumable = resumable
  ), class = "gb_stop_ctrl")
}


#' Control parameters for learner fitting and selection
#'
#' @description
#' Settings for how the base learners are fitted and selected at each
#' boosting iteration. The most important is `step_size`, which sets the
#' trade-off between accuracy and the number of iterations; `tau` sets
#' adaptive smoothness relaxation (ASR). The other settings can be left at
#' their defaults in most applications.
#'
#' @param step_size Positive number, usually in (0, 1]. The global step
#'   size nu, which shrinks every update of the selected learner. Smaller
#'   values tend to give more accurate fits but need more iterations. ASR
#'   can be more sensitive to the step size than boosting without it. As in
#'   `mboost`, the learners fit the negative gradient at the prior weights
#'   only, so the step relative to a Fisher scoring step is `step_size`
#'   times the Fisher weights; for example, large Poisson counts with a log
#'   link need a step size of roughly `0.1 / mean(y)`. A step size that is
#'   too large makes boosting diverge: [gbgam()] then stops with an error
#'   as soon as the training deviance increases or the negative gradient is
#'   not finite; refit with a smaller value. Default `0.1`, as in `mboost`.
#' @param learners_to_fit_after_each_iter Integer vector, or NULL. Learners
#'   that are updated after the selected learner in every iteration, in the
#'   given order, each fitted to the negative gradient at the updated linear
#'   predictor and shrunk by `step_size`. Useful for learners that should
#'   follow every step rather than compete in selection, for example:
#'   * **Intercept refit** (as the `always` argument of `bamlss::boost()`):
#'     `learners_to_fit_after_each_iter = 1L` with
#'     `excluded_learner_indices = 1L` refits the intercept after every
#'     iteration without letting it compete. The intercept is fitted to the
#'     gradient after the main update, not to the gradient the selected
#'     learner was fitted to, so this differs from `intercept_padding`
#'     except for the Gaussian family with identity link (see
#'     `intercept_padding`).
#'   * **Random effects**: a random-effects learner can be updated after
#'     every selection step.
#'
#'   Default `NULL` (none).
#' @param initialise_intercept Logical. If `TRUE`, the intercept \eqn{a}
#'   starts at the maximum likelihood estimate of the intercept-only model,
#'   the link function of the weighted mean response (for a model without
#'   `offset`). If `FALSE`, boosting starts from \eqn{a = 0} (the linear
#'   predictor equals the `offset`) with
#'   \eqn{\widetilde{\operatorname{edf}} = 0}{edf~ = 0}, and \eqn{a} is
#'   estimated by boosting only; the null deviance still uses the
#'   intercept-only estimate. `FALSE` is not recommended: the smooth
#'   learners are centred to sum to zero and rely on a well-started
#'   intercept. It cannot be used with links for which a linear predictor
#'   equal to the `offset` is invalid, such as `inverse` (Gamma) and
#'   `1/mu^2` (inverse Gaussian) without an offset; [gbgam()] then stops
#'   with an error. Default `TRUE`.
#' @param selection_loss_type Character. The loss by which the learner of
#'   each iteration is selected. One of:
#'   \describe{
#'     \item{`"rss"`}{the weighted residual sum of squares of the negative
#'       gradient after each learner's provisional update, shrunk by its
#'       learner step size nu_k. This is the usual rule in component-wise
#'       boosting (e.g. `mboost`), and the one under which the calibration to
#'       sedf 1 gives unbiased selection.}
#'     \item{`"deviance"`}{the deviance after the provisional update. It
#'       equals `"rss"` for the Gaussian family with identity link. For
#'       other families the gradient is not on the scale of the linear
#'       predictor, so an unshrunk provisional update can overshoot and make
#'       the comparison meaningless; use it only with
#'       `use_step_size_for_selection = TRUE`.}
#'   }
#'   Default `"rss"`.
#' @param use_step_size_for_selection Logical. If `FALSE`, the learners are
#'   compared on provisional updates shrunk by their learner step size nu_k
#'   only, as in `mboost`. Under this setting the calibration to sedf 1 makes
#'   the expected drop in the weighted residual sum of squares the same for
#'   every learner when no covariate has an effect (Hofner et al., 2011). If
#'   `TRUE`, the provisional updates are also shrunk by `step_size`; the
#'   expected drop then differs between learners, and selection favours the
#'   less shrunk learners (e.g. unpenalised ones). Default `FALSE`.
#' @param use_learner_step_size_for_update Logical. If `TRUE`, the update of
#'   the selected learner is shrunk by `step_size` times its learner step
#'   size nu_k, and not by `step_size` alone. nu_k is chosen so that the
#'   selection edf of every learner,
#'   \eqn{\text{sedf}_k = 2\nu_k \text{tr}(\mathbf{A}_k) -
#'   \nu_k^2 \text{tr}(\mathbf{A}_k^2)}{
#'   sedf_k = 2 nu_k tr(A_k) - nu_k^2 tr(A_k^2)}, equals one, where
#'   \eqn{\mathbf{A}_k}{A_k} is the learner's hat matrix. nu_k is used in
#'   selection either way. Default `FALSE`.
#' @param knots Named list, or NULL. Knot positions passed to
#'   [mgcv::smoothCon()]; each name must be a covariate of a smooth.
#'   Default `NULL`: `mgcv`'s knots.
#' @param tau Non-negative number, or NULL. The rate of ASR. Each learner
#'   \eqn{k} accumulates the changes in
#'   \eqn{\widetilde{\operatorname{edf}}}{edf~} of the iterations that
#'   selected it, \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k}. Whenever
#'   `tau` times \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} exceeds the
#'   learner's current ledf (the trace of its hat matrix), its smoothing
#'   parameter \eqn{\lambda_k} is lowered so that they match;
#'   \eqn{\lambda_k} is never raised. A target above the learner's number of
#'   columns sets \eqn{\lambda_k = 0}. `tau = 0` turns ASR off, and
#'   `tau = 1` lets the ledf follow
#'   \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} one for one. A
#'   `tau > 0` together with `gb_edf_ctrl("none")` is ignored with a
#'   warning. Default `NULL`: `1` when EDF are tracked, else `0`.
#' @param preselected_learner_indices Integer vector, or NULL. For advanced
#'   users: the learner to update at each iteration, one per iteration, in
#'   place of the selection step. It must have at least `max_iter` elements
#'   (see [gb_stop_ctrl()]). Useful for replaying a selection path. Default
#'   `NULL` (selection by `selection_loss_type`).
#' @param preselected_rhos Numeric vector, or NULL. Only with
#'   `preselected_learner_indices`, and of the same length. For advanced
#'   users: the log smoothing parameter (`rho` = log(lambda)) given to the
#'   selected learner after each iteration, in place of the ASR update.
#'   Useful for replaying a smoothing-parameter path. Default `NULL` (ASR
#'   by `tau`).
#' @param excluded_learner_indices Integer vector, or NULL. For advanced
#'   users: learners that are never selected. They can still be updated
#'   through `preselected_learner_indices` or
#'   `learners_to_fit_after_each_iter`. Default `NULL` (none).
#' @param intercept_padding Logical. If `TRUE`, every learner is treated as
#'   padded with an unpenalised intercept column: each update fits both the
#'   intercept-only learner (learner 1) and the selected learner to the
#'   same negative gradient, and adds both fits, each shrunk by
#'   `step_size`, in one iteration. Because all learners are centred, the
#'   hat matrix of a padded learner is the sum of those of its two parts.
#'   The intercept then does not compete in selection, so learner 1 must
#'   not appear in `preselected_learner_indices` or
#'   `learners_to_fit_after_each_iter`. Selection by `"rss"` is unchanged
#'   (padding adds the same constant for every learner), while `"deviance"`
#'   compares the provisional updates including the intercept shift. The
#'   change of \eqn{\widetilde{\operatorname{edf}}}{edf~} in the iteration
#'   is split between the two parts: the intercept's share is credited to
#'   learner 1, the rest to the selected learner, where it drives ASR. The
#'   updates of `learners_to_fit_after_each_iter` are padded as well. For
#'   the Gaussian family with identity link and
#'   `initialise_intercept = TRUE` the intercept fit is zero in every
#'   iteration, so padding has no effect. Default `FALSE`: the intercept is
#'   a learner that competes in selection.
#'
#' @return A list of class `"gb_learner_ctrl"`.
#'
#' @references
#' Hofner, B., Hothorn, T., Kneib, T., & Schmid, M. (2011).
#' A framework for unbiased model selection based on boosting.
#' *Journal of Computational and Graphical Statistics, 20*(4), 956--971.
#' \href{https://doi.org/10.1198/jcgs.2011.09220}{doi:10.1198/jcgs.2011.09220}
#'
#' @seealso [gb_edf_ctrl()], [gb_stop_ctrl()]
#'
#' @examples
#' # Defaults: step size 0.1; ASR on when EDF are tracked (tau = NULL,
#' # resolved to 1 by gbgam())
#' gb_learner_ctrl()
#'
#' # A smaller step size (more iterations needed), e.g. for Poisson counts
#' # with a mean of about 50
#' gb_learner_ctrl(step_size = 0.002)
#'
#' # ASR turned off
#' gb_learner_ctrl(tau = 0)
#'
#' # The intercept refitted after every iteration, without competing
#' # (as bamlss::boost(always = ); not the same as intercept padding)
#' gb_learner_ctrl(learners_to_fit_after_each_iter = 1L,
#'                 excluded_learner_indices = 1L)
#'
#' # Padded learners: the intercept is fitted to the same gradient as the
#' # selected learner in every iteration
#' gb_learner_ctrl(intercept_padding = TRUE)
#'
#' @export
gb_learner_ctrl <- function(
    step_size = 0.1,
    learners_to_fit_after_each_iter = NULL,
    initialise_intercept = TRUE,
    selection_loss_type = c("rss", "deviance"),
    use_step_size_for_selection = FALSE,
    use_learner_step_size_for_update = FALSE,
    knots = NULL,
    tau = NULL,
    preselected_learner_indices = NULL,
    preselected_rhos = NULL,
    excluded_learner_indices = NULL,
    intercept_padding = FALSE) {

  # Validate inputs
  checkmate::assert_number(step_size, lower = 0, finite = TRUE,
                           .var.name = "step_size")
  if (step_size <= 0) {
    stop("step_size must be strictly positive.", call. = FALSE)
  }
  selection_loss_type <- match.arg(selection_loss_type)
  checkmate::assert_flag(initialise_intercept,
                         .var.name = "initialise_intercept")
  checkmate::assert_flag(use_step_size_for_selection,
                         .var.name = "use_step_size_for_selection")
  checkmate::assert_flag(use_learner_step_size_for_update,
                         .var.name = "use_learner_step_size_for_update")

  checkmate::assert_integerish(learners_to_fit_after_each_iter, lower = 1,
                               any.missing = FALSE, null.ok = TRUE,
                               .var.name = "learners_to_fit_after_each_iter")

  checkmate::assert_list(knots, null.ok = TRUE, .var.name = "knots")

  checkmate::assert_number(tau, lower = 0, finite = TRUE, null.ok = TRUE,
                           .var.name = "tau")

  checkmate::assert_integerish(preselected_learner_indices, lower = 1,
                               any.missing = FALSE, null.ok = TRUE,
                               .var.name = "preselected_learner_indices")

  checkmate::assert_numeric(preselected_rhos, finite = TRUE,
                            any.missing = FALSE, null.ok = TRUE,
                            .var.name = "preselected_rhos")

  checkmate::assert_integerish(excluded_learner_indices, lower = 1,
                               any.missing = FALSE, null.ok = TRUE,
                               .var.name = "excluded_learner_indices")

  checkmate::assert_flag(intercept_padding, .var.name = "intercept_padding")
  if (intercept_padding) {
    # The intercept (learner 1) is part of every padded update, never a
    # learner of its own
    idx_args <- list(
      learners_to_fit_after_each_iter = learners_to_fit_after_each_iter,
      preselected_learner_indices     = preselected_learner_indices)
    for (nm in names(idx_args)) {
      if (1 %in% idx_args[[nm]]) {
        stop("intercept_padding = TRUE fits the intercept (learner 1) in ",
             "every iteration; remove 1 from ", nm, ".", call. = FALSE)
      }
    }
  }

  if (!is.null(preselected_rhos)) {
    if (is.null(preselected_learner_indices)) {
      stop("supply preselected_learner_indices when supplying ",
           "preselected_rhos.", call. = FALSE)
    }
    if (length(preselected_learner_indices) != length(preselected_rhos)) {
      stop("preselected_learner_indices and preselected_rhos must have the ",
           "same length.", call. = FALSE)
    }
  }

  structure(list(
    step_size                        = step_size,
    initialise_intercept             = initialise_intercept,
    learners_to_fit_after_each_iter  = learners_to_fit_after_each_iter,
    selection_loss_type              = selection_loss_type,
    use_step_size_for_selection      = use_step_size_for_selection,
    use_learner_step_size_for_update = use_learner_step_size_for_update,
    knots                            = knots,
    tau                              = tau,
    preselected_learner_indices      = preselected_learner_indices,
    preselected_rhos                 = preselected_rhos,
    excluded_learner_indices         = excluded_learner_indices,
    intercept_padding                = intercept_padding
  ), class = "gb_learner_ctrl")
}


#' Control parameters for the parametric bootstrap
#'
#' @description
#' Settings of the parametric bootstrap of [bootstrap_gbgam()]: the number
#' of replicates, which steps of the fit each replicate redoes on its own,
#' and how the intervals are widened for smoothing bias. Pass the result to
#' [gb_uq_ctrl()] as `control`, or give the same arguments to
#' [bootstrap_gbgam()]. Supplying `max_iter` with a mode other than
#' `"free_stopping"` results in an error.
#'
#' @param n_samples Integer, at least 2. Number of bootstrap replicates.
#'   Default `200`.
#' @param mode Character. Which steps of the fit each replicate redoes on
#'   its own; the steps it does not redo are copied from the main fit. One
#'   of:
#'   \describe{
#'     \item{`"fixed"`}{replays the main fit's learner selection, smoothing
#'       parameters and stopping iteration. It needs no tracking of the
#'       effective degrees of freedom (EDF; see [gb_edf_ctrl()]), so it is
#'       the fastest.}
#'     \item{`"free_smoothing"`}{replays the main fit's learner selection
#'       and stopping iteration; each replicate adapts its own smoothing
#'       parameters by adaptive smoothness relaxation (ASR).}
#'     \item{`"free_selection"`}{replays the main fit's stopping iteration;
#'       each replicate selects its own learners and adapts its own
#'       smoothing parameters by ASR.}
#'     \item{`"free_stopping"`}{each replicate selects its own learners,
#'       adapts its own smoothing parameters by ASR and chooses its own
#'       stopping iteration by the main fit's stopping method.}
#'   }
#'   All modes but `"fixed"` need EDF tracking in the main fit;
#'   `"free_stopping"` also needs `stopping_method = "aic"` or `"true_f"`.
#'   Default `"fixed"`.
#' @param alpha Number in (0, 0.5], or NULL. The miscoverage at which the
#'   intervals are widened for smoothing bias (see Widening in
#'   [bootstrap_gbgam()]): intervals at level 1 - `alpha`. `NULL` leaves
#'   `$Vp` unwidened. Default `0.05`.
#' @param correction Character. How many widening factors are found:
#'   `"model"`, one for the whole model (Wood, 2020); `"smooth"`, one per
#'   smooth and one per parametric coefficient; `"learner"`, one per
#'   learner. `"smooth"` and `"learner"` apply the rule of `"model"` to each
#'   group; factors above 3 or below 0.9 give a warning that names the
#'   groups (see Widening in [bootstrap_gbgam()]). Not used with
#'   `alpha = NULL`. Default `"model"`.
#' @param phi Positive number, or NULL. The dispersion with which the
#'   replicates are simulated for the families with a free dispersion
#'   (Gaussian, Gamma and inverse Gaussian); the other families ignore it.
#'   Default `NULL`: the fitted dispersion `mod$scale`, Fletcher's (2012)
#'   estimate with \eqn{\widetilde{\operatorname{edf}}}{edf~} model degrees
#'   of freedom (see Dispersion in Details of [bootstrap_gbgam()]).
#' @param max_iter Integer, at least 1, or NULL. Only with
#'   `mode = "free_stopping"`. The iteration cap of each replicate; a
#'   replicate that reaches it before its stopping method stops it is counted
#'   as censored (`$boot_censored` of the bootstrap result). Default `NULL`:
#'   the main fit's `max_iter`.
#'
#' @return A list of class `"gb_param_boot_ctrl"`.
#'
#' @seealso [bootstrap_gbgam()], [gb_uq_ctrl()], [smooth_intervals()]
#'
#' @references
#' Wood, S. N. (2020).
#' Rejoinder on: Inference and computation with generalized additive models and
#' their extensions.
#' *TEST, 29*(2), 354--358.
#' \href{https://doi.org/10.1007/s11749-020-00716-0}{doi:10.1007/s11749-020-00716-0}
#'
#' @examples
#' gb_param_boot_ctrl()
#' gb_param_boot_ctrl(n_samples = 500, mode = "free_stopping",
#'                    correction = "smooth")
#' @export
gb_param_boot_ctrl <- function(
    n_samples  = 200L,
    mode       = c("fixed", "free_smoothing", "free_selection",
                   "free_stopping"),
    alpha      = 0.05,
    correction = c("model", "smooth", "learner"),
    phi        = NULL,
    max_iter   = NULL) {

  checkmate::assert_integerish(n_samples, lower = 2L, len = 1L,
                               any.missing = FALSE, .var.name = "n_samples")
  check_choice(mode, eval(formals()$mode), "mode")
  mode <- match.arg(mode)
  checkmate::assert_number(alpha, lower = 0, upper = 1, null.ok = TRUE,
                           .var.name = "alpha")
  if (!is.null(alpha) && alpha <= 0)
    stop("alpha must be positive, or NULL.", call. = FALSE)
  if (!is.null(alpha) && alpha > 0.5)  # most likely a level given as alpha
    stop("alpha = ", alpha, " is above 0.5; alpha is the miscoverage, ",
         "e.g. 0.05 for 95% intervals.", call. = FALSE)
  correction <- match.arg(correction)
  checkmate::assert_number(phi, lower = 0, finite = TRUE, null.ok = TRUE,
                           .var.name = "phi")
  if (!is.null(phi) && phi <= 0) stop("phi must be positive.", call. = FALSE)
  checkmate::assert_integerish(max_iter, lower = 1, len = 1,
                               any.missing = FALSE, null.ok = TRUE,
                               .var.name = "max_iter")
  if (!is.null(max_iter) && mode != "free_stopping")
    stop("max_iter is only used with mode = \"free_stopping\".",
         call. = FALSE)
  structure(
    list(n_samples = n_samples, mode = mode, alpha = alpha,
         correction = correction, phi = phi, max_iter = max_iter),
    class = "gb_param_boot_ctrl")
}


#' Control parameters for Bayesian backfitting
#'
#' @description
#' Settings of [bayesian_backfitting()]: the number of posterior draws and
#' of burn-in sweeps, the prior of the unpenalised learners, and the
#' undersmoothing of learners with a low edf. Pass the result to
#' [gb_uq_ctrl()] as `control`, or give the same arguments to
#' [bayesian_backfitting()].
#'
#' @param n_samples Integer, at least 2. Number of posterior draws kept
#'   after the burn-in. Default `1000`.
#' @param n_burnin Integer, at least 0. Number of burn-in sweeps, which are
#'   discarded. Default `200`.
#' @param linear Character. The prior of the unpenalised (linear) learners.
#'   One of:
#'   \describe{
#'     \item{`"free"`}{a flat prior (\eqn{\lambda_k = 0}); in each sweep
#'       they are also drawn together with the intercept \eqn{a} in one
#'       joint block (see Details of [bayesian_backfitting()]).}
#'     \item{`"kept"`}{each keeps the smoothing parameter found for it in
#'       Stage 1 (see Details of [bayesian_backfitting()]).}
#'   }
#'   Default `"free"`.
#' @param undersmooth Logical. If `TRUE`, the smoothing parameter of every
#'   penalised learner that boosting selected and whose edf at the
#'   posterior mode is below `edf_floor` is lowered until its edf equals
#'   `edf_floor` (see Details of [bayesian_backfitting()]). Default `FALSE`.
#' @param edf_floor Positive number. The edf to which `undersmooth` lifts
#'   learners; ignored with `undersmooth = FALSE`. Default `1`.
#' @param phi Positive number, or NULL. The dispersion of the posterior for
#'   the families with a free dispersion (Gaussian, Gamma and inverse
#'   Gaussian), e.g. a known one; the other families ignore it. Default
#'   `NULL`: the restricted maximum likelihood estimate at the posterior
#'   mode for the final smoothing parameters (see Dispersion in Details of
#'   [bayesian_backfitting()]).
#'
#' @return A list of class `"gb_bayesian_ctrl"`.
#'
#' @seealso [bayesian_backfitting()], [gb_uq_ctrl()],
#'   [smooth_intervals()]
#'
#' @examples
#' gb_bayesian_ctrl()
#' gb_bayesian_ctrl(n_samples = 2000, linear = "kept")
#' @export
gb_bayesian_ctrl <- function(
    n_samples   = 1000L,
    n_burnin    = 200L,
    linear      = c("free", "kept"),
    undersmooth = FALSE,
    edf_floor   = 1,
    phi         = NULL) {

  checkmate::assert_integerish(n_samples, lower = 2L, len = 1L,
                               any.missing = FALSE, .var.name = "n_samples")
  checkmate::assert_integerish(n_burnin, lower = 0L, len = 1L,
                               any.missing = FALSE, .var.name = "n_burnin")
  linear <- match.arg(linear)
  checkmate::assert_flag(undersmooth, .var.name = "undersmooth")
  checkmate::assert_number(edf_floor, lower = 0, finite = TRUE,
                           .var.name = "edf_floor")
  if (edf_floor <= 0) stop("edf_floor must be positive.", call. = FALSE)
  checkmate::assert_number(phi, lower = 0, finite = TRUE, null.ok = TRUE,
                           .var.name = "phi")
  if (!is.null(phi) && phi <= 0) stop("phi must be positive.", call. = FALSE)
  structure(
    list(n_samples = n_samples, n_burnin = n_burnin, linear = linear,
         undersmooth = undersmooth, edf_floor = edf_floor, phi = phi),
    class = "gb_bayesian_ctrl")
}


#' Control parameters for uncertainty quantification
#'
#' @description
#' Chooses the method that [gbgam()] runs after the fit to attach a
#' covariance matrix of the coefficients (`$Vp`), and holds that method's
#' settings. `$Vp` is used by [vcov()], [predict(se.fit =
#' TRUE)][mgcv::predict.gam], [mgcv::plot.gam] and [smooth_intervals()].
#'
#' @param method Character. One of:
#'   \describe{
#'     \item{`"none"`}{no uncertainty quantification.}
#'     \item{`"bootstrap"`}{the parametric bootstrap of [bootstrap_gbgam()],
#'       with settings from [gb_param_boot_ctrl()].}
#'     \item{`"bayes"`}{Bayesian backfitting ([bayesian_backfitting()]),
#'       which samples the posterior of a penalised model whose smoothing
#'       parameters are chosen so that its posterior mode is close to the
#'       boosted smooths, with settings from [gb_bayesian_ctrl()].}
#'   }
#'   Default `"none"`, or the method of `control` when only `control` is
#'   supplied.
#' @param control A list from [gb_param_boot_ctrl()] (for `"bootstrap"`) or
#'   [gb_bayesian_ctrl()] (for `"bayes"`), or NULL. The settings of the
#'   method. Default `NULL`: the settings in `...`.
#' @param ... Named settings passed to [gb_param_boot_ctrl()] or
#'   [gb_bayesian_ctrl()], whichever `method` names, in place of `control`.
#'
#' @return A list of class `"gb_uq_ctrl"` with the `method` and its
#'   `control` (`NULL` for `"none"`).
#'
#' @details
#' The settings are given either as `control` or in `...`, not both; with
#' neither, the method's defaults are used. When `method` is not given,
#' it is taken from the class of `control`. A `control` made for the other
#' method, an unnamed setting in `...`, a name in `...` that is not a
#' setting of the method's control function, and any settings with
#' `method = "none"` result in an error.
#'
#' @seealso [gb_param_boot_ctrl()], [gb_bayesian_ctrl()],
#'   [bootstrap_gbgam()], [bayesian_backfitting()], [smooth_intervals()]
#'
#' @examples
#' # Default: no uncertainty quantification
#' gb_uq_ctrl()
#'
#' # Parametric bootstrap with the default settings
#' gb_uq_ctrl("bootstrap")
#'
#' # Bootstrap settings in ...
#' gb_uq_ctrl("bootstrap", n_samples = 300, mode = "free_selection")
#'
#' # Bootstrap settings as a control; the method is taken from it
#' gb_uq_ctrl(control = gb_param_boot_ctrl(mode = "free_stopping",
#'                                                correction = "smooth"))
#'
#' # Bayesian backfitting
#' gb_uq_ctrl("bayes", n_samples = 2000, n_burnin = 500)
#'
#' @export
gb_uq_ctrl <- function(method = c("none", "bootstrap", "bayes"),
                       control = NULL, ...) {
  ctors <- c(bootstrap = "gb_param_boot_ctrl", bayes = "gb_bayesian_ctrl")
  if (missing(method) && !is.null(control)) {
    hit <- vapply(ctors, function(cls) inherits(control, cls), logical(1))
    if (!any(hit))
      stop("control must be created with gb_param_boot_ctrl() or ",
           "gb_bayesian_ctrl().", call. = FALSE)
    method <- names(ctors)[hit]
  }
  method <- match.arg(method)
  dots   <- list(...)

  if (method == "none") {
    if (missing(method) && length(dots) > 0L)
      stop("settings in ... need a method: use ",
           "gb_uq_ctrl(\"bootstrap\", ...) or ",
           "gb_uq_ctrl(\"bayes\", ...).", call. = FALSE)
    if (!is.null(control) || length(dots) > 0L)
      stop("method = \"none\" takes no settings (control or ...).",
           call. = FALSE)
    return(structure(list(method = "none", control = NULL),
                     class = "gb_uq_ctrl"))
  }

  ctor <- ctors[[method]]
  if (!is.null(control)) {
    if (length(dots) > 0L)
      stop("supply the settings either as control or in ..., not both.",
           call. = FALSE)
    if (!inherits(control, ctor))
      stop("control must be created with ", ctor, "() for method = \"",
           method, "\".", call. = FALSE)
  } else {
    nm <- names(dots)
    if (length(dots) > 0L && (is.null(nm) || any(nm == "")))
      stop("the settings in ... must be named.", call. = FALSE)
    bad <- setdiff(nm, names(formals(ctor)))
    if (length(bad) > 0L)
      stop(bad[1L], " is not a setting of ", ctor, "() (method = \"",
           method, "\").", call. = FALSE)
    control <- do.call(ctor, dots)
  }
  structure(list(method = method, control = control),
            class = "gb_uq_ctrl")
}

#' The uncertainty quantification settings given to gbgam()
#'
#' Accepts what `gbgam(uq_control = )` takes: a list of arguments of
#' [gb_uq_ctrl()] (or its result), or a control from
#' [gb_param_boot_ctrl()] or [gb_bayesian_ctrl()], alone or as the only
#' element of a list.
#'
#' @param x The `uq_control` argument of [gbgam()].
#' @return A list of class `"gb_uq_ctrl"`.
#' @noRd
as_uq_ctrl <- function(x) {
  ctrl_classes <- c("gb_param_boot_ctrl", "gb_bayesian_ctrl")
  if (inherits(x, ctrl_classes)) return(gb_uq_ctrl(control = x))
  if (!is.list(x))
    stop("uq_control must be a list, as returned by ",
         "gb_uq_ctrl(), gb_param_boot_ctrl() or gb_bayesian_ctrl().",
         call. = FALSE)
  has <- vapply(ctrl_classes, function(cls)
    any(vapply(x, inherits, NA, what = cls)), NA)
  if (all(has))
    stop("uq_control holds both a gb_param_boot_ctrl() and a ",
         "gb_bayesian_ctrl(), but a model holds one covariance ($Vp). ",
         "Give one of them, and run the other method afterwards with ",
         "bootstrap_gbgam() or bayesian_backfitting().", call. = FALSE)
  if (length(x) == 1L && any(has))
    return(gb_uq_ctrl(control = x[[1L]]))
  do.call("gb_uq_ctrl", x)
}

# The bootstrap modes as the sets of steps each replicate redoes on its own.
boot_modes <- list(
  fixed          = character(),
  free_smoothing = "lambda",
  free_selection = c("lambda", "selection"),
  free_stopping  = c("lambda", "selection", "stopping")
)

#' Check that the main fit supports a bootstrap mode
#'
#' The modes other than `"fixed"` need EDF tracking (ASR); `"free_stopping"`
#' also needs stopping by AIC or `true_f`. Called by `gbgam()` before the
#' fit, so that an incompatible bootstrap does not fail after it, and by
#' `bootstrap_gbgam()`.
#'
#' @param mode A bootstrap mode (see [gb_param_boot_ctrl()]).
#' @param stopping_method,edf_compute_method Methods of the main fit.
#' @noRd
check_boot_requirements <- function(mode, stopping_method,
                                    edf_compute_method) {
  if (mode != "fixed" && identical(edf_compute_method, "none")) {
    stop("mode = \"", mode, "\" requires EDF tracking in the main fit: ",
         "fit it with gb_edf_ctrl(\"hutchpp\") or gb_edf_ctrl(\"exact\").",
         call. = FALSE)
  }
  if (mode == "free_stopping" && identical(stopping_method, "cv")) {
    stop("mode = \"free_stopping\" is not available with ",
         "stopping_method = \"cv\": the replicates would need their own ",
         "resampling fits. Use mode = \"free_selection\".", call. = FALSE)
  }
  if (mode == "free_stopping" && identical(stopping_method, "none")) {
    stop("mode = \"free_stopping\" requires the main fit to stop by ",
         "stopping_method = \"aic\" or \"true_f\".", call. = FALSE)
  }
  invisible(NULL)
}
