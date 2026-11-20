#' Intervals for the intercept plus each smooth
#'
#' @description
#' Intervals for the intercept `a` plus each smooth at each evaluation
#' point, as numbers, from the covariance attached by [bootstrap_gbgam()] or
#' [bayesian_backfitting()] (or by [gbgam()] with an `uq_control`). This
#' function is optional: [plot()][plot.gbgam()] and
#' [predict(se.fit = TRUE)][predict.gbgam()] use that covariance without it.
#' It is for what they do not give:
#'
#' - bootstrap intervals at a level other than the bootstrap's own `alpha`
#'   (`$Vp` is widened for one level only; this function finds the widening
#'   again at its own `alpha`), or with another `correction`;
#' - Bayesian backfitting intervals centred on the posterior mean, which
#'   [plot()][plot.gbgam()] and [predict()][predict.gbgam()] do not use;
#' - the limits of the intervals in matrices, to transform them, tabulate
#'   them or draw them in a figure of your own.
#'
#' The coverage of the intervals holds across the function (Marra & Wood,
#' 2012): on average over the evaluation points, about 1 - `alpha` of the
#' intervals contain the true effect, while at a single point the coverage
#' can be lower (or higher). Each interval includes the uncertainty of the
#' intercept and of no other term.
#'
#' @param mod A fitted model of class `"gbgam"`, as returned by [gbgam()],
#'   with `$Vp`.
#' @param alpha Number in (0, 0.5]. The miscoverage: the intervals have
#'   level 1 - `alpha` across the function. Default `0.05` (95%
#'   intervals).
#' @param centre Character. Where the intervals are centred: `"boosted"`,
#'   on the boosting estimate `mod$coefficients`; `"posterior"`, on the
#'   posterior mean `mod$coefficients_bayes` (Bayesian backfitting only).
#'   Default `"boosted"`.
#' @param correction Character, or NULL. How the intervals of a bootstrap
#'   fit are widened for smoothing bias: `"none"`, no widening; `"model"`,
#'   `"smooth"` or `"learner"`, as `correction` in [gb_param_boot_ctrl()],
#'   with the widening factors found at `alpha`. A Bayesian backfitting fit
#'   takes only `NULL` or `"none"`. Default `NULL`: the bootstrap's own
#'   `correction`, or `"none"` if it was run with `alpha = NULL` or the fit
#'   is from Bayesian backfitting.
#' @param newdata Data frame, or NULL. The data at which the smooths are
#'   evaluated. Default `NULL`: the data the model was fitted to.
#'
#' @return
#' A list with
#' \describe{
#'   \item{`fit`, `se`, `lower`, `upper`}{Matrices with one column per
#'     smooth, named by the smooth's label, and one row per row of
#'     `newdata`: the centre (intercept included), the standard error
#'     (widened, if a widening is used) and the interval's limits. Without
#'     `newdata`, one row per element of `fitted(mod)`, with `NA` in the
#'     rows excluded by `na.action = na.exclude`.}
#'   \item{`alpha`, `centre`}{The settings used.}
#'   \item{`factor`}{The widening factor(s) used, as `$widening_factor` of
#'     [bootstrap_gbgam()]; `NULL` without a widening.}
#' }
#'
#' @details
#' ## Intervals
#'
#' The centre and the standard error of smooth \eqn{j} are computed from
#' the columns of the intercept and of smooth \eqn{j} in the model matrix
#' (`type = "lpmatrix"`), for every smooth, `by`-variable smooths included.
#' The interval is the centre plus or minus the normal quantile
#' \eqn{\Phi^{-1}(1 - \frac{\alpha}{2})}{Phi^-1(1 - alpha/2)} times the
#' standard error, with \eqn{\Phi^{-1}}{Phi^-1} the standard normal quantile
#' function. For
#' models whose only parametric term is the intercept and whose smooths are
#' centred, the intervals with `centre = "boosted"` equal those of
#' `predict(type = "iterms", se.fit = TRUE)` shifted by the intercept, and
#' those drawn by `plot(mod, seWithMean = TRUE, se = qnorm(1 - alpha / 2),
#' shift = coef(mod)[1])`, whenever the covariance used is `$Vp`. This holds
#' for Bayesian backfitting, and for a bootstrap fit with `correction = NULL`
#' at the bootstrap's own `alpha` (at any `alpha` if the bootstrap was run with
#' `alpha = NULL`).
#'
#' ## Standard errors and widening
#'
#' For a Bayesian backfitting fit the standard errors come from `$Vp`. For
#' a bootstrap fit they come from the replicates stored by
#' [bootstrap_gbgam()]: with `correction = "none"` from the replicates'
#' plain covariance, otherwise widened by factors found again at `alpha`
#' (see Widening in [bootstrap_gbgam()]), so intervals at several levels
#' are each widened for their own level. At the bootstrap's own `alpha`
#' and `correction`, the factors are `$widening_factor`.
#'
#' @seealso [bootstrap_gbgam()], [bayesian_backfitting()],
#'   [gb_uq_ctrl()]
#'
#' @references
#' Marra, G., & Wood, S. N. (2012).
#' Coverage properties of confidence intervals for generalized additive model
#' components.
#' *Scandinavian Journal of Statistics, 39*(1), 53--74.
#' \href{https://doi.org/10.1111/j.1467-9469.2011.00760.x}{doi:10.1111/j.1467-9469.2011.00760.x}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x   <- runif(150)
#' dat <- data.frame(x = x, y = sin(2 * pi * x) + rnorm(150, sd = 0.3))
#' mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = dat)
#' mod <- bootstrap_gbgam(mod, n_samples = 50)
#'
#' # plot() draws intervals without smooth_intervals(). With these settings
#' # they are 95% intervals for the intercept plus the smooth
#' a <- coef(mod)[["(Intercept)"]]
#' plot(mod, seWithMean = TRUE, se = qnorm(0.975), shift = a)
#'
#' # smooth_intervals() adds another level: 80%, with its own widening
#' grid <- data.frame(x = seq(min(x), max(x), length.out = 200))
#' ci   <- smooth_intervals(mod, alpha = 0.2, newdata = grid)
#' matlines(grid$x, cbind(ci$lower, ci$upper), lty = 2, col = "firebrick")
#' legend("topright", legend = c("95%", "80%"), lty = 2,
#'        col = c("black", "firebrick"), bty = "n")
#'
#' # The limits as numbers
#' head(cbind(ci$lower, ci$upper))
#' }
#' @export
smooth_intervals <- function(mod, alpha = 0.05,
                             centre = c("boosted", "posterior"),
                             correction = NULL, newdata = NULL) {
  checkmate::assert_class(mod, "gbgam")
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  if (alpha > 0.5)  # most likely a level given as alpha
    stop("smooth_intervals: alpha = ", alpha, " is above 0.5; alpha is ",
         "the miscoverage, e.g. 0.05 for 95% intervals.", call. = FALSE)
  if (alpha <= 0)
    stop("smooth_intervals: alpha must be positive.", call. = FALSE)
  centre <- match.arg(centre)
  if (is.null(mod$Vp))
    stop("smooth_intervals: the model has no $Vp; run bootstrap_gbgam() or ",
         "bayesian_backfitting() first.", call. = FALSE)
  if (length(mod$smooth) == 0L)
    stop("smooth_intervals: the model has no smooth terms.", call. = FALSE)
  ic      <- mod$uq_control
  is_boot <- identical(ic$method, "bootstrap")
  if (!is.null(correction)) {
    checkmate::assert_choice(correction,
                             c("none", "model", "smooth", "learner"),
                             .var.name = "correction")
    if (!is_boot && correction != "none")
      stop("smooth_intervals: correction = \"", correction, "\" needs a ",
           "bootstrap fit (bootstrap_gbgam()).", call. = FALSE)
  }
  if (centre == "posterior" && is.null(mod$coefficients_bayes))
    stop("smooth_intervals: centre = \"posterior\" needs Bayesian ",
         "backfitting ($coefficients_bayes).", call. = FALSE)

  corr <- correction %||%
    if (is_boot && !is.null(ic$control$alpha)) ic$control$correction else
      "none"
  fac    <- NULL
  m2     <- mod
  if (is_boot && corr == "none") {
    m2$Vp <- stats::cov(mod$boot_coef)  # the replicates' plain covariance
  } else if (corr != "none") {
    inf   <- widen_covariance(mod, alpha, corr)
    fac   <- inf$factor
    m2$Vp <- inf$Vp
  }
  b  <- if (centre == "posterior") mod$coefficients_bayes else
    mod$coefficients
  X  <- if (is.null(newdata)) stats::predict(mod, type = "lpmatrix") else
    stats::predict(mod, newdata = newdata, type = "lpmatrix")
  i0 <- which(names(b) == "(Intercept)")
  labs <- vapply(mod$smooth, `[[`, "", "label")
  # a + f_j from the intercept's and smooth j's columns only, for any smooth
  fit <- se <- matrix(NA_real_, nrow(X), length(labs),
                      dimnames = list(NULL, labs))
  for (j in seq_along(labs)) {
    idx      <- c(i0, mod$smooth[[j]]$first.para:mod$smooth[[j]]$last.para)
    Xj       <- X[, idx, drop = FALSE]
    fit[, j] <- Xj %*% b[idx]
    se[, j]  <- sqrt(pmax(0, rowSums((Xj %*% m2$Vp[idx, idx]) * Xj)))
  }
  z    <- stats::qnorm(1 - alpha / 2)
  list(fit = fit, se = se, lower = fit - z * se, upper = fit + z * se,
       alpha = alpha, centre = centre, factor = fac)
}
