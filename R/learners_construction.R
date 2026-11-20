#' Initialise base learners for boosted generalised additive models
#'
#' @description
#' Builds the base learners of a GAM formula, as [gbgam()] does: the
#' intercept learner, one unpenalised learner per centred column of the other
#' parametric terms, and for each [mgcv::t2()] smooth, once its basis columns
#' are centred and those without variance removed, one penalised learner per
#' penalty on at least two columns and one unpenalised learner per other
#' column (see Details). All penalised learners start with the same selection
#' edf (sedf, \eqn{2 \operatorname{tr}(\mathbf{A}) -
#' \operatorname{tr}(\mathbf{A}^2)}{2 tr(A) - tr(A^2)}), so that learner
#' selection is approximately unbiased (Hofner et al., 2011).
#' Only `t2()` smooths are supported: `s()`, `te()` and
#' `ti()` terms result in an error.
#'
#' @param formula Formula. A GAM formula as in [mgcv::formula.gam], such as
#'   `y ~ t2(x1) + x2 + t2(x3, x4)`. An intercept is required: formulae with
#'   `~ 0 + ...` or `~ -1 + ...` result in an error. The `sp` argument of
#'   [mgcv::t2()] sets the starting smoothness of a penalised learner while
#'   keeping learner selection unbiased:
#'   * `sp = NA` (the default of `t2()`): the learner starts at sedf = 1,
#'     with step size \eqn{\nu_k = 1}.
#'   * `sp > 0`: `sp` is the smoothing parameter lambda itself; unlike in
#'     `mgcv`, the penalty is not rescaled (`S.scale = 1`). It must not give
#'     sedf below one.
#'   * `sp < 0`: `-sp` is the sedf of the learner without its step size; it
#'     must lie strictly between 1 and \eqn{r}
#'     (\eqn{1 < -\mathrm{sp} < r}{1 < -sp < r}), where \eqn{r} is the number
#'     of identifiable directions of the centred, weighted penalised group of
#'     columns (its number of columns unless some are removed, see Details).
#'     The step size \eqn{\nu_k} then brings the sedf of the scaled learner
#'     back to one.
#'
#'   The values of `sp` are matched to the penalties of a `t2()` term in the
#'   order of its penalty matrices. An `sp` given for a penalty on a single
#'   column is ignored, because that column is fitted as an unpenalised
#'   learner.
#' @param data Data frame or list. The variables of `formula`. A list may
#'   hold matrices, for example for smooths with matrix arguments. All
#'   elements have the same number of observations: rows for matrices,
#'   length for vectors. Default `list()`: the variables are taken from the
#'   environment of `formula`.
#' @param weights Numeric vector, or NULL. Prior weights, one per
#'   observation, used for the weighted centring, the Demmler-Reinsch basis
#'   of the penalised learners and the starting smoothness. Rows removed by
#'   the `NA` handling of the model frame are removed from `weights` as well.
#'   Default `NULL`: all weights are one.
#' @param knots Named list, or NULL. Knot positions passed to
#'   [mgcv::smoothCon()]; each name must be a covariate of a smooth.
#'   Default `NULL`: `mgcv`'s knots.
#'
#' @return
#'   A named list with seven elements:
#'   \describe{
#'     \item{`nsdf`}{Number of parametric (non-smooth) coefficients, including
#'       the intercept, as in [`mgcv`][mgcv::mgcv-package].}
#'     \item{`smooths`}{The list of `mgcv` smooth objects from
#'       [mgcv::smoothCon()], with `first.para` and `last.para` set. Smooths
#'       without a by-variable or matrix arguments get the class
#'       `"gbgam_centred"`, their centring constants in `$gbgam_cen` and a
#'       reference to the `gbgam` namespace in `$gbgam_ns`, so that `mgcv`
#'       predicts from the centred basis (see
#'       [Predict.matrix.gbgam_centred()]). Smooths with a by-variable and
#'       functional (matrix-argument) smooths are not centred in their
#'       prediction basis.}
#'     \item{`learners`}{A flat list of base learners in this order: the
#'       intercept, the other parametric columns, then the smooths in
#'       formula order (for each, its penalty groups in the order of its
#'       penalty matrices, then its unpenalised columns).}
#'     \item{`coef_names`}{A character vector of coefficient names, as in
#'      [`mgcv`][mgcv::mgcv-package].}
#'     \item{`ind0var`}{An integer vector of the positions in the full
#'       coefficient vector (parametric terms included) of the columns
#'       without variance that were dropped.}
#'     \item{`cen_const`}{A numeric vector of the weighted column means used
#'       for the sum-to-zero constraints, in the same order as
#'       `coef_names`. The intercept entry is zero.}
#'     \item{`cen_absorb`}{A logical vector in the same order as
#'       `coef_names`: `TRUE` for the columns whose centring shift is absorbed
#'       into the `"(Intercept)"` coefficient
#'       (parametric columns and smooths with a by-variable or matrix
#'       arguments), `FALSE` for the columns of `"gbgam_centred"` smooths.}
#'     }
#'
#' @details
#' The intercept is one unpenalised learner. Every other column of the
#' parametric [model matrix][stats::model.matrix] (e.g. each contrast column
#' of a factor) is centred on its prior-weighted mean and becomes one
#' unpenalised learner. The `t2()` terms are prepared in these steps:
#' 1. [mgcv::smoothCon()] builds their basis matrices.
#' 1. Each column is centred on its prior-weighted mean (the sum-to-zero
#'    constraint), and columns without variance are removed, such as the
#'    constant column of the null space (the direction of the intercept).
#'    For a factor `by` variable the variance is taken over the rows of the
#'    smooth's level, so the level's constant column is removed.
#' 1. The columns are split into one group per penalty and one column per
#'    unpenalised direction (Wood et al., 2013), so that each penalty is
#'    selected and relaxed on its own during boosting. Each unpenalised
#'    column, and each penalty group of a single column (e.g. the `"rr"`
#'    group of `t2(x, z, k = c(3, 3), bs = "cr")`), becomes one unpenalised
#'    learner.
#' 1. Each penalty group of at least two columns becomes one penalised
#'    learner, with the Demmler and Reinsch (1975)
#'    reparametrisation (see also Wood, 2017), computed from a pivoted QR
#'    decomposition of the weighted group and an eigendecomposition. If the
#'    group is rank deficient or nearly so, a truncated singular value
#'    decomposition is used instead, which removes the directions that the
#'    rows with positive weight do not identify (e.g. the aliased direction
#'    of centred random-effect dummies, or the columns of a factor level
#'    whose weights are all zero). A group with fewer than two identifiable
#'    directions results in an error.
#' 1. The smoothing parameter lambda = exp(`rho`) and the step size
#'    \eqn{\nu_k} of each penalised learner are set so that its sedf is one,
#'    where the sedf of a learner with hat matrix \eqn{\mathbf{A}_k}{A_k} is
#'    \eqn{2 \nu_k \operatorname{tr}(\mathbf{A}_k) - \nu_k^2
#'    \operatorname{tr}(\mathbf{A}_k^2)}{2 nu_k tr(A_k) - nu_k^2 tr(A_k^2)}.
#'    Then all learners give the same expected reduction of the residual
#'    sum of squares when none has a true effect.
#'
#' For smooths without a by-variable or matrix arguments, the centred basis
#' is also used by `mgcv`'s prediction and plotting methods, so their
#' term-wise predictions are the centred effects: the prior-weighted mean of
#' each effect over the training data is zero.
#'
#' ## Fields of a learner object
#' Every learner (class `"learner"`) holds:
#' * `var_label`: the learner's name;
#' * `beta`: its coefficients (on the Demmler-Reinsch basis `design` for
#'   penalised learners);
#' * `edf`: its accumulated edf
#'   \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} (the increases of the
#'   model edf credited to it while boosting);
#' * `eta`: its contribution to the linear predictor;
#' * `design`: its design matrix (for penalised learners the Demmler-Reinsch
#'   basis V = Z A, with Z the centred penalised group of columns).
#'
#' Unpenalised learners (class `c("unpenalised", "learner")`) also hold
#' `xtWx`, the weighted squared norm of their single column; those of
#' smooths also hold `term`, `label` and `by` of their smooth.
#'
#' Penalised learners (class `c("demmler_reinsch", "learner")`) also hold
#' `term`, `label` and `by` of their smooth, and:
#' * `A`: the Demmler-Reinsch map, the matrix that turns Z into the basis V
#'   (not the hat matrix);
#' * `D`: the eigenvalues of the penalty on the basis V, decreasing;
#' * `rho`: log(lambda); `rho0`: its value at construction;
#' * `lambda_at_zero`: `TRUE` once lambda has been set to zero (`rho` is then
#'   a large negative number);
#' * `nu_k`: the learner-specific step size; `nu_k0`: its value at
#'   construction.
#'
#' @seealso [gbgam()].
#'
#' @references
#' Demmler, A., & Reinsch, C. (1975).
#' Oscillation matrices with spline smoothing.
#' *Numerische Mathematik, 24*(5), 375--382.
#' \href{https://doi.org/10.1007/BF01437406}{doi:10.1007/BF01437406}
#'
#' Hofner, B., Hothorn, T., Kneib, T., & Schmid, M. (2011).
#' A framework for unbiased model selection based on boosting.
#' *Journal of Computational and Graphical Statistics, 20*(4), 956--971.
#' \href{https://doi.org/10.1198/jcgs.2011.09220}{doi:10.1198/jcgs.2011.09220}
#'
#' Wood, S. N. (2017).
#' *Generalized additive models: An introduction with R* (2nd ed.). Chapman and
#' Hall/CRC.
#'
#' Wood, S. N., Scheipl, F., & Faraway, J. J. (2013).
#' Straightforward intermediate rank tensor product smoothing in mixed models.
#' *Statistics and Computing, 23*(3), 341--360.
#' \href{https://doi.org/10.1007/s11222-012-9314-z}{doi:10.1007/s11222-012-9314-z}
#'
#' @examples
#' set.seed(1)
#' n <- 80
#' dat <- data.frame(
#'   y   = rnorm(n),
#'   x   = runif(n),
#'   z   = runif(n),
#'   lin = rnorm(n)
#' )
#'
#' # Univariate smooth + linear term
#' lb <- new_learners(
#'   y ~ lin + t2(x, k = 6, bs = "cr"),
#'   data    = dat,
#'   weights = rep(1, n)
#' )
#' names(lb) # "nsdf", "smooths", "coef_names", "learners", "ind0var",
#'           # "cen_const", "cen_absorb"
#' length(lb$learners) # intercept, lin, then the learners of t2(x)
#'
#' # Bivariate tensor-product interaction
#' lb2 <- new_learners(
#'   y ~ t2(x, z, k = c(4, 4), bs = c("cr", "cr")),
#'   data    = dat,
#'   weights = rep(1, n)
#' )
#' sapply(lb2$learners, function(l) l$var_label)
#'
#' @export
new_learners <- function(formula, data = list(), weights = NULL,
                         knots = NULL) {
  check_formula(formula)
  gp  <- mgcv::interpret.gam(formula)
  # Pass the evaluated weights into the model frame, so that rows dropped by
  # the NA handling are dropped from the weights as well
  mf_call <- make_mf_call(
    formula            = gp$fake.formula,
    call               = match.call(expand.dots = FALSE),
    values             = list(weights = weights),
    drop.unused.levels = FALSE
  )
  mf <- eval(mf_call, parent.frame())

  pmf_call          <- mf_call
  pmf_call$formula  <- gp$pf
  pmf               <- eval(pmf_call, parent.frame())
  pterms            <- attr(pmf, "terms")

  weights <- stats::model.weights(mf) %||% rep(1, nrow(mf))
  checkmate::assert_numeric(weights, lower = 0, finite = TRUE,
                            .var.name = "weights")
  smooths <- construct_smooths(mf, gp$smooth.spec, knots)
  new_learners_internal(
    smooths = smooths, data = mf, weights = weights,
    pterms = pterms, ahat = 0, edf0 = 0, knots = knots
  )
}

#' An internal function for creating learners given the smooth objects from
#' `mgcv`
#'
#' @param smooths list of smooth objects from [mgcv::smoothCon()].
#' @param data a `data.frame` or model frame containing the covariates.
#' @param weights numeric vector of observation weights (prior weights).
#' @param pterms `terms` object for the parametric part of the formula.
#' @param ahat numeric scalar; initial linear predictor value (e.g. the MLE
#'   intercept), absorbed into the intercept learner's `$beta`.
#' @param edf0 numeric scalar;
#'   \eqn{\widetilde{\operatorname{edf}}^{[0]}}{edf~^[0]} assigned to the
#'   intercept learner.
#' @param knots named list of knot positions passed through to
#'   [mgcv::smoothCon()].
#' @param clamp_sp logical; clamp a fixed `sp` that would give a penalised
#'   learner sedf < 1 under these weights to sedf = 1 instead of stopping
#'   (used for resampling fits, whose smaller total weight lowers the
#'   smoothing parameter at sedf = 1).
#' @noRd
new_learners_internal <- function(smooths, data, weights, pterms, ahat, edf0,
                                  knots, clamp_sp = FALSE) {
  # Centre the parametric model matrix; one base learner per column
  pX         <- stats::model.matrix(pterms, data)
  nsdf       <- ncol(pX) # For mgcv compatibility
  coef_names <- colnames(pX)
  lin_obj    <- new_linear_learners(pX, ahat = ahat, edf0 = edf0,
                                    weights = weights)
  lrs        <- lin_obj$learners
  cen_const  <- lin_obj$cen_const

  first.para <- ncol(pX) + 1L
  ind0var    <- integer(0)
  # Parametric centring shifts are absorbed into the intercept coefficient
  cen_absorb <- rep(TRUE, ncol(pX))
  # Centre, drop-zero-var, and split its columns into separate base learners
  for(j in seq_along(smooths)){
    sm <- smooths[[j]]

    obj <- smooth_to_learners(sm, weights = weights, clamp_sp = clamp_sp)
    lrs <- append(lrs, obj$learners)

    # Save objects for coef back-transform. Smooths centred in their
    # prediction basis keep their centring out of the intercept.
    cen_const    <- append(cen_const, obj$cen_const)
    centred      <- centre_in_basis(sm, data)
    cen_absorb   <- c(cen_absorb, rep(!centred, ncol(sm$X)))
    smooths[[j]] <- set_smooth_centring(sm, obj$cen_const, centred)
    # Positions of the dropped columns in the full coefficient vector
    model_level_ind0var     <- first.para - 1L + obj$ind0var
    ind0var                 <- append(ind0var, model_level_ind0var)
    coef_names              <- c(coef_names,
                                 paste0(sm$label, ".", seq_len(ncol(sm$X))))
    smooths[[j]]$first.para <- first.para
    smooths[[j]]$last.para  <- first.para + ncol(sm$X) - 1L
    first.para              <- first.para + ncol(sm$X)
  }

  list(
    nsdf       = nsdf,
    smooths    = smooths,
    coef_names = coef_names,
    learners   = lrs,
    ind0var    = ind0var,    # Zero-variance columns
    cen_const  = cen_const,  # Weighted means (0 for the intercept)
    cen_absorb = cen_absorb  # Centring shifts absorbed into the intercept
  )
}


#' Should a smooth be centred in its prediction basis?
#'
#' The weighted centring of a smooth's basis columns is carried by the smooth
#' itself (see `set_smooth_centring()`) only when `mgcv::PredictMat()` would
#' apply it in the right place, i.e. when no by-variable multiplication or
#' summation convention follows `mgcv::Predict.matrix()`:
#' * With a by-variable, the package centres the by-multiplied columns.
#'   A numeric-by smooth \eqn{z f(x)} is identifiable without centring, so
#'   its curve \eqn{f} is reported uncentred, as in `mgcv`. For a factor-by
#'   smooth, the level-specific mean is either identifiable (no factor main
#'   effect) or aliased with the factor main effect, and centring over all
#'   rows would not remove it from the level curve either.
#' * With matrix arguments (linear functional terms), the basis rows are
#'   summed after `mgcv::Predict.matrix()`.
#'
#' For these smooths the centring shift stays absorbed into the intercept
#' coefficient.
#'
#' @param sm An `mgcv` smooth object.
#' @param data The model frame the smooth was constructed from.
#' @return `TRUE` if the centring constants are subtracted in the basis.
#' @noRd
centre_in_basis <- function(sm, data) {
  if (!identical(sm$by, "NA")) return(FALSE)
  is_mat <- vapply(sm$term, function(v) is.matrix(data[[v]]), logical(1L))
  !any(is_mat)
}


#' Attach the centring of a smooth to the smooth object
#'
#' For a centred smooth, the weighted column means `cen_const` of its basis
#' are stored in `$gbgam_cen` and the class `"gbgam_centred"` is prepended,
#' so that `Predict.matrix.gbgam_centred()` (and hence `mgcv::PredictMat()`,
#' `predict.gam()` and `plot.gam()`) returns the centred basis used during
#' boosting. The attribute `"nCons"` (number of identifiability constraints,
#' set by `mgcv::smoothCon()` only when it absorbs constraints) is 1 for
#' centred smooths and 0 otherwise, so that `plot.gam(seWithMean = TRUE)` and
#' `predict.gam(type = "iterms")` work.
#'
#' A centred smooth also stores the `gbgam` namespace in `$gbgam_ns`. A
#' namespace is serialised by reference, so `readRDS()` of a saved model
#' calls `loadNamespace("gbgam")`, which registers
#' `Predict.matrix.gbgam_centred()` also in a session where only `mgcv` is
#' attached. Without it, `mgcv::Predict.matrix()` would silently fall back to
#' the uncentred basis there.
#'
#' @param sm An `mgcv` smooth object.
#' @param cen_const Numeric vector of weighted column means of `sm$X`.
#' @param centred Logical; result of `centre_in_basis()`.
#' @return `sm`, modified.
#' @noRd
set_smooth_centring <- function(sm, cen_const, centred) {
  if (centred) {
    sm$gbgam_cen <- unname(cen_const)
    sm$gbgam_ns  <- topenv(environment())
    class(sm)    <- c("gbgam_centred", class(sm))
  }
  attr(sm, "nCons") <- as.integer(centred)
  sm
}


#' Prediction matrix of a centred `gbgam` smooth
#'
#' Method for [mgcv::Predict.matrix()], used by [mgcv::PredictMat()] for
#' smooths of class `"gbgam_centred"`. It evaluates the uncentred basis with
#' the next method (e.g. that of `"t2.smooth"`) and subtracts the weighted
#' column means of the training basis stored in `object$gbgam_cen`, so that
#' `predict()` and `plot()` use the same sum-to-zero centred basis as the
#' boosting learners. The constants are not subtracted when a derivative of
#' the basis is requested (`deriv` set on the smooth or one of its margins),
#' because the derivative of a constant is zero.
#'
#' `gbgam` must therefore be installed wherever a saved model is used with
#' `predict()` or `plot()`. Centred smooths keep a reference to the `gbgam`
#' namespace (`$gbgam_ns`), so `readRDS()` of a saved model loads `gbgam` and
#' registers this method automatically, also when only `mgcv` is attached.
#'
#' @param object A smooth of class `"gbgam_centred"`, as stored in `$smooth`
#'   of a [gbgam()] fit.
#' @param data A data frame or list of covariate values.
#' @return The centred prediction matrix.
#' @keywords internal
#' @exportS3Method mgcv::Predict.matrix
Predict.matrix.gbgam_centred <- function(object, data) {
  X <- NextMethod()
  if (smooth_has_deriv(object)) return(X)
  X - rep(object$gbgam_cen, each = nrow(X))
}


#' Whether a derivative of a smooth's basis is requested
#'
#' @param sm An `mgcv` smooth object.
#' @return `TRUE` if `sm$deriv` or the `deriv` of any margin is non-zero.
#' @noRd
smooth_has_deriv <- function(sm) {
  is_deriv <- function(d) !is.null(d) && any(d != 0)
  margin_deriv <- vapply(sm$margin %||% list(),
                         function(m) is_deriv(m$deriv), logical(1L))
  is_deriv(sm$deriv) || any(margin_deriv)
}


#' Initialise learners from parametric terms (intercept required)
#'
#' @param pX model.matrix of parametric terms. Must include the
#'   "(Intercept)" column as the first column.
#' @param ahat initial linear predictor values, typically \eqn{g(\bar{y})},
#'   where \eqn{g} is the link function.
#' @param edf0 \eqn{\widetilde{\operatorname{edf}}^{[0]}}{edf~^[0]}
#'   (typically 1, after initialising the intercept).
#' @param weights internal weights used by the model (i.e. `prior.weights`).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`learners`}{A list of base learner objects of class
#'       `unpenalised`, one per column of `pX`.}
#'     \item{`cen_const`}{A named numeric vector of weighted column means
#'       used for centring, with the intercept entry set to zero.}
#'   }
#' @noRd
new_linear_learners <- function(pX, ahat, edf0, weights){
  stopifnot(is.matrix(pX))
  stopifnot(colnames(pX)[1]=="(Intercept)")

  # Centre pX
  cen_const <- apply(pX, 2, FUN=stats::weighted.mean, weights, simplify = TRUE)
  cen_const[1] <- 0 # Don't centre the intercept
  cen_pX <- t(t(pX)-cen_const)

  # Create learners
  learners <- vector("list", ncol(pX))
  for (j in seq_len(ncol(pX))) {
    vname <- colnames(cen_pX)[j]
    learners[[j]] <- new_unpenalised_learner(
      var_label = vname,
      x = cen_pX[,j],
      weights = weights
    )

    # Absorb offset into the intercept
    if(vname=="(Intercept)"){
      learners[[j]]$beta <- ahat
      learners[[j]]$edf <- edf0
    }
  }
  return(list(learners=learners, cen_const=cen_const))
}



#' Initialise learners from `mgcv` smooth objects
#'
#' @param sm one of the elements generated from [mgcv::smoothCon] using a
#'   [mgcv::t2] constructor.
#' @param weights vector of internal weights (`prior.weights`)
#' @param tol numeric; variance threshold below which a basis column is
#'   considered zero-variance and dropped. Default `1e-8`.
#' @param clamp_sp logical; clamp a fixed `sp` giving sedf < 1 to sedf = 1
#'   instead of stopping (see `new_learners_internal()`).
#'
#' @return list of `learners` (penalised or not), `cen_const` (the weighted
#'   column means used for the sum-to-zero constraints) and `ind0var` (the
#'   columns of `sm$X` with zero variance, which are dropped).
#' @noRd
smooth_to_learners <- function(sm, weights, tol = 1e-8, clamp_sp = FALSE) {
  stopifnot(inherits(sm, "t2.smooth"))
  if (nrow(sm$X) < ncol(sm$X)) {
    stop("The basis of ", sm$label, " has more columns than rows.")
  }
  # Centre
  w_means <- apply(sm$X, 2, stats::weighted.mean, w = weights,
                   simplify = TRUE)
  cenX    <- t(t(sm$X) - w_means)

  # Drop zero-variance columns
  if (sm$by == "NA" || is.null(sm$by.level)) {
    # No by-variable, or a numeric/matrix by-variable: use all rows
    in_level <- rep(TRUE, nrow(sm$X))
  } else {
    # Factor by: sm$X[i, ] == 0 if observation i is not in sm$by.level, so
    # assess the variance only over the rows belonging to this level
    in_level <- rowSums(abs(sm$X) >= 1e-10) > 0
  }
  is0var  <- apply(cenX[in_level, , drop = FALSE], 2,
                   \(x) stats::var(x) < tol,
                   simplify = TRUE)
  ind0var <- if (any(is0var)) which(is0var) else NULL

  # Partition into learners (based on penalty)
  design_mats <- partition_t2_basis_matrix(sm, cenX = cenX, ind0var = ind0var)
  # A one-column penalised block cannot reach sedf = 1 as a ridge learner,
  # so it is fitted as an unpenalised learner
  is_pen    <- sapply(design_mats, \(X) ncol(X) >= 2)
  pen_index <- attr(design_mats, "pen_index")

  bls <- list()
  for (l in seq_along(design_mats)) {
    X       <- design_mats[[l]]
    t2label <- names(design_mats)[[l]]

    if (is_pen[[l]]) {
      # sp of the sm$S penalty this block came from (not its position among
      # the blocks kept as penalised learners)
      sp       <- sm$sp[pen_index[l]] %||% NA
      bls[[l]] <- new_dr_learner(X, sm = sm, weights = weights,
                                 clamp_sp = clamp_sp,
                                 t2label = t2label, sp = sp)
    } else {
      bls[[l]] <- new_nullspace_learner(X, sm = sm, weights = weights,
                                        t2label = t2label)
    }
  }

  return(list(learners = bls, cen_const = w_means, ind0var = ind0var))
}




# Partition design matrices -----------------------------------------------


#' Partition Tensor-Product Basis Matrix by Penalty
#'
#' @description
#' Decomposes a tensor-product smooth's basis matrix into independent
#' ridge-penalised and unpenalised components via diagonal decomposition.
#'
#'
#' @param sm An `mgcv` smooth object of class `"t2.smooth"`, as returned by
#'   [mgcv::smoothCon()].
#' @param cenX Centred basis matrix.
#' @param ind0var Indices of zero-variance columns to drop.
#'
#' @return A list of basis matrices (penalised blocks, then single
#'   unpenalised columns) with
#'   names indicating the penalty/margin dimension (penalised blocks keep the
#'   names of `sm$S`, the unpenalised columns are named `"n1"`, `"n2"`, ...).
#'   The attribute `"pen_index"` gives, for each matrix, the index of the
#'   `sm$S` penalty it came from (`NA` for unpenalised columns).
#'
#' @details
#'   This function enables scale-invariant tensor-product interactions
#'   following Wood et al. (2013).
#'   Each component can be fitted independently during boosting,
#'   allowing adaptive smoothness relaxation along each dimension.
#'
#' @noRd
partition_t2_basis_matrix <- function(sm, cenX, ind0var) {

  checkmate::assert_class(sm, "t2.smooth")

  # Identify penalised and nullspace column groups using the ORIGINAL penalty
  # structure (before dropping zero-variance columns), so that nullspace
  # detection is not affected by dropping.
  # t2 penalties are 0/1 diagonal blocks, so diag(S) != 0 marks the
  # penalised columns
  ind_sep_penalised <- lapply(sm$S, \(S) which(abs(diag(S)) > 1e-14))
  pen_index         <- seq_along(ind_sep_penalised)
  ind_penalised     <- do.call(c, ind_sep_penalised)
  col_is_pen        <- seq_len(ncol(sm$X)) %in% ind_penalised

  # Nullspace (unpenalised) columns, dropping zero-variance ones
  if (any(!col_is_pen)) {
    ind_unpen <- which(!col_is_pen)
    ind_sep_unpenalised <- as.list(ind_unpen[!(ind_unpen %in% ind0var)])
  } else {
    ind_sep_unpenalised <- list()
  }

  # Drop zero-variance columns from penalised groups and remove groups that
  # become empty (can happen with functional predictors where some basis
  # columns are constant after the linear-functional integration).
  if (length(ind0var) > 0L) {
    ind_sep_penalised <- lapply(ind_sep_penalised,
                                \(ind) ind[!(ind %in% ind0var)])
    keep              <- lengths(ind_sep_penalised) > 0L
    ind_sep_penalised <- ind_sep_penalised[keep]
    pen_index         <- pen_index[keep]
  }

  # Distinct labels for the unpenalised (nullspace) columns
  names(ind_sep_unpenalised) <- sprintf("n%d", seq_along(ind_sep_unpenalised))
  ind_sep <- c(ind_sep_penalised, ind_sep_unpenalised)
  Xs      <- lapply(ind_sep, \(ind) cenX[, ind, drop = FALSE])
  attr(Xs, "pen_index") <- c(
    pen_index, rep(NA_integer_, length(ind_sep_unpenalised))
  )

  return(Xs)
}


# Initiate learners from design matrices/formula --------------------------

#' Unpenalised learner for a null-space column x of smooth sm
#'
#' Not merged with new_unpenalised_learner(): its fields come in another
#' order (term, label, by first) and xtWx is computed as a cross-product, so
#' a merged constructor would change the learner objects and their xtWx.
#' @noRd
new_nullspace_learner <- function(x, sm, weights, t2label = NULL) {
  if (is.matrix(x) && !(length(x) %in% dim(x))) {
    stop("new_nullspace_learner does not support matrices with more than ",
         "one vector.")
  }
  x       <- as.vector(x)
  n_total <- length(x)
  if (is.null(t2label) || !nzchar(t2label)) t2label <- "nullspace"
  learner <- new_learner(var_label = paste0(sm$label, t2label),
                         eta = rep(0, n_total))
  learner$term  <- sm$term
  learner$label <- sm$label
  learner$by    <- sm$by
  Wx            <- weights * x

  learner$xtWx   <- as.vector(crossprod(x, Wx))
  learner$design <- matrix(x, ncol = 1)
  class(learner) <- c("unpenalised", "learner")
  return(learner)
}

#' Demmler-Reinsch map and penalty of a ridge-penalised block
#'
#' With the pivoted QR decomposition \eqn{\mathbf{W}^{1/2}\mathbf{Z} =
#' \mathbf{Q}\mathbf{R}} and the spectral decomposition
#' \eqn{\mathbf{R}^{-\top}\mathbf{R}^{-1} = \mathbf{U}\mathbf{D}\mathbf{U}^\top}
#' (the penalty is the identity), the Demmler-Reinsch basis is
#' \eqn{\mathbf{V} = \mathbf{Z}\mathbf{R}^{-1}\mathbf{U}}, with
#' \eqn{\mathbf{V}^\top\mathbf{W}\mathbf{V} = \mathbf{I}} and penalty
#' \eqn{\mathbf{D}}. The eigenvalues are \eqn{d_j = 1/\sigma_j^2}, where
#' \eqn{\sigma_j} are the singular values of \eqn{\mathbf{W}^{1/2}\mathbf{Z}}.
#' When the QR decomposition is rank deficient or
#' \eqn{\sigma_{\min}/\sigma_{\max} \le} `tol`, the basis comes from a
#' truncated singular value decomposition instead, which keeps the directions
#' with \eqn{\sigma_j > } `tol` \eqn{\sigma_{\max}} (e.g. it drops the aliased
#' direction of centred random-effect dummies). Both give the same basis when
#' every \eqn{\sigma_j} passes this threshold.
#'
#' @param Z Numeric matrix. The centred penalised block.
#' @param weights Numeric vector. Prior weights, one per row of `Z`.
#' @param tol Numeric scalar. Relative singular value threshold.
#' @return A list with `A`, the map from `Z` to the basis (`V = Z %*% A`),
#'   and `D`, the decreasing penalty eigenvalues.
#' @noRd
demmler_reinsch_basis <- function(Z, weights, tol = 1e-7) {
  WZ  <- sqrt(weights) * Z
  p   <- ncol(Z)
  qrz <- qr(WZ, tol = tol)
  if (qrz$rank == p) {
    Ri <- backsolve(qr.R(qrz), diag(p))
    e  <- eigen(crossprod(Ri), symmetric = TRUE)
    if (e$values[p] / e$values[1L] > tol^2) {
      A <- matrix(0, p, p)
      A[qrz$pivot, ] <- Ri %*% e$vectors
      return(list(A = A, D = e$values))
    }
  }
  sv  <- svd(WZ)
  r   <- sum(sv$d > tol * sv$d[1L])
  ord <- rev(seq_len(r))
  list(A = sv$v[, ord, drop = FALSE] %*% diag(1 / sv$d[ord], r),
       D = 1 / sv$d[ord]^2)
}

#' Demmler-Reinsch learner for the ridge-penalised block Z of smooth sm
#'
#' sm is an element of smoothCon(); t2label names the block of the t2() term.
#' clamp_sp = TRUE: a fixed sp giving sedf < 1 under these weights is clamped
#' to sedf = 1 instead of rejected (resampling fits of a valid main fit, whose
#' smaller total weight lowers the rho at sedf = 1).
#' @noRd
new_dr_learner <- function(Z, sm, weights, sp, t2label = NULL,
                           clamp_sp = FALSE) {
  if ((!is.na(sp)) && sp == 0) {
    stop("Smooth learners without smoothing penalty are not allowed.")
  }
  n_total <- nrow(Z)
  learner <- new_learner(var_label = paste0(sm$label, t2label %||% ""),
                         eta = rep(0, n_total))
  learner$term  <- sm$term
  learner$label <- sm$label
  learner$by    <- sm$by

  # Demmler-Reinsch basis for the ridge penalty (S = I on t2 blocks)
  dr <- demmler_reinsch_basis(Z, weights)
  r  <- length(dr$D)
  if (r < 2L) {
    stop("The penalised block \"", learner$var_label, "\" has rank ", r,
         " after centring; at least 2 identifiable directions are needed ",
         "for a penalised learner with sedf = 1.", call. = FALSE)
  }
  # map Z -> V (the Demmler-Reinsch map), not the hat matrix
  learner$A      <- A <- dr$A
  learner$D      <- D <- dr$D
  learner$design <- Z %*% A
  learner$beta   <- rep(0, r)
  # rho at sedf = 1: the largest rho allowed
  rho_sedf1 <- init_rho(D, target_edf = 1, label = learner$var_label)

  learner$rho <- if (is.na(sp)) {
    # No sp: sedf = 1
    rho_sedf1
  } else if (sp > 0) {
    # sp is lambda
    rho_proposed <- log(sp)
    if (clamp_sp) rho_proposed <- min(rho_proposed, rho_sedf1)
    if (rho_proposed > rho_sedf1) {
      stop("Smoothing parameter sp=", sp, " results in learner sedf < 1, ",
           "which isn't supported by the method.")
    }
    rho_proposed
  } else if (sp < 0) {
    # -sp is the target sedf
    target_edf <- -sp
    if (target_edf <= 1) {
      stop("sedf proposals must be larger than one (sp < -1) for unbiased ",
           "learner selection to work.")
    }
    if (target_edf >= r) {
      stop("sp = ", sp, " requests sedf >= ", r, ", the rank of the ",
           "centred penalised basis of ", learner$var_label, "; use -", r,
           " < sp < -1.", call. = FALSE)
    }
    init_rho(D, target_edf = -sp, label = learner$var_label)
  }
  learner$lambda_at_zero <- FALSE        # TRUE once rho is rho_lambda_zero
  learner$rho0           <- learner$rho  # Restored on bootstrap reset
  class(learner) <- c("demmler_reinsch", class(learner))

  learner$nu_k  <- min(1, find_nu_k(learner, target_df = 1))
  learner$nu_k0 <- learner$nu_k  # Restored with rho0 on bootstrap reset

  return(learner)
}


#' Unpenalised learner for a parametric column x
#' @noRd
new_unpenalised_learner <- function(x, var_label, weights) {
  learner <- new_learner(var_label = var_label, eta = rep(0, length(x)))
  learner$design <- as.matrix(x)
  x2             <- as.vector(learner$design)^2
  learner$xtWx   <- sum(x2 * weights)
  return(structure(learner, class = c("unpenalised", "learner")))
}


#' General base learner constructor
#' @noRd
new_learner <- function(var_label, beta = 0, edf = 0, eta = NULL) {
  res <- list(
    var_label = var_label,
    beta = beta,
    edf = edf,
    eta = eta
  ) |>
    structure(class = "learner")
  return(res)
}


#' Recover `mgcv`-scale coefficients from fitted learners
#'
#' @description
#' Turns the coefficients of the base learners into a coefficient vector on
#' the basis of `mgcv`, so that `mgcv`'s prediction, plotting and
#' diagnostic methods can use it. [gbgam()] computes the model's
#' `$coefficients` this way; all the inputs are fields of a fitted model.
#'
#' @param learners List. The base learners, `mod$learners` (or
#'   `new_learners()$learners`), with their coefficients in field `beta`.
#' @param coef_names Character vector. The coefficient names,
#'   `mod$coef_names`.
#' @param ind0var Integer vector. The positions of the dropped columns
#'   without variance, `mod$ind0var`.
#' @param cen_const Numeric vector. The centring constants (weighted column
#'   means), `mod$cen_const`, in the order of `coef_names`; the entry of
#'   `"(Intercept)"` is zero.
#' @param cen_absorb Logical vector, or NULL. Which columns have their
#'   centring shift absorbed into the intercept, `mod$cen_absorb`, in the
#'   order of `coef_names`. Default `NULL`: every column, which matches the
#'   uncentred basis of [mgcv::smoothCon()].
#'
#' @return A named numeric vector of coefficients on the `mgcv` prediction
#'   basis (the [mgcv::smoothCon()] basis, centred for `"gbgam_centred"`
#'   smooths), named by `coef_names`: the parametric coefficients first,
#'   then the smooth coefficients term by term.
#'
#' @details
#' The function undoes the three changes made when the learners were built
#' (see [new_learners()]):
#' 1. The coefficients of each penalised learner are multiplied by its
#'    Demmler-Reinsch map (field `A`), which gives the coefficients of the
#'    centred columns of its group. When those columns are rank deficient,
#'    this gives the coefficient vector of smallest norm. The coefficient of
#'    an unpenalised learner is used as it is.
#' 1. Zeros are put back at the positions of the dropped columns without
#'    variance (`ind0var`), so that the result matches the full
#'    [mgcv::smoothCon()] basis.
#' 1. Each centred column is its original column minus its centring
#'    constant, so the fit from the centred columns equals the fit from the
#'    original columns minus the sum of the constants times the
#'    coefficients. For the columns flagged in `cen_absorb` (parametric
#'    columns and smooths with a by-variable or matrix arguments), this sum
#'    is subtracted from the intercept, and they predict from the original
#'    basis. Smooths of class `"gbgam_centred"` predict from the centred
#'    basis directly (see [Predict.matrix.gbgam_centred()]), so their shift
#'    stays out of the intercept and their term-wise predictions are
#'    centred.
#'
#' ## Fields of the learners
#'
#' The function reads only the class of each learner, `beta`, and for
#' penalised learners `A` (or, for learners without `A`, the factors `R` and
#' `U` of the same map). It writes no field. `rho`, `nu_k`, `edf` and `eta` do
#' not enter, so the result depends on the coefficients alone. See
#' [new_learners()] for what each field holds.
#'
#' The coefficients are joined in the order of `learners`, which has to be
#' the order of `mod$learners`, and matched to `coef_names` with `ind0var`
#' left out. So a learner's `beta` can be changed, as in the example below,
#' but learners cannot be removed or reordered.
#'
#' @seealso [new_learners()] for the fields of a learner.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = runif(100))
#' dat$y <- sin(2 * pi * dat$x) + rnorm(100, sd = 0.3)
#' mod <- gbgam(y ~ t2(x, k = 8, bs = "cr"), data = dat, verbose = FALSE)
#' b <- learners_to_mgcv_coef(mod$learners, mod$coef_names, mod$ind0var,
#'                            mod$cen_const, mod$cen_absorb)
#' all.equal(b, coef(mod))
#'
#' # The coefficients give the fit on the mgcv prediction basis
#' X <- predict(mod, type = "lpmatrix")
#' all.equal(as.vector(X %*% b), as.vector(predict(mod)))
#'
#' # The contribution of one learner: set its coefficients to zero
#' lrs <- mod$learners
#' k   <- which(vapply(lrs, inherits, TRUE, "demmler_reinsch"))[1]
#' lrs[[k]]$beta[] <- 0
#' b0  <- learners_to_mgcv_coef(lrs, mod$coef_names, mod$ind0var,
#'                              mod$cen_const, mod$cen_absorb)
#' all.equal(drop(X %*% (b - b0)), mod$learners[[k]]$eta,
#'           check.attributes = FALSE)
#' @noRd
learners_to_mgcv_coef <- function(learners, coef_names, ind0var, cen_const,
                                  cen_absorb = NULL) {

  # Reverse Demmler Reinsch reparametrisation per learner, then concatenate
  undone_DR_betas <- unlist(lapply(learners, function(lr) {
    if (inherits(lr, "demmler_reinsch")) {
      if (!is.null(lr$A)) {
        as.vector(lr$A %*% lr$beta)
      } else {
        # Learners without the map A keep its factors R and U
        as.vector(backsolve(lr$R, lr$U %*% lr$beta))
      }
    } else {
      as.vector(lr$beta)
    }
  }))

  # Insert zeros for the dropped (zero-variance) columns
  if(length(ind0var)>0 ){
    checkmate::assert_integerish(ind0var)
    expected_len <- length(coef_names) - length(ind0var)
    if (length(undone_DR_betas) != expected_len) {
      stop(sprintf(
        paste("learners_to_mgcv_coef: %d learner betas but %d non-zero-var",
              "coefficient slots. "),
        length(undone_DR_betas), expected_len
      ), "This is a bug in partition_t2_basis_matrix -- please report it.",
      call. = FALSE)
    }
    coefficients <- rep(0, length(coef_names))
    coefficients[-ind0var] <- undone_DR_betas
  }else{
    coefficients <- undone_DR_betas
  }
  names(coefficients) <- coef_names

  # Shift intercept
  # E.g. intercept + (x - xbar) beta = (intercept - xbar beta) + x beta
  # Columns of "gbgam_centred" smooths keep their centring in the basis.
  cen_absorb     <- cen_absorb %||% rep(TRUE, length(coefficients))
  stopifnot(length(cen_absorb) == length(coefficients))
  centring_shift <- sum((cen_const * coefficients)[cen_absorb])
  coefficients["(Intercept)"] <- coefficients["(Intercept)"]-centring_shift
  return(coefficients)
}
