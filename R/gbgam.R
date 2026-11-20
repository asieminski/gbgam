#' @importFrom stats gaussian binomial poisson Gamma inverse.gaussian
NULL

#' @importFrom stats setNames
NULL

#' @importFrom utils flush.console
NULL

#' @importFrom stats nobs weights coef fitted qnorm sd uniroot vcov
NULL


#' Fast and simple gradient-boosted generalised additive models
#'
#' @description
#' Fits an exponential-family generalised additive model (GAM) by
#' component-wise gradient boosting, which selects the terms and their
#' smoothness from the data. It speeds up convergence by adaptive smoothness
#' relaxation (ASR), chooses the stopping iteration automatically (by AIC by
#' default) and, on request, attaches a covariance matrix of the
#' coefficients. Smooth terms are built with [mgcv::t2()], with any of the
#' marginal bases it accepts; `s()`, `te()` and `ti()` terms result in an
#' error, so an `mgcv` formula is used after rewriting them as `t2()` terms.
#' The fitted model inherits from class `"gam"` and works with `mgcv` methods
#' such as
#' [mgcv::predict.gam()] and [mgcv::plot.gam()]; [summary()][summary.gbgam()]
#' and [vcov()][vcov.gbgam()] have their own methods.
#'
#' @param formula Formula. A GAM formula as in [mgcv::formula.gam], with
#'   linear terms and smooth terms built with [mgcv::t2()] (see
#'   [new_learners()] for the terms that are not supported). An intercept is
#'   required. The `sp` argument of `t2()` sets the starting smoothness of
#'   the term's learners (see [new_learners()]).
#' @param family Family. A [stats::family] object, a family function or its
#'   name: [stats::gaussian()], [stats::binomial()], [stats::poisson()],
#'   [stats::Gamma()] or [stats::inverse.gaussian()]. Quasi-likelihood
#'   families, such as [stats::quasipoisson()], are not supported. Default
#'   `gaussian()`.
#' @param data Data frame or list. The variables of `formula`. A list may
#'   hold matrices, for example for smooths with matrix arguments. All
#'   elements have the same number of observations: rows for matrices,
#'   length for vectors. Default `list()`: the variables are taken from the
#'   environment of `formula`.
#' @param weights Numeric vector, or NULL. Non-negative prior weights, as a
#'   vector or as an unquoted expression that may refer to columns of `data`
#'   (e.g. `weights = w` or `weights = w * k`). As in [stats::glm()] and
#'   [mgcv::gam()], a column of `data` takes precedence over a variable of
#'   the same name in the calling environment. The weights are case
#'   (frequency) weights: an integer weight \eqn{w_i} counts as \eqn{w_i}
#'   copies of observation \eqn{i}. Their sum acts as the number of
#'   observations in the dispersion estimate and the log-likelihood, so
#'   rescaling them changes AIC stopping; rescale weights that are not
#'   frequency weights (e.g. normalised survey weights) to sum to the number
#'   of observations. Unlike in `mgcv`, the weights also apply to the
#'   sum-to-zero constraints. Default `NULL`: all weights are one.
#' @param subset Vector, or NULL. The observations to use, as in
#'   [stats::model.frame()]. Default `NULL`: all observations.
#' @param na.action Function. What to do with `NA` values, passed to
#'   [stats::model.frame()]. Default `getOption("na.action")`, usually
#'   [stats::na.omit()].
#' @param offset Numeric vector, or NULL. Known offsets added to the linear
#'   predictor; like `weights`, it may refer to columns of `data`. Default
#'   `NULL`: no offset.
#' @param drop.unused.levels Logical. Drop unused factor levels, passed to
#'   [stats::model.frame()]. Default `FALSE`.
#' @param learner_control List. Settings from [gb_learner_ctrl()]: how the
#'   learners are built, selected and updated, including the global step
#'   size \eqn{\nu} (`step_size`; a smaller value makes divergence less
#'   likely) and ASR (`tau`). Default `list()`: the defaults of
#'   [gb_learner_ctrl()].
#' @param stop_control List. Settings from [gb_stop_ctrl()]: the stopping
#'   method (`"aic"`; `"cv"` for a resampling estimate set up with
#'   [gb_resampling()]; `"true_f"`; or `"none"`), the number of iterations
#'   without improvement before stopping (the patience, `iter_no_improv`
#'   and `iter_no_improv_rel`), the limits
#'   `max_iter` and `max_runtime_mins`, and whether the fit can be continued
#'   with [extend_gbgam()] (`resumable`). Default `list()`: the defaults of
#'   [gb_stop_ctrl()].
#' @param edf_control List. Settings from [gb_edf_ctrl()]: how the model
#'   edf \eqn{\widetilde{\operatorname{edf}}}{edf~} is computed: by the
#'   `"hutchpp"` estimator (with `z` probe vectors), `"exact"` for
#'   small data sets, or not at all (`"none"`). Default `list()`: the
#'   defaults of [gb_edf_ctrl()].
#' @param uq_control List. The uncertainty quantification run after
#'   the fit: settings from [gb_uq_ctrl()], or the settings of one
#'   method directly, from [gb_param_boot_ctrl()] (parametric bootstrap) or
#'   [gb_bayesian_ctrl()] (Bayesian backfitting); the method is then taken
#'   from their class. Giving both results in an error: a model holds one
#'   covariance `$Vp`, so run the other method afterwards with
#'   [bootstrap_gbgam()] or [bayesian_backfitting()]. Default `list()`: none.
#' @param save_fitted_values Logical. Keep the fitted means of every
#'   iteration in `$fitted_values`. Default `FALSE`.
#' @param verbose Logical. Show a progress display of the boosting loop and
#'   progress bars of the uncertainty quantification. Default `TRUE`. The
#'   display gives the iteration, the elapsed time, the model
#'   \eqn{\widetilde{\operatorname{edf}}}{edf~} and its largest increment
#'   (when it is tracked) and the current and the best stopping loss (when a
#'   loss is tracked). A console that supports cursor movement also shows a
#'   table of the penalised learners: their selections, their
#'   \eqn{\widetilde{\operatorname{edf}}}{edf~} and their signed share of the
#'   change of the stopping loss (negative when their selections raised the
#'   loss). The RStudio console shows one line, overwritten in place, and
#'   other output (e.g. a log file) gets one line at most every 5 seconds.
#'
#' @returns A fitted model of class `"gbgam"`, inheriting from
#'   `c("gam", "glm", "lm")`. Besides the fields used by
#'   [`mgcv`][mgcv::mgcv-package] methods (`coefficients`, `fitted.values`,
#'   `linear.predictors`, `smooth`, `model`, `prior.weights`, ...), it
#'   contains:
#'   \describe{
#'     \item{`mstop`}{the selected (best) iteration; the model is the one at
#'       this iteration.}
#'     \item{`metrics`}{per-iteration records of the whole run:
#'       `stopping_conditions` (named logical: `no_improvements`,
#'       `mstop_reached`, `time_limit_exceeded`; `mstop_reached` is `TRUE`
#'       when the run hit `max_iter` and is not related to `mstop`),
#'       `k_opts` (index of the learner selected in each iteration), `rhos`
#'       (its `rho` = log(lambda) after the iteration), `post_rhos` (the
#'       same for `learners_to_fit_after_each_iter`) and `ls` (the stopping
#'       loss, starting with iteration 0; the resampling estimate for
#'       `"cv"`).}
#'     \item{`edfs`}{the model edf
#'       \eqn{\widetilde{\operatorname{edf}}^{[m]}}{edf~^[m]} after each
#'       iteration, starting with iteration 0 (only with EDF tracking).}
#'     \item{`edf.total`}{the model edf at `mstop`, the sum of the learners'
#'       edfs (`edf` and `edf1` spread them over the coefficients). `NA`, as
#'       are `edf`, `edf1`, `rank` and `df.residual`, when
#'       `edf_compute_method = "none"`.}
#'     \item{`learners`}{the base learners, with their coefficients `beta`,
#'       accumulated edf `edf` and `rho` = log(lambda) (see
#'       [new_learners()] for all their fields).}
#'     \item{`weights`}{the working weights \eqn{w \mu'(\eta)^2 / V(\mu)}
#'       at the fitted values, as in [stats::glm()] and [mgcv::gam()].}
#'     \item{`user_weights`}{the `weights` supplied by the user (ones by
#'       default); `prior.weights` multiplies them by the number of trials
#'       for a `cbind(successes, failures)` binomial response.}
#'     \item{`scale`, `sig2`}{the estimated dispersion phi (Fletcher, 2012);
#'       `sig2` only for families with a free dispersion.}
#'     \item{`intercept0`}{the MLE of the intercept-only model.}
#'     \item{`converged`}{`TRUE` when loss-based stopping ended the run
#'       (no improvement within the patience), `FALSE` when it
#'       ended at `max_iter` or the time limit.}
#'     \item{`cv`}{for `stopping_method = "cv"`: the resampling settings,
#'       the estimate `risk` for each iteration (starting with iteration 0)
#'       and its components (see the section on the fitted model's `$cv` in
#'       [gb_resampling()]).}
#'     \item{`resume`}{the state of the last iteration, used by
#'       [extend_gbgam()] (`NULL` with `resumable = FALSE`).}
#'     \item{`Vp`}{the covariance matrix of the coefficients attached by
#'       `uq_control`: the (widened) bootstrap covariance, or the
#'       posterior covariance from Bayesian backfitting. See
#'       [bootstrap_gbgam()] and [bayesian_backfitting()] for the further
#'       fields, such as `boot_mstop`, `coefficients_bayes`, `bb_info` and
#'       `n_samples_done` (the replicates or draws completed).}
#'     \item{`uq_control`}{the settings of the uncertainty
#'       quantification that attached `Vp`, as returned by
#'       [gb_uq_ctrl()]; `gb_uq_ctrl("none")` when the model
#'       has no `Vp` (none was requested, or the time limit left too few
#'       replicates or draws).}
#'     \item{`fitted_values`}{the fitted means of every iteration, one
#'       column each (only with `save_fitted_values = TRUE`).}
#'     \item{`elapsed_secs`}{the run time of the boosting fit.}
#'   }
#'
#' @details
#'
#' ## Overview
#'
#' `gbgam()` combines the formulae and smooth constructors of
#' [`mgcv`][mgcv::mgcv-package] with component-wise gradient boosting. It is
#' designed as:
#' * an alternative to `mgcv` when the number of
#'   candidate smooth terms makes standard GAM fitting too slow;
#' * a higher-level interface for users of the `mboost` package, with simpler
#'   model setup, automatic choice of the stopping iteration, support for
#'   most `mgcv` [smooth terms][mgcv::smooth.terms] (including tensor-product
#'   interactions) and fast convergence for complex smooth effects;
#' * an entry point to GAMs for users who know boosting, who can encode
#'   additivity and smoothness as modelling assumptions.
#'
#' For gradient boosting without these features, but with a wider choice of
#' loss functions, see the `mboost` package (Hofner et al., 2014).
#'
#' The steps of a fit are:
#' 1. the smooth terms are built by [mgcv::smoothCon()] from the
#'    [t2()][mgcv::t2] terms of the formula;
#' 1. the columns of each smooth get sum-to-zero constraints, the columns
#'    without variance (within the level, for a factor `by` variable) are
#'    removed, and the others are split into one group per penalty and one
#'    column per unpenalised direction; each group of at least two columns
#'    gets a Demmler-Reinsch reparametrisation (columns orthonormal in the
#'    prior-weighted inner product, and a diagonal penalty) and is one
#'    penalised base learner, and
#'    each other column is one unpenalised base learner;
#' 1. the learners are fitted one at a time, each iteration adding the
#'    learner that reduces the loss most;
#' 1. the model edf is tracked (by default with the NA-Hutch++ trace
#'    estimator), and drives both ASR and the choice of the stopping
#'    iteration;
#' 1. on request, a parametric bootstrap or Bayesian backfitting attaches a
#'    covariance matrix `Vp` of the coefficients.
#'
#' ## The model
#'
#' `gbgam()` fits a generalised additive model of the form
#' \deqn{g(\mathbb{E}[Y_i \mid x_i]) = \eta_i =
#'       a + f_1(x_i) + \cdots + f_J(x_i)}{
#'       g(E[Y_i | x_i]) = eta_i = a + f_1(x_i) + ... + f_J(x_i)}
#' where `g` is the link function, `a` the intercept, and each `f_j` a
#' penalised basis expansion from a [mgcv::t2()] term.
#'
#' ### Base learners
#'
#' The construction of Wood et al. (2013) gives smoothing splines whose
#' ridge penalties act on separate groups of basis columns. The basis columns
#' are centred on their prior-weighted mean (the sum-to-zero constraint), and
#' columns without variance (within the level, for a factor `by` variable),
#' such as the constant column of the null space, are removed. Each group of
#' at least two penalised columns then becomes one penalised base learner,
#' and each other column (an unpenalised column, or the only column of a
#' penalty) one unpenalised base learner; see [new_learners()].
#'
#' For smooths without a `by` variable or matrix arguments, `mgcv`'s
#' prediction and plotting methods
#' use the same centred basis, so [mgcv::plot.gam()] and
#' `predict(type = "terms")` show the centred effects \eqn{f_j} (with
#' prior-weighted mean zero over the training data), and their standard
#' errors leave out the uncertainty of the intercept. Smooths with
#' a `by` variable and linear functional (matrix-argument) smooths are shown
#' uncentred, with their centring absorbed into the intercept. Saved models
#' need `gbgam` installed for `predict()`; `readRDS()` loads it automatically.
#'
#' The Demmler-Reinsch reparametrisation of each penalised learner makes a
#' change of its smoothing parameter cheap, which ASR needs.
#'
#' ### Component-wise boosting
#'
#' The model is fitted by component-wise gradient boosting
#' (Buhlmann & Yu, 2003), which at each iteration `m`:
#'
#' 1. computes the pseudo-residuals, the gradient of the log-likelihood with
#'  respect to the linear predictor;
#' 1. fits every base learner to the pseudo-residuals;
#' 1. selects the learner \eqn{\hat{k}} with the smallest weighted residual
#'  sum of squares of the pseudo-residuals (`selection_loss_type = "rss"`);
#' 1. adds its fit, shrunk by the global step size \eqn{\nu} (`step_size`),
#'  to the linear predictor (see [gb_learner_ctrl()]).
#'
#' ### Unbiased learner selection
#'
#' To make learner selection approximately unbiased when no learner has a
#' true effect, all penalised learners start at the same selection edf,
#' sedf = 1, where sedf is twice the trace of the learner's hat matrix
#' \eqn{\mathbf{A}_k} minus the trace of its square. Then every learner
#' gives the same expected reduction of the (weighted) residual sum of
#' squares when there is no signal (Hofner et al., 2011). The `sp` argument
#' of `t2()` changes this start (see [new_learners()]).
#'
#' ### Adaptive smoothness relaxation (ASR)
#'
#' With fixed smoothing parameters, boosting can converge extremely slowly
#' for wiggly effects. ASR relaxes them during the fit. Each learner `k`
#' accumulates the model edf credited to it,
#' \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} (see below). Whenever
#' `tau` times \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} exceeds the
#' learner's current leverage edf, ledf = \eqn{\operatorname{tr}(
#' \mathbf{A}_k)}{tr(A_k)}, its smoothing parameter \eqn{\lambda_k} is
#' lowered so that they match; \eqn{\lambda_k} is never raised, and it is
#' set to zero when the target exceeds the number of columns of the
#' learner. `tau` is 1 by default (see [gb_learner_ctrl()]). This widens
#' the model space of the learners that are selected often. To keep
#' selection unbiased after each relaxation, the learner's fit is scaled by
#' a learner-specific step size \eqn{\nu_k} in the selection step, chosen
#' so that the sedf of the scaled hat matrix \eqn{\nu_k \mathbf{A}_k}{nu_k
#' A_k} stays at one.
#'
#' ## Effective degrees of freedom
#'
#' Tracking the model edf is central to both ASR and stopping. `gbgam()`
#' uses the definition of Efron (2004): the trace of the derivative of the
#' fitted means with respect to the response,
#' \eqn{\widetilde{\operatorname{edf}}^{[m]} = \operatorname{tr}(d
#' \tilde{\boldsymbol{\mu}}^{[m]} / d \boldsymbol{y})}{edf~^[m] = tr(d mu~^[m]
#' / dy)}, evaluated at \eqn{\boldsymbol{y} =
#' \tilde{\boldsymbol{\mu}}^{[m]}}{y = mu~^[m]} (for non-canonical links,
#' see Where the Jacobian is evaluated). This derivative (the Jacobian) is
#' updated in every iteration by a recursive formula.
#'
#' The definition approximates the expected prediction error for
#' exponential-family likelihoods, but it treats the learner selections and
#' smoothing parameters as fixed in advance rather than chosen from the
#' data. It therefore usually underestimates the edf, because it leaves out
#' the *search degrees of freedom*. For ASR this downward bias is acceptable
#' (it limits overfitting), but it can make AIC stop late. A resampling
#' estimate of the prediction error accounts for the search, at a higher
#' cost: `gb_stop_ctrl("cv", resampling = gb_resampling())` (see Stopping).
#'
#' ### Assigning edf to learners
#'
#' The accumulated edf of learner `k`,
#' \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k}, is the sum of the
#' increases of the model edf over the iterations that selected it. For link
#' functions other than the identity this is an approximation: each update
#' of the linear predictor also changes how the earlier updates depend on
#' the response, and so the edf of the learners selected earlier.
#'
#' ### Where the Jacobian is evaluated
#'
#' For canonical link functions, the derivatives have no term in the
#' residuals, so evaluating them at the observed \eqn{\boldsymbol{y}}{y} or
#' at \eqn{\boldsymbol{y} = \tilde{\boldsymbol{\mu}}^{[m]}}{y = mu~^[m]}
#' gives the same result. For non-canonical link functions the derivatives
#' of the pseudo-residuals and of the starting intercept have terms
#' proportional to the residuals,
#' \eqn{\boldsymbol{y} - \tilde{\boldsymbol{\mu}}^{[m]}}{y - mu~^[m]} and
#' \eqn{\boldsymbol{y} - \tilde{\boldsymbol{\mu}}^{[0]}}{y - mu~^[0]};
#' `gbgam` evaluates them at the observed \eqn{\boldsymbol{y}}{y} rather
#' than at \eqn{\boldsymbol{y} = \tilde{\boldsymbol{\mu}}^{[m]}}{y =
#' mu~^[m]}, where they would vanish.
#' \eqn{\widetilde{\operatorname{edf}}^{[m]}}{edf~^[m]} then no longer
#' corresponds exactly to Efron's (2004) edf and AIC.
#'
#' ### Cost
#'
#' Updating the full \eqn{n \times n} Jacobian costs
#' \eqn{\mathcal{O}(n^2 p_{\hat{k}})}{O(n^2 p_k)} per iteration, where
#' \eqn{p_{\hat{k}}}{p_k} is the number of columns of the selected learner.
#' Instead, `gbgam()` keeps the product of the Jacobian with a fixed random
#' \eqn{n \times z} matrix (a sketch) and applies the non-adaptive Hutch++
#' estimator (Meyer et al., 2021) to it, at
#' \eqn{\mathcal{O}(n z p_{\hat{k}} + n z^2)}{O(n z p_k + n z^2)} cost per
#' iteration (see [gb_edf_ctrl()]).
#'
#' ## Stopping
#'
#' By default, the model edf gives the AIC, \eqn{-2 \log L + c \,
#' \widetilde{\operatorname{edf}}^{[m]}}{-2 log L + c edf~^[m]} with
#' \eqn{c} = `AIC_k` (2 by default), which chooses the stopping iteration
#' before the fit starts to overfit. The log-likelihood includes the
#' dispersion phi, estimated again in every iteration by Fletcher's (2012)
#' bias-corrected estimator. `AIC_k = log(nrow(data))` gives the BIC
#' instead. `iter_no_improv` and `iter_no_improv_rel` in [gb_stop_ctrl()]
#' set how many iterations without improvement of the AIC are run before
#' stopping: by default `max(50, ceiling(0.1 * m*))`, with m* the best
#' iteration so far.
#'
#' `gb_stop_ctrl("cv", resampling = gb_resampling(...))` stops by a
#' resampling estimate of the prediction error instead (k-fold
#' cross-validation, bootstrap, .632+ bootstrap or subsampling): the
#' resampling fits run alongside the main fit, so no number of iterations
#' has to be chosen in advance, and EDF tracking is then needed only for
#' ASR. The estimate is stored in `$cv`.
#'
#' `max_iter` limits the number of iterations and `max_runtime_mins` the
#' run time, e.g. for quick experiments. For simulation studies, the true
#' linear predictor `true_f` can be supplied, so that the fit stops at the
#' smallest mean squared error (the oracle stopping rule) and the noise of
#' choosing the stopping iteration is removed.
#'
#' `mstop` is the best iteration of the run, which typically goes on
#' longer: until the patience runs out without improvement, `max_iter` or
#' the time limit. By default (`gb_stop_ctrl(resumable = TRUE)`) the
#' fitted model keeps the state of its last iteration in `$resume`, so that
#' a run that ended too early can be continued with [extend_gbgam()]
#' instead of being refitted. The state includes the EDF-tracking matrices
#' (\eqn{n \times z} for `"hutchpp"`, \eqn{n \times n} for `"exact"`), which
#' can make the fitted object several times larger; use `resumable = FALSE`
#' when many fits are kept or saved and none will be extended.
#'
#' ## Tensor-product interactions
#'
#' `gbgam()` fits smooth interactions with [t2()][mgcv::t2]. Unlike
#' isotropic smooths, [t2()][mgcv::t2] is scale-invariant: it does not
#' assume that the covariates are measured on the same scale (Wood et al.,
#' 2013).
#'
#' Tensor products have several penalties, and one edf target cannot
#' identify all of their smoothing parameters. In the decomposition of Wood
#' et al. (2013), each penalty of the tensor product acts on its own group
#' of columns, so each penalty gives a separate base learner. The learners
#' of one interaction share the basis evaluation but are selected and
#' relaxed independently, so the fit finds the smoothness of the
#' interaction in each direction without assumptions about how the
#' penalties are scaled relative to each other.
#'
#' ## Uncertainty quantification
#'
#' After the fit, `gbgam()` can attach a covariance matrix `Vp` of the
#' coefficients, chosen with `uq_control` (by default none). It is
#' used by [stats::vcov()], [predict(se.fit = TRUE)][mgcv::predict.gam()],
#' [mgcv::plot.gam()] and [smooth_intervals()].
#'
#' * `gb_uq_ctrl("bootstrap")` (or `gb_param_boot_ctrl()`) runs a
#'   parametric bootstrap ([bootstrap_gbgam()]): responses are simulated
#'   from the fitted model and refitted. By default (mode `"fixed"`) each
#'   replicate replays the main fit's selection and smoothing parameters;
#'   other modes let replicates redo these steps (see
#'   [gb_param_boot_ctrl()]). The covariance of the replicates'
#'   coefficients is widened for smoothing bias, by default by one factor
#'   for the whole model (Wood, 2020).
#' * `gb_uq_ctrl("bayes")` (or `gb_bayesian_ctrl()`) runs Bayesian
#'   backfitting ([bayesian_backfitting()]; Hastie & Tibshirani, 2000),
#'   which samples the posterior of a penalised model whose smoothing
#'   parameters are chosen so that its posterior mode is close to the
#'   boosted smooths. As for the bootstrap, the intervals have their
#'   coverage across the function, not at each point (see
#'   [smooth_intervals()]).
#'
#' For long fits, a safer workflow is to fit without uncertainty
#' quantification and call
#' [bootstrap_gbgam()] or [bayesian_backfitting()] afterwards.
#'
#' ## Errors
#'
#' A step size that is too large makes boosting diverge. `gbgam()` then
#' stops with an error of class `"gbgam_diverged"` as soon as the training
#' deviance increases (or the negative gradient becomes non-finite); its
#' field `iter` gives that iteration. Refit with a smaller `step_size`.
#'
#' Settings of `uq_control` that the fit cannot support (e.g. a
#' bootstrap mode that needs EDF tracking, with `edf_compute_method =
#' "none"`) result in an error before the fit starts. If the uncertainty
#' quantification fails after the fit, the error has class
#' `"gbgam_uq_error"` and keeps the fitted model in its field `fit`.
#'
#' @seealso [gb_learner_ctrl()], [gb_stop_ctrl()], [gb_resampling()],
#'   [gb_edf_ctrl()], [gb_uq_ctrl()], [gb_param_boot_ctrl()] and
#'   [gb_bayesian_ctrl()] for the settings; [extend_gbgam()] to continue a
#'   fit; [bootstrap_gbgam()], [bayesian_backfitting()] and
#'   [smooth_intervals()] for uncertainty quantification; [new_learners()]
#'   for the base learners.
#'
#' @references
#' Buhlmann, P., & Yu, B. (2003).
#' Boosting with the L2 loss: Regression and classification.
#' *Journal of the American Statistical Association, 98*(462), 324--339.
#' \href{https://doi.org/10.1198/016214503000125}{doi:10.1198/016214503000125}
#'
#' Efron, B. (2004).
#' The estimation of prediction error: Covariance penalties and
#' cross-validation.
#' *Journal of the American Statistical Association, 99*(467), 619--632.
#' \href{https://doi.org/10.1198/016214504000000692}{doi:10.1198/016214504000000692}
#'
#' Fletcher, D. J. (2012).
#' Estimating overdispersion when fitting a generalized linear model to sparse
#' data.
#' *Biometrika, 99*(1), 230--237.
#' \href{https://doi.org/10.1093/biomet/asr083}{doi:10.1093/biomet/asr083}
#'
#' Hastie, T., & Tibshirani, R. (2000).
#' Bayesian backfitting (with comments and a rejoinder by the authors).
#' *Statistical Science, 15*(3), 196--223.
#' \href{https://doi.org/10.1214/ss/1009212815}{doi:10.1214/ss/1009212815}
#'
#' Hofner, B., Hothorn, T., Kneib, T., & Schmid, M. (2011).
#' A framework for unbiased model selection based on boosting.
#' *Journal of Computational and Graphical Statistics, 20*(4), 956--971.
#' \href{https://doi.org/10.1198/jcgs.2011.09220}{doi:10.1198/jcgs.2011.09220}
#'
#' Hofner, B., Mayr, A., Robinzonov, N., & Schmid, M. (2014).
#' Model-based boosting in R: A hands-on tutorial using the R package mboost.
#' *Computational Statistics, 29*(1--2), 3--35.
#' \href{https://doi.org/10.1007/s00180-012-0382-5}{doi:10.1007/s00180-012-0382-5}
#'
#' Meyer, R. A., Musco, C., Musco, C., & Woodruff, D. P. (2021).
#' Hutch++: Optimal stochastic trace estimation.
#' In *Symposium on Simplicity in Algorithms (SOSA)* (pp. 142--155). Society
#' for Industrial and Applied Mathematics.
#' \href{https://doi.org/10.1137/1.9781611976496.16}{doi:10.1137/1.9781611976496.16}
#'
#' Wood, S. N. (2020).
#' Rejoinder on: Inference and computation with generalized additive models and
#' their extensions.
#' *TEST, 29*(2), 354--358.
#' \href{https://doi.org/10.1007/s11749-020-00716-0}{doi:10.1007/s11749-020-00716-0}
#'
#' Wood, S. N., Scheipl, F., & Faraway, J. J. (2013).
#' Straightforward intermediate rank tensor product smoothing in mixed models.
#' *Statistics and Computing, 23*(3), 341--360.
#' \href{https://doi.org/10.1007/s11222-012-9314-z}{doi:10.1007/s11222-012-9314-z}
#'
#' @examples
#' \donttest{
#' ## --- Minimal Gaussian example (mgcv-style setup) ------------------------
#' ## gbgam uses Hutch++ EDF, AIC-based stopping and ASR by default.
#' ## Uncertainty quantification is added on request with
#' ## uq_control = gb_uq_ctrl("bootstrap") (or "bayes", or
#' ## directly gb_param_boot_ctrl() or gb_bayesian_ctrl()), or with
#' ## bootstrap_gbgam() or bayesian_backfitting() after the fit.
#' set.seed(1)
#' n  <- 200
#' x  <- runif(n)
#' f  <- sin(2 * pi * x)
#' y  <- rnorm(n, f, sd = 0.3)
#' dat <- data.frame(x = x, y = y)
#'
#' ## Drop-in replacement for mgcv::gam in a simple smooth model
#' m <- gbgam(
#'   y ~ t2(x, k = 10, bs = "cr"),
#'   data = dat
#' )
#' summary(m)
#' plot(m, residuals = TRUE)
#'
#' ## --- Using mgcViz for visualisation (if available) ----------------------
#' ## mgcViz requires Vp to be set; attach a quick bootstrap first.
#' m <- bootstrap_gbgam(m, n_samples = 20)
#' if (requireNamespace("mgcViz", quietly = TRUE)) {
#'   library(mgcViz)
#'   b <- getViz(m)
#'   o <- plot( sm(b, 1) )
#'   o + l_fitLine(colour = "red") +
#'     l_rug(mapping = aes(x = x, y = y), alpha = 0.8) +
#'     l_ciLine(mul = 5, colour = "blue", linetype = 2) +
#'     l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
#' }
#' ## ---- Plain component-wise boosting baseline (no ASR, no EDF, no CI) -----
#' ## Without EDF tracking there is no ASR (tau defaults to 0) and no AIC, so
#' ## a fixed number of iterations gives plain component-wise L2 boosting,
#' ## close to mboost::gamboost() with centred base learners and an explicit
#' ## intercept learner, but using gbgam's learner construction.
#' m_fixed <- gbgam(
#'   y ~ t2(x, k = 10, bs = "cr"),
#'   data         = dat,
#'   stop_control = gb_stop_ctrl("none", max_iter = 100),
#'   edf_control  = gb_edf_ctrl("none")
#' )
#' ## See also gb_learner_ctrl(learners_to_fit_after_each_iter = ) for refitting
#' ## the intercept after every iteration instead of letting it compete, and
#' ## gb_learner_ctrl(intercept_padding = TRUE) for intercept-padded learners.
#' print(m_fixed)
#'
#' ## ---- Warm-starting learners via sp to accelerate convergence ----------
#' ## The sp argument of t2() has a special meaning in gbgam: sp = -5 starts
#' ## the learner at sedf 5 instead of 1 (see ?new_learners). This helps for
#' ## smooth terms that are expected to be complex.
#' set.seed(2)
#' f2 <- sin(2 * pi * x) + 0.4 * sin(6 * pi * x)  # more complex: needs more EDF
#' y2 <- rnorm(n, f2, sd = 0.3)
#' dat2 <- data.frame(x = x, y = y2)
#' # The learner begins closer to its final smoothing level, so ASR needs fewer
#' # iterations to relax it fully, and AIC-based stopping fires earlier.
#' m_fast <- gbgam(
#'   y ~ t2(x, k = 15, bs = "cr", sp = -5),
#'   data = dat2
#' )
#' cat("Warm start (sedf 5): mstop =", m_fast$mstop, "\n")
#'
#'
#' ## ---- Gaussian tensor-product: sin(x)*sin(3z) --------------------------
#' ## f(x,z) = sin(x)*sin(3z) is smooth but anisotropic. t2(x, z) penalises
#' ## curvature in each direction independently, giving scale-invariant
#' ## smoothing.
#' set.seed(3)
#' n2 <- 500
#' x2 <- runif(n2, 0, pi)
#' z2 <- runif(n2, 0, pi)
#' f_tp <- sin(x2) * sin(3*z2)
#' y_tp <- rnorm(n2, f_tp, sd = 0.3)
#' dat_tp <- data.frame(x = x2, z = z2, y = y_tp)
#'
#' ## Quick look at the data-generating surface on a grid
#' gx <- seq(0, pi, length.out = 40)
#' gz <- seq(0, pi, length.out = 40)
#' grd <- expand.grid(x = gx, z = gz)
#' f_grd <- matrix(sin(grd$x) * sin(3*grd$z),
#'                 nrow = length(gx),
#'                 ncol = length(gz))
#' contour(gx, gz, f_grd, xlab = "x", ylab = "z",
#'       col = hcl.colors(64, "viridis"), main = "True f(x,z) = sin(x)*sin(3z)")
#' points(x2, z2, pch = ".", cex = 2)   # observed locations
#'
#' ## Each of the three t2 penalties becomes its own learner (and each
#' ## unpenalised column another), so the smoothness in x and in z is
#' ## relaxed separately
#' m_tp <- gbgam(
#'   y ~ t2(x, z, k = c(8, 8)),
#'   data         = dat_tp,
#'   stop_control = gb_stop_ctrl("aic"),
#'   uq_control   = gb_uq_ctrl("bootstrap", n_samples = 20)
#' )
#' print(m_tp)
#' summary(m_tp)
#' plot(m_tp)
#'
#' ## ---- Poisson: mgcv::gamSim test case ---------------------------------
#' ## gamSim(eg = 1, dist = "poisson") provides a standard four-term benchmark
#' ## with f0(x0) + f1(x1) + f2(x2); x3 has no true effect (variable selection).
#' set.seed(6)
#' dat_pois <- gamSim(1, n = 400, dist = "poisson", scale = .1, verbose = FALSE)
#'
#' m_pois <- gbgam(
#'   y ~ t2(x0, bs = "cr", k = 10) +
#'     t2(x1, bs = "cr", k = 10) +
#'     t2(x2, bs = "cr", k = 10) +
#'     t2(x3, bs = "cr", k = 10),
#'   data       = dat_pois,
#'   family     = poisson(),
#'   # the bootstrap's own settings: the method is taken from their class
#'   uq_control = gb_param_boot_ctrl(n_samples = 20)
#' )
#' print(m_pois)
#' ## x3 has no effect: small share of the loss reduction (column %d_stop_loss)
#' summary(m_pois)
#' plot(m_pois, pages = 1)
#'
#' ## ---- Binomial: mgcv::gamSim test case --------------------------------
#' set.seed(7)
#' dat_bin <- gamSim(1, n = 400, dist = "binary", scale = .33, verbose = FALSE)
#'
#' m_bin <- gbgam(
#'   y ~ t2(x0, bs = "cr", k = 15) +
#'     t2(x1, bs = "cr", k = 15) +
#'     t2(x2, bs = "cr", k = 15) +
#'     t2(x3, bs = "cr", k = 15),
#'   data         = dat_bin,
#'   family       = binomial(),
#'   stop_control = gb_stop_ctrl("aic"),
#'   uq_control   = gb_uq_ctrl("bootstrap", n_samples = 10)
#' )
#' print(m_bin)
#' summary(m_bin)
#' plot(m_bin, pages = 1)
#' }
#' @export
gbgam <-
  function(formula,
           family = gaussian(),
           data = list(),
           weights = NULL,
           subset = NULL,
           na.action,
           offset = NULL,
           drop.unused.levels = FALSE,

           learner_control = list(),  # gb_learner_ctrl()
           stop_control    = list(),  # gb_stop_ctrl()
           edf_control     = list(),  # gb_edf_ctrl()
           uq_control      = list(),  # gb_uq_ctrl()

           save_fitted_values = FALSE,
           verbose = TRUE) {
  # Start the clock first, so that the run time covers the whole call
  start_time <- Sys.time()


# Model frame -------------------------------------------------------------

  # Extract gbgam calls
  cl <- match.call() # call needed in gam object for update to work
  cl_no_dots <- match.call(expand.dots=FALSE)

  # Evaluate weights, subset and offset in `data` (enclosed by the caller,
  # like glm/gam) only when the expression refers to a column of `data`,
  # e.g. `weights = w` or `weights = w * k`. Otherwise, or if that
  # evaluation fails (e.g. `weights = d$w` forwarded by a wrapper, where `d`
  # is not visible from the caller), force the promise, so that local
  # variables of wrappers, `...` forwarding and base-named user vectors
  # (e.g. `scale`) resolve in their own environment.
  caller <- parent.frame()
  if (!missing(weights)) {
    weights <- eval_in_data(substitute(weights), function() weights,
                            data, caller)
  }
  if (!missing(subset)) {
    subset <- eval_in_data(substitute(subset), function() subset,
                           data, caller)
  }
  if (!missing(offset)) {
    offset <- eval_in_data(substitute(offset), function() offset,
                           data, caller)
  }

  # Create a model.frame for all data, evaluated in the caller's environment
  gp <- mgcv::interpret.gam(formula) # interpret the formula
  # Pass weights, subset and offset as evaluated values to avoid
  # parent.frame() scoping issues; model.frame() still applies subset and
  # na.action to them.
  mf_call <- make_mf_call(formula=gp$fake.formula,
                          call=cl_no_dots,
                          values=list(weights = weights, subset = subset,
                                      offset = offset),
                          drop.unused.levels=drop.unused.levels)
  mf <- eval(mf_call, caller) # the model frame now contains all the data
  n <- nrow(mf)
  if (n < 2)
    stop("gbgam: fewer than 2 observations without NA.", call. = FALSE)
  # Extract parametric terms (pterms) from parametric model.frame (pmf)
  pmf_call <- mf_call
  pmf_call$formula <- gp$pf # parametric part of gp$fake.formula
  pmf <- eval(pmf_call, caller)


  # Prior weights supplied by the user (ones by default)
  user_weights <- stats::model.weights(mf) %||% rep(1, n)

  # Initialise offset
  offset <- stats::model.offset(mf) %||% rep(0, n)

# Control arguments -------------------------------------------------------

  # edf_control first: its edf_compute_method sets the tau default
  edf_control <- do.call("gb_edf_ctrl", edf_control)

  # tau = NULL resolves to 1 (ASR on) when EDF is computed, 0 otherwise.
  # A tau > 0 with EDF disabled is ignored with a warning.
  no_edf <- edf_control$edf_compute_method == "none"
  if (is.null(learner_control$tau)) {
    learner_control$tau <- if (!no_edf) 1 else 0
  } else if (isTRUE(learner_control$tau > 0) && no_edf) {
    warning("gbgam: tau = ", learner_control$tau, " (adaptive smoothness ",
            "relaxation) has no effect with edf_compute_method = \"none\"; ",
            "tau is set to 0.", call. = FALSE)
    learner_control$tau <- 0
  }

  learner_control   <- do.call("gb_learner_ctrl",   learner_control)
  if (learner_control$step_size > 1) {
    warning("gbgam: step_size > 1 is not recommended.", call. = FALSE)
  }
  stop_control      <- do.call("gb_stop_ctrl",      stop_control)
  uq_control <- as_uq_ctrl(uq_control)

  check_true_f_length(stop_control$true_f, n)

  # Response and family
  Y <- stats::model.response(mf)
  family <- new_family(family = family)

# Checks ------------------------------------------------------------------

  check_formula(formula)
  # Pass the weights as supplied (NULL when absent), not the defaulted
  # `user_weights`, so the fractional-binomial warning can fire.
  check_family_inputs(family, model.resp = Y,
                      weights = stats::model.weights(mf))
  checkmate::assert_flag(save_fitted_values)
  checkmate::assert_flag(verbose)
  checkmate::assert_numeric(weights, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(offset, finite = TRUE)
  check_control_args_compatible(
    learner_control,
    stop_control = stop_control,
    edf_control  = edf_control,
    uq_control   = uq_control)

# Fit ---------------------------------------------------------------------

  # Initialise the model object
  mod <- list()
  mod$family <- family
  # Store the original call and additional control parameters
  mod$call <- cl
  mod$verbose <- verbose

  # append with numeric vectors y and prior.weights
  mod <- c(mod, init_numeric_response_with_weights(Y, family = family,
                                                   weights = user_weights))

  # Extract info from model.frames
  mod$model <- mf
  mod$na.action <- attr(mf, "na.action")
  mod$terms <- attr(mf, "terms")
  mod$pterms <- attr(pmf, "terms")
  mod$offset <- offset
  # User-supplied weights. `$weights` is set to the working weights in
  # finalise_gbgam() (glm/mgcv convention).
  mod$user_weights <- user_weights
  # Save levels and contrasts of parametric factors only, as mgcv does.
  # Factors that appear only inside smooths are handled by the smooths.
  # contrasts is NULL when there is no parametric factor: predict.gam()
  # passes it to model.matrix(), which rejects an unnamed empty list.
  pX            <- stats::model.matrix(mod$pterms, mf)
  mod$contrasts <- attr(pX, "contrasts")
  mod$xlevels   <- stats::.getXlevels(mod$pterms, mf)
  # Map parametric columns to terms; predict.gam(type = "terms") needs it
  mod$assign    <- attr(pX, "assign")


  # Initialise linear predictor eta. intercept0 (the MLE of the
  # intercept-only model) is always kept for the null deviance; the boosting
  # start a0 is the MLE only when initialise_intercept = TRUE, and 0 otherwise.
  mod$intercept0 <- find_intercept_mle(
    y = mod$y,
    family = family,
    prior.weights = mod$prior.weights,
    offset = mod$offset
  )
  a0 <- if (learner_control$initialise_intercept) mod$intercept0 else 0
  mod$linear.predictors <- a0 + offset
  if (!learner_control$initialise_intercept &&
      is.function(family$valideta) &&
      !isTRUE(family$valideta(mod$linear.predictors))) {
    stop("gbgam: initialise_intercept = FALSE starts boosting at eta = ",
         "offset, which is outside the valid range of the \"", family$link,
         "\" link; use initialise_intercept = TRUE.", call. = FALSE)
  }

  # Initialise learners
  # Mgcv smooth terms
  mod$smooth <- construct_smooths(mf,
                                  mgcv_terms = gp$smooth.spec,
                                  knots = learner_control$knots)
  # Learners (natural decomposition, centring)
  compute_edf <- edf_control$edf_compute_method!="none"
  learner_bundle <- new_learners_internal(smooths=mod$smooth,
    data=mf, # Same as data but cleaned up
    pterms=mod$pterms,
    weights = mod$prior.weights,
    ahat = a0,
    edf0 = as.numeric(compute_edf && learner_control$initialise_intercept),
    knots = learner_control$knots)
  mod$nsdf <- learner_bundle$nsdf
  mod$smooth <- learner_bundle$smooths
  mod$learners <- learner_bundle$learners
  mod$coef_names <- learner_bundle$coef_names
  mod$ind0var <- learner_bundle$ind0var
  mod$cen_const <- learner_bundle$cen_const
  mod$cen_absorb <- learner_bundle$cen_absorb

  check_learner_index_bounds(
    learner_control$preselected_learner_indices,
    learner_control$excluded_learner_indices,
    after      = learner_control$learners_to_fit_after_each_iter,
    n_learners = length(mod$learners)
  )
  if (learner_control$intercept_padding) {
    check_intercept_padding(mod$learners)
  }
  # Data-driven selection needs at least one candidate learner
  if (is.null(learner_control$preselected_learner_indices)) {
    check_selection_candidates(
      n_learners        = length(mod$learners),
      excluded          = learner_control$excluded_learner_indices,
      intercept_padding = learner_control$intercept_padding
    )
  }

  # sedf_k = 2 tr(A_k) - tr(A_k^2) that every learner starts at
  mod$target_tr2AminAA <- 1
  # $uq_control describes the $Vp the model carries: none until the
  # uncertainty quantification attaches one (update.gbgam() reruns what
  # it records)
  mod$uq_control <- gb_uq_ctrl("none")
  # Starting model EDF (the intercept, if initialised); init_derivs()
  # replaces it by tr(dmu0/dy) when EDF is tracked
  mod$edf.total <- if(learner_control$initialise_intercept) 1 else 0

  # Combine all objects (must precede init_derivs so that edf_compute_method,
  # z, cweights, and initialise_intercept are available on mod)
  mod <- c(mod, learner_control, edf_control, stop_control)

  # Initialise EDF-tracking derivatives (no-op when edf_compute_method = "none")
  mod <- init_derivs(mod, n)



  # Iterate
  final_model <- if (stop_control$stopping_method == "cv") {
    fit_gbgam_cv(mod, start_time, save_fitted_values, verbose,
                 resumable = stop_control$resumable)
  } else {
    fit_gbgam(mod, start_time, save_fitted_values, verbose,
              resumable = stop_control$resumable)
  }
  final_model <- finalise_gbgam(final_model, formula, gp)

# Uncertainty quantification ----------------------------------------------

  elapsed_mins <- as.numeric(Sys.time()-start_time, units="mins")
  # Clamp at 0: boosting may have used up the whole budget, in which case the
  # samplers warn and return the fitted model unchanged.
  remaining_mins <- if (is.null(stop_control$max_runtime_mins)) NULL else
    max(0, stop_control$max_runtime_mins - elapsed_mins -
             finalise_margin_mins)
  # An error here keeps the fit: it is attached to the condition as `$fit`
  final_model <- tryCatch(
    apply_uq(final_model,
             uq_control,
             verbose = verbose,
             max_runtime_mins = remaining_mins),
    error = function(e) stop_uq_failed(e, final_model,
                                       uq_control$method))

  return(final_model)
}


#' Attach the requested uncertainty quantification to a fitted model
#'
#' @param mod A finalised `gbgam` object.
#' @param uq_control Output of [gb_uq_ctrl()].
#' @param verbose Logical; show progress.
#' @param max_runtime_mins Remaining time budget, or `NULL`.
#' @noRd
apply_uq <- function(mod, uq_control, verbose, max_runtime_mins) {
  ctrl <- uq_control$control
  switch(
    uq_control$method,
    none = mod,
    bootstrap = bootstrap_gbgam(
      mod,
      n_samples        = ctrl$n_samples,
      mode             = ctrl$mode,
      alpha            = ctrl$alpha,
      correction       = ctrl$correction,
      phi              = ctrl$phi,
      max_iter         = ctrl$max_iter,
      verbose          = verbose,
      max_runtime_mins = max_runtime_mins
    ),
    bayes = bayesian_backfitting(
      mod,
      n_samples        = ctrl$n_samples,
      n_burnin         = ctrl$n_burnin,
      linear           = ctrl$linear,
      undersmooth      = ctrl$undersmooth,
      edf_floor        = ctrl$edf_floor,
      phi              = ctrl$phi,
      verbose          = verbose,
      max_runtime_mins = max_runtime_mins
    )
  )
}


#' Stop when the uncertainty quantification fails after the main fit,
#' keeping the fit
#'
#' Raises an error of class `gbgam_uq_error` whose `$fit` field holds
#' the fitted model without uncertainty quantification, so that a long
#' fit is not lost.
#'
#' @param e The condition raised by the uncertainty quantification.
#' @param fit The fitted model.
#' @param method Uncertainty quantification method.
#' @noRd
stop_uq_failed <- function(e, fit, method) {
  msg <- paste0(
    "Uncertainty quantification (method = \"", method,
    "\") failed after the main fit: ",
    conditionMessage(e), "\nThe fitted model is kept in the `fit` field of ",
    "this error (class \"gbgam_uq_error\"); catch it with ",
    "tryCatch() and run bootstrap_gbgam() or bayesian_backfitting() on it.")
  stop(structure(
    class = c("gbgam_uq_error", "error", "condition"),
    list(message = msg, call = NULL, fit = fit, parent = e)
  ))
}

#' Initialise a numeric response vector with corresponding model weights
#'
#' @param model.resp output from [stats::model.response()]
#' @param family a [stats::family] object
#' @param weights output from [stats::model.weights()] or a vector of ones.
#'
#' @returns a named list with a numeric response vectors `y` and `prior.weights`
#' @noRd
init_numeric_response_with_weights <- function(model.resp, family, weights) {
  # Binomial formula with cbind(successes, failures)~.
  if(is.matrix(model.resp)&&ncol(model.resp)==2){
    n_trials <- model.resp[,1]+model.resp[,2]
    # Rows with zero trials carry zero prior weight; set y = 0 for them, as
    # binomial()$initialize does, to avoid NaN
    y <- ifelse(n_trials == 0, 0, model.resp[, 1] / n_trials)
    prior.weights <- weights*n_trials
    # All other formulae
  }else{
    y <- model.resp
    if (family$family == "binomial") {
      if (is.factor(y)) y <- as.numeric(y == levels(y)[2])
      if (is.logical(y)) y <- as.numeric(y)
    }
    prior.weights <- weights
  }
  list(y=y, prior.weights=prior.weights)
}

#' Evaluate a model.frame argument, in `data` when it refers to a column
#'
#' @param expr The unevaluated argument.
#' @param fallback A function forcing the argument's promise.
#' @param data The `data` argument of [gbgam()].
#' @param env The calling environment of [gbgam()].
#' @return The value of the argument.
#' @noRd
eval_in_data <- function(expr, fallback, data, env) {
  if (!any(all.vars(expr) %in% names(data))) return(fallback())
  tryCatch(eval(expr, data, env), error = function(e) fallback())
}

#' Create a call to create a model.frame data object
#'
#' Constructs a call to [stats::model.frame()] forwarding `data` and
#' `na.action` from the parent call, and `weights`, `subset` and `offset` as
#' evaluated values.
#'
#' @param formula fake.formula object from [mgcv::interpret.gam()], which
#'   contains all relevant variables
#' @param call a `match.call(expand.dots=FALSE)` object
#' @param values named list of the evaluated `weights`, `subset` and
#'   `offset` (`NULL` entries are left out); injected directly into the
#'   `model.frame` call to avoid a parent-frame scoping issue with
#'   unevaluated symbols (e.g. an unevaluated `weights` symbol resolving to
#'   `base::weights`, or `..1` forwarded by a wrapper).
#' @param drop.unused.levels logical flag passed to [stats::model.frame()].
#' @noRd
make_mf_call <- function(formula, call, values, drop.unused.levels) {
  args <- c(as.list(call)[c("data", "na.action")], values)
  args <- args[!vapply(args, is.null, logical(1))]

  as.call(c(
    list(quote(stats::model.frame), formula = formula),
    args,
    list(drop.unused.levels = drop.unused.levels)
  ))
}
