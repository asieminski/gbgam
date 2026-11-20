#' Bayesian backfitting for a fitted `gbgam` model
#'
#' @description
#' Attaches a posterior covariance matrix `$Vp` to a fitted [gbgam()] model.
#' The smoothing parameters of a penalised additive model are chosen so
#' that its posterior mode is close to the boosted smooths (Stage 1), and
#' its posterior is then sampled by Bayesian backfitting (Hastie &
#' Tibshirani, 2000), with one sampling block per smooth (Stage 2). Two
#' options change the Stage-1 smoothing parameters before sampling:
#' `linear = "free"` sets those of the linear learners to zero, and
#' `undersmooth = TRUE` lowers those of learners with a low edf. Either
#' moves the posterior mode away from the boosted fit.
#' [plot()][plot.gbgam()] and [predict(se.fit = TRUE)][predict.gbgam()] then
#' use `$Vp`, with no further step; [smooth_intervals()] gives intervals
#' centred on the posterior mean, and the limits as numbers.
#'
#' @inheritParams gb_bayesian_ctrl
#' @param mod A fitted model of class `"gbgam"`, as returned by [gbgam()].
#' @param verbose Logical. Show a progress bar. Default `FALSE`.
#' @param max_runtime_mins Non-negative number, or NULL. Time limit in
#'   minutes for the whole call: Stage 1, the posterior mode and the
#'   sampling sweeps all count against it, and the sampler gets the time
#'   that is left. When the limit is reached, the posterior is computed from
#'   the draws made so far (with a warning; their number is in
#'   `$n_samples_done`), and with fewer than 2 draws after the burn-in `mod`
#'   is returned unchanged (with a warning). Default `NULL` (no limit).
#'
#' @return
#' The original `mod` object with additional fields:
#' \describe{
#'   \item{`$Vp`}{Posterior covariance matrix of `mod$coefficients` (the
#'     scale of `mgcv`'s coefficients), intercept included. Dimension `p x p`
#'     where `p = length(mod$coefficients)`.}
#'   \item{`$coefficients_bayes`}{Posterior mean of the coefficients, on the
#'     same scale. It differs from the boosting estimate `mod$coefficients`
#'     (see Linear learners); [smooth_intervals()] can centre intervals on
#'     either.}
#'   \item{`$beta_samples`}{List with one `p_k x draws` matrix per learner:
#'     the posterior draws of each learner's coefficients on its own scale,
#'     one column per draw. This is the transpose of the orientation of
#'     `$boot_learner_coef` from [bootstrap_gbgam()].}
#'   \item{`$mh_accept_rate`}{Acceptance rate of each learner's draws after
#'     the burn-in: all 1 for the Gibbs sampler (Gaussian family, identity
#'     link). For Metropolis-Hastings, the rate of the block of the learner's
#'     smooth (of the intercept's own block for the intercept), and `NA` for
#'     learners held at zero; the rate of the joint linear block is in
#'     `$bb_info$linear_accept`.}
#'   \item{`$bb_info`}{A list describing the posterior:
#'     \describe{
#'       \item{`stage1`}{how Stage 1 ended: `"converged"` or `"stalled"`
#'         (see Stage 1), or `"no_learners"` when boosting selected no
#'         learner;}
#'       \item{`stage1_sweeps`}{the number of Stage-1 sweeps, not counting
#'         the fit at the starting smoothing parameters;}
#'       \item{`discrepancy`}{the discrepancy left at the end of Stage 1
#'         (see Stage 1);}
#'       \item{`phi`}{the dispersion used: `phi`, or its estimate (see
#'         Dispersion);}
#'       \item{`linear`}{the setting `linear`;}
#'       \item{`linear_accept`}{the acceptance rate of the joint linear
#'         block, `NA` without one;}
#'       \item{`aliased`}{the learners removed as exactly aliased (see Linear
#'         learners);}
#'       \item{`undersmoothed`}{`NULL`, or a data frame of the learners
#'         lifted to `edf_floor`: `term`, `learner`, `rho_from` and `rho_to`
#'         (`rho` = log(lambda) before and after), `edf_from` and `edf_to`;}
#'       \item{`lambda`}{a data frame of every learner's `term`, `learner`,
#'         `status` (`"finite"`, `"zero"` for lambda_k = 0, or `"inf"` for
#'         lambda_k = infinity) and `rho` (`NA` unless finite).}
#'     }}
#'   \item{`$n_samples_done`}{The number of posterior draws kept after the
#'     burn-in: `n_samples`, or fewer when `max_runtime_mins` stopped the
#'     sampler. [bootstrap_gbgam()] records its number of replicates in the
#'     same field.}
#'   \item{`$uq_control`}{`gb_uq_ctrl("bayes")` with the
#'     settings that were run (`n_samples` as requested).}
#' }
#' Uncertainty quantification already attached to `mod` (by
#' [bootstrap_gbgam()] or an earlier call) is replaced. If the time limit
#' leaves fewer than 2 draws, `mod` is returned unchanged, with its earlier
#' results if it had any, and a warning. Either way, `$uq_control`
#' describes the `$Vp` that the returned model carries.
#'
#' @details
#' ## The model
#'
#' The learners of each model term form one smooth: all learners of a
#' `t2()` term, or one parametric covariate. Each learner \eqn{k} has its
#' own smoothing parameter \eqn{\lambda_k} and, on its Demmler-Reinsch
#' basis, the prior
#' \deqn{\boldsymbol{\beta}_k \sim \mathcal{N}\{\mathbf{0},
#'   \phi(\lambda_k \mathbf{D}_k)^{-1}\},}{beta_k ~ N(0, phi (lambda_k
#'   D_k)^(-1)),}
#' where \eqn{\mathbf{D}_k}{D_k} holds the penalty eigenvalues. An
#' unpenalised column is scaled to unit weighted norm and given
#' \eqn{\mathbf{D}_k = 1}{D_k = 1}.
#' A learner with \eqn{\lambda_k = \infty} is removed (held at zero), and
#' one with \eqn{\lambda_k = 0} has a flat prior, as has the intercept \eqn{a}.
#' An intercept that boosting holds at zero (excluded with
#' `excluded_learner_indices`, not in `learners_to_fit_after_each_iter`, and
#' with `initialise_intercept = FALSE`; see [gb_learner_ctrl()]) stays at
#' zero in the posterior.
#'
#' ## Stage 1 (MAP approximation): smoothing parameters for the boosted fit
#'
#' The smoothing parameters are chosen so that the posterior mode is close
#' to the boosted smooths. Closeness is measured by the discrepancy: the sum
#' over smooths of the squared distance, weighted by the prior weights,
#' between the smooth at the posterior mode and the boosted smooth (the
#' intercept is not included). Learners that boosting never selected stay
#' removed. Each of the other \eqn{K} learners that is penalised starts at
#' the \eqn{\lambda_k} whose trace equals its accumulated edf in the boosted
#' fit, \eqn{\operatorname{tr}(\mathbf{F}_k) = \sum_i 1 / (1 + \lambda_k
#' d_{ki}) = \widetilde{\operatorname{edf}}_k}{tr(F_k) = sum_i 1 / (1 +
#' lambda_k d_ki) = edf~_k}, with \eqn{d_{ki}}{d_ki} its penalty eigenvalues;
#' at \eqn{\lambda_k = 0} if \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k}
#' is at least its number of columns, and at its \eqn{\lambda_k} in the
#' boosted fit if \eqn{\widetilde{\operatorname{edf}}_k}{edf~_k} is not
#' positive (always so with `edf_compute_method = "none"`). An unpenalised
#' learner, a single column
#' \eqn{\mathbf{x}}{x} with \eqn{\mathbf{D}_k = 1}{D_k = 1}, starts at the
#' \eqn{\lambda_k} for which one penalised update of it alone, from the
#' boosted fit, returns its boosted coefficient
#' \eqn{\tilde\beta_k}{beta~_k},
#' \eqn{\lambda_k = \mathbf{x}^\top \mathbf{W}_z (\mathbf{z} -
#' \tilde{\boldsymbol\eta}) / \tilde\beta_k}{lambda_k = x' W_z (z - eta~) /
#' beta~_k}, or at \eqn{\lambda_k = 0} if this is not positive; here
#' \eqn{\mathbf{z}}{z} and \eqn{\mathbf{W}_z}{W_z} are the working response
#' and weights at the boosted fit \eqn{\tilde{\boldsymbol\eta}}{eta~}. The
#' penalised model is first fitted at these smoothing parameters, by
#' backfitting from the boosted fit until the fit converges; Stage 1 then
#' runs sweeps over the \eqn{K} learners, one learner at a time.
#'
#' For learner \eqn{k} of smooth \eqn{j}, a sweep first proposes a new
#' \eqn{\lambda_k}. With the other smooths fixed, the coefficients of smooth
#' \eqn{j} after one penalised update are a function of \eqn{\lambda_k};
#' the proposal is the \eqn{\lambda_k} (0 and \eqn{\infty} included) whose
#' update is closest to the boosted smooth \eqn{j}, found on a grid of
#' log(lambda_k) and refined by a one-dimensional search. A proposal that
#' would change smooth \eqn{j} is checked: the penalised
#' model is refitted at it by backfitting until the fit converges, and the
#' proposal is kept if the discrepancy of that fit does not increase.
#' Otherwise the step in log(lambda_k) is halved and the model refitted, up
#' to four times; if no step is kept, \eqn{\lambda_k} stays as it was. So
#' the whole model is fitted to convergence after each update of a single
#' \eqn{\lambda_k}.
#'
#' Stage 1 ends when a sweep changes neither the intercept nor any smooth
#' by more than a small tolerance, and fails if that takes more than 200
#' sweeps. `$bb_info$stage1` is `"converged"` if the last sweep kept every
#' proposal it tried and `"stalled"` if it rejected one or had none for a
#' learner (a singular system skips the learner in that sweep). Both are used
#' in the same way. No tolerance is applied to the discrepancy that
#' is left, which is reported in `$bb_info$discrepancy`.
#'
#' ## Linear learners
#'
#' With `linear = "free"`, every unpenalised learner gets
#' \eqn{\lambda_k = 0}, whether boosting selected it, shrank it or never
#' selected it. Only learners excluded with `excluded_learner_indices` (and
#' never fitted), and learners whose column is zero on every row with
#' positive weight, stay out; so does an intercept held at zero (see The
#' model), which is then not in the joint linear block. The linear learners
#' are all unpenalised learners: the columns of parametric terms and the
#' null-space columns of every smooth (for a `t2(x, z)` term, its learners
#' `n1`, `n2` and `n3`).
#' After the smooths, the intercept and all linear learners are drawn as
#' one joint linear block in each sweep (and fitted jointly when the
#' posterior mode is computed), so that the sampler moves freely along
#' collinear linear parts. This joint update of all unpenalised effects
#' follows the modified backfitting of Buja et al. (1989), but the smooths
#' are not restricted to the complement of the linear space: a null-space
#' column is also drawn with the other learners of its smooth, so it is
#' updated twice per sweep. Under concurvity, boosting splits a linear
#' trend shared by correlated covariates along its path; freeing the linear
#' learners lets the posterior keep the trade-off between those linear
#' parts, as `mgcv` does. With `linear = "kept"`, each linear learner keeps
#' its Stage-1 smoothing parameter and is drawn only with its own smooth
#' (no linear block).
#'
#' Columns with a flat prior (the intercept and the freed linear learners)
#' can be exactly aliased, for example the level constants of a
#' `by = factor` smooth with the intercept; the posterior would then be
#' improper. As `mgcv` does for identifiability, such columns are found by a
#' pivoted QR decomposition and their learners removed; they are listed in
#' `$bb_info$aliased`. Removing them does not change the space of fitted
#' functions.
#'
#' Freeing the linear learners also gives up the shrinkage that boosting
#' applied to them, so the posterior mode and mean move away from the
#' boosted fit, and the intervals centred on the boosted fit (`predict()`)
#' and on the posterior mean ([smooth_intervals()] with
#' `centre = "posterior"`) differ.
#'
#' ## Dispersion
#'
#' The dispersion \eqn{\phi} is 1 for Poisson and binomial responses.
#' Otherwise it is `phi` when given (e.g. a known dispersion), else the
#' restricted maximum likelihood (REML) estimate at the fixed smoothing
#' parameters (Wood, 2011, eq. 4, approach (a)), at the posterior mode for
#' the final smoothing parameters (after `linear` and any undersmoothing).
#' It maximises the Laplace approximate restricted likelihood in
#' \eqn{\phi}, the root of
#' \deqn{\frac{D_p}{2\phi^2} + l_s'(\phi) + \frac{M_p}{2\phi} = 0,
#'   \qquad D_p = \sum_i d_i(y_i, \mu_i) + \sum_k \lambda_k
#'   \boldsymbol\beta_k^\top \mathbf{D}_k \boldsymbol\beta_k,}{
#'   D_p / (2 phi^2) + l_s'(phi) + M_p / (2 phi) = 0,
#'   D_p = sum_i d_i(y_i, mu_i) + sum_k lambda_k beta_k' D_k beta_k,}
#' with \eqn{d_i} the deviance contribution of observation \eqn{i} (prior
#' weight included), the penalty summed over the learners with a finite
#' \eqn{\lambda_k} (unpenalised columns scaled to unit weighted norm),
#' \eqn{l_s} the saturated log-likelihood ([mgcv::fix.family.ls()]) and
#' \eqn{M_p} the dimension of the null space of the penalty: the number of
#' coefficients with a flat prior, the intercept (unless held at zero) and
#' the columns of every learner with \eqn{\lambda_k = 0}, such as the freed
#' linear learners (removed learners are not counted). For the Gaussian and
#' inverse Gaussian families
#' \eqn{l_s'(\phi) = -n_+ / (2\phi)}{l_s'(phi) = -n_+ / (2 phi)}, so
#' \deqn{\hat\phi = \frac{D_p}{n_+ - M_p},}{phi = D_p / (n_+ - M_p),}
#' with \eqn{n_+} the number of observations with non-zero weight; for the
#' Gamma family the root is found numerically. The estimate needs neither
#' the edf of the model nor a cross-product of the model matrix, and costs
#' \eqn{\mathcal{O}(n + P)}{O(n + P)}. The restricted likelihood is exact
#' for Gaussian responses with the identity link and a Laplace
#' approximation otherwise. This is a different estimate from `mod$scale`,
#' the dispersion of the boosted fit (Fletcher, 2012; see [gbgam()]). The
#' smoothing parameters and \eqn{\phi} are held fixed: the posterior is
#' conditional on them, as `mgcv`'s `Vp` is conditional on the estimated
#' smoothing parameters.
#'
#' ## Undersmoothing
#'
#' A learner shrunk to a low edf has a posterior concentrated near zero.
#' With `undersmooth = TRUE`, every penalised learner that boosting
#' selected, with a finite \eqn{\lambda_k} and an edf below `edf_floor`
#' at the posterior mode, has its \eqn{\lambda_k} lowered until its edf
#' equals `edf_floor`; the posterior mode is then recomputed. The edf of a
#' learner is the sum of the diagonal elements, in its own columns, of the
#' shrinkage matrix \eqn{\mathbf{F}_j = (\mathbf{V}_j^\top \mathbf{W}_z
#' \mathbf{V}_j + \mathbf{S}_j)^{-1} \mathbf{V}_j^\top \mathbf{W}_z
#' \mathbf{V}_j}{F_j = (V_j' W_z V_j + S_j)^(-1) V_j' W_z V_j} of its smooth
#' \eqn{j}, with the other smooths fixed at the mode; for a
#' Demmler-Reinsch learner on its own, the sum of its shrinkage factors.
#' The learners of one smooth are lifted in turn, each at the smoothing
#' parameters of those already lifted. Learners at or above the floor,
#' removed learners and learners with a flat prior are left as they are.
#'
#' ## Stage 2: Bayesian backfitting
#'
#' Stage 2 starts with the steps above: it applies `linear`, computes the
#' posterior mode, undersmooths (and recomputes the mode), and estimates
#' \eqn{\phi} at the final mode, so their failures are Stage-2 errors (see
#' Errors). The posterior at the
#' resulting smoothing parameters is then sampled with one block for the
#' intercept, one per smooth, and the linear block with `linear = "free"`.
#' Gaussian responses
#' with the identity link use a Gibbs sampler (each block drawn from its
#' full conditional); all other families and links use
#' Metropolis-Hastings, with a proposal from one Fisher-scoring (IWLS) step
#' at the current state and the reverse proposal at the candidate.
#'
#' ## Intervals
#'
#' `$Vp` includes the intercept, so intervals can include its uncertainty,
#' with coverage across the function (Marra & Wood, 2012).
#' `plot(mod, seWithMean = TRUE)` draws such intervals and
#' `predict(mod, type = "iterms", se.fit = TRUE)` computes them, for centred
#' smooths in models whose only parametric term is the intercept;
#' [smooth_intervals()] gives them for any smooth. `plot()` and `predict()`
#' centre the intervals on the boosting estimate `mod$coefficients`;
#' [smooth_intervals()] can centre them on the posterior mean instead.
#'
#' ## Cost
#'
#' No step forms the full weighted cross-product of the model matrix; each
#' block uses only its own columns. For the Gaussian family with the
#' identity link, each block's weighted cross-product and Cholesky factor
#' are computed once, so a sampling sweep costs
#' \eqn{\mathcal{O}(nP + \sum_j p_j^2)}{O(nP + sum_j p_j^2)}, for \eqn{n}
#' observations, \eqn{p_j} coefficients in block \eqn{j} and
#' \eqn{P = \sum_j p_j}{P = sum_j p_j}. For other families,
#' Metropolis-Hastings builds a forward and a reverse proposal per block,
#' each with its own weighted cross-product and Cholesky factor, so a
#' sampling sweep costs
#' \eqn{\mathcal{O}(n \sum_j p_j^2 + \sum_j p_j^3)}{O(n sum_j p_j^2 +
#' sum_j p_j^3)}. In Stage 1, every proposal or halved step
#' that is checked costs a penalised fit: backfitting until the fit
#' converges (at most 1000 backfitting sweeps). After the fit at the
#' starting smoothing parameters, a sweep therefore runs at most \eqn{5K}
#' penalised fits. Their number, the number
#' of Stage-1 sweeps and the backfitting sweeps each fit needs all depend
#' on the data; backfitting needs more sweeps when smooths are concurve.
#'
#' ## Errors
#'
#' If Stage 1 fails, an error of class `"map_approx_failed"` is signalled;
#' if Stage 2 fails, an error of class `"bb_failed"`. Both carry a `$status`
#' naming the failure. Stage 1:
#' \describe{
#'   \item{`"not_converged"`}{no convergence within the allowed sweeps;}
#'   \item{`"no_mode"`}{the penalised fit at the starting smoothing
#'     parameters did not converge within 1000 backfitting sweeps, its
#'     working quantities failed, or no part of a backfitting sweep could be
#'     kept;}
#'   \item{`"search_failed"`}{the search over log(lambda_k) for one learner
#'     failed;}
#'   \item{`"nonfinite"`}{non-finite working quantities;}
#'   \item{`"not_pd"`}{a matrix that is not positive definite, also after
#'     the perturbation of the backfitting sweeps.}
#' }
#' Stage 2:
#' \describe{
#'   \item{`"no_mode"`}{no posterior mode;}
#'   \item{`"aliased"`}{a penalised learner whose \eqn{\lambda_k} is 0 (a
#'     flat prior) is aliased with other columns with a flat prior;}
#'   \item{`"no_residual_df"`}{no residual degrees of freedom for
#'     \eqn{\phi}: \eqn{n_+ \le M_p}{n_+ <= M_p} (see Dispersion);}
#'   \item{`"zero_dispersion"`}{an estimate of \eqn{\phi} that is not
#'     positive and finite (e.g. a zero penalised deviance);}
#'   \item{`"floor_unreachable"`}{with `undersmooth = TRUE`, a learner whose
#'     edf stays below `edf_floor` however small its \eqn{\lambda_k}, for
#'     example a floor above its number of columns;}
#'   \item{`"floor_no_root"`}{no crossing of the floor found between two
#'     values of log(lambda_k);}
#'   \item{`"nonfinite"`, `"not_pd"`}{as in Stage 1;}
#'   \item{`"nan"`}{a Metropolis-Hastings acceptance ratio that is not a
#'     number.}
#' }
#'
#' @seealso [gb_bayesian_ctrl()], [gb_uq_ctrl()],
#'   [smooth_intervals()], [gbgam()]
#'
#' @references
#' Buja, A., Hastie, T., & Tibshirani, R. (1989).
#' Linear smoothers and additive models.
#' *The Annals of Statistics, 17*(2), 453--510.
#' \href{https://doi.org/10.1214/aos/1176347115}{doi:10.1214/aos/1176347115}
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
#' Marra, G., & Wood, S. N. (2012).
#' Coverage properties of confidence intervals for generalized additive model
#' components.
#' *Scandinavian Journal of Statistics, 39*(1), 53--74.
#' \href{https://doi.org/10.1111/j.1467-9469.2011.00760.x}{doi:10.1111/j.1467-9469.2011.00760.x}
#'
#' Wood, S. N. (2011).
#' Fast stable restricted maximum likelihood and marginal likelihood
#' estimation of semiparametric generalized linear models.
#' *Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology), 73*(1), 3--36.
#' \href{https://doi.org/10.1111/j.1467-9868.2010.00749.x}{doi:10.1111/j.1467-9868.2010.00749.x}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' x <- sort(runif(120))
#' f <- sin(2 * pi * x)
#' y <- f + rnorm(120, sd = 0.3)
#' dat <- data.frame(x = x, y = y)
#'
#' mod <- gbgam(
#'   y ~ t2(x, k = 8, bs = "cr"),
#'   data = dat
#' )
#'
#' mod_bb <- bayesian_backfitting(mod, n_samples = 500, n_burnin = 100)
#'
#' # plot() now draws intervals. With these settings they are 95% intervals
#' # for the intercept plus the smooth, centred on the boosting estimate
#' a <- coef(mod_bb)[["(Intercept)"]]
#' plot(mod_bb, seWithMean = TRUE, se = qnorm(0.975), shift = a)
#' lines(x, f, col = "firebrick")
#'
#' # smooth_intervals() centres them on the posterior mean instead
#' ci <- smooth_intervals(mod_bb, centre = "posterior")
#' matlines(x, cbind(ci$lower, ci$upper), lty = 2, col = "steelblue")
#' legend("topright", lty = c(2, 2, 1),
#'        col = c("black", "steelblue", "firebrick"),
#'        legend = c("boosting estimate", "posterior mean", "true function"),
#'        bty = "n")
#' }
#' @export
bayesian_backfitting <- function(mod,
                                  n_samples        = 1000L,
                                  n_burnin         = 200L,
                                  linear           = c("free", "kept"),
                                  undersmooth      = FALSE,
                                  edf_floor        = 1,
                                  phi              = NULL,
                                  verbose          = FALSE,
                                  max_runtime_mins = NULL) {

  start <- proc.time()[["elapsed"]]  # the time limit covers the whole call
  checkmate::assert_class(mod, "gbgam")
  ctrl <- gb_bayesian_ctrl(n_samples = n_samples, n_burnin = n_burnin,
                           linear = linear, undersmooth = undersmooth,
                           edf_floor = edf_floor, phi = phi)
  linear <- ctrl$linear
  check_verbose(verbose)
  check_runtime_mins(max_runtime_mins)
  mod_in <- mod  # returned if the time limit leaves too few draws
  mod    <- remove_uq_fields(mod)

  post <- boost_to_backfit(mod, n_samples = n_samples, n_burnin = n_burnin,
                           linear = linear, undersmooth = undersmooth,
                           edf_floor = edf_floor, phi = phi,
                           verbose = verbose,
                           max_runtime_mins = max_runtime_mins, start = start)
  n_total   <- n_samples + n_burnin
  n_sampled <- max(0L, post$n_done - n_burnin)
  if (post$n_done < n_total) {
    if (n_sampled < 2L) {
      warning(sprintf(paste0(
        "bayesian_backfitting: time limit reached with fewer than 2 ",
        "posterior samples (%d / %d sweeps). Returning the model ",
        "unchanged, without posterior results."),
        post$n_done, n_total), call. = FALSE)
      return(mod_in)
    }
    warning(sprintf(
      "bayesian_backfitting: time limit reached after %d / %d sampling sweeps.",
      n_sampled, n_samples), call. = FALSE)
  }

  mod <- attach_posterior_draws(mod, post$beta_samples)
  # Gibbs: every draw accepted; 1 also for learners held at zero
  mod$mh_accept_rate <- if (is_gaussian_identity(mod$family))
    rep(1, length(mod$learners)) else post$accept
  mod$bb_info <- post$info
  mod$n_samples_done    <- as.integer(n_sampled)
  mod$uq_control <- gb_uq_ctrl("bayes", control = ctrl)
  mod
}


#' Posterior draws on the `mgcv` scale
#'
#' Maps each draw of the learner coefficients to the `mgcv` coefficients with
#' learners_to_mgcv_coef() and attaches their covariance `$Vp`, mean
#' `$coefficients_bayes` and the draws `$beta_samples`.
#' @param beta_samples List of `p_k x T` draws, one element per learner.
#' @noRd
attach_posterior_draws <- function(mod, beta_samples) {

  learners   <- mod$learners
  coef_names <- mod$coef_names
  ind0var    <- mod$ind0var
  cen_const  <- mod$cen_const
  K          <- length(learners)
  n_samples  <- ncol(beta_samples[[1L]])

  # Recover mgcv-scale coefficients for each sample
  p          <- length(mod$coefficients)
  coef_mat   <- matrix(NA_real_, nrow = n_samples, ncol = p,
                       dimnames = list(NULL, names(mod$coefficients)))

  # Temporarily swap in sample betas to reuse learners_to_mgcv_coef()
  for (s in seq_len(n_samples)) {
    learners_s <- learners
    for (k in seq_len(K)) {
      learners_s[[k]]$beta <- beta_samples[[k]][, s]
    }
    coef_mat[s, ] <- learners_to_mgcv_coef(
      learners   = learners_s,
      coef_names = coef_names,
      ind0var    = ind0var,
      cen_const  = cen_const,
      cen_absorb = mod$cen_absorb
    )
  }

  mod$Vp                 <- stats::cov(coef_mat)
  # Posterior mean on the mgcv scale; generally differs from the boosting
  # estimate mod$coefficients, on which predict() centres the intervals
  mod$coefficients_bayes <- colMeans(coef_mat)
  mod$beta_samples       <- beta_samples
  mod
}
