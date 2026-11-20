#' Mean, dmu/deta and variance at the linear predictor eta
#'
#' @param eta numeric vector; the linear predictor.
#' @param family a [stats::family] object (`linkinv`, `mu.eta`, `variance`).
#' @return list(mu, mu_eta, Vmu).
#' @noRd
family_eta_mu_derivs <- function(eta, family) {
  mu <- family$linkinv(eta)
  list(mu     = mu,
       mu_eta = family$mu.eta(eta),                # dmu/deta
       Vmu    = family$variance(mu))               # V(mu)
}

#' Generate random Rademacher matrix
#'
#' @param n A single integer specifying the number of rows.
#' @param p A single integer specifying the number of columns.
#'
#' @returns
#' An `n` by `p` matrix with entries randomly sampled from -1 and 1.
#' @noRd
rrademacher <- function(n, p){
  G <- sample(c(-1,1), n*p, replace = TRUE) |>
    matrix(n, p)
  return(G)
}

# UPDATE DERIVATIVES ------------------------------------------------------

#' Update derivatives
#'
#' @param mod A model object.
#' @param eta_new A numeric vector of the new linear predictor values.
#' @param k_opt A single integer specifying the learner index.
#' @param step_size The step size actually applied to the selected learner's
#'   update, i.e. `step_size * nu_k` when
#'   `use_learner_step_size_for_update = TRUE`, and `step_size` otherwise.
#' @param int_step_size Step size of the padded intercept update (intercept
#'   padding), or `NULL` (default) when the intercept is not padded. The
#'   intercept learner (learner 1) is then differentiated with respect to the
#'   same gradient \eqn{W du/dy}, and both increments are added to
#'   \eqn{d\eta/dy} in this iteration.
#'
#' @returns
#' The updated derivatives list from the model. Under intercept padding it
#' also holds `dedf_int`, the intercept increment's share of the change of
#' \eqn{\widetilde{\operatorname{edf}}}{edf~}: the increment
#' \eqn{\mathrm{diag}(h'(\eta^{[m+1]})) \mathbf{1} \Delta_{\mathrm{int}}}
#' of \eqn{d\tilde{\mu}/dy}: its exact trace, or the difference of the Hutch++
#' estimates with and without it. It is credited to learner 1, and is `NULL`
#' otherwise.
#' @noRd
update_derivatives <- function(mod, eta_new, k_opt, step_size = mod$step_size,
                               int_step_size = NULL){
  # Extract variables
  derivs <- mod$derivs
  fam    <- mod$family
  # d(u_m)/d(y) must use mu_m (pre-update) because u_m = omega(mu_m)*(y-mu_m).
  mu_old <- fam$linkinv(mod$linear.predictors)
  lopt   <- mod$learners[[k_opt]]

  # Compute derivatives. Wdudy is W %*% du/dy (prior weights folded in).
  pad    <- !is.null(int_step_size)  # intercept padding
  Wdudy  <- du_dy(derivs$dmudy, mu = mu_old, y = mod$y, family = fam,
                  GGG = derivs$GGG, weights = mod$prior.weights)
  Delta_dgammady <- delta_dgamma_dy(lopt, Wdudy, step_size = step_size)
  if (pad) {
    # Padded intercept: step * 1^T (W du/dy) / sum(w), from the same Wdudy
    lint      <- mod$learners[[1L]]
    Delta_int <- delta_dgamma_dy(lint, Wdudy, step_size = int_step_size)
  }
  rm(Wdudy)  # n x n without sketching: free it early
  detady <- derivs$detady + lopt$design %*% Delta_dgammady
  if (pad) detady <- detady + lint$design %*% Delta_int
  # d(mu_{m+1})/d(eta_{m+1}) uses the post-update eta.
  dmudy  <- dmu_dy(eta_new, detady, fam)

  # Save
  derivs$detady <- detady
  derivs$dmudy  <- dmudy

  derivs$dedf_int <- NULL
  if (pad) {
    derivs$dedf_int <- padded_intercept_edf(derivs, dmudy, eta_new, fam,
                                            lint, Delta_int)
  }

  return(derivs)
}

#' Share of the padded intercept in the change of edf~
#'
#' edf~ with minus without the intercept's increment
#' diag(h'(eta_new)) 1 Delta_int of dmu~/dy. It is credited to learner 1, so
#' that the ASR update of k_opt only sees its own share.
#' @param lint The intercept learner; `Delta_int` its increment of
#'   d(gamma_1)/dy.
#' @noRd
padded_intercept_edf <- function(derivs, dmudy, eta_new, fam, lint,
                                 Delta_int) {
  hp <- if (identical(fam$link, "identity")) 1 else
    fam$mu.eta(as.vector(eta_new))
  a_int <- hp * as.vector(lint$design)
  if (is.null(derivs$GGG)) {
    # Exact: the trace of the rank-one increment a_int %*% Delta_int
    sum(a_int * as.vector(Delta_int))
  } else {
    # Hutch++ is not linear in its argument, so take the difference of
    # the two estimates (as a separate intercept update would)
    Gs   <- derivs$Gs %||% split_mat_cols(derivs$GGG, derivs$zs)
    E_wo <- dmudy - outer(a_int, as.vector(Delta_int))
    hutchpp(Gs, split_mat_cols(dmudy, derivs$zs)) -
      hutchpp(Gs, split_mat_cols(E_wo, derivs$zs))
  }
}

#' dmu~/dy = diag(h'(eta)) deta/dy
#'
#' For the identity link h' = 1, so the same matrix is returned (no copy).
#' @noRd
dmu_dy <- function(eta, detady, family) {
  if (identical(family$link, "identity")) return(detady)
  family$mu.eta(as.vector(eta)) * detady
}

#' Compute the prior-weighted derivative of the pseudo-residual vector \eqn{u}
#'
#' @param dmudy A Jacobian matrix where
#'   `dmudy[i,j]` \eqn{= d \tilde{\mu}_i/d y_j} (or its sketch
#'   `dmudy %*% GGG`).
#' @param mu A numeric vector of the fitted mean values.
#' @param y A numeric vector of the observed responses.
#' @param family A family object.
#' @param GGG Optional. A random sub-Gaussian matrix for sketching or `NULL`.
#' @param weights Optional prior weights \eqn{w}. If supplied, the rows are
#'   scaled by \eqn{w}, i.e. \eqn{W\,du/dy} is returned.
#'
#' @returns
#' A matrix representing the (prior-weighted) derivative of the residual
#' vector with respect to `y` (sketched by `GGG` if supplied).
#' @noRd
du_dy <- function(dmudy, mu, y, family, GGG, weights = NULL){
  # u_i = omega(mu_i)*(y_i - mu_i)
  # du_i/dy_j = omega'_i*(y_i-mu_i)*dmudy[i,j]
  #             + omega_i*(delta_{ij} - dmudy[i,j])
  #           = (omega'_i*(y_i-mu_i) - omega_i)*dmudy[i,j] + omega_i*delta_{ij}
  # (omega' is the derivative d omega / d mu)
  # Row scalings are length-n vectors, so folding the prior weights in here
  # avoids a separate n x z (or n x n) copy.
  w    <- if (is.null(weights)) 1 else weights
  omu  <- omega(mu, family)
  domu <- domega_dmu(mu, family)
  dudy <- (w * (domu * (y - mu) - omu)) * dmudy
  if (is.null(GGG)) {
    # Add w * omega to the diagonal via linear indexing (no `diag<-` copy)
    n         <- length(mu)
    idx       <- seq.int(1L, by = n + 1L, length.out = n)
    dudy[idx] <- dudy[idx] + w * omu
  } else {
    dudy <- dudy + (w * omu) * GGG
  }
  return(dudy)
}

#' Derivative of omega with respect to mu
#'
#' @description
#' Computes the derivative of \eqn{\omega(\mu) = (V(\mu) \cdot g'(\mu))^{-1}}
#' with respect to \eqn{\mu}, where \eqn{V(\mu)} is the variance function and
#' \eqn{g'(\mu)} is the derivative of the link function.
#'
#' \deqn{\frac{d\omega}{d\mu} = -\frac{V'(\mu)\,g'(\mu) + V(\mu)\,g''(\mu)}{(V(\mu)\,g'(\mu))^2}}
#'
#' @param mu A numeric vector of means.
#' @param family A family object with `variance`, `linkfun`, `mu.eta`,
#'  `dvar`, and `d2link` methods (see [mgcv::fix.family.link()] and
#'  [mgcv::fix.family.var()]).
#' @noRd
domega_dmu <- function(mu, family){
  V <- family$variance(mu)
  eta <- family$linkfun(mu)
  dg <- 1/family$mu.eta(eta)
  dV <- family$dvar(mu)
  d2g <- family$d2link(mu)
  return(-(dV*dg + V*d2g) / ((V*dg)^2))
}

#' Increment of d(gamma_k)/dy for the selected learner k
#'
#' step_size * F_k X_k^T W du/dy, where `Wdudy` already includes the prior
#' weights W (see du_dy()) and F_k is the diagonal shrinkage of the
#' learner (1 / (1 + lambda D) for Demmler-Reinsch learners, 1 / x^T W x
#' for unpenalised ones).
#' @noRd
delta_dgamma_dy <- function(learner, Wdudy, step_size = 1) {
  step_F <- if (inherits(learner, "demmler_reinsch")) {
    step_size / (1 + learner$D * exp(learner$rho))
  } else {
    # A zero-norm column (no weight on its rows) has a zero fit (fit_learner())
    xtWx <- as.vector(learner$xtWx)
    if (xtWx > 0) step_size / xtWx else 0
  }
  step_F * crossprod(learner$design, Wdudy)
}

# COMPUTE EDF -------------------------------------------------------------

#' Model edf~: tr(dmu~/dy), exact or by Hutch++
#' @noRd
model_edf <- function(mod, derivs) {
  if (mod$edf_compute_method == "exact"){
    sum(diag(derivs$dmudy))
  } else {
    Es <- split_mat_cols(derivs$dmudy, derivs$zs)
    # The probe blocks are constant, so they are split once in init_derivs()
    Gs <- derivs$Gs %||% split_mat_cols(derivs$GGG, derivs$zs)
    hutchpp(Gs, Es)
  }
}



#' Split matrix columns into a list of sub-matrices
#'
#' @description
#' Splits the columns of a matrix into a list of sub-matrices according to
#' the column counts specified in `zs`.
#'
#' @param mat A matrix to split.
#' @param zs An integer vector specifying the number of columns in each split.
#'
#' @returns A list of matrices, each containing the corresponding columns.
#' @noRd
split_mat_cols <- function(mat, zs) {
  checkmate::assert_matrix(mat)
  # seq_len() so that an empty block gives no columns, not a reversed range
  offsets <- c(0, cumsum(zs)[-length(zs)])
  zs_ind  <- Map(function(o, k) o + seq_len(k), offsets, zs)
  lapply(zs_ind, \(z) mat[, z, drop = FALSE])
}

#' NA-Hutch++ estimate of tr(d mu~/dy)
#'
#' @description
#' Implements a right-sketch-only variant of the NA-Hutch++ stochastic trace
#' estimator for the trace of the Jacobian of the fitted means with respect
#' to the response, written here
#' \eqn{J = d\tilde{\mu}/dy}{J = d mu~/dy}. Given sketches
#' \eqn{E_j = J G_j}, the trace is split into a low-rank component and a
#' Hutchinson correction:
#' \deqn{\mathrm{tr}(J) \approx \mathrm{tr}(\tilde J) +
#'   \frac{1}{z_3} \mathrm{tr}\{G_3^\top (J - \tilde J) G_3\},}
#' where \eqn{z_3} = `ncol(Gs[[3]])` and
#' \eqn{\tilde J = E_1 (G_2^\top E_1)^{+} E_2^\top
#'   = J G_1 (G_2^\top J G_1)^{+} G_2^\top J^\top}.
#' This equals the NA-Hutch++ low-rank term of Meyer et al. (2021) only when
#' \eqn{J} is symmetric, because the left sketch \eqn{G_2^\top J} is not
#' available from forward-mode sketches and is replaced by
#' \eqn{(J G_2)^\top}. The estimator stays unbiased for any \eqn{J}, because
#' \eqn{G_3} is independent of \eqn{G_1} and \eqn{G_2}.
#'
#' @param Gs A list of 3 sketching matrices.
#' @param Es A list of 3 sketched matrices \eqn{E_j = J G_j}.
#'
#' @returns A scalar estimate of \eqn{\mathrm{tr}(J)}.
#'
#' @references
#' Meyer, R. A., Musco, C., Musco, C., & Woodruff, D. P. (2021).
#' Hutch++: Optimal stochastic trace estimation.
#' In *Symposium on Simplicity in Algorithms (SOSA)* (pp. 142--155). Society
#' for Industrial and Applied Mathematics.
#' \href{https://doi.org/10.1137/1.9781611976496.16}{doi:10.1137/1.9781611976496.16}
#'
#' @noRd
hutchpp <- function(Gs, Es) {
  core   <- MASS::ginv(crossprod(Gs[[2]], Es[[1]]))      # (G2^T E1)^+
  tr_low <- tr_prod(core, crossprod(Es[[2]], Es[[1]]))   # tr(J~)

  # Hutchinson estimate of tr(J - J~) with the probes G3
  z3        <- ncol(Gs[[3]])
  G3tE1     <- crossprod(Gs[[3]], Es[[1]])
  G3tE2     <- crossprod(Gs[[3]], Es[[2]])
  hutch_J   <- tr_xprod(Gs[[3]], Es[[3]]) / z3
  # tr(G3^T J~ G3) = tr(core (G3^T E2)^T (G3^T E1))
  hutch_low <- tr_prod(core, crossprod(G3tE2, G3tE1)) / z3

  tr_low + (hutch_J - hutch_low)
}
