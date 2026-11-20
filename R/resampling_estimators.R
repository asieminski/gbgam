# Bootstrap estimators of the prediction error (Efron & Tibshirani, 1997),
# written for errors already computed elsewhere and vectorised, so that they
# apply to a whole path of errors at once (e.g. one per boosting iteration).
#
# Notation of Efron and Tibshirani (1997): err is the apparent (training)
# error, Err(1) the leave-one-out bootstrap (out-of-bag) error and gamma the
# no-information error.


#' The .632 bootstrap estimate of the prediction error
#'
#' Err(.632) = 0.368 err + 0.632 Err(1) (Efron & Tibshirani, 1997, eq. 24).
#'
#' @param apparent Apparent error err (numeric vector).
#' @param oob Out-of-bag error Err(1), of the same length.
#' @return The estimate, of the same length.
#' @noRd
err_632 <- function(apparent, oob) {
  0.368 * apparent + 0.632 * oob
}


#' The .632+ bootstrap estimate of the prediction error
#'
#' Efron and Tibshirani (1997), eqs. 28-32. The out-of-bag error is truncated
#' at the no-information error, Err(1)' = min(Err(1), gamma), so that the
#' relative overfitting rate R' lies in \[0, 1\] (eq. 31), and the estimate is
#' Err(.632) plus a correction that grows with R' (eq. 32). Where Err(1) is at
#' most gamma this equals (1 - w) err + w Err(1) with w = 0.632 /
#' (1 - 0.368 R') (eq. 29). Implementations differ where Err(1) > gamma > err:
#' `pec` and `peperr` return gamma there and `ipred` does not truncate R, while
#' eq. 32, used here, returns 0.632 Err(1) + 0.368 gamma (R' = 1).
#'
#' @param apparent Apparent error err (numeric vector).
#' @param oob Out-of-bag error Err(1), of the same length.
#' @param no_information No-information error gamma, of the same length.
#' @return A list with the estimate `risk`, the relative overfitting rate
#'   `relative_overfitting` (R') and the weight `weight` (w) on Err(1).
#' @references
#' Efron, B., & Tibshirani, R. (1997). Improvements on cross-validation: The
#' .632+ bootstrap method. *Journal of the American Statistical Association,
#' 92*(438), 548--560.
#' \href{https://doi.org/10.1080/01621459.1997.10474007}{doi:10.1080/01621459.1997.10474007}
#' @noRd
err_632plus <- function(apparent, oob, no_information) {
  # eq. 31: Err(1)' = min(Err(1), gamma), and R' = (Err(1)' - err) /
  # (gamma - err) where Err(1)' and gamma exceed err, else 0
  oob_truncated <- pmin(oob, no_information)
  overfits      <- oob_truncated > apparent & no_information > apparent
  relative_overfitting <- ifelse(
    overfits, (oob_truncated - apparent) / (no_information - apparent), 0)
  # eq. 32: Err(.632) plus the correction for overfitting
  risk <- err_632(apparent, oob) +
    (oob_truncated - apparent) * 0.368 * 0.632 * relative_overfitting /
      (1 - 0.368 * relative_overfitting)
  # eq. 29: the weight on Err(1)
  weight <- 0.632 / (1 - 0.368 * relative_overfitting)
  list(risk = risk, relative_overfitting = relative_overfitting,
       weight = weight)
}


#' No-information error of the deviance
#'
#' The deviance per unit weight if responses and fitted means were paired at
#' random: gamma = sum_i sum_j w_i w_j d(y_i, mu_j) / (sum_i w_i)^2, the
#' estimate of Efron and Tibshirani (1997, eq. 26), with the unit deviance d
#' of `family` and prior weights w. The unit deviance is in the q class of
#' Efron (2004, eq. 3.1): d(y, mu) = q(mu) + q'(mu) (y - mu) - q(y) for a
#' concave q, with -q'(mu) / 2 the natural parameter. Expanding this
#' definition gives, for any c, d(y, mu) = d(y, c) + d(c, mu) +
#' (y - c) (q'(mu) - q'(c)). With c = ybar, the weighted mean of y, the last
#' term vanishes in the weighted sum over i, because
#' sum_i w_i (y_i - ybar) = 0, so the double sum equals
#' sum_i w_i d(y_i, ybar) / sum_i w_i + sum_j w_j d(ybar, mu_j) / sum_j w_j
#' and takes \eqn{\mathcal{O}(n)}{O(n)} operations.
#'
#' @param y Response.
#' @param mu Fitted means.
#' @param w Prior weights.
#' @param family A family object with `dev.resids()`.
#' @return The no-information error (scalar).
#' @references
#' Efron, B. (2004). The estimation of prediction error: Covariance penalties
#' and cross-validation. *Journal of the American Statistical Association,
#' 99*(467), 619--632.
#' \href{https://doi.org/10.1198/016214504000000692}{doi:10.1198/016214504000000692}
#'
#' Efron, B., & Tibshirani, R. (1997). Improvements on cross-validation: The
#' .632+ bootstrap method. *Journal of the American Statistical Association,
#' 92*(438), 548--560.
#' \href{https://doi.org/10.1080/01621459.1997.10474007}{doi:10.1080/01621459.1997.10474007}
#' @noRd
no_information_deviance <- function(y, mu, w, family) {
  y_bar <- sum(w * y) / sum(w)
  y_to_y_bar  <- sum(family$dev.resids(y, y_bar, w))
  y_bar_to_mu <- sum(family$dev.resids(rep(y_bar, length(mu)), mu, w))
  (y_to_y_bar + y_bar_to_mu) / sum(w)
}
