#' Update precision for observations (k = 1)
#'
#' Full conditional \eqn{\phi_{\epsilon,1} \sim \mathrm{Gamma}(a + n/2,\; b + \sum r^2 / 2)}
#' where residuals \eqn{r = z - [1, U_{T1:T2}] \beta_{1,\cdot}}.
#'
#' @param z Numeric vector of observed values for years \eqn{T_1{:}T_2}.
#' @param U_sub Numeric matrix \eqn{(T_2-T_1+1) \times H}, rows \eqn{U_{T_1:T_2,\cdot}}.
#' @param beta_row Numeric length-\eqn{(H+1)} vector: intercept + loadings for \eqn{k=1}.
#' @param a Shape hyperparameter.
#' @param b Rate hyperparameter.
#'
#' @return A single \code{rgamma} draw for \eqn{\phi_{\epsilon,1}}.
#' @keywords internal
#' @noRd
update_phiep1 <- function(z, U_sub, beta_row, a, b) {
    resid <- z - cbind(1, U_sub) %*% beta_row
    shape <- length(z) / 2 + a
    rate <- b + sum(resid^2) / 2
    rgamma(1, shape, rate)
}

#' Update precision for model series (k >= 2)
#'
#' Full conditional \eqn{\phi_{\epsilon,k} \sim \mathrm{Gamma}(a + T/2,\; b + \sum r^2 / 2)}.
#'
#' @param Yk Numeric vector length \eqn{T}: model k series.
#' @param U Numeric matrix \eqn{T \times H}: latent factors.
#' @param beta_row Numeric length-\eqn{(H+1)} vector: intercept + loadings for model k.
#' @param a Shape hyperparameter.
#' @param b Rate hyperparameter.
#'
#' @return A single \code{rgamma} draw for \eqn{\phi_{\epsilon,k}}.
#' @keywords internal
#' @noRd
update_phiepk <- function(Yk, U, beta_row, a, b) {
    resid <- Yk - cbind(1, U) %*% beta_row
    shape <- length(Yk) / 2 + a
    rate <- b + sum(resid^2) / 2
    rgamma(1, shape, rate)
}

#' Update global shrinkage deltas \eqn{\delta_\ell} for \eqn{\tau_h = \prod_{\ell \le h} \delta_\ell}
#'
#' Implements the gamma full-conditionals for \eqn{\delta_1} and \eqn{\delta_\ell} (\eqn{\ell \ge 2})
#' following the multiplicative gamma process shrinkage prior.
#'
#' @param U Numeric matrix \eqn{T \times H}: latent factors.
#' @param X Numeric matrix \eqn{T \times J}: spline basis.
#' @param gamma Numeric matrix \eqn{J \times H}: spline coefficients.
#' @param phiet Numeric length-\eqn{H} vector of factor precisions \eqn{\phi_{\eta h}}.
#' @param del Numeric length-\eqn{H} current values of \eqn{\delta_\ell}.
#' @param a1 Shape hyperparameter for \eqn{\delta_1}.
#' @param a2 Shape hyperparameter for \eqn{\delta_{\ell\ge2}}.
#'
#' @return Numeric length-\eqn{H} updated \eqn{\delta} vector.
#' @keywords internal
#' @noRd
update_del <- function(U, X, gamma, phiet, del, a1, a2) {
    H <- ncol(U)
    T <- nrow(U)
    del_out <- del

    # Precompute RSS per factor: || U[,h] - X %*% gamma[,h] ||^2
    rss <- vapply(1:H, function(h) {
        r <- U[, h] - X %*% gamma[, h]
        sum(r * r)
    }, numeric(1))

    # δ1 | ·
    shape <- a1 + T * H / 2
    rate <- 1 + sum(vapply(1:H, function(h) {
        idx <- if (h >= 2) 2:h else integer(0)  # product of δ2…δh (empty product → 1)
        prod(del_out[idx]) * phiet[h] * rss[h] / 2
    }, numeric(1)))
    del_out[1] <- rgamma(1, shape, rate)

    # δℓ | · for ℓ = 2,…,H
    if (H > 1) {
        for (ell in 2:H) {
            shape <- T * (H - ell + 1) / 2 + a2
            rate <- 1 + sum(vapply(ell:H, function(h) {
                idx <- setdiff(seq_len(h), ell)  # δ1…δh except δℓ
                prod(del_out[idx]) * phiet[h] * rss[h] / 2
            }, numeric(1)))
            del_out[ell] <- rgamma(1, shape, rate)
        }
    }
    del_out
}
