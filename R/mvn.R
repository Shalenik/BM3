#' Draw from N(A^{-1} b, A^{-1}) via Cholesky
#'
#' Generates an n-dimensional draw from a multivariate normal with mean \eqn{A^{-1} b}
#' and covariance \eqn{A^{-1}}, using a Cholesky factorization of \eqn{A}.
#'
#' @details
#' This uses the identity:
#' \deqn{y = A^{-1} b + A^{-1/2} \, \varepsilon, \qquad \varepsilon \sim N(0, I_n)}
#' implemented by solving with the Cholesky factor \eqn{C} such that \eqn{A = C^\top C}.
#'
#' @param A Numeric \eqn{n \times n} positive-definite precision matrix.
#' @param b Numeric length-\eqn{n} vector (right-hand side).
#' @param n Integer dimension (\eqn{n}); must match \code{nrow(A)} and \code{length(b)}.
#'
#' @return Numeric length-\eqn{n} vector sampled from \eqn{N(A^{-1} b, A^{-1})}.
#' @keywords internal
#' @noRd
mvn <- function(A, b, n) {
    # C is upper-triangular Cholesky factor of A: A = t(C) %*% C
    C <- chol(A)
    # CI = C^{-1}
    CI <- solve(C)
    # Sample: y = C^{-1} ( C^{-T} b + z ),  z ~ N(0, I)
    y <- CI %*% (t(CI) %*% b + rnorm(n))
    return(y)
}
