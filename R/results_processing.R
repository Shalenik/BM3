#' Compute Posterior Mean Curve
#'
#' Returns the posterior mean curve from a BFACT object.
#' @param BFACT A BFACT object with a $ysim element.
#' @return Numeric vector of posterior means.
#' @export
posterior_mean_curve <- function(BFACT) {
    stopifnot(!is.null(BFACT$ysim))
    Y <- as.matrix(BFACT$ysim[floor(nrow(BFACT$ysim) / 2):nrow(BFACT$ysim), 2:ncol(BFACT$ysim), drop = FALSE])
    colMeans(Y, na.rm = TRUE)
}

#' Sampled Posterior Curves
#'
#' Returns a sample of posterior curves from a BFACT object.
#' @param BFACT A BFACT object with a $ysim element.
#' @return Numeric matrix of posterior curves.
#' @export
sample_posterior <- function(BFACT, ncurves = NULL) {
    stopifnot(!is.null(BFACT$ysim))
    if (is.null(ncurves)) {
        ncurves <- nrow(BFACT$ysim) / 2
    }
    indices <- sample(floor(nrow(BFACT$ysim) / 2):nrow(BFACT$ysim), ncurves)
    Y <- as.matrix(BFACT$ysim[indices, 2:ncol(BFACT$ysim), drop = FALSE])
    Y
}
