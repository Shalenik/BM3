#' Compute Posterior Mean Curve
#'
#' Returns the posterior mean curve from a BM3 object.
#' @param BM3 A BM3 object with a $ysim element.
#' @return Numeric vector of posterior means.
#' @export
posterior_mean_curve <- function(BM3) {
    stopifnot(!is.null(BM3$ysim))
    Y <- as.matrix(BM3$ysim[floor(nrow(BM3$ysim) / 2):nrow(BM3$ysim), 2:ncol(BM3$ysim), drop = FALSE])
    colMeans(Y, na.rm = TRUE)
}

#' Sampled Posterior Curves
#'
#' Returns a sample of posterior curves from a BM3 object.
#' @param BM3 A BM3 object with a $ysim element.
#' @return Numeric matrix of posterior curves.
#' @export
sample_posterior <- function(BM3, ncurves = NULL) {
    stopifnot(!is.null(BM3$ysim))
    if (is.null(ncurves)) {
        ncurves <- nrow(BM3$ysim) / 2
    }
    indices <- sample(floor(nrow(BM3$ysim) / 2):nrow(BM3$ysim), ncurves)
    Y <- as.matrix(BM3$ysim[indices, 2:ncol(BM3$ysim), drop = FALSE])
    Y
}
