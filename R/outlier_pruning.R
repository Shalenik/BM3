#' RMSE Threshold Function (Mean + k SD)
#'
#' Computes a threshold for outlier detection as the mean plus
#' `dial` standard deviations.
#'
#' @param rmse_vals A numeric vector of RMSE values.
#' @param dial Numeric multiplier for the standard deviation. Default is 2.5.
#'
#' @return A numeric threshold value.
#' @export
mean_sd_threshold <- function(rmse_vals, dial = 2.5) {
    mean(rmse_vals, na.rm = TRUE) + dial * stats::sd(rmse_vals, na.rm = TRUE)
}

#' Detect Outlier Climate Models Using RMSE Thresholding
#'
#' Iteratively removes climate model columns whose RMSE from the ensemble mean
#' exceeds a user-defined threshold. Assumes input is a matrix or data frame
#' with time series in rows and models in columns.
#'
#' @param data A numeric matrix or data frame with model outputs.
#' @param start_index Integer, starting row index (e.g., corresponding to 2015).
#' @param end_index Integer, ending row index (e.g., corresponding to 2100).
#' @param threshold_fn A function taking a numeric vector of RMSE values and
#'   returning an RMSE threshold. Default is [mean_sd_threshold()].
#' @param threshold_dial Numeric multiplier passed to `threshold_fn`. Default is 2.5.
#'
#' @return A character vector of column names identified as outliers.
#' @export
detect_outliers <- function(data, start_index, end_index,
                            threshold_fn = mean_sd_threshold,
                            threshold_dial = 2.5) {
    if (is.null(dim(data))) {
        stop("`data` must be a matrix or data frame.")
    }

    data <- as.matrix(data)
    n_rows <- nrow(data)
    n_cols <- ncol(data)

    if (n_rows == 0 || n_cols == 0) {
        return(character(0))
    }

    if (start_index < 1 || end_index > n_rows || start_index > end_index) {
        stop("`start_index` and `end_index` must be within data row bounds.")
    }

    outliers <- character(0)
    data <- data[start_index:end_index, , drop = FALSE]
    flag <- FALSE

    rmse <- function(x, y) {
        sqrt(mean((x - y)^2, na.rm = TRUE))
    }

    while (!flag) {
        ensemble_mean <- rowMeans(data, na.rm = TRUE)
        rmse_values <- apply(data, 2, function(model) rmse(model, ensemble_mean))
        max_rmse_index <- which.max(rmse_values)
        max_rmse_value <- rmse_values[max_rmse_index]
        threshold <- threshold_fn(rmse_values, threshold_dial)

        if (max_rmse_value > threshold) {
            outliers <- c(outliers, colnames(data)[max_rmse_index])
            data <- data[, -max_rmse_index, drop = FALSE]
        } else {
            flag <- TRUE
        }

        if (ncol(data) == 0) {
            flag <- TRUE
        }
    }

    unique(outliers)
}