#' Find Elbow Point in a Curve (Perpendicular Distance Method)
#'
#' Identifies the "elbow" point in a metric-vs-H curve using the perpendicular distance method.
#' Useful for model selection when fitting multiple H values (e.g., identifying the optimal H
#' by finding where RMSE or variance improvements plateau).
#'
#' @param x Numeric vector of x-values (typically H values).
#' @param y Numeric vector of y-values (typically RMSE or variance).
#' @param minimize Logical. If TRUE (default), assumes lower y is better (e.g., RMSE).
#'   If FALSE, assumes higher y is better (e.g., likelihood).
#'
#' @return Numeric. The x-value (H) at the identified elbow point.
#'
#' @details
#' The elbow method finds where the curve transitions from steep to flat. The algorithm:
#' 1. Identifies the decreasing region of the curve
#' 2. Normalizes x and y to [0,1] scale
#' 3. Draws a line from the start to the end of the decreasing region
#' 4. Calculates perpendicular distance from each point to this line
#' 5. Returns the x-value (H) with maximum perpendicular distance
#'
#' This is particularly useful for model selection when the optimal H is unclear
#' from the metric alone.
#'
#' @examples
#' \dontrun{
#' # Find elbow in RMSE curve
#' H_values <- 2:7
#' rmse_values <- c(0.85, 0.42, 0.38, 0.36, 0.35, 0.35)
#' elbow_h <- find_elbow(H_values, rmse_values, minimize = TRUE)
#' print(elbow_h) # Likely returns H around 3-4
#' }
#'
#' @export
find_elbow <- function(x, y, minimize = TRUE) {
    if (length(x) < 3) {
        return(x[1])
    }

    # If maximizing, flip the curve
    if (!minimize) {
        y <- -y
    }

    # Skip any initial increase - find where decrease starts
    start_idx <- 1
    for (i in 1:(length(y) - 1)) {
        if (y[i + 1] < y[i]) {
            start_idx <- i
            break
        }
    }

    # Find the minimum from the decreasing region onwards
    y_subset <- y[start_idx:length(y)]
    min_idx_rel <- which.min(y_subset)
    min_idx <- start_idx + min_idx_rel

    # Only analyze points from start of decrease to minimum
    x_mono <- x[start_idx:min_idx]
    y_mono <- y[start_idx:min_idx]

    if (length(x_mono) < 3) {
        return(x_mono[1])
    }

    # Normalize x and y to 0-1 scale for fair distance calculation
    x_norm <- (x_mono - min(x_mono)) / (max(x_mono) - min(x_mono))
    y_norm <- (y_mono - min(y_mono)) / (max(y_mono) - min(y_mono))

    # Line from first to last point in the decreasing region
    x1 <- x_norm[1]
    y1 <- y_norm[1]
    x2 <- x_norm[length(x_norm)]
    y2 <- y_norm[length(y_norm)]

    # Calculate perpendicular distance from each point to the line
    distances <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
        sqrt((y2 - y1)^2 + (x2 - x1)^2)

    # The elbow is the point with maximum perpendicular distance
    max_dist <- max(distances)
    candidates <- which(distances >= max_dist)
    elbow_idx <- candidates[1] # Take the first (earliest H) among close candidates

    return(x_mono[elbow_idx])
}
