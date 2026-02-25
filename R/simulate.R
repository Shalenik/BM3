#' Simulate Climate Data
#'
#' Generate synthetic climate model outputs and observations for testing.
#' Creates a latent factor model with H latent factors, K climate models,
#' and realistic temporal structure.
#'
#' @param H_true Integer. True number of latent factors.
#' @param K Integer. Number of climate models (default 20).
#' @param years_start Integer. Start year (default 1850).
#' @param years_end Integer. End year (default 2100).
#' @param obs_end Integer. Last year of observations (default 2022).
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#'   - Y: matrix of climate model outputs (n_time × K)
#'   - z: vector of observations (n_time)
#'   - years: vector of years
#'   - y_obs_full: full observations including future
#'   - K: number of models
#'   - H_true: true H value
#'
#' @details
#' The generative model:
#' Y_kt ~ N(mu_k + U %*% beta_k(t), tau2_k)
#' z_t ~ N(mu + U %*% beta(t), sigma2)
#'
#' where U is a basis matrix (n_time × H) with smooth temporal patterns.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_data(H_true = 2, K = 20, seed = 123)
#' dim(sim$Y) # [1] 251 20 (years 1850-2100)
#' length(sim$z) # [1] 173 (years 1850-2022)
#' }
#'
#' @export
simulate_data <- function(H_true = 2, K = 20,
                          years_start = 1850, years_end = 2100,
                          obs_end = 2022, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    # Setup time grid
    years <- seq(years_start, years_end)
    obs_years <- seq(years_start, obs_end)
    n_time <- length(years)
    n_obs <- length(obs_years)

    message(sprintf("Simulating data: H=%d, K=%d, years %d-%d", H_true, K, years_start, years_end))

    # Create latent factors U with smooth temporal structure
    # Use B-spline or polynomial basis
    t_scaled <- (1:n_time) / n_time # [0, 1]
    U <- matrix(0, nrow = n_time, ncol = H_true)

    for (h in 1:H_true) {
        # Smooth sinusoidal patterns
        U[, h] <- sin(h * pi * t_scaled) + 0.5 * cos(h * pi * t_scaled / 2)
    }
    U <- scale(U) # Standardize

    # True time-varying weights for observations
    beta_true <- matrix(0, nrow = H_true, ncol = n_time)
    for (h in 1:H_true) {
        # Linear trend + noise
        beta_true[h, ] <- (h * t_scaled) + 0.1 * rnorm(n_time)
    }

    # Generate climate model outputs
    mu_Y <- rnorm(K, mean = 10, sd = 5) # Model-specific means
    tau2_Y <- rep(0.5, K) # Model-specific variances

    Y <- matrix(0, nrow = n_time, ncol = K)
    for (k in 1:K) {
        Y[, k] <- mu_Y[k] + U %*% (beta_true * (0.8 + 0.2 * rnorm(H_true))) +
            rnorm(n_time, 0, sqrt(tau2_Y[k]))
    }

    # Generate true observations
    mu_obs <- 15
    sigma2_obs <- 0.3
    z_true <- mu_obs + U %*% beta_true + rnorm(n_time, 0, sqrt(sigma2_obs))

    # Mask future observations (only keep past)
    z <- c(z_true[1:n_obs], rep(NA, n_time - n_obs))

    result <- list(
        Y = Y,
        z = z_true[1:n_obs], # Only observed part
        z_full = z_true, # Full including future
        years = years,
        years_obs = obs_years,
        K = K,
        H_true = H_true,
        U_true = U,
        beta_true = beta_true,
        mu_true = mu_obs,
        sigma2_true = sigma2_obs
    )

    message(sprintf("  Generated: Y is %d × %d, z is length %d", n_time, K, n_obs))

    return(result)
}
