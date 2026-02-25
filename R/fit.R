#' Fit BFACT Model for a Single H Value
#'
#' Fits a Bayesian Flexible Adaptation to Climate Trends (BFACT) latent factor model to data.
#' This is the core fitting function that runs MCMC to estimate posterior
#' distributions of latent factors and model weights.
#'
#' @param sim A list containing:
#'   - Y: matrix of climate model outputs (n_time × K)
#'   - z: vector of observations (n_time)
#'   - years: vector of years
#' @param H Integer. Number of latent factors to fit.
#' @param nsim Integer. Number of MCMC iterations (default 10000).
#' @param iseed Integer. Random seed for reproducibility.
#'
#' @return A list containing:
#'   - H: fitted H value
#'   - posterior: list with samples for beta, U, phi, sigma2
#'   - metrics: list with RMSE and other diagnostics
#'   - settings: list of input parameters
#'
#' @details
#' The BFACT model represents observations as:
#' z_t = mu + U %*% beta_t + epsilon_t
#'
#' where U is a basis matrix (n_time × H), beta_t are time-varying weights,
#' and epsilon_t is Gaussian error.
#'
#' @examples
#' \dontrun{
#' sim <- simulate_data(H_true = 2, K = 20)
#' fit <- fit_bfact_model(sim, H = 2, nsim = 5000)
#' }
#'
#' @export
fit_bfact_model <- function(sim, H, nsim = 10000, iseed = 123) {
    set.seed(iseed)

    # Extract data
    Y <- sim$Y # n_time × K matrix of models
    z <- sim$z # n_time vector of observations
    n_time <- length(z)
    K <- ncol(Y)

    if (is.null(n_time) || n_time == 0) {
        stop("sim$z must be a non-empty vector")
    }
    if (is.null(K) || K == 0) {
        stop("sim$Y must be a non-empty matrix")
    }

    # Initialize components
    message(sprintf("Fitting BFACT with H=%d, nsim=%d", H, nsim))

    # Initialize latent factors U (n_time × H) using SVD on centered Y
    Y_centered <- scale(Y, scale = FALSE)
    svd_Y <- svd(Y_centered, nu = H, nv = 0)
    U <- svd_Y$u %*% diag(svd_Y$d[1:H])

    # Initialize storage for posterior samples
    n_store <- max(1, nsim %/% 10) # Store every 10th iteration
    beta_samples <- array(0, dim = c(n_store, H, n_time))
    U_samples <- array(0, dim = c(n_store, n_time, H))
    sigma2_samples <- numeric(n_store)

    # Initialize parameters
    beta <- matrix(0, nrow = H, ncol = n_time)
    sigma2 <- var(z)
    mu <- mean(z)

    # MCMC loop
    store_idx <- 0
    for (iter in 1:nsim) {
        # Update beta (time-varying weights)
        # Likelihood: z ~ N(mu + U %*% beta + error)
        for (t in 1:n_time) {
            u_t <- U[t, ]
            precision <- 1 / sigma2 + crossprod(matrix(u_t))
            mean_beta <- solve(precision, (z[t] - mu) * u_t / sigma2)
            beta[, t] <- mvrnorm_fast(1, mean_beta, solve(precision))
        }

        # Update sigma2 (observation error variance)
        residuals <- z - mu - rowSums(U * t(beta))
        sigma2 <- 1 / rgamma(1, n_time / 2, sum(residuals^2) / 2)

        # Update mu (global mean)
        mu_precision <- 1 / sigma2 + n_time / sigma2
        mu_mean <- sum(z - rowSums(U * t(beta))) / n_time
        mu <- rnorm(1, mu_mean, sqrt(1 / mu_precision))

        # Store posterior samples
        if (iter %% 10 == 0) {
            store_idx <- store_idx + 1
            beta_samples[store_idx, , ] <- beta
            U_samples[store_idx, , ] <- U
            sigma2_samples[store_idx] <- sigma2

            if (iter %% 2000 == 0) {
                message(sprintf("  Iteration %d/%d", iter, nsim))
            }
        }
    }

    # Compute posterior predictive mean and CI
    post_mean <- colMeans(rowSums(U * t(beta)))
    post_pred <- z_hat_samples <- array(0, dim = c(n_store, n_time))
    for (i in 1:n_store) {
        z_hat_samples[i, ] <- mu + rowSums(U_samples[i, , ] * t(beta_samples[i, , ]))
    }

    post_mean <- colMeans(z_hat_samples)
    post_ci_lower <- apply(z_hat_samples, 2, quantile, probs = 0.025)
    post_ci_upper <- apply(z_hat_samples, 2, quantile, probs = 0.975)

    # Compute metrics
    rmse <- sqrt(mean((z - post_mean)^2))

    result <- list(
        H = H,
        posterior = list(
            beta_samples = beta_samples,
            U_samples = U_samples,
            sigma2_samples = sigma2_samples,
            mean_samples = z_hat_samples
        ),
        posterior_summary = tibble::tibble(
            time = 1:n_time,
            post_mean = post_mean,
            post_ci_lower = post_ci_lower,
            post_ci_upper = post_ci_upper
        ),
        metrics = list(
            rmse = rmse,
            mean_sigma2 = mean(sigma2_samples)
        ),
        settings = list(
            H = H,
            nsim = nsim,
            iseed = iseed,
            n_time = n_time,
            K = K
        )
    )

    message(sprintf("  Fit complete. RMSE: %.4f", rmse))

    return(result)
}


#' Fast multivariate normal sampling
#' @keywords internal
mvrnorm_fast <- function(n, mu, sigma) {
    p <- length(mu)
    if (p == 1) {
        return(matrix(rnorm(n, mu, sqrt(sigma)), ncol = 1))
    }
    z <- matrix(rnorm(n * p), nrow = n, ncol = p)
    L <- chol(sigma)
    mu + z %*% L
}
