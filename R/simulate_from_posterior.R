#' Simulate Data from a Posterior Draw of a Fitted BM3 Object
#'
#' Generates a simulated dataset from a single posterior draw of a fitted BM3 object.
#' @param fit A fitted BM3 object (as returned by BM3()).
#' @param s Integer index of the posterior draw to use (default: random draw).
#' @param years_start First year (default: 1850).
#' @param years_end Last year (default: 2100).
#' @param T2 Last observed time index (default: 173).
#' @param T1 First observed time index (default: 1).
#' @param J Number of spline basis functions (default: 6).
#' @param OBS_COL Which column of Y is the observed series (default: 1).
#' @param seed Random seed for reproducibility (default: NULL).
#' @return A list with simulated data and metadata.
#' @export
simulate_from_posterior <- function(
  fit,
  s = NULL,
  years_start = 1850,
  years_end = 2100,
  T2 = 173,
  T1 = 1,
  J = 6,
  OBS_COL = 1,
  seed = NULL
) {
    if (!is.null(seed)) set.seed(seed)
    beta_all <- fit$beta_samples
    u_all <- fit$U_samples
    phi_all <- fit$phiep_samples
    S <- dim(beta_all)[1]
    K <- dim(beta_all)[2]
    T_post <- dim(u_all)[2]
    H <- dim(u_all)[3]
    years <- seq.int(years_start, years_end)
    T <- length(years)
    if (T_post != T) stop("Time length mismatch: posterior T=", T_post, " but years imply T=", T)
    if (T2 < 1 || T2 > T) stop("Invalid T2: ", T2, " (must be between 1 and T=", T, ")")
    if (OBS_COL < 1 || OBS_COL > K) stop("Invalid OBS_COL: ", OBS_COL, " (K=", K, ")")
    if (is.null(s)) s <- sample.int(S, 1)
    beta_s <- beta_all[s, , , drop = TRUE]
    u_s <- u_all[s, , , drop = TRUE]
    phi_s <- phi_all[s, , drop = TRUE]
    mu <- matrix(0, nrow = T, ncol = K)
    for (k in seq_len(K)) {
        mu[, k] <- beta_s[k, 1]
        for (h in seq_len(H)) {
            mu[, k] <- mu[, k] + beta_s[k, h + 1] * u_s[, h]
        }
    }
    Y <- matrix(NA_real_, nrow = T, ncol = K)
    for (k in seq_len(K)) {
        sd_k <- sqrt(1 / phi_s[k])
        Y[, k] <- rnorm(T, mean = mu[, k], sd = sd_k)
    }
    y_obs_full <- as.numeric(Y[, OBS_COL])
    y_masked <- y_obs_full
    if (T2 < T) y_masked[(T2 + 1):T] <- NA_real_
    z <- y_obs_full[T1:T2]
    sim <- list(
        years      = years,
        T          = T,
        T1         = T1,
        T2         = T2,
        J          = J,
        H_true     = H,
        K          = K,
        Y          = Y,
        mu         = mu,
        z          = z,
        y_obs_full = y_obs_full,
        y_masked   = y_masked
    )
    return(sim)
}
