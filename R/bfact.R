#' Bayesian Latent Factor Model for Climate Projections (BFACT)
#'
#' Runs a Gibbs sampler to estimate latent factor trajectories U, spline
#' coefficients gamma, loadings beta, observation/model precisions phiep,
#' and factor precisions phiet.
#'
#' @param Y Numeric matrix T x m of model outputs (each column a model series).
#' @param z Numeric vector of observed anomalies for years T1:T2.
#' @param T Integer total number of time points.
#' @param T1,T2 Integers giving the observation window indices (1 ≤ T1 ≤ T2 ≤ T).
#' @param H Integer number of latent factors.
#' @param iseed Integer RNG seed for reproducibility.
#' @param J Integer number of spline basis functions (natural splines via splines::ns).
#' @param nsim Integer number of posterior draws.
#' @param mga Prior precision (ridge) for spline coefficients gamma. Default 0.1.
#' @param mbe Prior precision (ridge) for loadings beta. Default 0.1.
#' @param aep Gamma shape for observation/model precisions phiep. Default 0.1.
#' @param bep Gamma rate  for observation/model precisions phiep. Default 0.1.
#' @param aet Gamma shape for factor precisions phiet. Default 0.1.
#' @param bet Gamma rate  for factor precisions phiet. Default 0.1.
#' @param a1  Gamma shape for shrinkage delta_1. Default 3.
#' @param a2  Gamma shape for shrinkage delta_{>=2}. Default 3.
#'
#' @return A list with elements: X, ysim, U_samples, beta_samples, gamma_samples
#' @importFrom splines ns
#' @importFrom stats rnorm dnorm rgamma
#' @export
BFACT <- function(Y, z, T, T1, T2, H, iseed, J, nsim,
                  mga = 0.1, mbe = 0.1,
                  aep = 0.1, bep = 0.1,
                  aet = 0.1, bet = 0.1,
                  a1 = 3, a2 = 3) {
    base::set.seed(iseed)

    m <- ncol(Y)
    M <- m + 1L

    # --- Design blocks & caches ---
    X <- splines::ns(1:T, J) # spline basis (T x J)
    SXX <- t(X) %*% X # cached crossprod for gamma updates
    E0 <- diag(T)[T1:T2, , drop = FALSE] # selector for obs window (n_obs x T)
    E0TE0 <- t(E0) %*% E0 # diagonal with ones in obs rows

    # --- State ---
    U <- matrix(0, T, H)
    gamma <- matrix(0, J, H)
    beta <- matrix(0, M, H + 1L)
    phiet <- rep(1, H)
    phiep <- rep(1, M)
    del <- rep(1, H)
    ysim <- matrix(NA_real_, nsim, T + 1L)

    # --- Storage ---
    gamma_samples <- array(NA_real_, dim = c(nsim, J, H))
    beta_samples <- array(NA_real_, dim = c(nsim, M, H + 1L))
    U_samples <- array(NA_real_, dim = c(nsim, T, H))
    phiep_samples <- matrix(NA_real_, nsim, M)
    beta0.sim <- numeric(nsim)

    # --- Gibbs iterations ---
    for (isim in 1:nsim) {
        if (isim %% 500 == 0) cat("Iteration", isim, "\n")

        # --- Update gamma ---
        for (h in 1:H) {
            c1 <- phiet[h] * prod(del[1:h])
            A <- diag(mga, J) + c1 * SXX
            b <- c1 * t(X) %*% U[, h]
            gamma[, h] <- mvn(A, b, J)
        }
        # gamma should be 6x2
        gamma_samples[isim, , ] <- gamma

        # --- Update beta ---
        for (k in 1:M) {
            Uk <- if (k == 1) cbind(1, U[T1:T2, ]) else cbind(1, U)
            A <- diag(mbe, H + 1) + phiep[k] * t(Uk) %*% Uk
            b <- phiep[k] * t(Uk) %*% (if (k == 1) z else Y[, k - 1])
            beta[k, ] <- mvn(A, b, H + 1)
        }
        # beta should be 16x3
        beta_samples[isim, , ] <- beta

        # --- Update U ---
        for (h in 1:H) {
            # Precision matrix Aprior: A = c1*I
            c1 <- phiet[h] * prod(del[1:h])
            A <- diag(1, T) * c1

            # Add contributions to A from CMIP and Observed
            for (k in 2:M) {
                A <- A + phiep[k] * beta[k, h + 1]^2 * diag(1, T)
            }
            A <- A + phiep[1] * beta[1, h + 1]^2 * E0TE0

            # Contribution from observed data (z)
            # Columns of U excluding the h-th factor (if H > 1)
            U_obs_other <- if (H > 1) U[T1:T2, -h, drop = FALSE] else matrix(0, T2 - T1 + 1, 0)
            # Corresponding beta weights, excluding intercept and h-th factor
            beta_obs_other <- if (H > 1) beta[1, -c(1, h + 1)] else numeric(0)
            # Compute residual: z minus fixed intercept beta[1,1] and contributions from other factors
            #   ytil_z ≈ z - (intercept + all other factor contributions)
            ytil_z <- z - beta[1, 1] - U_obs_other %*% beta_obs_other

            # Start b with prior term: c1 * X * gamma
            b <- c1 * X %*% gamma[, h]
            # Add contribution from observations:
            #   t(E0) %*% residuals = transpose design matrix times residuals
            b <- b + phiep[1] * beta[1, h + 1] * t(E0) %*% ytil_z


            # Contribution from climate models (Y)
            for (k in 2:M) {
                U_mod_other <- if (H > 1) U[, -h, drop = FALSE] else matrix(0, T, 0)
                beta_mod_other <- if (H > 1) beta[k, -c(1, h + 1)] else numeric(0)
                # Compute residual for model k: observed model minus expected
                ytil_Y <- Y[, k - 1] - beta[k, 1] - U_mod_other %*% beta_mod_other

                b <- b + phiep[k] * beta[k, h + 1] * ytil_Y
            }

            # Sample the h-th latent factor column from its conditional posterior:
            #   Normal(A^{-1} * b, A^{-1})
            U[, h] <- mvn(A, b, T)
        }

        # U should be 251x2
        U_samples[isim, , ] <- U

        # --- Update phiet ---
        for (h in 1:H) {
            resid <- U[, h] - X %*% gamma[, h]
            scale <- prod(del[1:h])
            shape <- aet + T / 2
            rate <- sum((resid^2) * scale) / 2 + bet
            phiet[h] <- rgamma(1, shape = shape, rate = rate)
        }

        # --- Update del ---
        del <- update_del(U, X, gamma, phiet, del, a1, a2)

        # --- Update phiep ---
        phiep[1] <- update_phiep1(z, U[T1:T2, ], beta[1, ], aep, bep)
        for (k in 2:M) {
            phiep[k] <- update_phiepk(Y[, k - 1], U, beta[k, ], aep, bep)
        }
        phiep_samples[isim, ] <- phiep

        # --- Generate ysim ---
        ysim[isim, ] <- c(phiep[1], beta[1, 1] + X %*% gamma %*% beta[1, -1])
        beta0.sim[isim] <- beta[1, 1]

        # --- Deviance ---
        loglik <- sum(dnorm(z,
            mean = cbind(1, U[T1:T2, ]) %*% beta[1, ],
            sd = sqrt(1 / phiep[1]), log = TRUE
        ))
        for (k in 2:M) {
            loglik <- loglik + sum(dnorm(Y[, k - 1],
                mean = cbind(1, U) %*% beta[k, ],
                sd = sqrt(1 / phiep[k]), log = TRUE
            ))
        }
    }

    return(list(
        X = X,
        H = H,
        nsim = nsim,
        ysim = ysim,
        U_samples = U_samples,
        beta_samples = beta_samples,
        gamma_samples = gamma_samples,
        phiep_samples = phiep_samples
    ))
}
