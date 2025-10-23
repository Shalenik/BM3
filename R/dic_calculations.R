# DIC ----
#' Compute Deviance Information Criterion (DIC)
#'
#' Calculates the Deviance Information Criterion and its components using posterior
#' samples from a Bayesian model fit.
#' Called by BM3 function during the MCMC sampling process.
#'
#' @param X Spline basis matrix (T x J)
#' @param Y Matrix of model outputs (T x m)
#' @param z Observational vector of length (T2 - T1 + 1)
#' @param gammas Array of gamma samples (nsim x J x H)
#' @param betas Array of beta samples (nsim x M x H+1)
#' @param phieps Matrix of precision samples (nsim x M)
#' @param Us Array of latent factor samples (nsim x T x H)
#' @param deviances Vector of deviance values (length nsim)
#' @param T1 Starting time index for observational window
#' @param T2 Ending time index for observational window
#'
#' @return A list containing DIC, effective number of parameters, mean deviance,
#'         and deviance at posterior means
dic.calculations <- function(X, Y, z, gammas, betas, phieps, Us, deviances, T1, T2) {
    nsim <- length(deviances)
    M <- ncol(phieps)
    T <- nrow(Y)
    noburn <- (nsim / 2):nsim # Second half of samples

    # Mean deviance from posterior
    mean_deviance <- mean(deviances[noburn])

    # Posterior means of parameters
    mean_gamma <- apply(gammas[noburn, , ], c(2, 3), mean)
    mean_beta <- apply(betas[noburn, , ], c(2, 3), mean)
    mean_phiep <- apply(phieps[noburn, ], 2, mean)
    mean_U <- apply(Us[noburn, , ], c(2, 3), mean)

    # Reconstruct U1 and UK using mean_U
    U1_mean <- cbind(1, mean_U[T1:T2, ])
    UK_mean <- cbind(1, mean_U)

    # Log-likelihood at posterior mean
    loglik_obs <- sum(dnorm(z,
                            mean = U1_mean %*% mean_beta[1, ],
                            sd = sqrt(1 / mean_phiep[1]),
                            log = TRUE
    ))

    loglik_models <- sum(sapply(2:M, function(k) {
        sum(dnorm(Y[, k - 1],
                  mean = UK_mean %*% mean_beta[k, ],
                  sd = sqrt(1 / mean_phiep[k]),
                  log = TRUE
        ))
    }))

    dev_post_mean <- -2 * (loglik_obs + loglik_models)
    p_D <- mean_deviance - dev_post_mean
    DIC <- mean_deviance + p_D

    return(list(
        DIC = DIC,
        p_D = p_D,
        mean_deviance = mean_deviance,
        deviance_posterior_mean = dev_post_mean
    ))
}
