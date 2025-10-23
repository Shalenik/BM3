# tests/testthat/test-dic_calculations.R
test_that("dic.calculations returns correct structure and numeric values", {
    # Minimal synthetic setup consistent with the function's expectations
    set.seed(123)

    # Dimensions
    T  <- 6L
    H  <- 2L
    J  <- 2L
    M  <- 3L  # 1 obs + 2 model series
    T1 <- 2L
    T2 <- 5L
    n_obs <- T2 - T1 + 1L

    # Posterior-mean objects we'll replicate across draws
    mean_U     <- cbind(seq(0.1, 0.6, length.out = T),
                        seq(-0.2, 0.3, length.out = T))       # T x H
    mean_beta  <- rbind(                                   # M x (H+1)
        c(0.5,  1.0, -0.3),  # k=1 obs
        c(0.2, -0.5,  0.8),  # k=2 model
        c(-0.1, 0.4,  0.2)   # k=3 model
    )
    mean_phiep <- c(4.0, 1.0, 2.0)  # precisions (=> variances 0.25, 1, 0.5)

    # Build design matrices at the posterior means
    U1_mean <- cbind(1, mean_U[T1:T2, , drop = FALSE])  # n_obs x (H+1)
    UK_mean <- cbind(1, mean_U)                         # T x (H+1)

    # Generate data from the same posterior means
    mu_obs <- as.vector(U1_mean %*% t(mean_beta[1L, , drop = FALSE]))  # n_obs x 1 -> vector
    z <- rnorm(n_obs, mean = mu_obs, sd = sqrt(1 / mean_phiep[1L]))

    B_models <- t(mean_beta[2:M, , drop = FALSE])                      # (H+1) x (M-1)
    MU <- UK_mean %*% B_models                                         # T x (M-1)
    Y <- matrix(NA_real_, nrow = T, ncol = M - 1L)
    for (j in 1:(M - 1L)) {
        Y[, j] <- rnorm(T, mean = MU[, j], sd = sqrt(1 / mean_phiep[j + 1L]))
    }

    # Arrays across draws (identical across iterations for deterministic test)
    nsim <- 4L
    gammas <- array(0, dim = c(nsim, J, H))                 # not used by function
    betas  <- array(rep(mean_beta, each = nsim),
                    dim = c(nsim, M, H + 1L))
    phieps <- matrix(rep(mean_phiep, each = nsim), nrow = nsim, ncol = M)
    Us     <- array(rep(mean_U, each = nsim),
                    dim = c(nsim, T, H))

    # Helper to compute deviance at posterior means (same closed-form as function)
    dev_post_mean_ref <- local({
        sig2_obs <- 1 / mean_phiep[1L]
        rss_obs  <- sum((z - mu_obs)^2)
        loglik_obs <- -0.5 * (length(z) * log(2 * pi * sig2_obs) + rss_obs / sig2_obs)

        sig2_mod <- 1 / mean_phiep[2:M]
        rss_term <- sum((Y - MU)^2 / rep(sig2_mod, each = T))
        const_term <- T * sum(log(2 * pi * sig2_mod))
        loglik_models <- -0.5 * (const_term + rss_term)

        -2 * (loglik_obs + loglik_models)
    })

    # Case 1: identical deviances across draws => p_D ~ 0 and DIC == mean deviance
    deviances <- rep(dev_post_mean_ref, nsim)

    res <- dic.calculations(
        X = matrix(1, T, J),
        Y = Y,
        z = z,
        gammas = gammas,
        betas = betas,
        phieps = phieps,
        Us = Us,
        deviances = deviances,
        T1 = T1,
        T2 = T2
    )

    expect_type(res, "list")
    expect_named(res, c("DIC", "p_D", "mean_deviance", "deviance_posterior_mean"))
    expect_true(all(sapply(res, is.numeric)))
    expect_equal(res$deviance_posterior_mean, dev_post_mean_ref, tolerance = 1e-8)
    expect_equal(res$mean_deviance, dev_post_mean_ref, tolerance = 1e-8)
    expect_equal(res$p_D, 0, tolerance = 1e-8)
    expect_equal(res$DIC, res$mean_deviance, tolerance = 1e-8)
})

test_that("dic.calculations uses back half (burn-in) and handles odd nsim", {
    set.seed(456)

    # Reuse a simple deterministic shape; only deviances matter for this test
    T <- 4L; H <- 2L; J <- 2L; M <- 2L; T1 <- 2L; T2 <- 3L
    mean_U <- cbind(
        seq(0.0, 0.3, length.out = T),
        seq(0.1, 0.4, length.out = T)
    )                             # T x H
    mean_beta <- rbind(
        c(0.1,  0.2, -0.1),  # k=1 obs: (intercept, load1, load2)
        c(-0.2, 0.4,  0.3)   # k=2 model
    )
    mean_phiep <- c(2, 3)        # precisions

    U1_mean <- cbind(1, mean_U[T1:T2, , drop = FALSE])  # (T2-T1+1) x 3
    UK_mean <- cbind(1, mean_U)                         # T x 3
    # obs vector (length 2)
    z <- as.vector(U1_mean %*% t(mean_beta[1, , drop = FALSE])) + c(0.01, -0.02)
    # model series: Y is T x (M-1) = T x 1
    Y <- matrix(as.vector(UK_mean %*% t(mean_beta[2, , drop = FALSE])) + 0.03, ncol = 1L)

    nsim <- 4L
    gammas <- array(0, dim = c(nsim, J, H))                       # nsim x 2 x 2
    betas  <- array(rep(mean_beta, each = nsim), dim = c(nsim, M, H + 1L))
    phieps <- matrix(rep(mean_phiep, each = nsim), ncol = M)
    Us     <- array(rep(mean_U, each = nsim), dim = c(nsim, T, H))

    # Compute deviance at posterior means once
    dev_post_mean_ref <- local({
        # obs
        sig2o <- 1 / mean_phiep[1]; mu_obs <- as.vector(U1_mean %*% mean_beta[1, ])
        lko <- -0.5 * (length(z) * log(2 * pi * sig2o) + sum((z - mu_obs)^2) / sig2o)
        # model
        mu_m <- as.vector(UK_mean %*% mean_beta[2, ]); sig2m <- 1 / mean_phiep[2]
        lkm <- -0.5 * (length(mu_m) * log(2 * pi * sig2m) + sum((Y[, 1] - mu_m)^2) / sig2m)
        -2 * (lko + lkm)
    })

    # Back-half check for even nsim=4: take deviances = c(100, 100, 200, 200)
    # mean_deviance should be 200 (average of last two draws)
    deviances <- c(100, 100, 200, 200)

    res <- dic.calculations(
        X = matrix(1, T, J),
        Y = Y, z = z,
        gammas = gammas, betas = betas, phieps = phieps, Us = Us,
        deviances = deviances, T1 = T1, T2 = T2
    )
    expected_mean_dev <- mean(deviances[as.integer((nsim / 2):nsim)])
    expect_equal(res$mean_deviance, expected_mean_dev)

    # Odd nsim=3 → floor(3/2)+1 = 2 .. 3 => mean of c(150, 300) = 225
    nsim2 <- 3L
    gammas2 <- gammas[1:nsim2, , , drop = FALSE]
    betas2  <- betas [1:nsim2, , , drop = FALSE]
    phieps2 <- phieps[1:nsim2, , drop = FALSE]
    Us2     <- Us    [1:nsim2, , , drop = FALSE]
    deviances2 <- c(50, 150, 300)

    res2 <- dic.calculations(
        X = matrix(1, T, J),
        Y = Y, z = z,
        gammas = gammas2, betas = betas2, phieps = phieps2, Us = Us2,
        deviances = deviances2, T1 = T1, T2 = T2
    )
    expected_mean_dev2 <- mean(deviances2[as.integer((nsim2 / 2):nsim2)])
    expect_equal(res2$mean_deviance, expected_mean_dev2)
    expect_true(is.finite(res2$DIC))
})

test_that("dic.calculations errors on dimension mismatch between Y and phieps", {
    # M must equal 1 + ncol(Y)
    T <- 3L; H <- 2L; J <- 2L
    Y <- matrix(0, nrow = T, ncol = 2L)        # M-1 = 2 => M should be 3
    z <- c(0, 0)
    T1 <- 1L; T2 <- 2L
    nsim <- 2L

    # phieps with wrong M (here M=4)
    phieps <- matrix(1, nrow = nsim, ncol = 4L)
    gammas <- array(0, dim = c(nsim, J, H))
    betas  <- array(0, dim = c(nsim, 4L, H + 1L))
    Us     <- array(0, dim = c(nsim, T, H))
    deviances <- c(10, 12)

    expect_error(
        dic.calculations(
            X = matrix(1, T, J), Y = Y, z = z,
            gammas = gammas, betas = betas, phieps = phieps, Us = Us,
            deviances = deviances, T1 = T1, T2 = T2
        )
    )
})
