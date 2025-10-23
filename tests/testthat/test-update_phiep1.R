test_that("update_phiep1 returns positive scalar and reacts to residual scale", {
    set.seed(123)
    H <- 2
    U_sub <- matrix(rnorm(10 * H), 10, H)
    beta  <- c(0.5, -0.2, 0.3)
    z     <- as.vector(cbind(1, U_sub) %*% beta + rnorm(10, sd = 0.5))

    val1 <- update_phiep1(z, U_sub, beta, a = 1, b = 1)
    expect_gt(val1, 0)
    expect_length(val1, 1)

    # Larger residuals -> larger rate -> typically higher precision draws
    z_noisier <- z + rnorm(10, sd = 2)
    val2 <- update_phiep1(z_noisier, U_sub, beta, a = 1, b = 1)
    expect_gt(val2, 0)
})
