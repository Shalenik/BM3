test_that("update_phiepk returns positive scalar", {
    set.seed(456)
    T <- 12; H <- 2
    U <- matrix(rnorm(T * H), T, H)
    beta <- c(0.2, 0.3, -0.4)
    Yk <- as.vector(cbind(1, U) %*% beta + rnorm(T, sd = 1.1))

    val <- update_phiepk(Yk, U, beta, a = 1, b = 1)
    expect_gt(val, 0)
    expect_length(val, 1)
})
