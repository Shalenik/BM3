test_that("update_del returns H-length vector with positive entries", {
    set.seed(789)
    T <- 10; J <- 3; H <- 2
    X <- matrix(rnorm(T * J), T, J)
    gamma <- matrix(rnorm(J * H), J, H)
    U <- matrix(rnorm(T * H), T, H)
    phiet <- runif(H, 0.5, 2)
    del <- rep(1, H)

    out <- update_del(U, X, gamma, phiet, del, a1 = 2, a2 = 3)
    expect_type(out, "double")
    expect_length(out, H)
    expect_true(all(out > 0))
})
