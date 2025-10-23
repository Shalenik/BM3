test_that("mvn returns a length-n numeric and respects mean/cov form", {
    set.seed(1)
    n <- 4
    A <- crossprod(matrix(rnorm(n^2), n, n)) + diag(1e-6, n)  # p.d.
    b <- rnorm(n)

    y <- mvn(A, b, n)

    expect_type(y, "double")
    expect_length(y, n)

    # Monte Carlo mean should be near A^{-1} b if we average many draws
    set.seed(2)
    B <- replicate(2000, mvn(A, b, n))
    mu_hat <- rowMeans(B)
    mu_true <- solve(A, b)
    expect_equal(mu_hat, mu_true, tolerance = 0.1)
})
