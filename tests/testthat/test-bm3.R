# tests/testthat/test-bm3.R
test_that("BM3 runs and returns objects with expected dimensions (smoke test)", {
    skip_on_cran()
    set.seed(42)

    T <- 30; H <- 2; J <- 6; m <- 18
    T1 <- 5; T2 <- 25; nsim <- 20

    Y <- matrix(stats::rnorm(T * m), T, m)
    z <- stats::rnorm(T2 - T1 + 1)

    # Use defaults (mga=mbe=aep=aet=bep=bet=0.1, a1=a2=3)
    res <- BM3(Y, z, T, T1, T2, H, iseed = 99, J, nsim)

    expect_true(is.matrix(res$X))
    expect_equal(dim(res$U_samples),     c(nsim, T, H))
    expect_equal(dim(res$beta_samples),  c(nsim, m + 1, H + 1))
    expect_equal(dim(res$gamma_samples), c(nsim, J, H))
    expect_equal(dim(res$phiep_samples), c(nsim, m + 1))
    expect_length(res$deviances, nsim)
    expect_true(is.list(res$DIC))
})
