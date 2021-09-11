test_that("sample_gen() generates simulation data", {
  n <- 10
  p <- 2
  mu <- c(0, 0)
  sigma <- matrix(c(1,0,0,1), nrow=2)
  labels <- c(1,2)
  hpara_func <- list(
    lambda_func = function(p) rnorm(p, 0, 1),
    sigma_func = function(p) rchisq(p, 2) + 1
  )
  res <- sample_gen(n, p, mu, sigma, labels, hpara_func)
  expect_identical(length(sample_gen(n, p, mu, sigma, labels, hpara_func)), length(res))
})
