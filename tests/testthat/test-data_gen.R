test_that("data_gen() generates simulation data", {
  n <- 10
  p <- 2
  centers <- 2
  mu <- c(0, 0)
  omega <- matrix(c(1, 0, 0, 1), nrow = 2)
  labels <- c(1, 2)
  size <- 5
  hparam_func <- list(
    lambda_func = function(p) stats::rnorm(p, 0, 1),
    sigma_func = function(p) stats::rchisq(p, 2) + 1
  )
  
  res <- data_gen(n, p, centers, mu, omega, labels, size, hparam_func)
  expect_identical(length(data_gen(n, p, centers, mu, omega, labels, size, hparam_func)), length(res))
})
