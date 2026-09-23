## tests for mvn_compare_means()/mvn_compare_cov() (R/validate.R)

test_that("mvn_compare_means gives the right RMSE, errors on shape mismatch", {
  out <- mvn_compare_means(true_mu = c(0, 0), est_mu = c(1, 1))
  expect_s3_class(out, "data.frame")
  expect_named(out, c("true", "est", "rmse"))
  expect_equal(out$rmse[1], 1) # sqrt(mean((0-1)^2, (0-1)^2)) = 1
  expect_equal(nrow(out), 2)
  expect_true(all(out$rmse == out$rmse[1])) # rmse repeated on every row

  expect_error(
    mvn_compare_means(c(0, 0, 0), c(1, 1)),
    "same number of elements"
  )

  ## matrices work too (flattened element-wise)
  m_true <- matrix(0, 2, 2)
  m_est <- matrix(c(1, 1, -1, -1), 2, 2)
  out2 <- mvn_compare_means(m_true, m_est)
  expect_equal(out2$rmse[1], 1)
})


test_that("mvn_compare_cov uses only the upper triangle; errors on mismatch", {
  true_Sigma <- diag(1, 3)
  est_Sigma <- diag(1, 3)
  est_Sigma[1, 2] <- est_Sigma[2, 1] <- 0.6 # off-diagonal difference
  ## lower-triangle-only differences should NOT affect the result
  ## upper.tri(..., diag=TRUE) is compared
  est_Sigma_lower_only <- diag(1, 3)
  ## asymmetric on purpose, lower triangle only
  est_Sigma_lower_only[2, 1] <- 0.6

  out <- mvn_compare_cov(true_Sigma, est_Sigma)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("true", "est", "rmse"))
  expect_equal(nrow(out), 6) # upper.tri(3x3, diag=TRUE) has 6 elements
  expect_true(out$rmse[1] > 0)

  out_lower_only <- mvn_compare_cov(true_Sigma, est_Sigma_lower_only)
  ## upper triangle identical, so RMSE should be 0
  expect_equal(out_lower_only$rmse[1], 0)

  expect_error(
    mvn_compare_cov(diag(1, 2), diag(1, 3)),
    "same dimensions"
  )
})
