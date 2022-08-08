test_that("Checking high-order design via non-parametric approach", {
  X <- X50[, 1:10]
  w0 <- rep(1/10, 10)
  X_moments <- estimate_moments(X)

  lambda <- c(1, 4, 10, 20)
  load("MVSK_non_parametric.RData")

  L_MVSK_result <- solve_MVSK_non_parametric(lambda, X_moments, method = "L-MVSK", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(L_MVSK_result$w, w_check, tolerance = 1e-4)

  Q_MVSK_result <- solve_MVSK_non_parametric(lambda, X_moments, method = "Q-MVSK", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(Q_MVSK_result$w, w_check, tolerance = 1e-4)

  DC_result <- solve_MVSK_non_parametric(lambda, X_moments, method = "DC", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(DC_result$w, w_check, tolerance = 1e-4)

  PGD_result <- solve_MVSK_non_parametric(lambda, X_moments, method = "PGD", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(PGD_result$w, w_check, tolerance = 1e-4)

  RFPA_result <- solve_MVSK_non_parametric(lambda, X_moments, method = "RFPA", maxiter = 10000, ftol = 1e-10, wtol = 1e-10, tau = 100)
  expect_equal(RFPA_result$w, w_check, tolerance = 1e-4)
})
