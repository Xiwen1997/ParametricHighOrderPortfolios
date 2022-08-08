test_that("Checking calculating moments from the parameters of a ghMST distribution", {
  X_parameters <- fit_ghMST(X50, nu_lb = 8.5 + runif(1))
  X_moments <- Compute_equivalent_moments(X_parameters)

  w <-  runif(50)
  w <- w/sum(w)

  phi1_skewt <- as.numeric(t(w) %*% as.vector(X_parameters$mu + X_parameters$a$a11 * X_parameters$gamma))
  phi1_np <- as.numeric(t(w) %*% X_moments$mu)
  expect_equal(phi1_skewt, phi1_np)

  phi2_skewt <- as.numeric(X_parameters$a$a21 * t(w) %*% X_parameters$scatter %*% w + X_parameters$a$a22 * ((t(w) %*% X_parameters$gamma) ** 2))
  phi2_np <- as.numeric(t(w) %*% X_moments$Sgm %*% w)
  expect_equal(phi2_skewt, phi2_np)

  phi3_skewt <- as.numeric(X_parameters$a$a31 * ((t(w) %*% X_parameters$gamma) ** 3) + X_parameters$a$a32 * (t(w) %*% X_parameters$gamma) * ( t(w) %*% (X_parameters$scatter) %*% w) )
  phi3_np <- as.numeric(t(w) %*% X_moments$Phi_mat %*% (w %x% w))
  expect_equal(phi3_skewt, phi3_np)

  phi4_skewt <- as.numeric(X_parameters$a$a41 * ((t(w) %*% X_parameters$gamma) ** 4) + X_parameters$a$a42 * (t(w) %*% (X_parameters$scatter) %*% w) * ((t(w) %*% X_parameters$gamma) ** 2) + X_parameters$a$a43 * ((t(w) %*% (X_parameters$scatter) %*% w) ** 2))
  phi4_np <- as.numeric(t(w) %*% X_moments$Psi_mat %*% (w %x% w %x% w))
  expect_equal(phi4_skewt, phi4_np)
})
