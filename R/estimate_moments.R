#' @title Estimate first four moment parameters of multivariate observations
#'
#' @description Estimate first four moments of multivariate observations, namely,
#' mean vector, covariance matrix, coskewness matrix, and cokurtosis matrix.
#'
#' @author Rui Zhou and Daniel P. Palomar
#'
#' @references
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms,"
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <doi:10.1109/TSP.2021.3051369>.
#'
#' @param X Data matrix.
#' @param adjust_magnitude Boolean indicating whether to adjust the order of magnitude of parameters.
#'                         NOTE: It is specially designed for \code{\link{design_MVSKtilting_portfolio}} function.
#'
#' @return A list containing the following elements:
#' \item{\code{mu}}{Mean vector.}
#' \item{\code{Sgm}}{Covariance matrix.}
#' \item{\code{Phi_mat}}{Co-skewness matrix.}
#' \item{\code{Psi_mat}}{Co-kurtosis matrix.}
#' \item{\code{Phi}}{Co-skewness matrix in vector form (collecting only the unique elements).}
#' \item{\code{Psi}}{Co-kurtosis matrix in vector form (collecting only the unique elements).}
#' \item{\code{Phi_shred}}{Partition on \code{Phi} (see reference).}
#' \item{\code{Psi_shred}}{Partition on \code{Psi} (see reference).}
#'
#'
#' @examples
#'
#' library(ParametricHighOrderPortfolios)
#' data(X50)
#'
#' X_moments <- estimate_moments(X50[, 1:10])
#'
#'
#' @importFrom PerformanceAnalytics M3.MM
#' @importFrom PerformanceAnalytics M4.MM
#'
#' @importFrom stats cov
#' @export
estimate_moments <- function(X, adjust_magnitude = FALSE) {
  M3.mat2vec <- get("M3.mat2vec", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  M3.vec2mat <- get("M3.vec2mat", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  M4.vec2mat <- get("M4.vec2mat", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)


  N <- ncol(X)
  mu  <- colMeans(X)
  Sgm <- cov(X)
  Phi <- M3.MM(X, as.mat = FALSE)
  Psi <- M4.MM(X, as.mat = FALSE)

  if (adjust_magnitude) {
    d <- abs(eval_portfolio_moments(w = rep(1/N, N), X_moments = list(mu = mu, Sgm = Sgm, Phi = Phi, Psi = Psi)))
    mu  <- mu  / d[1]
    Sgm <- Sgm / d[2]
    Phi <- Phi / d[3]
    Psi <- Psi / d[4]
  }

  Phi_mat <- M3.vec2mat(Phi, N)
  Phi_shred <- lapply(1:N, function(i) Phi_mat[, (1:N)+N*(i-1)])

  Psi_mat <- M4.vec2mat(Psi, N)
  Psi_shred <- list()
  for (i in 1:N) {
    tmp <- Psi_mat[, (1:N^2)+N^2*(i-1)]
    Psi_shred[[i]] <- M3.mat2vec(tmp)
    gc()
  }
  gc()

  return(list(mu = mu, Sgm = Sgm, Phi_mat = Phi_mat, Psi_mat = Psi_mat, Phi = Phi, Psi = Psi, Phi_shred = Phi_shred, Psi_shred = Psi_shred))
}
