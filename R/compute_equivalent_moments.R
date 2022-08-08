#' @title Computing equivalent sample moments given the parameters of a ghMST distribution
#'
#' @description Computing equivalent sample moments given the parameters of a ghMST distribution
#'
#' @author Xiwen Wang, Daniel P. Palomar
#'
#' @references
#' X. Wang, R. Zhou, J. Ying, and D. P. Palomar, "Efficient and Scalable High-Order Portfolios Design via Parametric Skew-t Distribution,"
#' Available in arXiv, 2022. <https://arxiv.org/pdf/2206.02412.pdf>.
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms,"
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <doi:10.1109/TSP.2021.3051369>.
#' @param parameters List of fitted parameters, including location vector, skewness vector, scatter matrix, and the degree of freedom..
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
#' @examples
#' library(ParametricHighOrderPortfolios)
#' data(X50)
#'
#' parameters <- fit_ghMST(X50)
#' # decide moment weights
#' equivalent_sample_moments <- Compute_equivalent_moments(parameters)
#'
#' @import Rcpp
#' @useDynLib ParametricHighOrderPortfolios
#' @export
Compute_equivalent_moments <- function(parameters) {
  N <- length(parameters$mu)
  equivalent_moments <- list()
  equivalent_moments$mu <- as.vector(parameters$mu + parameters$a$a11 * parameters$gamma)
  equivalent_moments$Sgm <- parameters$a$a21 * parameters$scatter + parameters$a$a22 * parameters$gamma %*% t(parameters$gamma)

  equivalent_moments$Phi_mat <- compute_Phi_hat(N, parameters)

  equivalent_moments$Psi_mat <- compute_Psi_hat(N, parameters)

  M3.mat2vec <- get("M3.mat2vec", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)
  M4.mat2vec <- get("M4.mat2vec", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)

  equivalent_moments$Phi <- M3.mat2vec(equivalent_moments$Phi_mat)
  equivalent_moments$Psi <- M4.mat2vec(equivalent_moments$Psi_mat)

  equivalent_moments$Phi_shred <- lapply(1:N, function(i) equivalent_moments$Phi_mat[, (1:N)+N*(i-1)])

  equivalent_moments$Psi_shred <- list()
  for (i in 1:N) {
    tmp <- equivalent_moments$Psi_mat[, (1:N^2)+N^2*(i-1)]
    equivalent_moments$Psi_shred[[i]] <- M3.mat2vec(tmp)
  }

  return(equivalent_moments)
}
