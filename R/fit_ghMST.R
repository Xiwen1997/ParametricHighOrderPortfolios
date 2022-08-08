#' @title Estimate the parameters of ghMST distribution from multivariate observations
#'
#' @description Using the package fitHeavyTail to estimate the parameters of ghMST distribution from multivariate observations, namely,
#' location vector (mu), skewness vector (gamma), scatter matrix (scatter), degree of freedom (nu), parameters a,
#' and the Cholesky decomposition of the scatter matrix (chol_Sigma).
#'
#' @author Xiwen Wang, Rui Zhou and Daniel P. Palomar
#'
#' @references
#' Aas, Kjersti and Ingrid Hobæk Haff. "The generalized hyperbolic skew student’st-distribution,"
#' Journal of financial econometrics, pp. 275-309, 2006.
#'
#' @param X Data matrix.
#' @param X Data matrix containing the multivariate time series (each column is one time series).
#' @param initial List of initial values of the parameters for the iterative estimation method.
#'                Possible elements include:
#'                \itemize{\item{\code{nu}: default is \code{4},}
#'                         \item{\code{mu}: default is the data sample mean,}
#'                         \item{\code{gamma}: default is the sample skewness vector,}
#'                         \item{\code{scatter}: default follows from the scaled sample covariance matrix,}
#'                         }
#' @param nu_lb Minimum value for the degree of freedom to maintain the existence of high-order moments (default is \code{9}).
#' @param max_iter Integer indicating the maximum number of iterations for the iterative estimation
#'                 method (default is \code{100}).
#' @param ptol Positive number indicating the relative tolerance for the change of the variables
#'             to determine convergence of the iterative method (default is \code{1e-3}).
#' @param ftol Positive number indicating the relative tolerance for the change of the log-likelihood
#'             value to determine convergence of the iterative method (default is \code{Inf}, so it is
#'             not active). Note that using this argument might have a computational cost as a convergence
#'             criterion due to the computation of the log-likelihood (especially when \code{X} is high-dimensional).
#' @param PXEM Logical value indicating whether to use the parameter expansion (PX) EM method to accelerating the convergence.
#' @param return_iterates Logical value indicating whether to record the values of the parameters (and possibly the
#'                        log-likelihood if \code{ftol < Inf}) at each iteration (default is \code{FALSE}).
#' @param verbose Logical value indicating whether to allow the function to print messages (default is \code{FALSE}).
#'
#' @return A list containing the following elements:
#'         \item{\code{mu}}{Location vector estimate (not the mean).}
#'         \item{\code{gamma}}{Skewness vector estimate.}
#'         \item{\code{scatter}}{Scatter matrix estimate.}
#'         \item{\code{nu}}{Degrees of freedom estimate.}
#'         \item{\code{chol_Sigma}}{Choleski decomposition of the Scatter matrix estimate.}
#'         \item{\code{a}}{A list of coefficients useful for later computation}
#'
#' @examples
#' library(ParametricHighOrderPortfolios)
#' data("X50")
#' ghMST_parameters <- fit_ghMST(X50)
#'
#' @import fitHeavyTail
#' @export
fit_ghMST <- function(X, initial = NULL, nu_lb = 9, max_iter = 100, ptol = 1e-3, ftol = Inf,
                     PXEM = TRUE, return_iterates = FALSE, verbose = FALSE) {

  options(nu_min = nu_lb)
  parameters <- fitHeavyTail::fit_mvst(X, max_iter = max_iter, ptol = ptol, ftol = ftol,
                                       PXEM = PXEM, return_iterates = return_iterates, verbose = verbose)
  return_paramters <- list()
  return_paramters$mu    <- parameters$mu
  return_paramters$nu    <- parameters$nu
  return_paramters$gamma <- parameters$gamma
  return_paramters$scatter <- parameters$scatter
  return_paramters$chol_Sigma <- chol(return_paramters$scatter)

  ## Compute a given paramters of skew-t model
  compute_a <- function(parameters) {
    a11 <- (parameters$nu)/(parameters$nu-2)
    a21 <- (parameters$nu)/(parameters$nu-2)
    a22 <- (2 * parameters$nu ** 2)/((parameters$nu - 2) ** 2 * (parameters$nu - 4))
    a31 <- 16 * (parameters$nu ** 3)/(((parameters$nu - 2) ** 3) * (parameters$nu - 4) * (parameters$nu - 6))
    a32 <- 6 * (parameters$nu ** 2)/(((parameters$nu -2) ** 2) * (parameters$nu - 4))
    a41 <- (12 * parameters$nu + 120) * (parameters$nu ** 4)/ (((parameters$nu - 2) ** 4) * (parameters$nu - 4) * (parameters$nu - 6) * (parameters$nu - 8))
    a42 <- (2 * parameters$nu + 4) * (parameters$nu ** 3)/(((parameters$nu - 2) ** 3) * (parameters$nu - 4) * (parameters$nu - 6)) * 6
    a43 <- (parameters$nu ** 2)/((parameters$nu - 2) * (parameters$nu - 4)) * 3
    return(list("a11" = a11, "a21" = a21, "a22" = a22,
                "a31" = a31, "a32" = a32,
                "a41" = a41, "a42" = a42, "a43" = a43))
  }

  return_paramters$a <- compute_a(return_paramters)
  return(return_paramters)
}
