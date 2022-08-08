#' @title Design MVSK portfolio via nloptr based on the parameters of generalized hyperbolic skew-t distribution
#'
#' @description Design MVSK portfolio via nloptr based on the parameters of generalized hyperbolic skew-t distribution:
#'
#' @author Xiwen Wang and Daniel P. Palomar
#'
#' @param lambda Numerical vector of length 4 indicating the weights of first four moments.
#' @param parameters List of fitted parameters, including location vector, skewness vector, scatter matrix, and the degree of freedom..
#' @param w0 Numerical vector indicating the initial value of portfolio weights.
#'
#' @return A list containing the following elements:
#' \item{\code{w_optimal}}{Optimal portfolio vector.}
#' \item{\code{obj_optimal}}{Optimal objective value.}
#' \item{\code{cpu_time}}{Computational time.}
#'
#' @import nloptr
#' @export
solve_MVSK_nloptr_skewt <- function(w0, lambda, parameters) {
  start_time <- proc.time()[3]
  lambda_max <- max(lambda)
  lambda <- lambda/lambda_max
  N <- length(w0)
  eval_g_eq <- function(w, lambda, parameters) {
    constr <- c(sum(w)-1)
    grad <- c(rep(1, N))
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  obj_value_skewt <- function(w, lambda, parameters) {
    obj <- 0
    wT_mu <- t(w) %*% parameters$mu
    wT_gamma <- t(w) %*% parameters$gamma
    wT_Sigma_w <- sum((parameters$chol_Sigma %*% w) ** 2)
    obj <- obj - lambda[1] * (wT_mu + parameters$a$a11 * wT_gamma)
    obj <- obj + lambda[2] * (parameters$a$a21 * wT_Sigma_w  + parameters$a$a22 * ((wT_gamma)**2) )
    obj <- obj - lambda[3] * (parameters$a$a31 * ((wT_gamma) ** 3) + parameters$a$a32 * (wT_gamma) * ( wT_Sigma_w) )
    obj <- obj + lambda[4] * (parameters$a$a41 * ((wT_gamma) ** 4) + parameters$a$a42 * (wT_Sigma_w) * ((wT_gamma) ** 2) + parameters$a$a43 * ((wT_Sigma_w) ** 2))
    return(as.numeric(obj))
  }

  nabla_f_skewt <- function(w, lambda, parameters) {
    wT_gamma <- t(w) %*% parameters$gamma
    wT_Sigma_w <- sum((parameters$chol_Sigma %*% w) ** 2)
    nf <- - lambda[1] * (parameters$mu + parameters$a$a11 * parameters$gamma) + lambda[2] * (2 * parameters$a$a21 * parameters$scatter %*% w + parameters$a$a22 * 2 * as.numeric(wT_gamma) * as.matrix(parameters$gamma) )
    nf <- nf - lambda[3] * (parameters$a$a31 * 3 * (as.numeric(wT_gamma) ** 2) * as.matrix(parameters$gamma) + parameters$a$a32 * (as.numeric(t(w) %*% (parameters$scatter) %*% w) * as.matrix(parameters$gamma) + 2 * as.numeric(wT_gamma) * parameters$scatter %*% w ))
    nf <- nf + lambda[4] * (4 * parameters$a$a41 * as.numeric((wT_gamma) ** 3) *as.matrix(parameters$gamma) + parameters$a$a42 * (2 * as.numeric((wT_gamma) ** 2) * parameters$scatter %*% w + 2 * as.numeric(wT_Sigma_w) * as.numeric((wT_gamma)) *as.matrix(parameters$gamma)) + 4 * parameters$a$a43 * (as.numeric(wT_Sigma_w)) * (parameters$scatter) %*% w   )
    return(nf)
  }

  local_opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                      "xtol_rel"=1.0e-14,"maxeval"=2500)

  res <- nloptr(x0 = c(w0),
                eval_f=obj_value_skewt,
                eval_grad_f=nabla_f_skewt,
                lb = c(rep(0,N)),
                ub = c(rep(Inf,N)),
                eval_g_eq = eval_g_eq,
                opts = list("algorithm" = "NLOPT_LD_SLSQP",
                            maxeval = 2500, "xtol_rel"=1.0e-14,
                            "local_opts" = local_opts ),
                lambda = lambda,
                parameters = parameters
  )
  w_optimal <- res$solution
  cpu_time <- as.numeric(proc.time()[3] - start_time)
  obj_optimal <- obj_value_skewt(w_optimal, lambda, parameters)

  return(list("w_optimal"   = w_optimal,
              "obj_optimal" = obj_optimal * lambda_max,
              "cpu_time"    = cpu_time))
}
