#' @title Design MVSK portfolio via nloptr based on sample moments
#'
#' @description Design MVSK portfolio via nloptr based on sample moments
#'
#' @author Xiwen Wang and Daniel P. Palomar
#'
#' @param lambda Numerical vector of length 4 indicating the weights of first four moments.
#' @param X_moments List of moment parameters, see \code{\link{estimate_moments}}.
#' @param w0 Numerical vector indicating the initial value of portfolio weights.
#'
#' @return A list containing the following elements:
#' \item{\code{w_optimal}}{Optimal portfolio vector.}
#' \item{\code{obj_optimal}}{Optimal objective value.}
#' \item{\code{cpu_time}}{Computational time.}
#'
#' @import nloptr
#' @export
solve_MVSK_nloptr_non_parametric <- function(lambda, X_moments, w0){
  N <- length(w0)
  start_time <- proc.time()[3]
  obj_value <- function(w, lambda, X_moments) {
    obj <- 0
    obj <- obj - lambda[1] * as.numeric(t(w) %*% X_moments$mu)
    obj <- obj + lambda[2] * as.numeric(t(w) %*% X_moments$Sgm %*% w)
    obj <- obj - lambda[3] * as.numeric(t(w) %*% X_moments$Phi_mat %*% (w %x% w))
    obj <- obj + lambda[4] * as.numeric(t(w) %*% X_moments$Psi_mat %*% (w %x% w %x% w))
    return(obj)
  }

  nabla_f <- function(w, lambda, X_moments) {
    N <- length(w)
    grads <- rep(0, N)
    grads <- grads - lambda[1] * as.vector(X_moments$mu)
    grads <- grads + lambda[2] * 2 * as.vector(X_moments$Sgm %*% w)
    grads <- grads - lambda[3] * 3 * as.vector(X_moments$Phi_mat %*% (w %x% w))
    grads <- grads + lambda[4] * 4 * as.vector(X_moments$Psi_mat %*% (w %x% w %x% w))
    return(grads)
  }

  eval_g_eq <- function(w, lambda, X_moments) {
    constr <- c(sum(w)-1)
    grad <- c(rep(1, N))
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  local_opts <- list( "algorithm" = "NLOPT_LD_SLSQP",
                      "xtol_rel"=1.0e-14,"maxeval"=2500)

  res <- nloptr(x0 = c(w0),
                eval_f=obj_value,
                eval_grad_f=nabla_f,
                lb = c(rep(0,N)),
                ub = c(rep(Inf,N)),
                eval_g_eq = eval_g_eq,
                opts = list("algorithm" = "NLOPT_LD_SLSQP",    maxeval = 2500, "xtol_rel"=1.0e-14,
                            "local_opts" = local_opts ),
                lambda = lambda,
                X_moments = X_moments
  )
  w_optimal <- res$solution
  cpu_time <- as.numeric(proc.time()[3] - start_time)
  obj_optimal <- obj_value(w_optimal, lambda, X_moments)

  return(list("w_optimal"   = w_optimal,
              "obj_optimal" = obj_optimal,
              "cpu_time"    = cpu_time))
}
