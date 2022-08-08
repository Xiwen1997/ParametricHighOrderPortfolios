#' @title Solve high-order portfolio without shorting based on non-parametric representation of the high-order moments
#'
#' @description Solve high-order portfolio without shorting based on non-parametric representation of the high-order moments
#' (i.e., mean, variance, skewness, and kurtosis):
#' \preformatted{
#'   minimize     - lmd1*(w'*mu) + lmd2*(w'*Sigma*w)
#'                - lmd3*(w'*Phi*w*w) + lmd4*(w'*Psi*w*w*w)
#'   subject to   w>=0, sum(w) == 1,
#' }
#'
#' @author Xiwen, Rui Zhou and Daniel P. Palomar
#'
#' @references
#' R. Zhou and D. P. Palomar, "Solving High-Order Portfolios via Successive Convex Approximation Algorithms,"
#' in \emph{IEEE Transactions on Signal Processing}, vol. 69, pp. 892-904, 2021.
#' <doi:10.1109/TSP.2021.3051369>.
#'
#' X. Wang, R. Zhou, J. Ying, and D. P. Palomar, "Efficient and Scalable High-Order Portfolios Design via Parametric Skew-t Distribution,"
#' Available in arXiv, 2022. <https://arxiv.org/pdf/2206.02412.pdf>.
#'
#' @param lmd Numerical vector of length 4 indicating the weights of first four moments.
#' @param X_moments List of moment parameters, see \code{\link{estimate_moments}}.
#' @param w_init Numerical vector indicating the initial value of portfolio weights.
#' @param method String indicating the algorithm method, must be one of: "Q-MVSK", "L-MVSK", "DC", "PGD", "RFPA", "SQUAREM".
#' @param tau_w Number (>= 0) guaranteeing the strong convexity of approximating function.
#' @param initial_eta initial eta for the PGD and RFPA methods
#' @param beta scaling factor for the PGD and RFPA methods
#' @param tau hyper-parameter tau for the PGD and SQUAREM methods
#' @param gamma Number (0 < gamma <= 1) indicating the initial value of gamma for the Q-MVSK method.
#' @param zeta Number (0 < zeta < 1) indicating the diminishing parameter of gamma for the Q-MVSK method.
#' @param maxiter Positive integer setting the maximum iteration.
#' @param ftol Positive number setting the convergence criterion of function objective.
#' @param wtol Positive number setting the convergence criterion of portfolio weights.
#' @param stopval Number setting the stop value of objective.
#'
#' @return A list containing the following elements:
#' \item{\code{w}}{Optimal portfolio vector.}
#' \item{\code{cpu_time_vs_iterations}}{Time usage over iterations.}
#' \item{\code{objfun_vs_iterations}}{Objective function over iterations.}
#' \item{\code{iterations}}{Iterations index.}
#' \item{\code{convergence}}{Boolean flag to indicate whether or not the optimization converged.}
#' \item{\code{moments}}{Moments of portfolio return at optimal portfolio weights.}
#'
#' @examples
#'
#' library(ParametricHighOrderPortfolios)
#' data(X50)
#'
#' # estimate moments
#' X_moments <- estimate_moments(X50[, 1:10])
#'
#' # decide moment weights
#' lmd <- c(1,4,10,20)
#'
#' # portfolio optimization
#' sol <- solve_MVSK_non_parametric(lmd, X_moments, method = "Q-MVSK")
#'
#'
#' @importFrom utils tail
#' @import PerformanceAnalytics
#' @export
solve_MVSK_non_parametric <- function(lmd, X_moments,
                                  w_init = rep(1/length(X_moments$mu), length(X_moments$mu)),
                                  method = c("Q-MVSK", "L-MVSK", "DC", "PGD", "RFPA", "SQUAREM"),
                                  initial_eta = 10, beta = 0.5, tau = 1e6,
                                  tau_w = 0, gamma = 1, zeta = 1e-8, maxiter = 1e3, ftol = 1e-6, wtol = 1e-6, stopval = -Inf) {
  method <- match.arg(method)
  derportm3 <- get("derportm3", envir = asNamespace("PerformanceAnalytics"), inherits = FALSE)

  # prep
  N <- length(X_moments$mu)
  Amat <- t(rbind(matrix(1, 1, N), diag(N)))
  bvec <- cbind(c(1, rep(0, N)))

  fun_eval <- function(w) {
    return_list <- list()

    if (method == "Q-MVSK") {
      return_list$H3 <- 6 * sapply(X_moments$Phi_shred, function(x) x%*%w)
      return_list$H4 <- 4 * sapply(X_moments$Psi_shred, function(x) derportm3(w, x))
      return_list$H34 <- - lmd[3] * return_list$H3 + lmd[4] * return_list$H4
      return_list$jac <- rbind("grad1" = X_moments$mu, "grad2" = 2 * c(X_moments$Sgm %*% w), "grad3" = (1/2) * c(return_list$H3 %*% w), "grad4" = (1/3) * c(return_list$H4 %*% w))

    } else {
      w_kron_w <- kronecker(w, w)
      return_list$jac <- rbind("grad1" = X_moments$mu, "grad2" = 2 * c(X_moments$Sgm %*% w), "grad3" = 3 * c(X_moments$Phi_mat %*% w_kron_w), "grad4" = 4 * c(X_moments$Psi_mat %*% kronecker(w_kron_w, w)))
    }
    return_list$obj <- sum(lmd * as.vector(return_list$jac %*% w) / c(-1, 2, -3, 4))

    return(return_list)
  }

  # initialization
  start_time <- proc.time()[3]

  if (method == "L-MVSK") {
    rho <- lmd[3]*.maxEigHsnS(S = X_moments$Phi, N = N, func = "max") + lmd[4]*.maxEigHsnK(K = X_moments$Psi, N = N, func = "max")
  }
  if (method == "DC") {
    rho <- 2 * lmd[2] * norm(X_moments$Sgm, "I") + lmd[3]*.maxEigHsnS(S = X_moments$Phi, N = N, func = "sum") + lmd[4]*.maxEigHsnK(K = X_moments$Psi, N = N, func = "sum")
  }

  if (method == "PGD" || method == "RFPA" || method == "SQUAREM") {
    # define the PGD update function
    PGD_update <- function(w, eta, g) {
      r <- w - eta * g
      if(sum(r-min(r)) <= 1) {
        gamma_tilde <- (1/N) * (-1+sum(r))
        w <- pmax(0, r - gamma_tilde)
      } else {
        # otherwise we need to find out the position
        r_vec <- sort(r)
        inv_cum_r_vec <- rev(cumsum(rev(r_vec)))

        n_down <- 0
        n_up <- N
        n <- round((n_down + n_up)/2)

        stop_sign <- FALSE
        while(n_up - n_down > 1) {
          fn <- inv_cum_r_vec[n] - (N-n+1) * r_vec[n]
          if(fn < 1) {
            n_up <- n
          } else if(fn > 1) {
            n_down <- n
          }
          n <- round((n_down + n_up)/2)
        }
        gamma_tilde <- (sum(inv_cum_r_vec[n_up]) -1)/(N-n_up+1)
        w <- pmax(0, r - gamma_tilde)
      }
      return(w)
    }

    # get the gradients from jacobians
    get_gradient <- function(jac) {
      return(-(lmd[1]*jac[1, ] - lmd[2]*jac[2, ] + lmd[3]*jac[3, ] - lmd[4]*jac[4, ]))
    }
  }



  wk <- w_init
  cpu_time <- c(0)
  objs  <- c()
  fun_k <- fun_eval(wk)
  objs <- c(objs, fun_k$obj)

  for (iter in 1:maxiter) {
    # record previous w
    w_old <- wk

    ## construct QP approximation problem (the symbol and scale is adjusted to match the format of solver quadprog::solve.QP)
    switch(method,
           "Q-MVSK" = {
             H_ncvx <- .apprxHessian(fun_k$H34)
             Qk <- 2*lmd[2]*X_moments$Sgm + H_ncvx + diag(tau_w, N)
             qk <- lmd[1]*X_moments$mu + lmd[3]*fun_k$jac[3, ] - lmd[4]*fun_k$jac[4, ] + H_ncvx%*%wk + tau_w*wk

             # solve the QP problem
             w_hat <- quadprog::solve.QP(Dmat = Qk, dvec = cbind(qk), Amat = Amat, bvec = bvec, meq = 1)$solution

             # update w
             wk <- wk + gamma * (w_hat - wk)
             gamma <- gamma * (1 - zeta * gamma)
           },
           "L-MVSK" = {
             Qk <- 2*lmd[2]*X_moments$Sgm + diag(rho, N)
             qk <- lmd[1]*X_moments$mu + lmd[3]*fun_k$jac[3, ] - lmd[4]*fun_k$jac[4, ] + rho*wk
             wk <- quadprog::solve.QP(Dmat = Qk, dvec = cbind(qk), Amat = Amat, bvec = bvec, meq = 1)$solution
           },
           "DC" = {
             Qk <- diag(rho, N)
             qk <- rho*wk + lmd[1]*fun_k$jac[1, ] - lmd[2]*fun_k$jac[2, ] + lmd[3]*fun_k$jac[3, ] - lmd[4]*fun_k$jac[4, ]
             wk <- quadprog::solve.QP(Dmat = Qk, dvec = cbind(qk), Amat = Amat, bvec = bvec, meq = 1)$solution
           },
           "PGD" ={
             current_obj <- objs[length(objs)]
             # compute the gradient
             eta <- initial_eta
             gk <- get_gradient(fun_k$jac)

             wk_next <- PGD_update(wk, eta, gk)
             fun_k_next <- fun_eval(wk_next)
             next_obj <- fun_k_next$obj

             # backtracking line search
             while(next_obj > current_obj + t(gk) %*% (wk_next - wk) + ((1/(2*eta)) * sum((wk-wk_next)**2)) ){
               eta <- eta * beta
               wk_next <- PGD_update(wk, eta, gk)
               fun_k_next <- fun_eval(wk_next)
               next_obj <- fun_k_next$obj
             }
             wk <- wk_next
           },
           "RFPA" = {
             current_obj <- objs[length(objs)]

             # Try SQUAREM acceleration
             gk <- get_gradient(fun_k$jac)

             # One iterate of update
             wk1 <- PGD_update(wk, 1/tau, gk)
             fun_k1 <- fun_eval(wk1)

             # Another iterate of update
             wk2 <- PGD_update(wk1, 1/tau, get_gradient(fun_k1$jac) )

             r <- wk1 - wk
             v <- (wk2 - wk1) - r
             if(max(abs(v)) == 0) {
               alpha <- -1
             } else {
               alpha <- - (norm(r, "2"))/(norm(v, "2"))
             }
             if(as.numeric(t(r) %*% v) < 0)
             {
               alpha <- max(alpha, (as.numeric(t(r) %*% r))/(as.numeric(t(r) %*% v)))
             }


             wkt <- wk - 2 * alpha * r + alpha * alpha * v
             fun_kt <- fun_eval(wkt)

             wk_next <- PGD_update(wkt, 1/tau, rep(0, N) )
             fun_k_next <- fun_eval(wk_next)

             next_obj <- fun_k_next$obj

             # If we need PGD to leave the current point
             if(next_obj > current_obj) {
               eta <- initial_eta

               wk_next <- PGD_update(wk, eta, gk)
               fun_k_next <- fun_eval(wk_next)
               next_obj <- fun_k_next$obj

               # backtracking line search
               while(next_obj > current_obj + t(gk) %*% (wk_next - wk) + ((1/(2*eta)) * sum((wk-wk_next)**2)) ){
                 eta <- eta * beta
                 wk_next <- PGD_update(wk, eta, gk)
                 fun_k_next <- fun_eval(wk_next)
                 next_obj <- fun_k_next$obj
               }
             }
             wk <- wk_next
           },
           "SQUAREM" = {
             current_obj <- objs[length(objs)]

             # Try SQUAREM acceleration
             gk <- get_gradient(fun_k$jac)

             # One iterate of update
             wk1 <- PGD_update(wk, 1/tau, gk)
             fun_k1 <- fun_eval(wk1)

             # Another iterate of update
             wk2 <- PGD_update(wk1, 1/tau, get_gradient(fun_k1$jac) )

             r <- wk1 - wk
             v <- (wk2 - wk1) - r
             if(max(abs(v)) == 0) {
               alpha <- -1
             } else {
               alpha <- min(-1, - (norm(r, "2"))/(norm(v, "2")))
             }

             wkt <- wk - 2 * alpha * r + alpha * alpha * v
             fun_kt <- fun_eval(wkt)

             wk_next <- PGD_update(wkt, 1/tau, get_gradient(fun_kt$jac) )
             fun_k_next <- fun_eval(wk_next)

             next_obj <- fun_k_next$obj

             # backtracking line search on alpha
             while(next_obj > current_obj) {
               alpha <- 0.5 * (alpha - 1)
               wkt <- wk - 2 * alpha * r + alpha * alpha * v
               fun_kt <- fun_eval(wkt)

               wk_next <- PGD_update(wkt, 1/tau, get_gradient(fun_kt$jac) )
               fun_k_next <- fun_eval(wk_next)

               next_obj <- fun_k_next$obj
             }
             wk <- wk_next
           },
           stop("Method unknown")
    )

    wk[which(wk < 0)] <- 0
    wk <- wk/sum(wk)

    # recording...
    cpu_time <- c(cpu_time, proc.time()[3] - start_time)

    if(method == "PGD" || method == "RFPA") {
      fun_k <- fun_k_next
    } else {
      fun_k <- fun_eval(wk)
    }

    objs <- c(objs, fun_k$obj)

    # termination criterion
    # has_w_converged <- all(abs(wk - w_old) <= .5 * wtol )
    # has_w_converged <- norm(wk - w_old, "2") <= wtol * norm(w_old, "2")
    has_w_converged <- all(abs(wk - w_old) <= .5 * wtol * (abs(wk) + abs(w_old)))
    has_f_converged <- abs(diff(tail(objs, 2))) <= .5 * ftol * sum(abs(tail(objs, 2)))
    has_cross_stopval <- tail(objs, 1) <= stopval

    if (has_w_converged || has_f_converged || has_cross_stopval) break
  }

  return(list(
    "w"                      = wk,
    "cpu_time_vs_iterations" = cpu_time,
    "objfun_vs_iterations"   = objs,
    "iterations"             = 0:iter,
    "convergence"            = !(iter == maxiter),
    "moments"                = as.vector(fun_k$jac %*% wk) / c(1, 2, 3, 4)
  ))

}

