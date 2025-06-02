#' cmgnd: Function for Clustering using Constrained Mixtures of Generalized Normal Distributions
#'
#' @description
#' Fits univariate constrained mixture of generalized normal distribution models
#' by imposing mixture partitions. Models are estimated by the ECM algorithm initialized by k-means.
#'
#' @param x A numeric vector of observations.
#' @param K An integer specifying the number of mixture components to fit.
#' Default is 2.
#' @param Cmu A binary vector indicating mixture components
#' for the location parameter. The k-th element is set to 1 if the
#' k-th mixture component belongs to the Cr partition, and 0 otherwise.
#' Default is \code{c(0,0)}, indicating no mixture partition with \code{K=2}.
#' @param Csigma A binary vector indicating mixture components
#' for the scale parameter. The k-th element is set to 1 if the k-th
#' mixture component belongs to the Cr partition, and 0 otherwise.
#' Default is \code{c(0,0)}, indicating no mixture partition with \code{K=2}.
#' @param Cnu A binary vector indicating mixture components
#' for the shape parameter. The k-th element is set to 1 if the k-th mixture component belongs to the Cr partition, and 0 otherwise.
#' Default is \code{c(0,0)}, indicating no mixture partition with \code{K=2}.
#' @param nstart An integer specifying the number of starting
#' points for the shape parameter. Default is 10.
#' @param theta A parameter matrix used to initialize the estimation
#' for the first starting point.
#' @param nustart A numeric vector containing the starting values for the shape parameter \code{nu}.
#' Default is \code{c(2,2)} for \code{K=2}.
#' @param nustartype A character string indicating whether the initialization of \code{nu} should be \code{"random"},
#' around the values in \code{nustart}, or \code{"exact"}, using the exact values in \code{nustart}.
#' @param gauss A logical value indicating if the algorithm should use the Gaussian distribution.
#' Default is \code{FALSE}.
#' @param laplace A logical value indicating if the algorithm should use the Laplace distribution.
#' Default is \code{FALSE}.
#' @param scale A logical value indicating whether the function should scale the data. Default is \code{TRUE}.
#' @param eps A numeric value specifying the tolerance level of the ECM algorithm.
#' @param maxit An integer specifying the maximum number of iterations.
#' @param verbose A logical value indicating whether to display running output. Default is \code{TRUE}.
#' @param sigbound A numeric vector of length two specifying the lower and upper bounds for resetting the sigma estimates.
#' Default value is \code{c(.01,5)}.
#' @param sr A character string specifying the type of convergence criterion to use.
#' The default is \code{"like"}, but \code{"parameter"} can be used for likelihood-based convergence.
#' @param eta A numeric value specifying the tolerance level for the likelihood-based convergence.
#' Default value is \code{.5}.
#' @details
#' The constrained mixture of generalized normal distributions (CMGND) model is an advanced statistical tool designed for
#' analyzing univariate data characterized by non-normal features such as asymmetry, multi-modality,
#' leptokurtosis, and heavy tails. This model extends the mixture of generalized normal
#' distributions (MGND) by incorporating constraints on the parameters, thereby reducing
#' the number of parameters to be estimated and improving model performance.
#' The CMGND model is defined by the following components:
#' \deqn{f(x|\theta) = \sum_{k=1}^{K} \pi_k f_k(x|\mu_k, \sigma_k, \nu_k)}
#' where:
#' \eqn{\pi_k} are the mixture weights, satisfying \eqn{0 < \pi_k < 1} and \eqn{\sum_{k=1}^{K} \pi_k = 1}.
#' \eqn{f_k(x|\mu_k, \sigma_k, \nu_k)} is the Generalized Normal Distribution for the k-th component with mean \eqn{\mu_k},
#' scale \eqn{\sigma_k}, and shape parameter \eqn{\nu_k}.
#'
#' The parameter space can be constrained by imposing equality constraints
#' such as \eqn{\mu_k = \mu_r}, \eqn{\sigma_k = \sigma_r}, and/or \eqn{\nu_k = \nu_r}
#' for all \eqn{k \in C_r}, where \eqn{C_r} is a partition of the set \eqn{\{1, 2, \ldots, K\}}.
#'
#' The \eqn{k \in C_r} partition for each parameter can be specified
#' by the binary vectors \code{Cmu}, \code{Csigma} and \code{Cnu}.
#' @returns
#' \item{\code{ll}}{The log-likelihood corresponding to the estimated model.}
#' \item{\code{nobs}}{Number of observations.}
#' \item{\code{parameters}}{Data frame of the estimated parameters.}
#' \item{\code{ic}}{Data frame of information criteria. AIC, BIC, HQIC and EDC are returned.}
#' \item{\code{res}}{Matrix of posterior probabilities or responsibilities.}
#' \item{\code{clus}}{Vector of group classifications.}
#' \item{\code{op_it}}{List containing three integers: \code{permstart} the
#' optimal starting value of the permutation of k-means solutions; \code{startnu} the
#' optimal starting value of the shape parameter; \code{iter} number of iterations.}
#' \item{\code{cputime}}{A numeric value indicating the cpu time employed.}
#' \item{\code{info}}{List containing a few of the original user inputs,
#' for use by other dedicated functions of the \code{cmgnd} class.}
#' @examples
#' \dontrun{
#' # Data simulation
#' pi <- c(0.5, 0.3, 0.2)
#' mu <- c(-5, 2, 7)
#' sigma <- c(1, 1, 2)
#' nu <- c(1, 2, 3)
#' n <- 500
#' set.seed(8824312)
#' x <- sim_cmgnd(n, pi, mu, sigma, nu)
#' # Unconstrained model estimation
#' Cmu <- c(0, 0, 0)
#' Csigma <- c(0, 0, 0)
#' Cnu <- c(0, 0, 0)
#' model_unc <- cmgnd(x$sim_data, nstart = 2, K = 3, Cmu, Csigma, Cnu)
#' model_unc$parameters
#' plot_cmgnd(x$sim_data, model_unc)
#' # Constrained model estimation with partition on the scale parameter
#' # Only the first two mixture components have common scale parameter
#' Csigma <- c(1, 1, 0)
#' model_con <- cmgnd(x$sim_data, nstart = 2, K = 3, Cmu, Csigma, Cnu)
#' model_con$parameters
#' plot_cmgnd(x$sim_data, model_con)
#' }
#' @references
#' Bazi, Y., Bruzzone, L., and Melgani, F. (2006). Image thresholding
#' based on the em algorithm and the generalized gaussian distribution.
#' Pattern Recognition, 40(2), pp 619–634.
#'
#' Wen, L., Qiu, Y., Wang, M., Yin, J., and Chen, P. (2022). Numerical characteristics and
#' parameter estimation of finite mixed generalized normal distribution. Communications in
#' Statistics - Simulation and Computation, 51(7), pp 3596–3620.
#'
#' Duttilo, P. (2024). Modelling finacial returns with mixture of
#' generalized normal distributions. Phd thesis, University "G.
#' d’Annunzio" of Chieti-Pescara, Pescara, IT, pp 1-166. Available
#' at \url{https://doi.org/10.48550/arXiv.2411.11847}
#'
#' Duttilo, P., Gattone, S.A., and Iannone, B. (2023). Mixtures
#' of generalized normal distributions and EGARCH models to
#' analyse returns and volatility of ESG and traditional investments.
#' AStA Adv Stat Anal, pp. 1-21
#' \url{https://doi.org/10.1007/s10182-023-00487-7}
#'
#' Duttilo, P. and Gattone, S. A. (2025). Enhancing parameter estimation in finite mixture
#' of generalized normal distributions. Computational Statistics, pp. 1-28
#' \url{https://doi.org/10.1007/s00180-025-01638-x}
#'
#' Duttilo, P., Kume A., and Gattone, S.A. (2023). Constrained Mixtures of Generalized
#' Normal Distributions. In: “SEAS IN Book of short papers 2023”, Pearson, pp.
#' 611-616, ISBN 9788891935618AAVV
#'
#' @export
cmgnd <- function(x, K = 2, Cmu = rep(0, K), Csigma = rep(0, K), Cnu = rep(0, K), nstart = 50,
                  theta = FALSE, nustart = rep(2, K), nustartype = "random", gauss = FALSE,
                  laplace = FALSE, scale = FALSE,eps = 10^-4, maxit = 999, verbose = TRUE,
                  sigbound = c(.1, 5), sr = "like",eta=0.5) {

  grid=1
  # scale
  if (scale == TRUE) {
    x <- scale(x)
  }

  # vector and matrix initialization
  N <- length(x)
  ll_optim <- -Inf
  mu_new <- matrix(NA, nrow = 1, ncol = K)
  sigma_new <- matrix(NA, nrow = 1, ncol = K)
  nu_new <- matrix(NA, nrow = 1, ncol = K)

  sigdata <- stats::sd(x)/((gamma(3/3)/gamma(1/3)))^0.5
  # starting point k-means
  set.seed(41895160)
  rob <- stats::kmeans(x, K, nstart = 10)
  cm1 <- rob$cluster
  mu_km <- tapply(x, cm1, mean)
  sigma_km <- tapply(x, cm1, stats::sd)
  n_perm <- 1
  # permutation
  p_mu <- RcppAlgos::permuteGeneral(mu_km, K)
  p_sigma <- RcppAlgos::permuteGeneral(sigma_km, K)
  if (K > 2) {
    n_perm <- nrow(p_mu)
  }
  # progress bar
  if (verbose == TRUE) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = n_perm,
      style = 3,
      width = n_perm, # Needed to avoid multiple printings
      char = "="
    )
    init <- numeric(n_perm)
    end <- numeric(n_perm)
  }

  for (n in 1:n_perm) {
    # print(n)
    if (verbose == TRUE) {
      init[n] <- Sys.time()
    }

    # input 1st cycle for
    ll_old <- -Inf
    mu_km <- p_mu[n, ]
    sigma_km <- p_sigma[n, ]

    # start 2nd cycle for
    for (sp in 1:nstart) {
      dif <- Inf
      it <- 0
      mu_new <- mu_km
      sigma_new <- sigma_km
      set.seed(sp)

      # mixture weights initialization
      set.seed(sp)
      res <- randgenu(N, K)
      nk <- apply(res, 2, sum)
      pi_new <- nk / N
      # shape parameter initialization
      if (gauss == TRUE) {
        for (k in 1:K) {
          nu_new[k] <- 2
        }
      } else if (laplace == TRUE) {
        for (k in 1:K) {
          nu_new[k] <- 1
        }
      } else if (nustartype == "random") {
        for (k in 1:K) {
          nu_new[k] <- stats::runif(1, 0.8, nustart[k] + 1)
        }
      } else {
        for (k in 1:K) {
          nu_new[k] <- nustart[k]
        }
      }
      nu_copy <- nu_new
      # initialization of the partition parameters
      for (k in 1:K) {
        if (Cmu[k] == 1) {
          mu_new[k] <- mean(mu_km[which(Cmu == 1)])
        }
        if (Csigma[k] == 1) {
          sigma_new[k] <- mean(sigma_km[which(Csigma == 1)])
        }
        if (Cnu[k] == 1) {
          nu_new[k] <- mean(nu_copy[which(Cnu == 1)])
        }
      }

      # 1st starting point with theta parameter matrix
      if (is.matrix(theta) & sp == 1) {
        pi_new <- theta[, 1]
        mu_new <- theta[, 2]
        sigma_new <- theta[, 3]
        nu_new <- theta[, 4]
      }

      # posterior probabilities initialization
      res <- responsibilities(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
      counter <- 0
      estold <- matrix(c(pi_new, mu_new, sigma_new,nu_new), nrow = K)

      # start cycle while

      check=0
      Ks=length(which(Cnu == 1))
      fds=100
      fd=rep(100,K)
      nu_check=rep(0,K)
      ind_cr_nu <- which(Cnu == 0)
      ind_s <- which(Cnu==1)
      ind_cr_mu <- which(Cmu == 0)
      ind_cr_sigma <- which(Csigma == 0)
      nuj=NA
      while (check<1) {
        it <- it + 1
        # update of the of the location parameter
        if(sr=="like"){
          ll_old_mu <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
        }
        if (any(Cmu == 1) & any(Cmu == 0)) {
          mu_update_j <- NA
          for (k in ind_cr_mu) {
            mu_update_j <- ls_mu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K)
          }
          mu_update_s <- ls_mu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cmu)[1]
          mu_new <- (mu_update_s * Cmu) + (mu_update_j * +(!Cmu))
        } else if (all(Cmu == 1)) {
          mu_new <- ls_mu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cmu)
        } else {
          for (k in 1:K) {
            mu_new[k] <- ls_mu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K)
          }
        }
        if (any(Csigma == 1) & any(Csigma == 0)) {
          sigma_update_j <- NA
          for (k in ind_cr_sigma) {
            sigma_update_j <- sigma_j(x, res[, k], mu_new[k], nu_new[k], sigdata, sigbound)
          }
          sigma_update_s <- ls_sigma_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, sigdata, sigbound, Csigma)[1]
          sigma_new <- (sigma_update_s * Csigma) + (sigma_update_j * +(!Csigma))
        } else if (all(Csigma == 1)) {
          sigma_update_s <- ls_sigma_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, sigdata, sigbound, Csigma)[1]
          sigma_new <- rep(sigma_update_s,K)
        } else {
          for (k in 1:K) {
            sigma_new[k] <- sigma_j(x, cbind(res[, k]), mu_new[k], nu_new[k], sigdata, sigbound)
          }
        }
        if(sr=="like"){
          ll_new_sigma <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
          if(sum(nu_check>0)==K){
            if((ll_new_sigma-ll_old_mu)<eps){check=1}
            if((ll_new_sigma-ll_old_mu)>eps){check=1}
          }
        }
        # update of the shape parameter
        if (gauss == TRUE) {
          for (k in 1:K) {
            nu_new[k] <- 2
            nu_check=rep(2,K)
          }
        } else if (laplace == TRUE) {
          for (k in 1:K) {
            nu_new[k] <- 1
          }
        } else {
          if (any(Cnu == 1) & any(Cnu == 0)) {
            #unconstrained nu
            nu_update_j <- NA
            for (k in ind_cr_nu) {
              if(sr=="like"){
                if(abs(fd[k])< eta){nu_check[k]=1}else{
                  nu_check[k]=0
                  nu_update_j <- ls_nu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, sr)
                  nuj<-nu_update_j
                }
                #derivative check
                fd[k]=firstderivative(x,res[,k],mu_new[k],sigma_new[k],nu_new[k])
              }else if(sr=="parameters"){
                nu_update_j <- ls_nu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, sr)
                nuj<-nu_update_j
              }
            }
            #constrained nu
            if(sr=="like"){
              if(abs(fds)<eta){nu_check[ind_s]=1}else{
                nu_check[ind_s]=0
                nu_update_s <- ls_nu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cnu,sr)
              }
            }else if(sr=="parameters"){
              nu_update_s <- ls_nu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cnu,sr)
            }
            #merge estimates

            nu_new <- (nu_update_s * Cnu) + (nuj * +(!Cnu))
            zs=res[, ind_s]
            mus=mu_new[ind_s]
            sigmas=sigma_new[ind_s]
            nus=nu_new[ind_s]
            fds<-firstderivative_s(x,Ks,zs,mus,sigmas,nus)
          } else if (all(Cnu == 1)) {
            if(sr=="like"){
              if(abs(fds)<eta){nu_check[ind_s]=1}else{
                nu_check[ind_s]=0
                nu_new <- ls_nu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cnu,sr)}
              zs=res[, ind_s]
              mus=mu_new[ind_s]
              sigmas=sigma_new[ind_s]
              nus=nu_new[ind_s]
              fds<-firstderivative_s(x,Ks,zs,mus,sigmas,nus)
            }else if(sr=="parameters"){
              nu_new <- ls_nu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cnu,sr)
            }
          } else {
            for (k in 1:K) {
              fd[k]=firstderivative(x,res[,k],mu_new[k],sigma_new[k],nu_new[k])
              if(sr=="like"){
                if(abs(fd[k])< eta){nu_check[k]=1}else{nu_check[k]=0}
                if(nu_check[k]==0){
                  nu_new[k] <- ls_nu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, sr)
                }
              }else if(sr=="parameters"){
                nu_new[k] <- ls_nu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, sr)
              }
            }
          }
        }
        # update of the responsibilities
        res <- responsibilities(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))

        # update of the weights
        pi_new <- rbind(colSums(res) / N)
        # log-likelihood
        ll_new <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
        #
        if (sr == "like"){
          dif <- abs(ll_new - ll_old)
          ll_old <- ll_new
        }else if(sr=="parameters"){
          estnew <- matrix(c(pi_new, mu_new, sigma_new,nu_new), nrow = K)
          dif <- sum((estnew - estold)^2)
          estold <- estnew
          if(dif<eps){check=1}
        }
        if (it > maxit) {
          #dif <- 0
          check=1
          nu_check=rep(1,K)
          print("max iteration")
          estnew <- matrix(c(pi_new, mu_new, sigma_new,nu_new), nrow = K)
          print(estnew)
        }
      } # end while cycle
      if (ll_new > ll_optim) {
        op <- list(pi = pi_new, mu = mu_new, sigma = sigma_new, nu = nu_new)
        ll_optim <- ll_new
        op_it <- list(permstart = n, nustart = sp, iter = it)
      }
    } # end 2nd cycle for

    # progress bar update
    if (verbose == TRUE) {
      end[n] <- Sys.time()
      utils::setTxtProgressBar(pb, n)
      time <- round(lubridate::seconds_to_period(sum(end - init)), 0)
      # Estimated remaining time based on the
      # mean time that took to run the previous iterations
      est <- n_perm * (mean(end[end != 0] - init[init != 0])) - time
      remainining <- round(lubridate::seconds_to_period(est), 0)
      cat(paste(
        " // Execution time:", time,
        " // Estimated time remaining:", remainining
      ), "")
    }
  } # end 1st cycle for

  # here collect results of all the 1st cycle starts
  # close progress bar
  if (verbose == TRUE) {
    close(pb)
  }

  # save
  mu_new <- op$mu
  sigma_new <- op$sigma
  nu_new <- op$nu
  pi_new <- op$pi
  parameters <- pframe(pi_new, mu_new, sigma_new, nu_new)
  res <- responsibilities(x, pi_new, mu_new, sigma_new, nu_new)

  # information criteria (ic): aic, bic, hqic and edc
  if (gauss == TRUE | laplace == TRUE) {
    npar <- rep(NA, 3)
    npar[1] <- K - 1
    if (all(Cmu == 0)) {
      npar[2] <- K
    } else {
      npar[2] <- sum(+(!Cmu)) + 1
    }
    if (all(Csigma == 0)) {
      npar[3] <- K
    } else {
      npar[3] <- sum(+(!Csigma)) + 1
    }
    npar <- sum(npar)
  } else {
    npar <- rep(NA, 4)
    npar[1] <- K - 1
    if (all(Cmu == 0)) {
      npar[2] <- K
    } else {
      npar[2] <- sum(+(!Cmu)) + 1
    }
    if (all(Csigma == 0)) {
      npar[3] <- K
    } else {
      npar[3] <- sum(+(!Csigma)) + 1
    }
    if (all(Cnu == 0)) {
      npar[4] <- K
    } else {
      npar[4] <- sum(+(!Cnu)) + 1
    }
    npar <- sum(npar)
  }
  aic <- -2 * ll_optim + 2 * npar
  bic <- -2 * ll_optim + npar * log(N)
  hqic <- -2 * ll_optim + 2 * npar * log(log(N))
  edc <- -2 * ll_optim + npar * 0.2 * sqrt(N)
  ic <- data.frame(aic, bic, hqic, edc)

  # cluster
  clus <- NA
  for (n in 1:N) {
    clus[n] <- which.max(res[n, ])
  }

  # info
  constraints <- data.frame(Cmu, Csigma, Cnu)
  if ((gauss == FALSE) & (laplace == FALSE)) {
    gnd <- TRUE
  } else {
    gnd <- FALSE
  }
  if (scale == TRUE) {
    info <- list(
      constraints = constraints, gnd = gnd, gauss = gauss, laplace = laplace, scalemeans = attr(x, "scaled:center"),
      scalesd = attr(x, "scaled:scale"), grid = grid)
  } else {
    info <- list(
      constraints = constraints, gnd = gnd, gauss = gauss, laplace = laplace,
      scalemeans = NULL, scalesd = NULL, grid = grid)
  }

  return(list(ll = ll_optim, nobs = N, npar = npar, parameters = parameters,
              ic = ic, res = res, clus = clus, op_it = op_it, cputime = time, info = info
  ))
}
