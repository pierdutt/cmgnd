.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
  #   ____  _    __  ____  _   _  ____
  #  / ___|| \\  / | / ___| |\\ | ||  _ \\
  # | |    | |\\/| || |  __ | \\| || | | |
  # | |___ | |  | || |__| || |\\ || |_| |
  # \\ ____||_|  |_| \\ ____||_| \\||____/
  #
  #   cmgnd: constrained mixture of generalised normal distributions
  #   Version 0.1.1
  ")
}


#' pframe
#'
#' pframe
#' @noRd
pframe <- function(pi, mu, sigma, nu) {
  parameters <- data.frame(
    col1 = c(pi),
    col2 = c(mu),
    col3 = c(sigma),
    col4 = c(nu)
  )
  colnames(parameters) <- c("pi", "mu", "sigma", "nu")
  return(parameters)
}

#' random partition of n objects with k groups
#'
#' random partition of n objects with k groups
#' @noRd
randgenu <- function(n, k) {
  # randomly generate a membership matrix U
  U <- matrix(stats::runif(n * k), nrow = n, ncol = k)
  U <- diag(1 / colSums(t(U))) %*% U
  ind <- max.col(U)
  U <- matrix(0, nrow = n, ncol = k)
  for (i in 1:n) {
    U[i, ind[i]] <- 1
  }
  return(U)
}


#' Update for Unconstrained Location Parameter without Adaptive Step-Size
#'
#' Updating equation for the unconstrained location parameter \eqn{\mu_k} without adaptive step-size.
#' @param x A numeric vector of observations.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @return \code{mu_j} returns a numeric value indicating the updating for the uncontrained shape parameter \eqn{\mu_k}.
#' @noRd
mu_j <- function(x, z, mu, nu) {
  num_comp1 <- sum(subset(z, x >= mu) * (subset(x, x >= mu) - mu)^(nu - 1))
  num_comp2 <- sum(subset(z, x < mu) * (mu - subset(x, x < mu))^(nu - 1))
  num <- num_comp1 - num_comp2
  den_comp1 <- sum(subset(z, x >= mu) * (subset(x, x >= mu) - mu)^(nu - 2) * (nu - 1))
  den_comp2 <- sum(subset(z, x < mu) * (mu - subset(x, x < mu))^(nu - 2) * (nu - 1))
  den <- den_comp1 + den_comp2
  mu_new <- mu +  (num / den)
  if (any(is.nan(mu_new))) {
    mu_new <- mu + stats::runif(1, -1, 1) * 1e-5
  }
  if (any(x == mu_new)) {
    mu_new <- mu_new + stats::runif(1, -1, 1) * 1e-5
  }
  return(mu_new)
}

#' Update for Unconstrained Location Parameter
#'
#' Updating equation for the unconstrained location parameter \eqn{\mu_k}.
#' @param x A numeric vector of observations.
#' @param pi_new A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu_new A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma_new A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu_new A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param res A matrix of posterior probabilities or responsibilities.
#' @param k An integer indicating the k-th mixture component.
#' @param K An integer specifying the number of mixture components to fit.
#' @return \code{ls_mu_j} returns a numeric value indicating the updating for the uncontrained location parameter \eqn{\mu_k}.
#' @noRd
ls_mu_j <- function(x, pi_new, mu_new, sigma_new, nu_new, res, k, K) {
  mu_temp <- mu_new
  ll_mu_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  mu_temp[k] <- mu_j(x, res[, k], mu_temp[k], nu_new[k])
  ll_mu_new <- log.likelihood(x, rbind(pi_new), rbind(mu_temp), rbind(sigma_new), rbind(nu_new))
  if (ll_mu_new>ll_mu_old) {
    mu_new[k] <- mu_temp[k]
  }
  return(mu_new[k])
}

#' Update for Constrained Location Parameters
#'
#' Updating equation for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}
#' @param x A numeric vector of observations.
#' @param K An integer specifying the number of mixture components to fit.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @return \code{mu_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
mu_s <- function(x, K, z, mu, nu) {
  xx <- x
  num <- den <- rep(NA, K)
  for (k in 1:K) {
    num[k] <- sum(subset(z[, k], xx >= mu[k]) * (subset(xx, xx >= mu[k]) - mu[k])^(nu[k] - 1)) - sum(subset(z[, k], xx < mu[k]) * (mu[k] - subset(xx, xx < mu[k]))^(nu[k] - 1))
    den[k] <- sum(subset(z[, k], xx >= mu[k]) * (subset(xx, xx >= mu[k]) - mu[k])^(nu[k] - 2) * (nu[k] - 1)) + sum(subset(z[, k], xx < mu[k]) * (mu[k] - subset(xx, xx < mu[k]))^(nu[k] - 2) * (nu[k] - 1))
  }
  r <- sum(num) / sum(den)
  mu_new <- rep((mu[1] + r), K)
  if (any(xx == mu_new[1])) {
    mu_new <- mu_new + stats::runif(1, -1, 1) * 1e-5
  }
  return(mu_new)
}

#' Update for Constrained Location Parameters
#'
#' Updating equation for constrained location parameters  \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @param x A numeric vector of observations.
#' @param pi_new A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu_new A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma_new A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu_new A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param K An integer specifying the number of mixture components to fit.
#' @param res A matrix of posterior probabilities or responsibilities.
#' @param Cmu A binary vector indicating mixture components
#' for the location parameter. The k-th element is set to 1 if the
#' k-th mixture component belongs to the Cr partition, and 0 otherwise.
#' Default is \code{c(0,0)}, indicating no mixture partition with \code{K=2}.
#' @return \code{ls_mu_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
ls_mu_s <- function(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cmu) {
  ll_mu_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  mu_old <- mu_new
  mu_new[which(Cmu == 1)] <- mu_s(
    x, length(which(Cmu == 1)), res[, which(Cmu == 1)],
    mu_new[which(Cmu == 1)], nu_new[which(Cmu == 1)])

  ll_mu_new <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))

  if (ll_mu_old > ll_mu_new) {
    mu_new <- mu_old
  }

  return(mu_new)
}

#' Update for Unconstrained Scale Parameter
#'
#' Updating equation for the unconstrained scale parameter \eqn{\sigma_k}.
#' @param x A numeric vector of observations.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param sigdata A numeric value representing the standard deviation of the data vector \code{x}.
#' @param sigbound A numeric vector of length two specifying the lower and upper bounds for resetting the sigma estimates.
#' Default value is \code{c(.01,5)}.
#' @return \code{sigma_j} returns a numeric value indicating the updating for the uncontrained scale parameter \eqn{\sigma_k}.
#' @noRd
sigma_j <- function(x, z, mu, nu, sigdata, sigbound) {

  num <- nu * sum(z * abs(x - mu)^nu)
  den <- sum(z)
  sigma_new <- (num / den)^(1 / nu)
  sigma_new <- max(sigbound[1] * sigdata, sigma_new)
  if (sigma_new > sigbound[2] * sigdata) {
    sigma_new <- sigdata
  }
  return(sigma_new)
}


#' Update for Constrained Scale Parameters with Adaptive Step-Size
#'
#' Updating equation for constrained scale parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r} with Adaptive Step-Size.
#' @param x A numeric vector of observations.
#' @param K An integer specifying the number of mixture components to fit.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param ss A numeric value specifying the fixed step-length.
#' @param sigdata A numeric value representing the standard deviation of the data vector \code{x}.
#' @param sigbound A numeric vector of length two specifying the lower and upper bounds for resetting the sigma estimates.
#' Default value is \code{c(.01,5)}.
#' @return \code{sigma_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
sigma_s <- function(x, K, z, mu, sigma, nu, ss, sigdata, sigbound) {
  sigma_new <- comp1 <- comp2 <- rep(NA, K)
  for (k in 1:K) {
    comp1[k] <- ((sum(z[, k]) * (-1 / sigma[k])) + (nu[k] / (sigma[k]^(nu[k] + 1))) * sum(z[, k] * abs(x - mu[k])^nu[k]))
    comp2[k] <- ((sum(z[, k]) * (1 / (sigma[k]^2))) + (nu[k] * (-nu[k] - 1) * (sigma[k]^(-nu[k] - 2)) * sum(z[, k] * abs(x - mu[k])^nu[k])))
  }
  r <- sum(comp1) / sum(comp2)
  for (k in 1:K) {
    sigma_new[k] <- sigma[1] - ss * r
    sigma_new[k] <- max(sigbound[1] * sigdata, sigma_new[k])
    if (sigma_new[k] > sigbound[2] * sigdata) {
      sigma_new[k] <- sigdata
    }
  }
  return(rbind(sigma_new))
}



#' Update for Constrained Scale Parameters
#'
#' Updating equation for constrained scale parameters  \eqn{\sigma_k = \sigma_r} for all \eqn{k \in C_r}.
#' @param x A numeric vector of observations.
#' @param pi_new A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu_new A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma_new A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu_new A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param K An integer specifying the number of mixture components to fit.
#' @param res A matrix of posterior probabilities or responsibilities.
#' @param grid A numeric vector specifying the fixed step-size for the scale parameter update. Default is \code{1}.
#' @param sigdata A numeric value representing the standard deviation of the data vector \code{x}.
#' @param sigbound A numeric vector of length two specifying the lower and upper bounds for resetting the sigma estimates.
#' Default value is \code{c(.01,5)}.
#' @param Csigma A binary vector indicating mixture components
#' for the scale parameter. The k-th element is set to 1 if the k-th
#' mixture component belongs to the Cr partition, and 0 otherwise.
#' @return \code{ls_sigma_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\sigma_k = \sigma_r} for all \eqn{k \in C_r}.
#' @noRd
ls_sigma_s <- function(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, sigdata, sigbound, Csigma) {
  ls_ll_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  sigma_old <- sigma_new
  G <- length(grid)
  lssigmasts <- matrix(NA, nrow = G, ncol = K)
  lsllsts <- NA

  for (i in 1:G) {
    lssigmasts[i, which(Csigma == 1)] <- sigma_s(x,
                                                 K = length(which(Csigma == 1)), res[, which(Csigma == 1)], mu_new[which(Csigma == 1)], sigma_new[which(Csigma == 1)],
                                                 nu_new[which(Csigma == 1)], grid[i], sigdata, sigbound
    )[1]


    lssigmasts[i, which(Csigma == 0)] <- sigma_new[which(Csigma == 0)]

    lsllsts[i] <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(lssigmasts[i, ]), rbind(nu_new))
  }


  sigma_new <- lssigmasts[which.max(lsllsts), which(Csigma == 1)]
  #print(sigma_new)
  sigma_new <- max(sigbound[1] * sigdata, sigma_new)

  if (sigma_new[1] > sigbound[2] * sigdata) {
    sigma_new <- sigdata
  }
  if (ls_ll_old > max(lsllsts)) {
    sigma_new <- sigma_old
  }
  #print(sigma_new)
  return(sigma_new)
}


#' Update for Unconstrained Shape Parameter with Adaptive Step-Size
#'
#' Updating equation for the unconstrained shape parameter \eqn{\nu_k} with and without adaptive step-size.
#' @param x A numeric vector of observations.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param sr A character string specifying the type of convergence criterion to use.
#' The default is \code{"parameter"}, but \code{"like"} can be used for likelihood-based convergence.
#' @return \code{nu_j} returns a numeric value indicating the updating for the uncontrained shape parameter \eqn{\nu_k}.
#' @noRd
nu_j <- function(x, z, mu, sigma, nu, sr) {
  indzero <- which((x - mu) == 0)
  A <- (1 / nu) * ((1 / nu) * digamma(1 / nu) + 1)
  B <- (abs((x - mu) / (sigma))^nu) * log(abs((x - mu) / (sigma)))
  C <- (-1 / nu^2) * (1 + (2 / nu) * digamma(1 / nu) + (1 / nu^2) * trigamma(1 / nu))
  D <- (abs((x - mu) / (sigma))^nu) * (log(abs((x - mu) / (sigma))))^2
  B[indzero] <- 0
  D[indzero] <- 0
  num <- (sum(z * A) - sum(z * B))
  den <- (sum(z * C) - sum(z * D))
  if (sr == "like"){
    nu_new <- nu - exp(-nu) * (num / den)
  }else{nu_new <- nu -  (num / den)}
  if ((is.na(nu_new))) {
    nu_new <- 2
  }
  if (nu_new < 0) {
    nu_new <- 2
  }
  return(nu_new)
}

#' firstderivative
#'
#' firstderivative
#' @noRd
firstderivative=function(x, z, mu, sigma, nu){
  indzero <- which((x - mu) == 0)
  A <- (1 / nu) * ((1 / nu) * digamma(1 / nu) + 1)
  B <- (abs((x - mu) / (sigma))^nu) * log(abs((x - mu) / (sigma)))
  B[indzero] <- 0
  num <- (sum(z * A) - sum(z * B))
  return(num)
}

#' Update for Unconstrained Shape Parameter
#'
#' Updating equation for the unconstrained shape parameter \eqn{\nu_k}.
#' @param x A numeric vector of observations.
#' @param pi_new A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu_new A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma_new A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu_new A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param res A matrix of posterior probabilities or responsibilities.
#' @param k An integer indicating the k-th mixture component.
#' @param K An integer specifying the number of mixture components to fit.
#' @param sr A character string specifying the type of convergence criterion to use.
#' The default is \code{"parameter"}, but \code{"like"} can be used for likelihood-based convergence.
#' @return \code{ls_nu_j} returns a numeric value indicating the updating for the uncontrained shape parameter \eqn{\nu_k}.
#' @noRd
ls_nu_j <- function(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, sr) {
  ll_nu_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  nu_old <- nu_new[k]
  nu_new[k] <- nu_j(x, res[,k], mu_new[k], sigma_new[k], nu_new[k], sr)
  ll_nu_new <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  if (ll_nu_old > ll_nu_new) {
    #print("negative dif")
    nu_new[k] <- nu_old
  }
  return(nu_new[k])
}


#' Update for Constrained Shape Parameters with Adaptive Step-Size
#'
#' Updating equation for constrained shape parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r} with Adaptive Step-Size.
#' @param x A numeric vector of observations.
#' @param K An integer specifying the number of mixture components to fit.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param sr A character string specifying the type of convergence criterion to use.
#' The default is \code{"parameter"}, but \code{"like"} can be used for likelihood-based convergence.
#' @return \code{nu_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
nu_s <- function(x, K, z, mu, sigma, nu,sr) {
  num <- den <- rep(NA, K)
  for (k in 1:K) {
    #print(k)
    num[k] <- (sum(z[, k] * (1 / nu[k]) * ((1 / nu[k]) * digamma(1 / nu[k]) + 1)) - sum(z[, k] * (abs((x - mu[k]) / (sigma[k]))^nu[k]) * log(abs((x - mu[k]) / (sigma[k])))))
    den[k] <- (sum(z[, k] * (-1 / nu[k]^2) * (1 + (2 / nu[k]) * digamma(1 / nu[k]) + (1 / nu[k]^2) * trigamma(1 / nu[k]))) - sum(z[, k] * (abs((x - mu[k]) / (sigma[k]))^nu[k]) * (log(abs((x - mu[k]) / (sigma[k]))))^2))
  }
  if(sr=="like"){
    nu_new <- matrix(nu[1] - exp(-nu[1]) * (sum(num) / sum(den)), ncol = K, nrow = 1)
  }else{

    nu_new <- matrix(nu[1] - (sum(num) / sum(den)), ncol = K, nrow = 1)

  }
  for (k in 1:K) {
    if (nu_new[k] < 0) {
      nu_new[k] <- 2
    }
  }
  return(nu_new)
}

#' firstderivative_s
#'
#' firstderivative_s
#' @noRd
firstderivative_s=function(x,K,z, mu, sigma, nu){
  num <- NA
  for (k in 1:K) {
    num[k] <- (sum(z[, k] * (1 / nu[k]) * ((1 / nu[k]) * digamma(1 / nu[k]) + 1)) - sum(z[, k] * (abs((x - mu[k]) / (sigma[k]))^nu[k]) * log(abs((x - mu[k]) / (sigma[k])))))
  }
  fd=sum(num)
  return(fd)
}

#' Update for Constrained Shape Parameters
#'
#' Updating equation for constrained shape parameters  \eqn{\nu_k = \nu_r} for all \eqn{k \in C_r}.
#' @param x A numeric vector of observations.
#' @param pi_new A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu_new A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma_new A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu_new A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param K An integer specifying the number of mixture components to fit.
#' @param res A matrix of posterior probabilities or responsibilities.
#' Default is \code{TRUE} to control the degeneracy of the shape parameter.
#' @param Cnu A binary vector indicating mixture components
#' for the shape parameter. The k-th element is set to 1 if the k-th mixture component belongs to the Cr partition, and 0 otherwise.
#' Default is \code{c(0,0)}, indicating no mixture partition with \code{K=2}.
#' @param sr A character string specifying the type of convergence criterion to use.
#' The default is \code{"parameter"}, but \code{"like"} can be used for likelihood-based convergence.
#' @return \code{ls_nu_s} returns a numeric value indicating the updating for constrained shape parameters \eqn{\nu_k = \nu_r} for all \eqn{k \in C_r}.
#' @noRd
ls_nu_s <- function(x, pi_new, mu_new, sigma_new, nu_new, K, res, Cnu, sr) {
  ll_nu_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  nu_old <- nu_new
  nu_new[ which(Cnu == 1)] <- nu_s(
    x, length(which(Cnu == 1)), res[, which(Cnu == 1)], mu_new[which(Cnu == 1)],
    sigma_new[which(Cnu == 1)], nu_new[which(Cnu == 1)],sr)
  ll_nu_new <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  if (ll_nu_old > ll_nu_new) {
    nu_new <- nu_old
  }
  return(nu_new)
}

#' Updating Equation for Posterior Probabilities
#'
#' Updating equation for posterior probabilities.
#' @param x A numeric vector of observations.
#' @param pi A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @return \code{responsibilities} returns a matrix of posterior probabilities.
#' @noRd
responsibilities <- function(x, pi, mu, sigma, nu) {
  N <- length(x)
  K <- length(pi)
  num <- matrix(NA, nrow = N, ncol = K)
  for (k in 1:K) {
    num[, k] <- pi[k] * dgnd(x, mu[k], sigma[k], nu[k])
  }
  # check for zeros
  nn <- c(num)
  nn[nn < 1e-10] <- 1e-10
  num <- matrix(nn, nrow = N, ncol = K)
  den <- rowSums(num)
  res <- num / den
  # check for nearly empty components
  rpi <- rbind(colSums(res) / N)
  if (any(rpi < 0.1)) {
    ant <- matrix(stats::runif(N * K), nrow = N)
    res <- diag(1 / rowSums(ant)) %*% ant
  }
  return(res)
}

#' Log-likelihood Function
#'
#' Log-likelihood function.
#' @param x A numeric vector of observations.
#' @param pi A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @return \code{loglikelihood} returns the log-likelihood value.
#' @export
#' @method log likelihood
log.likelihood <- function(x, pi, mu, sigma, nu) {
  N <- length(x)
  K <- length(pi)
  ll0 <- matrix(NA, nrow = N, ncol = K)
  for (k in 1:K) {
    ll0[, k] <- pi[, k] * dgnd(x, mu[, k], sigma[, k], nu[, k])
  }
  ll0[which(ll0==0)] <- 1e-51
  ll <- sum(log(rowSums(ll0)))
  return(ll)
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("..density.."))
}





