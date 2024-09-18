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

#' The Generalized Normal Distribution (GND)
#'
#' Density function for the GND with location parameter \code{mu},
#' scale parameter \code{sigma} and shape parameter \code{nu}.
#' @param x A numeric vector of observations.
#' @param mu A numeric value indicating the location parameter \eqn{\mu}.
#' @param sigma A numeric value indicating the scale parameter \eqn{\sigma}.
#' @param nu A numeric value indicating the shape parameter \eqn{\nu}.
#' @return \code{dgnd} returns the density.
#' @details
#' If \code{mu}, \code{sigma} and \code{nu} are not specified
#' they assume the default values of 0, 1 and 2, respectively.
#' The GND distribution has density
#' \deqn{ f_{GND}(x|\mu,\sigma,\nu)=\frac{\nu}{2\sigma\Gamma(1\mathbin{/}\nu)}\exp\Biggr\{-\Biggr|\frac{x-\mu}{\sigma}\Biggr|^\nu\Biggr\}.}
#' The shape parameter \eqn{\nu} controls both the peakedness and tail weights.
#' If \eqn{\nu=1} the GND reduces to the Laplace distribution and if \eqn{\nu=2}
#' it coincides with the normal distribution. It is noticed that \eqn{1<\nu<2}
#' yields an intermediate distribution between the normal and the Laplace distribution.
#' As limit cases, for \eqn{\nu\rightarrow\infty} the distribution tends to a uniform
#' distribution, while for \eqn{\nu\rightarrow0} it will be impulsive.
#' @references Nadarajah, S. (2005). A generalized normal distribution.
#' \emph{Journal of Applied Statistics}, \eqn{32(7):685–694}.
#' @export
dgnd <- function(x, mu = 0, sigma = 1, nu = 2) {
  comp1 <- (nu / (2 * sigma * (gamma(1 / nu))))
  comp2 <- exp(-abs((x - mu) / (sigma))^nu)
  density <- (comp1 * comp2)
  return(density)
}

#' Marginal Density Estimation for CMGND Models
#'
#' @description
#' This function estimates the marginal density for univariate constrained
#' mixture of generalized normal distribution (CMGND) models.
#'
#' @param x A numeric vector representing the observed data points.
#' @param parameters A matrix or data.frame containing the parameters of the CMGND model.
#' This can also be an object returned from the `cmgnd()` function, representing a previously estimated CMGND model.
#'
#' @return A vector of density estimates corresponding to the input data `x`.
#'
#' @details
#' The function computes the marginal density based on the provided parameters of the CMGND model.
#' It can handle both newly supplied parameters or those extracted from an existing CMGND model object.
#'
#' @seealso `cmgnd()` for estimating the model parameters.
#'
#' @export
dcmgnd <- function(x, parameters) {
  if (any(class(parameters) == "list")) {
    parameters <- parameters$parameters
  } else {
    parameters <- parameters
  }
  K <- nrow(parameters)
  den <- matrix(NA, ncol = K, nrow = length(x))
  for (k in 1:K) {
    den[, k] <- dgnd(x, parameters[k, 2], parameters[k, 3], parameters[k, 4]) * parameters[k, 1]
  }
  return(rowSums(den))
}

#' Update for Unconstrained Location Parameter without Adaptive Step-Size
#'
#' Updating equation for the unconstrained location parameter \eqn{\mu_k} without adaptive step-size.
#' @param x A numeric vector of observations.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param ss A number indicating the step-size.
#' @return \code{mu_j} returns a numeric value indicating the updating for the uncontrained shape parameter \eqn{\mu_k}.
#' @noRd
mu_j <- function(x, z, mu, nu, ss) {
  num_comp1 <- sum(subset(z, x >= mu) * (subset(x, x >= mu) - mu)^(nu - 1))
  num_comp2 <- sum(subset(z, x < mu) * (mu - subset(x, x < mu))^(nu - 1))
  num <- num_comp1 - num_comp2
  den_comp1 <- sum(subset(z, x >= mu) * (subset(x, x >= mu) - mu)^(nu - 2) * (nu - 1))
  den_comp2 <- sum(subset(z, x < mu) * (mu - subset(x, x < mu))^(nu - 2) * (nu - 1))
  den <- den_comp1 + den_comp2
  mu_new <- mu + ss * (num / den)
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
#' @param grid A numeric vector specifying the step-lengths that determine the largest increase in the \eqn{Q(\theta, \theta(m-1))}
#' function for the location parameter. Default is \code{c(0.1,0.2,0.5,0.8,1)}.
#' @param it An integer specifying the iteration number.
#' @return \code{ls_mu_j} returns a numeric value indicating the updating for the uncontrained location parameter \eqn{\mu_k}.
#' @noRd
ls_mu_j <- function(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, grid, it) {
  ls_ll_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  mu_old <- mu_new[k]
  G <- length(grid)
  lsmusts <- matrix(mu_new, nrow = G, ncol = K, byrow = T)
  lsllsts <- NA
  for (i in 1:G) {
    lsmusts[i, k] <- mu_j(x, res[, k], mu_new[k], nu_new[k], grid[i])
    lsllsts[i] <- log.likelihood(x, rbind(pi_new), rbind(lsmusts[i, ]), rbind(sigma_new), rbind(nu_new))
  }
  mu_new[k] <- lsmusts[which.max(lsllsts), k]
  if (ls_ll_old > max(lsllsts)) {
    mu_new[k] <- mu_old
  } # print("stay");print(c(it,k))
  return(mu_new[k])
}

#' Update for Constrained Location Parameters with Adaptive Step-Size
#'
#' Updating equation for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r} with Adaptive Step-Size.
#' @param x A numeric vector of observations.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param K An integer specifying the number of mixture components to fit.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param ss A numeric value specifying the step-length.
#' @return \code{mu_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
mu_s <- function(x, K, z, mu, nu, ss) {
  ss <- 1
  xx <- x
  num <- den <- rep(NA, K)
  for (k in 1:K) {
    num[k] <- sum(subset(z[, k], xx >= mu[k]) * (subset(xx, xx >= mu[k]) - mu[k])^(nu[k] - 1)) - sum(subset(z[, k], xx < mu[k]) * (mu[k] - subset(xx, xx < mu[k]))^(nu[k] - 1))
    den[k] <- sum(subset(z[, k], xx >= mu[k]) * (subset(xx, xx >= mu[k]) - mu[k])^(nu[k] - 2) * (nu[k] - 1)) + sum(subset(z[, k], xx < mu[k]) * (mu[k] - subset(xx, xx < mu[k]))^(nu[k] - 2) * (nu[k] - 1))
  }
  r <- sum(num) / sum(den)
  mu_new <- rep((mu[1] + ss * r), K)
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
#' @param grid A numeric vector specifying the step-lengths that determine the largest increase in the \eqn{Q(\theta, \theta(m-1))}
#' function for the location parameter. Default is \code{c(0.1,0.2,0.5,0.8,1)}.
#' @return \code{ls_mu_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
ls_mu_s <- function(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, Cmu) {
  grid <- 1
  ls_ll_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  mu_old <- mu_new
  G <- length(grid)
  lsmusts <- matrix(NA, nrow = G, ncol = K)
  lsllsts <- NA
  for (i in 1:G) {
    lsmusts[i, which(Cmu == 1)] <- mu_s(
      x, length(which(Cmu == 1)), res[, which(Cmu == 1)],
      mu_new[which(Cmu == 1)], nu_new[which(Cmu == 1)], grid[i]
    )
    lsmusts[i, which(Cmu == 0)] <- mu_new[which(Cmu == 0)]

    lsllsts[i] <- log.likelihood(x, rbind(pi_new), rbind(lsmusts[i, ]), rbind(sigma_new), rbind(nu_new))
  }


  mu_new <- lsmusts[which.max(lsllsts), which(Cmu == 1)]
  if (ls_ll_old > max(lsllsts)) {
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
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param K An integer specifying the number of mixture components to fit.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param ss A numeric value specifying the step-length.
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
#' @param grid A numeric vector specifying the step-lengths that determine the largest increase in the \eqn{Q(\theta, \theta(m-1))}
#' function for the scale parameter. Default is \code{c(0.1,0.2,0.5,0.8,1)}.
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
  # print(sigma_new)
  sigma_new <- max(sigbound[1] * sigdata, sigma_new)

  if (sigma_new[1] > sigbound[2] * sigdata) {
    sigma_new <- sigdata
  }

  if (ls_ll_old > max(lsllsts)) {
    sigma_new <- sigma_old
  }

  # print(sigma_new)
  return(sigma_new)
}


#' Update for Unconstrained Shape Parameter without Adaptive Step-Size
#'
#' Updating equation for the unconstrained shape parameter \eqn{\nu_k} without adaptive step-size.
#' @param x A numeric vector of observations.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param ss A number indicating the step-size.
#' @param inversenu A logical value indicating if the algorithm should use an adaptive step size for the shape parameter.
#' Default is \code{TRUE} to control the degeneracy of the shape parameter.
#' @return \code{nu_j} returns a numeric value indicating the updating for the uncontrained shape parameter \eqn{\nu_k}.
#' @noRd
nu_j <- function(x, z, mu, sigma, nu, ss, inversenu) {
  indzero <- which((x - mu) == 0)
  A <- (1 / nu) * ((1 / nu) * digamma(1 / nu) + 1)
  B <- (abs((x - mu) / (sigma))^nu) * log(abs((x - mu) / (sigma)))
  C <- (-1 / nu^2) * (1 + (2 / nu) * digamma(1 / nu) + (1 / nu^2) * trigamma(1 / nu))
  D <- (abs((x - mu) / (sigma))^nu) * (log(abs((x - mu) / (sigma))))^2
  B[indzero] <- 0
  D[indzero] <- 0
  num <- (sum(z * A) - sum(z * B))
  den <- (sum(z * C) - sum(z * D))
  nu_new <- nu - ss * (num / den)
  if (inversenu == TRUE) {
    nu_new <- nu - ss * nu^(-1) * (num / den)
  }
  if ((is.na(nu_new))) {
    nu_new <- 2
  }
  if (nu_new < 0) {
    nu_new <- 2
  }
  return(nu_new)
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
#' @param inversenu A logical value indicating if the algorithm should use an adaptive step size for the shape parameter.
#' Default is \code{TRUE} to control the degeneracy of the shape parameter.
#' @param grid A numeric vector specifying the step-lengths that determine the largest increase in the \eqn{Q(\theta, \theta(m-1))}
#' function for the shape parameter. Default is \code{c(0.1,0.2,0.5,0.8,1)}.
#' @param it An integer specifying the iteration number.
#' @return \code{ls_nu_j} returns a numeric value indicating the updating for the uncontrained shape parameter \eqn{\nu_k}.
#' @noRd
ls_nu_j <- function(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, inversenu, grid, it) {
  ls_ll_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  nu_old <- nu_new[k]
  G <- length(grid)
  lsnusts <- matrix(nu_new, nrow = G, ncol = K, byrow = T)
  lsllsts <- NA
  for (i in 1:G) {
    lsnusts[i, k] <- nu_j(x, res[, k], mu_new[k], sigma_new[k], nu_new[k], grid[i], inversenu)
    lsllsts[i] <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(lsnusts[i, ]))
  }
  nu_new[k] <- lsnusts[which.max(lsllsts), k]
  if (ls_ll_old > max(lsllsts)) {
    nu_new[k] <- nu_old
  }
  return(nu_new[k])
}

#' Update for Constrained Shape Parameters with Adaptive Step-Size
#'
#' Updating equation for constrained shape parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r} with Adaptive Step-Size.
#' @param x A numeric vector of observations.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#' @param K An integer specifying the number of mixture components to fit.
#' @param z A matrix of posterior probabilities or responsibilities.
#' @param ss A numeric value specifying the step-length.
#' @param inversenu A logical value indicating if the algorithm should use an adaptive step size for the shape parameter.
#' @return \code{nu_s} returns a numeric value indicating the updating for constrained location parameters \eqn{\mu_k = \mu_r} for all \eqn{k \in C_r}.
#' @noRd
nu_s <- function(x, K, z, mu, sigma, nu, ss, inversenu = FALSE) {
  num <- den <- rep(NA, K)
  for (k in 1:K) {
    num[k] <- (sum(z[, k] * (1 / nu[k]) * ((1 / nu[k]) * digamma(1 / nu[k]) + 1)) - sum(z[, k] * (abs((x - mu[k]) / (sigma[k]))^nu[k]) * log(abs((x - mu[k]) / (sigma[k])))))
    den[k] <- (sum(z[, k] * (-1 / nu[k]^2) * (1 + (2 / nu[k]) * digamma(1 / nu[k]) + (1 / nu[k]^2) * trigamma(1 / nu[k]))) - sum(z[, k] * (abs((x - mu[k]) / (sigma[k]))^nu[k]) * (log(abs((x - mu[k]) / (sigma[k]))))^2))
  }
  nu_new <- matrix(nu[1] - ss * (sum(num) / sum(den)), ncol = K, nrow = 1)
  if (inversenu == TRUE) {
    nu_new <- matrix(nu[1] - ss * nu[1]^(-1) * (sum(num) / sum(den)), ncol = K, nrow = 1)
  }
  for (k in 1:K) {
    if (nu_new[k] < 0) {
      nu_new[k] <- 2
    }
  }
  return(nu_new)
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
#' @param inversenu A logical value indicating if the algorithm should use an adaptive step size for the shape parameter.
#' Default is \code{TRUE} to control the degeneracy of the shape parameter.
#' @param grid A numeric vector specifying the step-lengths that determine the largest increase in the \eqn{Q(\theta, \theta(m-1))}
#' function for the shape parameter. Default is \code{c(0.1,0.2,0.5,0.8,1)}.
#' @return \code{ls_nu_s} returns a numeric value indicating the updating for constrained shape parameters \eqn{\nu_k = \nu_r} for all \eqn{k \in C_r}.
#' @noRd
ls_nu_s <- function(x, pi_new, mu_new, sigma_new, nu_new, K, res, inversenu, grid, Cnu) {
  ls_ll_old <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
  nu_old <- nu_new
  G <- length(grid)
  lsnusts <- matrix(NA, nrow = G, ncol = K)
  lsllsts <- NA
  for (i in 1:G) {
    lsnusts[i, which(Cnu == 1)] <- nu_s(
      x, length(which(Cnu == 1)), res[, which(Cnu == 1)], mu_new[which(Cnu == 1)],
      sigma_new[which(Cnu == 1)], nu_new[which(Cnu == 1)], grid[i], inversenu
    )
    lsnusts[i, which(Cnu == 0)] <- nu_new[which(Cnu == 0)]
    lsllsts[i] <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(lsnusts[i, ]))
  }
  nu_new <- lsnusts[which.max(lsllsts), which(Cnu == 1)]
  if (ls_ll_old > max(lsllsts)) {
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
#' @noRd
log.likelihood <- function(x, pi, mu, sigma, nu) {
  N <- length(x)
  K <- length(pi)
  ll0 <- matrix(NA, nrow = N, ncol = K)
  for (k in 1:K) {
    ll0[, k] <- pi[, k] * dgnd(x, mu[, k], sigma[, k], nu[, k])
  }
  ll0[which(ll0 == 0)] <- 1e-51
  ll <- sum(log(rowSums(ll0)))
  return(ll)
}



#' sim_cmgnd: Function to Simulate Univariate Constrained Mixtures of Generalized Normal Distributions
#'
#' @description
#' Simulate univariate constrained mixture of generalized normal distribution models.
#' Remeber to set the set.seed() before the function sim_cmgnd().
#'
#' @param n A numeric value indicating the total number of observations to simulate.
#' @param pi A numeric vector of the mixture weights \eqn{\pi_k}.
#' @param mu A numeric vector of the location parameter \eqn{\mu_k}.
#' @param sigma A numeric vector of the scale parameter \eqn{\sigma_k}.
#' @param nu A numeric vector of the shape parameter \eqn{\nu_k}.
#'
#' @returns
#' \item{\code{sim_data}}{The simulated data.}
#' \item{\code{sim_clus}}{The cluster indication of simulated data.}
#'
#' @export
sim_cmgnd <- function(n = 1000, pi = rep(0.5, 2), mu = c(1, 5), sigma = c(1, 1), nu = c(2, 2)) {
  K <- length(pi)
  n_lenght <- NA
  sim_clus <- sim_data <- list()
  for (k in 1:K) {
    n_lenght[k] <- round(n * pi[k])
    sim_data[[k]] <- gnorm::rgnorm(n_lenght[k], mu[k], sigma[k], nu[k])
    sim_clus[[k]] <- rep(k, n_lenght[k])
  }
  sim_data <- purrr::list_c(sim_data)
  sim_clus <- purrr::list_c(sim_clus)
  df <- data.frame(sim_data, sim_clus)
  df <- df[sample(nrow(df), n, replace = FALSE), ]
  return(list(sim_data = df$sim_data, sim_clus = df$sim_clus))
}


#' Plot Marginal and Mixture Component Densities of the CMGND Model
#'
#' @description
#' This function generates a plot displaying both the marginal density and individual mixture
#' component densities for univariate constrained mixture of generalized normal distribution (CMGND) models.
#' It visually represents how the different components of the mixture model contribute to the overall density.
#'
#' @param x A numeric vector representing the observed data points.
#' @param parameters A matrix or data.frame containing the parameters of the CMGND model.
#' Alternatively, this can be an object returned from the `cmgnd()` function, representing an estimated CMGND model.
#'
#' @return A plot illustrating the marginal density along with the densities of the individual mixture components for the given data `x`.
#'
#' @details
#' The function plots the overall (marginal) density curve for the CMGND model, as well as the density curves of each mixture component.
#' This visualization helps in understanding how each component contributes to the model and provides insights into the data distribution.
#'
#' @seealso `cmgnd()` for estimating the model parameters.
#'
#' @export
plot_cmgnd <- function(x, parameters) {
  if (any(class(parameters) == "list")) {
    parameters <- parameters$parameters
  } else {
    parameters <- parameters
  }
  n <- length(x)
  K <- nrow(parameters)
  ind1 <- components <- list()
  for (k in 1:K) {
    ind1[[k]] <- rep(k, n)
    components[[k]] <- gnorm::dgnorm(x, parameters[k, 2], parameters[k, 3], parameters[k, 4]) * parameters[k, 1]
  }
  ind1 <- rapply(ind1, as.character, how = "replace")
  ind1[[K + 1]] <- rep("Marginal", n)
  components[[K + 1]] <- dcmgnd(x, parameters)
  components <- purrr::list_c(components)
  ind1 <- purrr::list_c(ind1)
  df <- data.frame(x, components, ind1)
  df$ind1 <- factor(df$ind1, levels = c("Marginal", as.character(c(1:K))))
  plot.cmgnd <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = components, group = ind1)) +
    ggplot2::geom_line(ggplot2::aes(linetype = ind1, color = ind1)) +
    ggplot2::labs(x = "Data", y = "Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
  return(plot.cmgnd)
}

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
#' @param grid A numeric vector specifying the step-lengths that determine the largest increase in the \eqn{Q(\theta, \theta(m-1))}
#' function for the shape parameter. The default value is 1, but based on simulations, 0.2 is a suitable value.
#' @param inversenu A logical value indicating if the algorithm should use an adaptive step size for the shape parameter.
#' Default is \code{TRUE} to control the degeneracy of the shape parameter.
#' @param scale A logical value indicating whether the function should scale the data. Default is \code{TRUE}.
#' @param eps A numeric value specifying the tolerance level of the ECM algorithm.
#' @param maxit An integer specifying the maximum number of iterations.
#' @param verbose A logical value indicating whether to display running output. Default is \code{TRUE}.
#' @param sigbound A numeric vector of length two specifying the lower and upper bounds for resetting the sigma estimates.
#' Default value is \code{c(.01,5)}.
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
#' @references Bazi, Y., Bruzzone, L., and Melgani, F. (2006). Image thresholding based on the em algorithm and the generalized gaussian distribution. Pattern Recognition, 40(2):619–634.
#' @references Duttilo, P. (2024). Modelling finacial returns with mixture of generalized normal distributions. Phd thesis, University "G. d’Annunzio" of Chieti-Pescara, Pescara, IT, pp 1-166. Available
#' at \url{https://drive.google.com/file/d/16whH1O4pVeGu_VY2sN_SPDdbTZeE1jNI/view?usp=sharing}
#' @references Duttilo, P., Gattone, S.A., and Iannone, B. (2023). Mixtures of generalized normal
#' distributions and EGARCH models to analyse returns and volatility of ESG and traditional investments. AStA Adv Stat Anal, pp. 1-21
#' \url{https://doi.org/10.1007/s10182-023-00487-7}
#' @references Duttilo, P., Gattone, S.A., and Kume, A. (2023). Mixtures of generalized normal
#' distributions with constraints. In: “Programme and Abstracts 25th International
#' Conference on Computational Statistics (COMPSTAT 2023)”, IASC, pp. 21, ISBN
#' 9789073592414
#' @references Duttilo, P., Kume A., and Gattone, S.A. (2023). Constrained Mixtures of Generalized
#' Normal Distributions. In: “SEAS IN Book of short papers 2023”, Pearson, pp.
#' 611-616, ISBN 9788891935618AAVV
#' @references Wen, L., Qiu, Y., Wang, M., Yin, J., and Chen, P. (2022). Numerical characteristics and
#' parameter estimation of finite mixed generalized normal distribution. Communications in
#' Statistics - Simulation and Computation, 51(7):3596–3620
#' @export
cmgnd <- function(x, K = 2, Cmu = rep(0, K), Csigma = rep(0, K), Cnu = rep(0, K), nstart = 10,
                  theta = FALSE, nustart = rep(2, K), nustartype = "random", gauss = FALSE,
                  laplace = FALSE, grid = 1, inversenu = TRUE, scale = FALSE,
                  eps = 10^-4, maxit = 999, verbose = TRUE, sigbound = c(.01, 5)) {
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

  sigdata <- stats::sd(x)
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
      estold <- matrix(c(pi_new, mu_new, sigma_new), nrow = K)

      # start cycle while
      check <- rep(0, 2)
      while (sum(check) < 2) {
        it <- it + 1
        # update of the of the location parameter
        ls_ll_old_mu <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
        if (any(Cmu == 1) & any(Cmu == 0)) {
          ind_cr <- which(Cmu == 0)
          mu_update_j <- NA
          for (k in ind_cr) {
            mu_update_j <- ls_mu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, grid, it)
          }
          mu_update_s <- ls_mu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, Cmu)[1]
          mu_new <- (mu_update_s * Cmu) + (mu_update_j * +(!Cmu))
        } else if (all(Cmu == 1)) {
          mu_new <- ls_mu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, Cmu)
        } else {
          for (k in 1:K) {
            mu_new[k] <- ls_mu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, grid, it)
          }
        }
        ls_ll_new_mu <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
        if ((ls_ll_new_mu - ls_ll_old_mu) < eps) {
          check[1] <- 1
        }
        if ((ls_ll_new_mu - ls_ll_old_mu) > eps) {
          check[1] <- 0
        }
        # update of the scale parameter
        ls_ll_old_sigma <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
        if (any(Csigma == 1) & any(Csigma == 0)) {
          ind_cr <- which(Csigma == 0)
          sigma_update_j <- NA
          for (k in ind_cr) {
            sigma_update_j <- sigma_j(x, res[, k], mu_new[k], nu_new[k], sigdata, sigbound)
          }
          sigma_update_s <- ls_sigma_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, sigdata, sigbound, Csigma)[1]
          sigma_new <- (sigma_update_s * Csigma) + (sigma_update_j * +(!Csigma))
        } else if (all(Csigma == 1)) {
          sigma_update_s <- ls_sigma_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, grid, sigdata, sigbound, Csigma)[1]
          sigma_new <- rep(sigma_update_s, K)
        } else {
          for (k in 1:K) {
            sigma_new[k] <- sigma_j(x, cbind(res[, k]), mu_new[k], nu_new[k], sigdata, sigbound)
          }
        }
        ls_ll_new_sigma <- log.likelihood(x, rbind(pi_new), rbind(mu_new), rbind(sigma_new), rbind(nu_new))
        if ((ls_ll_new_sigma - ls_ll_old_sigma) < eps) {
          check[2] <- 1
        }
        if ((ls_ll_new_sigma - ls_ll_old_sigma) > eps) {
          check[2] <- 0
        }
        # update of the shape parameter
        if (gauss == TRUE) {
          for (k in 1:K) {
            nu_new[k] <- 2
          }
        } else if (laplace == TRUE) {
          for (k in 1:K) {
            nu_new[k] <- 1
          }
        } else {
          if (any(Cnu == 1) & any(Cnu == 0)) {
            ind_cr <- which(Cnu == 0)
            nu_update_j <- NA
            for (k in ind_cr) {
              nu_update_j <- ls_nu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, inversenu, grid, it)
            }
            nu_update_s <- ls_nu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, inversenu, grid, Cnu)[1]
            nu_new <- (nu_update_s * Cnu) + (nu_update_j * +(!Cnu))
          } else if (all(Cnu == 1)) {
            nu_new <- ls_nu_s(x, pi_new, mu_new, sigma_new, nu_new, K, res, inversenu, grid, Cnu)
          } else {
            for (k in 1:K) {
              nu_new[k] <- ls_nu_j(x, pi_new, mu_new, sigma_new, nu_new, res, k, K, inversenu, grid, it)
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
        if (it > maxit) {
          # print("max iteration")
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
  mu_new <- c(op$mu)
  sigma_new <- c(op$sigma)
  nu_new <- c(op$nu)
  pi_new <- c(op$pi)
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
      scalesd = attr(x, "scaled:scale"), grid = grid, inversenu = inversenu
    )
  } else {
    info <- list(
      constraints = constraints, gnd = gnd, gauss = gauss, laplace = laplace,
      scalemeans = NULL, scalesd = NULL, grid = grid, inversenu = inversenu
    )
  }

  return(list(
    ll = ll_optim, nobs = N, npar = npar, parameters = parameters,
    ic = ic, res = res, clus = clus, op_it = op_it, cputime = time, info = info
  ))
}
