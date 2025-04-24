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
#' @param model A character indicating the model type name.
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
plot_cmgnd <- function(x, parameters, model) {
  if (any(class(parameters) == "list")) {
    parameters <- parameters$parameters
  } else {
    parameters <- parameters
  }

  # Calculate the range of the x values
  x_range <- range(x)

  # Extend the range by 5 units in both directions
  x_extended_inf <- seq(x_range[1] - 5, x_range[1], length.out = 100)  # Extend lower bound
  x_extended_sup <- seq(x_range[2], x_range[2] + 5, length.out = 100)  # Extend upper bound

  # Combine the original x with the extended values
  x_full_range <- c(x_extended_inf, x, x_extended_sup)

  n <- length(x_full_range)
  K <- nrow(parameters)
  ind1 <- components <- list()

  # Compute the components (density functions) over the combined x values
  for (k in 1:K) {
    ind1[[k]] <- rep(k, n)
    components[[k]] <- gnorm::dgnorm(x_full_range, parameters[k, 2], parameters[k, 3], parameters[k, 4]) * parameters[k, 1]
  }

  ind1 <- rapply(ind1, as.character, how = "replace")
  ind1[[K+1]] <- rep("Marginal", n)
  components[[K+1]] <- dcmgnd(x_full_range, parameters)
  components <- purrr::list_c(components)
  ind1 <- purrr::list_c(ind1)

  # Now create the data frame with combined x values for plotting
  df <- data.frame(x = x_full_range, components, ind1)
  df$ind1 <- factor(df$ind1, levels = c("Marginal", as.character(c(1:K))))

  # Plot the density functions over the combined x values
  plot.cmgnd <- ggplot2::ggplot(data = df, ggplot2::aes(x = x, y = components, group = ind1)) +
    ggplot2::geom_line(ggplot2::aes(linetype = ind1, color = ind1)) +
    ggplot2::labs(x = "Data", y = "Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())+
    ggplot2::ggtitle(model)
  return(plot.cmgnd)
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
#' @param bins Number of bins. Defaults to 80.
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

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("..density.."))
}


hist_cmgnd <- function(x, parameters, bins = 80) {
  if (any(class(parameters) == "list")) {
    parameters <- parameters$parameters
  } else {
    parameters <- parameters
  }

  # Calculate the range of the x values
  x_range <- range(x)

  # Extend the range by 5 units in both directions
  x_extended_inf <- seq(x_range[1] - 5, x_range[1], length.out = 100)  # Extend lower bound
  x_extended_sup <- seq(x_range[2], x_range[2] + 5, length.out = 100)  # Extend upper bound

  # Combine the original x with the extended values
  x_full_range <- c(x_extended_inf, x, x_extended_sup)

  n <- length(x_full_range)
  K <- nrow(parameters)
  ind1 <- components <- list()

  # Compute the components (density functions) over the combined x values
  for (k in 1:K) {
    ind1[[k]] <- rep(k, n)
    components[[k]] <- gnorm::dgnorm(x_full_range, parameters[k, 2], parameters[k, 3], parameters[k, 4]) * parameters[k, 1]
  }

  ind1 <- rapply(ind1, as.character, how = "replace")
  ind1[[K+1]] <- rep("Marginal", n)
  components[[K+1]] <- dcmgnd(x_full_range, parameters)
  components <- purrr::list_c(components)
  ind1 <- purrr::list_c(ind1)

  # Now create the data frame with combined x values for plotting
  df <- data.frame(x = x_full_range, components, ind1)
  df$ind1 <- factor(df$ind1, levels = c("Marginal", as.character(c(1:K))))

  # Plot the density functions over the combined x values with a semi-transparent histogram
  plot.cmgnd <- ggplot2::ggplot() +
    ggplot2::geom_histogram(data = data.frame(x = x), ggplot2::aes(x = x, y = ..density..),
                            bins = bins, fill = "grey70", color = "grey50", alpha = 0.4) +
    ggplot2::geom_line(data = df, ggplot2::aes(x = x, y = components, linetype = ind1, color = ind1), size = 1) +
    ggplot2::labs(x = "Data", y = "Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())

  return(plot.cmgnd)
}


