% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{sim_cmgnd}
\alias{sim_cmgnd}
\title{sim_cmgnd: Function to Simulate Univariate Constrained Mixtures of Generalized Normal Distributions}
\usage{
sim_cmgnd(
  n = 1000,
  pi = rep(0.5, 2),
  mu = c(1, 5),
  sigma = c(1, 1),
  nu = c(2, 2)
)
}
\arguments{
\item{n}{A numeric value indicating the total number of observations to simulate.}

\item{pi}{A numeric vector of the mixture weights \eqn{\pi_k}.}

\item{mu}{A numeric vector of the location parameter \eqn{\mu_k}.}

\item{sigma}{A numeric vector of the scale parameter \eqn{\sigma_k}.}

\item{nu}{A numeric vector of the shape parameter \eqn{\nu_k}.}
}
\value{
\item{\code{sim_data}}{The simulated data.}
\item{\code{sim_clus}}{The cluster indication of simulated data.}
}
\description{
Simulate univariate constrained mixture of generalized normal distribution models.
Remeber to set the set.seed() before the function sim_cmgnd().
}
