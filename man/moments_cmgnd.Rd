% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{moments_cmgnd}
\alias{moments_cmgnd}
\title{Compute the First Four Moments of the CMGND Marginal Distribution}
\usage{
moments_cmgnd(parameters)
}
\arguments{
\item{parameters}{A matrix or data frame where each row corresponds to a component
of the mixture. Columns must be ordered as follows:
\describe{
  \item{1}{Mixing proportions \eqn{\pi_k}}
  \item{2}{Component means \eqn{\mu_k}}
  \item{3}{Scale parameters \eqn{\sigma_k}}
  \item{4}{Shape parameters \eqn{\nu_k}}
}}
}
\value{
A named list with the following elements:
\describe{
  \item{mean}{The marginal mean of the CMGND distribution}
  \item{var}{The marginal variance}
  \item{skew}{The marginal skewness}
  \item{kur}{The marginal kurtosis}
}
}
\description{
Computes the mean, variance, skewness, and kurtosis of the marginal distribution
of a univariate Constrained Mixture of Generalized Normal Distributions (CMGND) model,
given the model parameters.
}
\details{
The function assumes that the parameters define a valid CMGND model and uses
analytical expressions to compute the first four moments of the marginal distribution.
The shape parameter \eqn{\nu_k} governs the kurtosis of each component.
}
\seealso{
\code{\link{cmgnd}} for estimating CMGND model parameters.
}
