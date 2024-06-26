% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ou.lik.R
\name{ou.lik}
\alias{ou.lik}
\title{OU likelihood function}
\usage{
ou.lik(cache, p, method, c_t = NULL, c_a = NULL)
}
\arguments{
\item{cache}{A cache object created using epcTools.}

\item{p}{Numeric vector of starting parameters to be loaded corresponding to theta, sig2, and alpha. If NULL start.searcher function will be run instead.}

\item{method}{A string corresponding to the optimx::optimx algorithm to be run, recommended 'Nelder-Mead' for first search followed by 'BFGS'.}

\item{c_t}{A custom theta function under an EPC process, NULL for OU models.}

\item{c_a}{A custom alpha function under an EPC process, NULL for OU models.}
}
\value{
A table of parameters (theta, sig2, and alpha) with a likelihood, execution time, model convergence parameters from optimx::optimx (ignore), and AIC.
}
\description{
Performs a single maximum likelihood search under an Ornstein-Uhlenbeck process
}
