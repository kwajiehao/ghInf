% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centered_precgen2.R
\name{centered_precgen2}
\alias{centered_precgen2}
\title{centered_precgen2}
\usage{
centered_precgen2(i, j, flat_prior = TRUE, tau = 1, tau_a = 1,
  tau_b = 1, sigma_2 = 1)
}
\arguments{
\item{i}{number of nodes at level 1}

\item{j}{number of children nodes in level 2 per node at level 1}

\item{flat_prior}{determines whether to use the density with flat prior}

\item{tau}{variance of the root (level 0)}

\item{tau_a}{variance for parameters in level 1}

\item{tau_b}{variance for parameters in level 2}

\item{sigma_2}{variance of the observations}
}
\value{
list of non-zero entries of the precision matrix together with its corresponding indices
}
\description{
This function generates a precision matrix for a centered (every node in level 1 has the same number of children nodes)
2-level Gaussian hierarchical model (level 0 being the root, and level 3 being the observations) given the variances
for each level tau, tau_a, tau_b and sigma_2. It takes advantage of the sparse structure of the precision matrix for centered 2-level gaussian hierarchical
models and returns only the non-zero entries of the precision matrix together with its corresponding indices.
}
\details{
Assumptions:
  - variances are constant for parameters within the same level
  - all observations y_ij are equal to 0
  - the mean of the root parameter B is mu = 0
}
\examples{
i <- 2
j <- 3
centered_precgen2(i = i, j = j)
}
