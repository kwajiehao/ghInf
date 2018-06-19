#' precgen2
#'
#' This function generates a precision matrix for a general 2-level Gaussian hierarchical model
#' (level 0 being the root, and level 3 being the observations) given the variances for each level tau, tau_a, tau_b and
#' sigma_2. It takes advantage of the sparse structure of the precision matrix for centered 2-level gaussian hierarchical
#' models and returns only the non-zero entries of the precision matrix together with its corresponding indices.
#'
#' Assumptions:
#'   - variances are constant for parameters within the same level
#'   - all observations y_ij are equal to 0
#'   - the mean of the root parameter B is mu = 0
#'
#' @param i number of nodes at level 1
#' @param J vector specifying the number of children nodes in level 2 per node at level 1
#' @param flat_prior determines whether to use the density with flat prior
#' @param tau variance of the root (level 0)
#' @param tau_a variance for parameters in level 1
#' @param tau_b variance for parameters in level 2
#' @param sigma_2 variance of the observations
#'
#' @return list of non-zero entries of the precision matrix together with its corresponding indices
#' @export
#'
#' @examples
#' i <- 2
#' J <- c(1,3)
#' precgen2(i = i, J = J)
precgen2 <- function(i, J, flat_prior = TRUE, tau = 1, tau_a = 1, tau_b = 1, sigma_2 = 1){

  # this function takes advantage of the sparse structure of the precision matrix for any general 2-level gaussian
  # hierarchical model and returns only the non-zero entries of the precision matrix together with its corresponding indices

  # tau represents the variance of the root (level 0)
  # tau_a is the variance for parameters in level 1
  # tau_b is the variance for parameters in level 2
  # sigma_2 is the variance of the observations

  # defensive programming
  #     - to ensure that submitted vector is of the right length
  #     - to ensure numeric input
  if (length(J)!=i){
    stop("vector of children nodes is not the right length!")
  }

  if (!is.numeric(J)){
    stop("submitted vector is of wrong type")
  }

  # i is a number specifying the number of children
  # J is a vector specifying the number of children for each node in the first level

  # n is the total number of parameters
  j <- sum(J)
  n <- j + i + 1

  # store the variances in a list
  # possible extensions:
  #     - accounting for possibility that each node has a different variance
  var_params <- rep(0, (n+1))
  var_params[1:j] <- tau_b
  var_params[(j+1):(j+i)] <- tau_a
  var_params[(j+i+1)] <- tau
  var_params[n+1] <- sigma_2
  inv_var_params <- 1 / var_params

  # store indices in lists
  diag <- seq(1, n)
  num_nondiag <- j + i # number of non-diagional elements in half of the matrix
  nondiag_x <- rep(0, num_nondiag)
  nondiag_y <- nondiag_x

  # fill in the precision matrix according to derivations in two steps:

  # 1. diagonals
  diag_entries <- rep(0, n)

  for (m in 1:j){ # B_ij or level 2 terms
    diag_entries[m] <- inv_var_params[n+1] + inv_var_params[m]
  }

  for (m in 1:i){ # B_i or level 1 terms
    diag_entries[(m+j)] <- J[m] * inv_var_params[(m)] + inv_var_params[(m+j)]
  }

  if (flat_prior==FALSE){ # B or level 0 term
    diag_entries[n] <- inv_var_params[n] + sum(inv_var_params[(j+1):(j+i)])
  } else{
    diag_entries[n] <- sum(inv_var_params[(j+1):(j+i)])
  }

  # 2. off-diagonals
  nondiags <- rep(0, num_nondiag) # store non-diagonal entries
  count <- 0

  for (b in 1:i){
    for (c in 1:J[b]){
      nondiags[count + c] <- - inv_var_params[count + c]
      nondiag_x[count + c] <- count + c
      nondiag_y[count + c] <- j + b
    }

    count <- count + J[b]
    # print(count)
  }

  for (a in 1:i){
    nondiags[j + a] <- - inv_var_params[j + a]
    nondiag_x[j + a] <- j + a
    nondiag_y[j + a] <- j + i + 1
  }


  # vector of indices and non-zero entries
  index_i <- c(diag, nondiag_x, nondiag_y)
  index_j <- c(diag, nondiag_y, nondiag_x)
  non_zero_entries <- c(diag_entries, nondiags, nondiags)

  return(list(indices_i = index_i, indices_j = index_j, entries = non_zero_entries))

}
