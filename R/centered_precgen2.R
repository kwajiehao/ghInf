#' centered_precgen2
#'
#' This function generates a precision matrix for a centered (every node in level 1 has the same number of children nodes)
#' 2-level Gaussian hierarchical model (level 0 being the root, and level 3 being the observations) given the variances
#' for each level tau, tau_a, tau_b and sigma_2. It takes advantage of the sparse structure of the precision matrix for centered 2-level gaussian hierarchical
#' models and returns only the non-zero entries of the precision matrix together with its corresponding indices.
#'
#' Assumptions:
#'   - variances are constant for parameters within the same level
#'   - all observations y_ij are equal to 0
#'   - the mean of the root parameter B is mu = 0
#'
#' @param i number of nodes at level 1
#' @param j number of children nodes in level 2 per node at level 1
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
#' j <- 3
#' centered_precgen2(i = i, j = j)
centered_precgen2 <- function(i, j, flat_prior = TRUE, tau = 1, tau_a = 1, tau_b = 1, sigma_2 = 1){

  # k is the total number of nodes at level 2
  k <- i*j

  # store the variances in a list
  # possible extensions:
  #     - accounting for possibility that each node has a different variance
  n <- i + k + 2
  var_params <- rep(0, n)
  var_params[1] <- tau
  var_params[2:(2+i-1)] <- tau_a
  var_params[(2+i):(2+i+k-1)] <- tau_b
  var_params[n] <- sigma_2
  inv_var_params <- 1 / var_params

  # store indices in lists
  num_nondiag <- (k + i) # number of nondiagonal elements
  nondiag_i <- rep(0, num_nondiag)
  nondiag_j <- nondiag_i
  diag <- seq(1, (n-1)) # we store the indices of the (n-1) diagonal entries separately

  # fill in the precision matrix in two steps:
  # 1. diagonals

  diag_entries <- rep(0, (n-1)) # diagonal entries

  for (m in 1:k){ # B_ij or level 2 terms
    diag_entries[m] <- inv_var_params[n] + inv_var_params[(1+i+m)]
  }

  for (m in 1:i){ # B_i or level 1 terms
    diag_entries[m+k] <- inv_var_params[(m+1)] +
      sum(inv_var_params[(1 + i + ((m-1) * j) + 1):(1 + i + ((m-1) * j) + j)])
  }

  if (flat_prior==FALSE){ # B or level 0 term
    diag_entries[(n-1)] <- inv_var_params[1] + sum(inv_var_params[2:(1+i)])
  } else{
    diag_entries[(n-1)] <- sum(inv_var_params[2:(1+i)])
  }

  # 2. off-diagonals

  nondiags <- rep(0, num_nondiag) # nondiagonal entries

  for (m in 1:i){
    for (l in 1:j){ # according to derived formula for posterior inference
      nondiags[((m-1)*j + l)] <- - inv_var_params[(1+i+(l+(m-1)*j))]
      nondiag_i[((m-1)*j + l)] <- (m+k)
      nondiag_j[((m-1)*j + l)] <- ((m-1)*j + l)
    }
    nondiags[(k+m)] <- - inv_var_params[m+1]
    nondiag_i[(k+m)] <- (n-1)
    nondiag_j[(k+m)] <- (m+k)
  }

  # vector of indices and non-zero entries
  index_i <- c(diag, nondiag_i, nondiag_j)
  index_j <- c(diag, nondiag_j, nondiag_i)
  non_zero_entries <- c(diag_entries, nondiags, nondiags)

  return(list(indices_i = index_i, indices_j = index_j, entries = non_zero_entries))
}

