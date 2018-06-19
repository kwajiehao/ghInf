#' centered_precgen3
#'
#' This function generates a precision matrix for a centered (every node in level 1 has the same number of children nodes)
#' 3-level Gaussian hierarchical model (level 0 being the root, and level 4 being the observations) given the variances
#' for each level tau, tau_a, tau_b, tau_c, and sigma_2. It takes advantage of the sparse structure of the precision matrix for centered 2-level gaussian hierarchical
#' models and returns only the non-zero entries of the precision matrix together with its corresponding indices.
#'
#' Assumptions:
#'   - variances are constant for parameters within the same level
#'   - all observations y_ijk are equal to 0
#'   - the mean of the root parameter B is mu = 0
#'
#' @param i number of nodes at level 1
#' @param j number of children nodes in level 2 per node at level 1
#' @param k number of children nodes in level 3 per node at level 2
#' @param flat_prior determines whether to use the density with flat prior
#' @param tau variance of the root (level 0)
#' @param tau_a variance for parameters in level 1
#' @param tau_b variance for parameters in level 2
#' @param tau_c variance for parameters in level 3
#' @param sigma_2 variance of the observations
#'
#' @return list of non-zero entries of the precision matrix together with its corresponding indices
#' @export
#'
#' @examples
#' i <- 2
#' j <- 3
#' k <- 2
#' centered_precgen3(i = i, j = j, k = k)
centered_precgen3 <- function(i, j, k, flat_prior = TRUE, tau = 1, tau_a = 1, tau_b = 1, tau_c = 1, sigma_2 = 1){

  # this function takes advantage of the sparse structure of the precision matrix for 3-level gaussian hierarchical
  # models and returns only the non-zero entries of the precision matrix together with its corresponding indices

  # tau represents the variance of the root (level 0)
  # tau_a is the variance for parameters in level 1
  # tau_b is the variance for parameters in level 2
  # tau_c is the variance for paramters in level 3
  # sigma_2 is the variance of the observations

  # i represents the number of nodes at level 1
  # j represents the number of children nodes in level 2 per node at level 1
  # k represents the number of children nodes in level 3 per node at level 2

  # ijk is the total number of nodes at level 3
  ijk <- i*j*k

  # ij is the total number of nodes at level 2
  ij <- i*j

  # store the variances in a list
  # possible extensions:
  #     - accounting for possibility that each node has a different variance
  n <- 1 + i + ij + ijk + 1
  var_params <- rep(0, n)
  var_params[1] <- tau
  var_params[(1+1):(1+i)] <- tau_a
  var_params[(1+i+1):(1+i+ij)] <- tau_b
  var_params[(1+i+ij+1):(1+i+ij+ijk)] <- tau_c
  var_params[n] <- sigma_2
  inv_var_params <- 1 / var_params

  # store indices in lists
  num_nondiag <- (ijk + ij + i) # number of non-diagonal elements
  nondiag_i <- rep(0, num_nondiag)
  nondiag_j <- nondiag_i
  diag <- seq(1, (n-1)) # we store the indices of the (n-1) diagonal entries separately

  # fill in the precision matrix according to derivations in two steps:
  # 1. diagonals

  diag_entries <- rep(1000, (n-1))

  for (m in 1:ijk){ # B_ijk or level 3 terms
    diag_entries[m] <- inv_var_params[n] + inv_var_params[(1+i+ij+m)]
  }

  for (m in 1:ij){ # B_ij or level 2 terms
    diag_entries[m+ijk] <- inv_var_params[1+i+m] +
      sum(inv_var_params[(1 + i + ij + (m-1)*k + 1): (1 + i + ij + (m-1)*k + k)])
  }

  for (m in 1:i){ # B_i or level 1 terms
    diag_entries[m+ijk+ij] <- inv_var_params[(1+m)] +
      sum(inv_var_params[(1 + i + ((m-1) * j) + 1):(1 + i + ((m-1) * j) + j)])
  }

  if (flat_prior==FALSE){ # B or level 0 term
    diag_entries[(n-1)] <- inv_var_params[1] + sum(inv_var_params[2:(1+i)])
  } else{
    diag_entries[(n-1)] <- sum(inv_var_params[2:(1+i)])
  }

  # 2. off-diagonals

  nondiags <- rep(1001, num_nondiag) # nondiagonal entries

  for (m in 1:i){
    for (l in 1:j){
      for (a in 1:k){ # according to derived formula for posterior inference
        nondiags[((((m-1)*j) + (l-1))*k + a)] <- -inv_var_params[1+i+ij+((((m-1)*j) + (l-1))*k + a)]
        nondiag_i[((((m-1)*j) + (l-1))*k + a)] <- (ijk + (m-1)*j + l)
        nondiag_j[((((m-1)*j) + (l-1))*k + a)] <- ((((m-1)*j) + (l-1))*k + a)
      }
      nondiags[((m-1)*j + l + ijk)] <- -inv_var_params[1+i+((m-1)*j + l)]
      nondiag_i[((m-1)*j + l + ijk)] <- (m+ijk+ij)
      nondiag_j[((m-1)*j + l + ijk)] <- ((m-1)*j + l + ijk)
    }
    nondiags[(ijk + ij + m)] <- - inv_var_params[(m+1)]
    nondiag_i[(ijk + ij +m)] <- (n-1)
    nondiag_j[(ijk + ij +m)] <- (m+ijk+ij)
  }

  # vector of indices and non-zero entries
  index_i <- c(diag, nondiag_i, nondiag_j)
  index_j <- c(diag, nondiag_j, nondiag_i)
  non_zero_entries <- c(diag_entries, nondiags, nondiags)

  return(list(indices_i = index_i, indices_j = index_j, entries = non_zero_entries, inv = inv_var_params))
}
