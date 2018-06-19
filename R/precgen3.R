#' precgen3
#'
#' This function generates a precision matrix for a general 3-level Gaussian hierarchical model (level 0 being the root,
#'and level 4 being the observations) given the variances for each level tau, tau_a, tau_b, tau_c, and sigma_2. It takes
#' advantage of the sparse structure of the precision matrix for centered 2-level gaussian hierarchical
#' models and returns only the non-zero entries of the precision matrix together with its corresponding indices.
#'
#' Assumptions:
#'   - variances are constant for parameters within the same level
#'   - all observations y_ijk are equal to 0
#'   - the mean of the root parameter B is mu = 0
#'
#' @param i number of nodes at level 1
#' @param J vector specifying the number of children nodes in level 2 per node at level 1
#' @param K vector specifying the number of children nodes in level 3 per node at level 2
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
#' J <- c(1,3)
#' K <- c(1,2,3,4)
#' precgen3(i = i, J = J, K = K)
precgen3 <- function(i, J, K, flat_prior = TRUE, tau = 1, tau_a = 1, tau_b = 1, tau_c = 1, sigma_2 = 1){

  # this function takes advantage of the sparse structure of the precision matrix for any general 3-level gaussian
  # hierarchical model and returns only the non-zero entries of the precision matrix together with its corresponding indices

  # tau represents the variance of the root (level 0)
  # tau_a is the variance for parameters in level 1
  # tau_b is the variance for parameters in level 2
  # tau_c is the variance for parameters in level 3
  # sigma_2 is the variance of the observations

  # i is a number specifying the number of children
  # J is a vector specifying the number of children for each node in the first level
  # K is a vector specifying the number of children for each node in the second level

  # n is the total number of parameters
  j <- sum(J)
  k <- sum(K)
  n <- k + j + i + 1

  len_j <- length(J)
  len_k <- length(K)

  # defensive programming
  #     - to ensure that submitted vectors are of the right length
  #     - to ensure numeric input
  if (len_j !=i){
    stop("vector of children nodes is not the right length for level 2!")
  }

  if (len_k !=j){
    stop("vector of children nodes is not the right length for level 3!")
  }

  if (!is.numeric(J)){
    stop("submitted vector for level 2 is of wrong type")
  }

  if (!is.numeric(K)){
    stop("submitted vector for level 3 is of wrong type")
  }


  # store the variances in a list
  # possible extensions:
  #     - accounting for possibility that each node has a different variance
  var_params <- rep(0, (n+1))
  var_params[1:k] <- tau_c
  var_params[(k+1):(k+j)] <- tau_b
  var_params[(k+j+1):(k+j+i)] <- tau_a
  var_params[(k+j+i+1)] <- tau
  var_params[n+1] <- sigma_2
  inv_var_params <- 1 / var_params

  # store indices in lists
  diag <- seq(1, n)
  num_nondiag <- k + j + i # number of non-diag in half of the matrix
  nondiag_x <- rep(0, num_nondiag)
  nondiag_y <- nondiag_x

  # fill in the precision matrix according to derivations in two steps:

  # 1. diagonals
  diag_entries <- rep(0, n)

  for (m in 1:k){ # B_ijk or level 3 terms
    diag_entries[m] <- inv_var_params[n+1] + inv_var_params[m]
  }

  for (m in 1:j){ # B_ij or level 2 terms
    diag_entries[(m+k)] <- K[m] * inv_var_params[m] + inv_var_params[m+k]
  }

  for (m in 1:i){ # B_i or level 1 terms
    diag_entries[(m+k+j)] <- J[m] * inv_var_params[k+m] + inv_var_params[(m+k+j)]
  }

  if (flat_prior==FALSE){ # B or level 0 term
    diag_entries[n] <- inv_var_params[n] + sum(inv_var_params[(k+j+1):(k+j+i)])
  } else{
    diag_entries[n] <- sum(inv_var_params[(k+j+1):(k+j+i)])
  }

  # 2. off-diagonals
  nondiags <- rep(0, num_nondiag)
  count1 <- k
  count2 <- 0

  for (d in 1:j){
    # print(c(d, K[d]))
    for (e in 1:K[d]){
      nondiags[count2 + e] <- - inv_var_params[count2 + e]
      nondiag_x[count2 + e] <- count2 + e
      nondiag_y[count2 + e] <- k + d
    }

    count2 <- count2 + K[d]
    # print(count2)
  }


  for (b in 1:i){
    for (c in 1:J[b]){
      nondiags[count1 + c] <- - inv_var_params[count1 + c]
      nondiag_x[count1 + c] <- count1 + c
      nondiag_y[count1 + c] <- k + j + b
    }

    count1 <- count1 + J[b]
    # print(count1)
  }

  for (a in 1:i){
    nondiags[k + j + a] <- - inv_var_params[k + j + a]
    nondiag_x[k + j + a] <- k + j + a
    nondiag_y[k + j + a] <- k + j + i + 1
  }


  # vector of indices and non-zero entries
  index_i <- c(diag, nondiag_x, nondiag_y)
  index_j <- c(diag, nondiag_y, nondiag_x)
  non_zero_entries <- c(diag_entries, nondiags, nondiags)

  return(list(indices_i = index_i, indices_j = index_j, entries = non_zero_entries))

}
