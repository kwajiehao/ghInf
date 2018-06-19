#' centered_gibbs3
#'
#' Gibbs sampler for centered 3-level Gaussian hierarchical model according to derived full conditionals
#'
#' Assumptions:
#'   - variances are constant for parameters within the same level
#'   - all observations y_ijk are equal to 0
#'   - the mean of the root parameter B is mu = 0
#'   - assume a naive sampler where the variances are not updated
#'
#' @param i number of nodes at level 1
#' @param j number of children nodes in level 2 per node at level 1
#' @param k number of children nodes in level 3 per node at level 2
#' @param ndraws number of samples the user wants to have
#' @param burnin number of samples to throw away at the start as the gibbs sampler warms up
#' @param flat_prior determines whether to use the density with flat prior
#' @param tau variance of the root (level 0)
#' @param tau_a variance for parameters in level 1
#' @param tau_b variance for parameters in level 2
#' @param tau_c variance for parameters in level 3
#' @param sigma_2 variance of the observations
#'
#' @return list of means and the samples
#' @export
#'
#' @examples
#' i <- 2
#' j <- 3
#' k <- 2
#' ndraws <- 10000
#' burnin <- 1000
#' centered_gibbs3(i = i, j = j, k = k, ndraws = ndraws, burnin = burnin)
centered_gibbs3 <- function(i, j, k, ndraws, burnin, flat_prior = TRUE, tau = 1, tau_a = 1, tau_b = 1, tau_c = 1, sigma_2 = 1){

  # ndraws is the number of samples the user wants to have
  # burnin is the number of samples to throw away at the start as the gibbs sampler warms up

  # tau represents the variance of the root (level 0)
  # tau_a is the variance for parameters in level 1
  # tau_b is the variance for parameters in level 2
  # tau_c is the variance for parameters in level 3
  # sigma_2 is the variance of the observations

  # possible extensions:
  #     - note that for a truly independent sample to be obtained, only 1 in every m samples from the gibbs sampler
  #       should be retained where m is sufficiently large (since subsequent samples are highly correlated). this
  #       implementation does not take the correlation into account.


  # generate a matrix to store results and initialize at zero for all parameters
  #     - note that the labelling order is lexicographic: B_111, ..., B_ijk, B_11, ..., B_ij, B_1, ..., B_i, B
  n <- i*j*k + i*j + i + 1
  means <- matrix(0, n, (ndraws + burnin))
  results <- matrix(0, n, (ndraws + burnin))

  for (a in 2:(ndraws+burnin)){
    for (b in 1:i){
      for (c in 1:j){
        for (d in 1:k){
          # update B_ijk (level 3) according to previous values of B_ij
          means[((b-1)*j*k + (c-1)*k + d)] <- results[ (i*j*k + (b-1)*j + c) , (a-1)] * sigma_2 / (tau_c + sigma_2)
          results[((b-1)*j*k + (c-1)*k + d)] <- stats::rnorm(1, means[((b-1)*j*k + (c-1)*k + d)],
                                                      sqrt((tau_c*sigma_2)/(tau_c+sigma_2)))

        }
        # update of B_ij (level 2) according to new values of B_ijk and previous values of B_i
        means[(i*j*k + (b-1)*j + c),a] <- (tau_b * mean(results[((b-1)*j*k + (c-1)*k + 1) : ((b-1)*j*k + (c-1)*k + k), a])
                                           + (tau_c / k) * results[(i*j*k + i*j + b), (a-1)]) / (tau_b + (tau_c/k))
        results[(i*j*k + (b-1)*j + c),a] <- stats::rnorm(1, means[(i*j*k + (b-1)*j + c),a],
                                                  sqrt((tau_b*tau_c/k)/(tau_b+(tau_c/k))))
      }

      # update B_i (level 1) according to new values of B_ij and previous valeus of B
      means[(i*j*k + i*j + b),a] <- (tau_a * mean(results[(i*j*k + (b-1)*j + 1):(i*j*k + (b-1)*j + j),a])
                                     + results[n,(a-1)] * tau_b/j) / (tau_a + tau_b/j)
      results[(i*j*k + i*j + b),a] <-  stats::rnorm(1, means[(i*j*k + i*j + b),a] ,
                                             sqrt((tau_a * tau_b/j)/(tau_a + tau_b/j)))
    }


    # update of B (level 0) according to new values of B_i
    if (flat_prior == TRUE){
      means[n,a] <- mean(results[(i*j*k + i*j +1):(n-1),a])
      results[n,a] <- stats::rnorm(1, means[n,a], sqrt(tau_a/i))

    } else{
      means[n,a] <- mean(results[(i*j*k+i*j+1):(n-1),a]) * tau / (tau_a/i + tau)
      results[n,a] <- stats::rnorm(1, means[n,a], sqrt((tau_a * tau/i) / (tau_a/i + tau)))
    }

  }

  # return a list with the means and the samples
  list(means = means, samples = results)
}

