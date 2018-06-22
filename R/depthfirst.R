#' depthfirst
#'
#' This function returns the sequence which permutes the rows and columns of the precision matrix for a centered 2- or 3-level
#' gaussian hierarchical model which has been labelled lexicographically (that is, for a 2-level model with i nodes in the
#' first level and j child nodes per node in the first level so that the columns and rows of the precision matrix are
#' labelled B_11, ..., B_ij, B_1, ..., B_i, B) into a depth-first labelling (for example, for a 2-level model, B_11, B_12,
#' ..., B_1J, B_1, B_21, ..., B_2J, B_2, ...)
#'
#' @param i number of nodes at level 1
#' @param j number of children nodes in level 2 per node at level 1
#' @param k number of children nodes in level 3 per node at level 2
#' @param levels the number of levels in chosen model
#'
#' @return a sequence of indices which permutes the original matrix
#' @export
#'
#' @examples
#' i <- 2
#' j <- 3
#' k <- 2
#' levels <- 3
#' depthfirst(i = i, j = j, k = k, levels = levels)
depthfirst <- function(i, j, k = NA, levels){

  # levels represents the number of levels of the gaussian hierarchical model in question. for now, this function only
  # works with centered 2- and 3-level parametrizations
  if (levels > 3 | levels < 2){
    stop("number of levels chosen is not supported!")
  }

  # i represents the number of nodes at level 1
  # j represents the number of children nodes in level 2 per node at level 1
  # k represents the number of children nodes in level 3 per node at level 2


  # the trick is to keep each index and add numbers so that they become the right index
  if (levels == 2){

    n <- i*j + i + 1# total number of params
    o <- seq(1, n) # ascending order
    increment <- seq(0, (i-1))

    for (k in 1:i){
      o[(j * increment[k] + 1):(j * increment[k] + j)] <- o[(j * increment[k] + 1):(j * increment[k] + j)] + increment[k]
      o[((i*j) + k)] <- j*(increment[k]+1) + k
    }
  }

  if (levels == 3){

    n <- i*j*k + i*j + i + 1 # total number of params
    o <- seq(1, n) # ascending order
    increment <- seq(0, (i-1))

    for (m in 1:i){
      for (l in 1:j){
        for (h in 1:k){
          # print(paste0(m,l,h))
          # print((m-1)*j*k + (l-1)*k + h)
          o[(m-1)*j*k + (l-1)*k + h] <- (m-1)*(j*k + j + 1) + (l-1)*(k+1) + h

        }
        # print(i*j*k + (m-1)*j +l)
        o[i*j*k + (m-1)*j +l] <- (m-1)*(j*k + j + 1) + (l)*(k+1)
      }
      # print((i*j*k + i*j + m) + (j*k + j) + 1)
      o[i*j*k + i*j + m] <- (m-1)*(j*k + j + 1) + (j*k + j) + 1
    }
  }

  return(order(o)[seq(1,n)])
}

# # test to show how depthfirst works
# I <- 2
# J <- 3
# K <- 2
# n <- I*J*K + I*J + I + 1
# depthfirst_test <- centered_precgen3(i=I, j=J, k=K, flat_prior = TRUE, tau = 1, tau_a = 1, tau_b = 1, tau_c = 1, sigma_2 = 1)
#
# test_mat <- matrix(0, (I*J*K + I*J + I + 1), (I*J*K + I*J + I + 1))
# test_mat[cbind(depthfirst_test$indices_i, depthfirst_test$indices_j)] <- depthfirst_test$entries
#
# # without permutation, I use a lexicographic ordering
# cols <- rep('B', n)
#
# for (i in 1:I){
#   for (j in 1:J){
#     for (k in 1:K){
#       cols[ (i-1)*K*J + K*(j-1) + k] <- paste0('B', i, j, k)
#     }
#     cols[I*J*K + (i-1)*J + j] <- paste0('B', i, j)
#   }
#   cols[I*J*K + (I*J) + i] <- paste0('B', i)
# }
#
# rows <- cols
# colnames(test_mat) <- cols
# rownames(test_mat) <- rows
#
# og <- depthfirst(i = 2, j = 3, k = 2, levels = 3)
# og
# length(og)
#
# ## compare original matrix with permuted matrix
# test_mat
# test_mat[og,og]
