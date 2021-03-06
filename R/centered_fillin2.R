#' centered_fillin2
#'
#' This function returns the ratio of the number of zeros obtained from the cholesky decomposition of the
#' specified precision matrix over the optimal cholesky fill-in for a centered 2-level gaussian hierarchical model.
#' If the cholesky decomposition is optimal, the number returned is 1, and if it is sub-optimal the number returned
#' is greater than 1. The function also returns the time taken for the cholesky decomposition.
#'
#' A choice of 2 packages is provided: spam and Matrix.
#' A choice of 3 permute methods is provided: no permutation, random permutation, and depth-first permutation.
#'
#' @param method determines which linear algebra package to use
#' @param permute_method determines how the precision matrix is permuted before use with the chosen linear algebra package
#' @param i number of nodes at level 1
#' @param j number of children nodes in level 2 per node at level 1
#' @param permute determines whether or not to use default row permutation algorithm before matrix factorization to reduce fill-in
#' @param ... for other params used in centered_precgen2
#'
#' @return the fill-in ratio with respect to the optimal fill-in, and the time taken for the cholesky decomposition
#' @export
#'
#' @examples
#' i <- 2
#' j <- 3
#' centered_fillin2(method = 'matrix', permute_method = 'random', i = i, j = j)
centered_fillin2 <- function(method = 'both', permute_method = 'none', i, j, permute = TRUE,...){

  # sparse precision generation
  q <- centered_precgen2(i, j, ...)

  # permutation of sequence of row/columns according to user's choice of permute method
  o <- seq(1, (i*j + i + 1))

  if (permute_method=="none"){
    og <- o
  }

  if (permute_method=="random"){
    og <- sample(o)
  }

  if (permute_method=="depthfirst"){
    og <- depthfirst(i,j, levels = 2)
  }

  if (permute_method == "reverse"){
    og <- rev(o)
  }

  # number of variables ignoring the root variable
  d_minus_one <- i*j + i

  if (method == "spam"){
    # print(c(i,j))
    Q_spam <- spam::spam(0, (i*j+i+1), (i*j+i+1))
    Q_spam[cbind(q$indices_i, q$indices_j)] <- q$entries
    Q_spam <- Q_spam[og,og]

    # start time
    start <- Sys.time()

    # cholesky factor
    cholesky <- spam::chol(Q_spam, pivot = permute)

    # end time
    end <- Sys.time() - start

    # count zeros
    nonzero_count <- utils::tail(cholesky@rowpointers, 1) - 1 - cholesky@dimension[1]

    # return list of ratios and time taken
    return(list(diff = (nonzero_count / d_minus_one), time_taken = end, cky = cholesky, q=q))
  }

  if (method == "matrix"){
    # print(c(i,j))
    Q <- Matrix::sparseMatrix(q$indices_i, q$indices_j, x = q$entries)
    Q <- Q[og,og]

    # start time
    start <- Sys.time()

    # cholesky factor
    cholesky <- Matrix::chol(Q, pivot = permute)

    # end time
    end <- Sys.time() - start

    # count zeros
    nonzero_count <- utils::tail(cholesky@p,1) - cholesky@Dim[1]


    # return list of ratios and time taken
    return(list(diff = (nonzero_count / d_minus_one), time_taken = end, cky = cholesky))
  }

  if (method == "both"){
    # print(c(i,j))

    ## matrix method
    Q <- Matrix::sparseMatrix(q$indices_i, q$indices_j, x = q$entries)
    Q <- Q[og,og]

    # start time
    start <- Sys.time()

    # cholesky factor
    cholesky <- Matrix::chol(Q, pivot = permute)

    # end time
    end <- Sys.time() - start

    # count zeros
    nonzero_count <- utils::tail(cholesky@p,1) - cholesky@Dim[1]


    ## spam method
    Q_spam <- spam::spam(0, (i*j+i+1), (i*j+i+1))
    Q_spam[cbind(q$indices_i, q$indices_j)] <- q$entries
    Q_spam <- Q_spam[og,og]

    # start time
    start <- Sys.time()

    # cholesky factor
    cholesky_spam <- spam::chol(Q_spam, pivot = permute)

    # end time
    end2 <- Sys.time() - start

    # count zeros
    nonzero_count2 <- utils::tail(cholesky_spam@rowpointers, 1) - 1 - cholesky_spam@dimension[1]


    # return list of ratios and time taken
    return(list(diff_matrix = (nonzero_count / d_minus_one), diff_spam = (nonzero_count2 / d_minus_one),
                time_taken_matrix = end, time_taken_spam = end2))
  }

}
