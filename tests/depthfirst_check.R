#################### centered 2-level ####################
I <- 2
J <- 3
n <- (I*J + I + 1)
centered_test2 <- ghInf::centered_precgen2(i = I, j = J)

# without permutation, I use a lexicographic ordering
cols <- rep('B', n)

for (i in 1:I){
  for (j in 1:J){
    cols[(i-1)*J + j] <- paste0('B', i, j)
  }
  cols[(I*J) + i] <- paste0('B', i)
}

rows <- cols

# construct spam matrix
Q_spam <- spam::spam(0, n, n) # note that this is different for the centered functions
Q_spam[cbind(centered_test2$indices_i, centered_test2$indices_j)] <- centered_test2$entries
Q_spam_mat <- as.matrix(Q_spam)
colnames(Q_spam_mat) <- cols
rownames(Q_spam_mat) <- rows

# obtain the depth-first permuted matrix
og2 <- ghInf::depthfirst(i = i, j = j, levels = 2)
Q_spam_depthfirst <- Q_spam[og2, og2]
Q_depthfirst <- as.matrix(Q_spam_depthfirst)

# obtain the randomly permuted matrix
o <- seq(1, n)
og <- sample(o)
Q_spam_random <- Q_spam[og, og]
Q_spam_random_labelled <- Q_spam_mat[og, og]

# find the ordering used in cholesky decomposition
cholesky <- spam::chol(Q_spam_random, pivot = TRUE)
og3 <- spam::ordering(cholesky, inv = FALSE)

# compare with depth-first permutation
ordered_Q <-Q_spam_random_labelled[og3, og3]
sum(ordered_Q != Q_depthfirst)
ordered_Q




#################### centered 3-level ####################
I <- 2
J <- 3
K <- 3
n <- (I*J*K + I*J + I + 1)
centered_test3 <- ghInf::centered_precgen3(i = I, j = J, k = K)

# without permutation, I use a lexicographic ordering
cols <- rep('B', n)

for (i in 1:I){
  for (j in 1:J){
    for (k in 1:K){
      cols[ (i-1)*K*J + K*(j-1) + k] <- paste0('B', i, j, k)
    }
    cols[I*J*K + (i-1)*J + j] <- paste0('B', i, j)
  }
  cols[I*J*K + (I*J) + i] <- paste0('B', i)
}

rows <- cols

# construct spam matrix
Q_spam <- spam::spam(0, n, n) # note that this is different for the centered functions
Q_spam[cbind(centered_test3$indices_i, centered_test3$indices_j)] <- centered_test3$entries
Q_spam_mat <- as.matrix(Q_spam)
colnames(Q_spam_mat) <- cols
rownames(Q_spam_mat) <- rows
Q_spam_mat

# obtain the depth-first permuted matrix
og2 <- ghInf::depthfirst(i = I, j = J, k = K, levels = 3)
Q_spam_depthfirst <- Q_spam[og2, og2]
Q_depthfirst <- Q_spam_mat[og2, og2]
Q_depthfirst

# obtain the randomly permuted matrix
o <- seq(1, n)
og <- sample(o)
Q_spam_random <- Q_spam[og, og]
Q_spam_random_labelled <- Q_spam_mat[og, og]

# find the ordering used in cholesky decomposition
cholesky <- spam::chol(Q_spam_random, pivot = TRUE)
og3 <- spam::ordering(cholesky, inv = FALSE)

# compare with depth-first permutation
ordered_Q <- Q_spam_random_labelled[og3, og3]
sum(ordered_Q != Q_depthfirst)
ordered_Q
# indeed true depth-first
