# illustrates the labelling of different permutations for a 2-level and 3-level model
I <- 2
J <- 4

# without permutation, I use a lexicographic ordering
n <- I*J + I + 1
test <- centered_precgen2(i = I, j = J)
test_mat <- sparseMatrix(test$indices_i, test$indices_j, x = test$entries)
test_mat <- as.matrix(test_mat)

cols <- rep('B', n)

for (i in 1:I){
  for (j in 1:J){
    cols[(i-1)*J + j] <- paste0('B', i, j)
  }
  cols[(I*J) + i] <- paste0('B', i)
}

rows <- cols

colnames(test_mat) <- cols
rownames(test_mat) <- rows
print(test_mat)

# depth-first ordering
og <- depthfirst(i = I, j = J, levels = 2)
print(test_mat[og, og])

# random ordering
o <- seq(1, n)
og <- sample(o)
print(test_mat[og, og])

# reverse ordering
o <- seq(1, n)
og <- rev(o)
print(test_mat[og, og])
