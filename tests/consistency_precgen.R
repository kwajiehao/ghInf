## testing that the 2-level general precision matrix generator matches with the centered precision matrix generator

# centered
i <- 5
j <- 10
centered_test2 <- ghInf::centered_precgen2(i = i, j = j)

# general
J <- rep(10, 5)
test2 <- ghInf::precgen2(i = i, J = J)

# compare
Q_centered2 <- Matrix::sparseMatrix(centered_test2$indices_i, centered_test2$indices_j, x = centered_test2$entries)
Q2 <- Matrix::sparseMatrix(test2$indices_i, test2$indices_j, x = test2$entries)
sum(Q_centered2 != Q2) # equals to 0 so they are the same

## testing that the 3-level general precision matrix generator matches with the centered precision matrix generator

# centered
i <- 5
j <- 10
k <- 2
centered_test3 <- ghInf::centered_precgen3(i = i, j = j, k = k)

# general
J <- rep(10, 5)
K <- rep(2, sum(J))
test3 <- ghInf::precgen3(i = i, J = J, K = K)

# compare
Q_centered3 <- Matrix::sparseMatrix(centered_test3$indices_i, centered_test3$indices_j, x = centered_test3$entries)
Q3 <- Matrix::sparseMatrix(test3$indices_i, test3$indices_j, x = test3$entries)
sum(Q_centered3 != Q3) # equals to 0 so they are the same
