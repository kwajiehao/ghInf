# shows the user how to convert the output from the precision generation functions into a matrix
# for both the Matrix and spam package

# load libraries
library(Matrix)
library(spam)

i <- 5
J <- seq(1, i)
j <- sum(J)
test <- precgen2(i = i, J = J)

# creating a precision matrix using Matrix package
Q <- sparseMatrix(test$indices_i, test$indices_j, x = test$entries)
print(Q)

# creating a precision matrix using spam package
Q_spam <- spam(0, (i + j + 1), (i + j + 1)) # note that this is different for the centered functions
Q_spam[cbind(test$indices_i, test$indices_j)] <- test$entries
print(as.matrix(Q_spam))
