# compares performance in terms of fill-in ratio of nonzero elements for both spam and Matrix packages
# for 3 different permutations (for 2-level models)

# import libraries
library(Matrix)
library(spam)
library(spam64)
library(ggplot2)

##### matrix #####

# number of children at each level
level1 <- c(10, 50, 100)
l <- length(level1)

# generate a data frame to store results
test2_cols <- rep(level1, l)
test2_rows <- c()
for (i in 1:l){
  test2_rows <- c(test2_rows, rep(level1[i], l))
}
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "matrix", permute_method = "none",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with depth-first labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "matrix", permute_method = "depthfirst",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with reverse labelling (this takes some time)
results_reverse <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "matrix", permute_method = "reverse",
                             i = level1[i], j = level1[j])
    results_reverse[(i-1)*l + j] <- temp$diff
  }
}

print(results_reverse) # not optimal! with reverse labelling the ratio is very bad.

# with random labelling
results_random <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "matrix", permute_method = "random",
                             i = level1[i], j = level1[j])
    results_random[(i-1)*l + j] <- temp$diff
  }
}

print(results_random) # not optimal! the random labelling it is not as bad, but the ratio is still not great.
