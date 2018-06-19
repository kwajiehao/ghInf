# compares performance in terms of fill-in ratio of nonzero elements for both spam and Matrix packages
# for 3 different permutations (for 2-level models)

# import libraries
library(Matrix)
library(spam)
library(spam64)
library(ggplot2)

##### spam #####

# number of children at each level
level1 <- c(10, 50, 100)
l <- length(level1)

# results
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "spam", permute_method = "none",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with depth-first labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "spam", permute_method = "depthfirst",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with reverse labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "spam", permute_method = "reverse",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with random labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- centered_fillin2(method = "spam", permute_method = "random",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in
