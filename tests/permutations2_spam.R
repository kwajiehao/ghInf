# compares performance in terms of fill-in ratio of nonzero elements for both spam and Matrix packages
# for 3 different permutations (for 2-level models)

##### spam #####

# number of children at each level
level1 <- c(10, 50, 100)
l <- length(level1)

##########  results with pivot = TRUE, allowing package to permute rows and columns before factorization ##########
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "none",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with depth-first labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "depthfirst",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with reverse labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "reverse",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with random labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "random",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

####################  results with pivot = FALSE ####################
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "none",
                             i = level1[i], j = level1[j])
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with depth-first labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "depthfirst",
                             i = level1[i], j = level1[j], permute = FALSE)
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # optimal fill-in

# with reverse labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "reverse",
                             i = level1[i], j = level1[j], permute = FALSE)
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # not optimal!

# with random labelling
results <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "spam", permute_method = "random",
                             i = level1[i], j = level1[j], permute = FALSE)
    results[(i-1)*l + j] <- temp$diff
  }
}

print(results) # not optimal

# reverse and random labelling requires row and column permutation to obtain the optimal fill-in
