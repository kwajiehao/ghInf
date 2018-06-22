# plot the performance of the cholesky function in the Matrix package without applying permutation algorithms
# prior to matrix factorization. performance is measured in terms of fill-in relative to optimal fill-in
# for random and reverse permutations to the labelling of the precision matrix

##### matrix #####

# number of children at each level
level1 <- c(10, 50, 100)
l <- length(level1)

# generate a data frame to store results

# store values of i and j, and the number of parameters
test2_cols <- rep(level1, l)
test2_rows <- c()
for (i in 1:l){
  test2_rows <- c(test2_rows, rep(level1[i], l))
}

num_params <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    num_params[(i-1)*l + j] <- level1[i]*level1[j] + level1[i] + 1
  }
}

# with reverse labelling (this takes some time)
results_reverse <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    temp <- ghInf::centered_fillin2(method = "matrix", permute_method = "reverse",
                             i = level1[i], j = level1[j], permute = FALSE)
    results_reverse[(i-1)*l + j] <- temp$diff
  }
}

print(results_reverse) # not optimal! with reverse labelling the ratio is very bad.

# with random labelling
reps <- 30 # reduce this number since this might take a while
results <- rep(0, reps) # take the average to account for random effects

results_random <- rep(0, l^2)
for (i in 1:l){
  for (j in 1:l){
    print(paste0("level 1: ", level1[i], ", level 2: ", level1[j]))
    for (k in 1:reps){
      print(paste0("iteration", k))
      temp <- ghInf::centered_fillin2(method = "matrix", permute_method = "random",
                               i = level1[i], j = level1[j], permute = FALSE)
      results[k] <- temp$diff
    }

    results_random[(i-1)*l + j] <- mean(results)
  }
}

print(results_random) # not optimal! the random labelling it is not as bad, but the ratio is still not great.

########## Performance ##########

# calculate percentage of matrix filled in
perc_fill_reverse <- rep(0, l^2)
perc_fill_random <- perc_fill_reverse

for (i in 1:length(num_params)){
  perc_fill_reverse[i] <- results_reverse[i] * (num_params[i]-1)/num_params[i]^2
  perc_fill_random[i] <- results_random[i] * (num_params[i]-1)/num_params[i]^2
}

# plot performance against number of parameters
df_reverse <- data.frame(i = test2_rows, j = test2_cols, results = results_reverse, no_params = num_params, percentage_fill = perc_fill_reverse)
df_random <- data.frame(i = test2_rows, j = test2_cols, results = results_random, no_params = num_params, percentage_fill = perc_fill_random)

plot(df_reverse$no_params, df_reverse$results, type = 'p',
     main = 'ratio of fill-in to optimal fill-in (reverse permutation)', xlab = 'no. of params', ylab = 'ratio')
abline(lm(df_reverse$results ~ df_reverse$no_params), col = 'black')
plot(df_random$no_params, df_random$results, type = 'p',
     main = 'ratio of fill-in to optimal fill-in (random permutation)', xlab = 'no. of params', ylab = 'ratio')
abline(lm(df_random$results ~ df_random$no_params), col = 'black')

plot(df_reverse$no_params, df_reverse$results, type = 'p',
     main = 'ratio of fill-in to optimal fill-in', xlab = 'no. of params', ylab = 'ratio')
points(df_random$no_params, df_random$results, col = 'red', pch = "*")
abline(lm(df_reverse$results ~ df_reverse$no_params), col = 'black')
abline(lm(df_random$results ~ df_random$no_params), col = 'red')
legend(0, 5000, legend=c("reverse","random"), col=c("black","red"), pch=c("o","*"))


# plot percentage of cholesky factor filled in

plot(df_reverse$no_params, df_reverse$percentage_fill, type = 'p', main = 'percentage of entries that are non-zero (reverse permutation)',
     xlab = 'no. of params', ylab = 'percentage filled')
abline(lm(df_reverse$percentage_fill ~ (df_reverse$no_params^2)), col="red")
# scatter.smooth(df_reverse$no_params, df_reverse$percentage_fill, span = 0.9,
#              main = "percentage fill of cholesky factor against number of parameters (reverse permutation)",
#               xlab = "no. of params", ylab = "percentage fill", ylim = c(0, 0.5))

plot(df_random$no_params, df_random$percentage_fill, type = 'p', main = 'percentage of entries that are non-zero (random permutation)',
     xlab = 'no. of params', ylab = 'percentage filled')
abline(lm(df_random$percentage_fill ~ (df_random$no_params^2)), col="red")
# scatter.smooth(df_random$no_params, df_random$percentage_fill, span = 0.9,
#               main = "percentage fill of cholesky factor against number of parameters (random permutation)",
#               xlab = "no. of params", ylab = "percentage fill", ylim = c(0, 0.5))

plot(df_reverse$no_params, df_reverse$percentage_fill, type = 'p', main = 'percentage of entries that are non-zero of cholesky factor',
     xlab = 'no. of params', ylab = 'percentage filled', ylim = c(0, 0.5))
abline(lm(df_reverse$percentage_fill ~ (df_reverse$no_params^2)), col="black")
points(df_random$no_params, df_random$percentage_fill, pch = '*', col = "red")
abline(lm(df_random$percentage_fill ~ (df_random$no_params^2)), col = "red")
legend(0, 0.3, legend=c("reverse","random"), col=c("black","red"), pch=c("o","*"))

# ggplot plots
# gg <- ggplot2::ggplot(data = df_reverse, aes(x= no_params, y = results, col = "reverse")) + geom_point() +
#   geom_point(data = df_random, aes(x= no_params, y = results, col = "random")) +
#   geom_smooth(method = "lm", size = 0.5, se = FALSE) +
#   geom_smooth(data = df_random, method = "lm", color = "red", size = 0.5, se = FALSE) +
#   labs(title = "Ratio of fill-in to optimal-fill-in (Matrix package)", x = "number of parameters",
#        y = "ratio")
# gg + ggplot2::labs(color="permutation")
