# confirm the correctness of the precision generating mechanism by comparing samples
# from exact inference and the gibbs sampler. we assume a mean of 0 for the root parameter
# and compare samples of the root variable B.

# sample specifications
burn <- 100
size <- 10000

# generate samples from the 3-level gibbs sampler
i <- 10
j <- 11
k <- 12
results_gibbs <- ghInf::centered_gibbs3(i = i, j = j, k = k, ndraws = size, burnin = burn, flat_prior = TRUE,
                                 tau = 1, tau_a = 2, tau_b = 1, tau_c = 2, sigma_2 = 2)
samples_gibbs <- results_gibbs$samples[,(burn+1):(size+burn)]

# trace plot for 3-level gibbs samples
df <- data.frame(iterations = seq(1, size), B = samples_gibbs[nrow(samples_gibbs),])
trace3 <- plot(df$iterations, df$B, type = 'l', main = 'Trace plot of beta',
               xlab = 'iterations', ylab = 'beta')

# generate samples with the precision matrix
results_exact <- ghInf::centered_precgen3(i = i, j = j, k = k, flat_prior = TRUE, tau = 1, tau_a = 2, tau_b = 1,
                                   tau_c = 2, sigma_2 = 2)
Q_exact <- Matrix::sparseMatrix(results_exact$indices_i, results_exact$indices_j, x = results_exact$entries)
cov_exact <- solve(Q_exact)
samples_exact <- MASS::mvrnorm(n = size, mu = rep(0, dim(cov_exact)[1]), Sigma = cov_exact)
samples_exact <- t(samples_exact)

# qqplot of samples from gibbs vs exact sampling
df2 <- data.frame(gauss_B = sort(samples_exact[nrow(samples_exact),]),
                  B = sort(samples_gibbs[nrow(samples_gibbs),]))
qq2 <- plot(df2$gauss_B, df2$B, main = "QQ plot of beta samples from Gibbs vs
            Gaussian Simulation", xlab = "Gaussian Simulation", ylab = "Gibbs")
abline(a = 0, b = 1, col = 'red')
# 45 degree line matches closely with qqplot: samples are consistent
