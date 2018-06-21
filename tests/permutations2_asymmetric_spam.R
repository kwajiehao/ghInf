# compares performance in terms of fill-in ratio of nonzero elements for the spam package
# for general non-symmetric 2-level models. The goal is to see if the asymmetry of the tree
# structure, achieved by having a different number of children per node in a subsequent level,
# affects the fill-in.

# I try 3 different types of asymmetry structures:
#     - increasing order (1st node in level 2 has 1 child, 2nd has 2, and so on)
#     - pyramid order (if there are 5 nodes in level 2, the number of children per node is 1,2,3,2,1
#       respectively)
#     - extreme order (if there are 5 nodes in level 2, the number of children per node is 5,1,1,1,5
#       respectively)

# increasing order
i = 100
J = seq(1, i)

results_inc <- ghInf::fillin2(method="spam", i=i, J=J)
results_inc$ diff

# pyramid
i = 100
J = c(seq(1, floor(i/2)), seq(floor(i/2),1))
if (i%%2 == 1){J = c(seq(1, floor(i/2)), floor(i/2)+1, seq(floor(i/2),1))}

results_pyr <- ghInf::fillin2(method="spam", i=i, J=J)
results_pyr$diff

# extremes
i = 100
J = rep(1,i)
J[1] <- i
J[i] <- i

results_ext <- ghInf::fillin2(method="spam", i=i, J=J)
results_ext$diff

# all three tests show that optimal fill-in is returned

