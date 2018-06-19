## Tests
Images of generated plots are included in corresponding folders. The tests are to:

Check that the general precision matrix generation function works by comparing it with the centered parametrization generation function:
- consistency_precgen.R

Check that the precision generation function for centered parametrizations work by comparing samples drawn from exact inference with that from a gibbs sampler:
- precision_accuracy2.R
- precision_accuracy3.R


Check the performance of each package in obtaining the optimal fill-in for the cholesky factor when permuting the labels of the precision matrix:
- permutations2_Matrix.R
- permutations2_spam.R
- permutations2_asymmetric_Matrix.R
- permutations2_asymmetric_spam.R

Show the poor performance of the Matrix function when the labels are permuted in reverse and random order:
- fillinPerf_Matrix.R
