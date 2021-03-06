# ghInf

The goal of ghInf is to provide tools for inference for general Gaussian hierarchical models. Note that since this package is for personal use, all observations are assumed to be 0 for now, and the mean of the root parameter also assumed to be zero (subject to updating in the future). The functions can be found in the R folder. This includes functions which:

Generate the precision matrix for centered 2- or 3-level Gaussian hierarchical models:
- centered_precgen2.R  
- centered_precgen3.R

Generate the precision matrix for any 2- or 3-level Gaussian hierarchical models:
- precgen2.R  
- precgen3.R

Calculate the fill-in ratio of the calculated Cholesky factor relative to the optimal fill-in for both spam and Matrix packages, and for different label permutations:
- centered_fillin2.R
- centered_fillin3.R
- fillin2.R
- fillin3.R

A Gibbs sampler for centered 2- or 3-level Gaussian hierarchical models:
- centered_gibbs2.R
- centered_gibbs3.R

A helper function for permuting the label of the precision matrix into a depth-first order from a lexographical order:
- depthfirst.R

## Installation

You can install the released version of ghInf by using the following code.

```r
library(devtools)
install_github("kwajiehao/ghInf")
```

You may also use the githubinstall package.

```r
library(githubinstall)
githubinstall("ghInf")
```


## Tests

Test examples which show how the functions work can be found in the "tests" folder. 
