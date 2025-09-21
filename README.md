
Distributed Multivariate Spline Estimation for Longitudinal Image-on-Scalar Regression
================
Hyunjae Cho, Yaotian Wang, and Shan Yu

2025-05-03

The codes are to implement simulation study I & II. The algorithm fits global and HD methods through a parallel algorithm based on domain decomposition.

## Required packages

Our code requires the following packages:
  
  ```{r}
# required packages
require(Triangulation) # install_github("FIRST-Data-Lab/Triangulation")
require(BPST) # install_github("FIRST-Data-Lab/BPST")
require(mgcv)
require(MGLM)
require(parallel)
require(prodlim)
require(MASS)
require(Matrix)
require(Rcpp)
require(pracma)
```

## Code for simulation studies in the paper

Code example for the result of global method in Table 1 in the simulation studies: *global_simu1.R*

Code example for the result of HD method in Table 1 in the simulation studies: *HD_simu1.R*

Code example for the result of global method in Table 2 in the simulation studies: *global_simu2.R*
  
Code example for the result of HD method in Table 2 in the simulation studies:: *HD_simu2.R*
  
## Main functions and package structure
  
### Major functions
  
- *simul1.data.gen.R*: generates random datasets for simulation study I scenario.
- *simul2.data.gen.R*: generates random datasets for simulation study II scenario.
- *bivar.gc.fit.R*: fits only bivariate estimators (without reference-level trajectories term).
- *gc.fit.R*: fits estimators through global methods.
- *gc.fit.dc.R*: fits estimators through a parallel algorithm based on domain decomposition.

### Functions for model fitting

- *gc.fit.gcv.R* : finds optimal tuning parameters using GCV.
- *basis.tensor.R*: generates bivariate and trivariate spline basis matrix.
- *energy.tensor.R*: generates bivariate and trivariate energy functions.

### Functions for domain decomposition

- *local.fit.R* : fits a local estimators for each subregions.
- *local.fit.bivar.R* : fits a local bivariate estimators for each subregions.
- *ring.dc.R*: identifies the neighborhood of a triangle within a triangulation.
- *sampling.HC.R* : identifies the index of a triangle within a triangulation using Hilbert space filling-curve.
- *basis.tensor.local.R*: generates bivariate and trivariate spline basis matrix in each subregions.

