# Distributed Multivariate Spline Estimation for Longitudinal Image-on-Scalar Regression

**Authors:** Hyunjae Cho, Yaotian Wang, and Shan Yu

This repository contains the code and supplementary data to implement Simulation Study I & II as described in the paper. The algorithm fits global and High-Dimensional (HD) methods using a parallel computing approach based on domain decomposition.

---

## 1. Required Packages

Our code requires the following R packages. Please ensure they are installed before running the scripts:

* require(Triangulation) # install_github("FIRST-Data-Lab/Triangulation")
* require(BPST)          # install_github("FIRST-Data-Lab/BPST")
* require(TPST) #install_github("funstatpackages/TPST")
* require(gghilbertstrings) #install_github("cran/gghilbertstrings")
* require(mgcv)
* require(MGLM)
* require(parallel)
* require(prodlim)
* require(MASS)
* require(splines2)
* require(Matrix)
* require(Rcpp)
* require(pracma)
* require(grpnet)
* require(dplyr)
* require(ggplot2)
* require(colorRamps)

---

## 2. Example Codes

Example datasets can be generated as follows to demonstrate the proposed methods.

### Simulation Study I 
Generated using funcs/simul1.data.gen.R.
* Y: True response $Y$ matrix.
* X: Covariate $X$ matrix.
* M.vec: Block diagonal matrix containing the elements of the weight matrix $M$.
* V: Spatial coordinates $V$ matrix.
* bivar.true: True bivariate coefficient functions $\alpha(v)$ and $\beta(v)$.
* ind.inside: Spatial indices of the grid points located strictly inside the horseshoe boundary.

### Simulation Study II 
Generated using funcs/simul2.data.gen.R.
* Y: True response $Y$ matrix.
* X: Covariate $X$ matrix.
* M.vec: Block diagonal matrix containing the elements of the weight matrix $M$.
* V: Spatial coordinates $V$ matrix.
* tij: Vector of observed time points.
* trivar.true: True trivariate coefficient function $\alpha(v, t)$.
* bivar.alpha: True bivariate coefficient function $\alpha(v)$.
* bivar.beta: True bivariate coefficient function $\beta(v)$.
* ind.inside: Spatial indices of the grid points located strictly inside the horseshoe boundary.

---

## 3. Reproducible Code for Simulation Studies

The following scripts reproduce the results and figures presented in the paper. We recommend setting the working directory to the root of this folder before running.

* Table 1 (Simulation Study I):
    * global_simu1.R: Reproduces the results for the Global method.
    * HD_simu1.R: Reproduces the results for the HD method.
* Table 2 (Simulation Study II):
    * global_simu2.R: Reproduces the results for the Global method.
    * HD_simu2.R: Reproduces the results for the HD method.
* MISE Calculation & Visualization Examples:
    * example_code/simu1.example.R: Code example for the Global method in Simulation Study I.
    * example_code/simu2.example.R: Code example for the HD method in Simulation Study II.
* Model Checking & Prediction (Section 5):
    * example_code/comparision.example.R: Demonstrates the prediction performance and model evaluation.

---

## 4. Main Functions and Package Structure

All helper functions are located under the funcs/ directory.

### Major Functions
* simul1.data.gen.R: Generates random datasets for the Simulation Study I scenario.
* simul2.data.gen.R: Generates random datasets for the Simulation Study II scenario.
* bivar.gc.fit.R: Fits only bivariate estimators (excluding the reference-level trajectories term).
* gc.fit.R: Fits estimators using the global method.
* gc.fit.dc.R: Fits estimators using a parallel algorithm based on domain decomposition.

### Functions for Model Fitting
* gc.fit.cv.R: Finds optimal tuning parameters using Cross-Validation (CV).
* gc.fit.ho.R: Obtains SSE/MSE to evaluate model performance using training and testing datasets.
* basis.tensor.R: Generates bivariate and trivariate spline basis matrices.
* energy.tensor.R: Generates bivariate and trivariate energy functions.

### Functions for Domain Decomposition
* local.fit.R: Fits local estimators for each sub-region.
* local.fit.bivar.R: Fits local bivariate estimators for each sub-region.
* ring.dc.R: Identifies the neighborhood of a triangle within a triangulation mesh.
* sampling.HC.R: Identifies the index of a triangle within a triangulation using a Hilbert space-filling curve.
* basis.tensor.local.R: Generates bivariate and trivariate spline basis matrices for local sub-regions.

### Functions for Model Checking and Comparisons
* trivar.fit.R: Fits baseline trivariate coefficient functions.
* trivar.fit.cv.R: Finds optimal tuning parameters for the baseline trivariate model using Cross-Validation (CV).
* trivar.fit.ho.R: Calculates SSE/MSE for the baseline trivariate model using training and testing datasets.
* plot.fun.R: Visualizes the true and estimated bivariate coefficient functions over the spatial domain.
