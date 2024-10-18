# EMS for latent variable selection in M3PL model
## Introduction

The expectation model selection (EMS) algorithm and EM-based L1 (EML1) penalized method for latent variable selection in multidimensional 3-parameter logistic (M3PL) model are provided. Both methods adopt the Gauss-Hermite quadrature to numerically approximate the multidimensional integral. The EMS applies the Newton's method for parameter updating, while EML1 uses the coordinate descent algorithm. The two methods can be directly used for M2PL by setting the initial vaule of $c_j$s to zero.

The **codes** directory contains the following 4 files:

- M3plmEMS_algorithm.cpp includes the c++ implementations of EMS for M3PLM.
- M3plmEMS_fcn.R comprises R functions for EMS along with a simple test that users can directly invoke.
- M3plmEML1_algorithm.cpp includes the c++ implementations of EM-based L1 penalized (EML1) method for M3PLM.
- M3plmEML1_fcn.R comprises R functions for EML1 along with a simple test that users can directly invoke.

Since these implementations rely on Gauss-Hermite quadrature, please ensure that the R package 'mvQuad' is installed before running the code. Besides, the c++ implementations are based on R-packages 'Rcpp', 'RcppArmadillo' and 'RcppClock'. To run the examples.R, the R-packages 'magrittr' and 'MASS' are required.

## Citation

To cite these codes in publications, please use the following reference:

Shang, L., Xu, P. F., Shan, N., Tang, M. L, & Zheng, Q. Z. (2024). The improved EMS algorithm for latent variable selection in M3PL model. Applied Psychological Measurement. Accepted for publication.
