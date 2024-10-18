# EMS for latent variable selection in M3PL model
## Introduction

The generalized expectation model selection (GEMS) algorithm is proposed for dealing with the model selection problem in presence of missing data. For the latent variable selection in multidimension two-parameter logistic model (M2PLM), we present an efficient implementation of GEMS to find the optimal model (i.e., the structure of item-trait relationships) and the parameter estimates (including the item discrimination and difficulty parameters) under optimal model which results in the smallest BIC value. The GEMS for M2PLM is more computationally efficient than the EMS proposed by Xu et al. (2022).

The **codes** directory contains the following 4 files:

- M3plmEMS_algorithm.cpp includes the c++ implementations of EMS for M3PLM.
- M3plmEMS_fcn.R comprises R functions for EMS along with a simple test that users can directly invoke.
- M3plmEML1_algorithm.cpp includes the c++ implementations of EM-based L1 penalized (EML1) method for M3PLM.
- M3plmEML1_fcn.R comprises R functions for EML1 along with a simple test that users can directly invoke.

These implementations rely on Gauss-Hermite quadrature. Please ensure that the R package 'mvQuad' is installed before running the code. The c++ implementations are based on R-packages 'Rcpp', 'RcppArmadillo' and 'RcppClock'. To run the examples.R, the R-packages 'magrittr' and 'MASS' are required.

## Citation

To cite these codes in publications, please use the following reference:

Shang, L., Xu, P. F., Shan, N., Tang, M. L, & Zheng, Q. Z. (2024). The improved EMS algorithm for latent variable selection in M3PL model. Applied Psychological Measurement. Accepted for publication.
