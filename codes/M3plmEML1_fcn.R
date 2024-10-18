#' -----------------------------------------------------------------------------
#' title:  "Latent variable selection for M3pl model by EML1"
#' author: "Laixu3107"
#' date:   "2023.11.05"
#' -----------------------------------------------------------------------------

if(sys.nframe() == 0L){rm(list=ls()); gc()}

# ---- required packages ------------------------------------------------------
library(mvQuad)     # GH quadrature points and weights
Rcpp::sourceCpp("./M3plmEML1_algorithm.cpp")


# ---- some useful function ----------------------------------------------------
calcu_CR <- function(A_t, A_opt, fixed, col_swap=F){
  # Calculate CR for column swapping case.
  # 2022.07.02, Laixu3017.
  
  if(!is.null(fixed)){ # with constraints on A
    A_t   <- A_t[,-fixed]
    A_opt <- A_opt[,-fixed]
  }
  J <- ncol(A_t)
  K <- nrow(A_t)
  
  CR    <- 0
  permu <- 1:K
  opt_permu <- 0
  
  if(!col_swap){
    CR <- sum((A_opt!=0)==(A_t!=0))/K/J
  }
  else{ # without constraints on A
    
    all_permu <- gtools::permutations(K,K,1:K)
    
    for(i in 1:nrow(all_permu)){
      
      CR_permu <- sum((A_opt[all_permu[i,],]!=0)==(A_t!=0))/K/J
      
      if(CR_permu > CR){
        CR <- CR_permu
        opt_permu <- i
      }
    }
    
    permu <- all_permu[opt_permu,]
  }
  
  return(list(CR=CR, permu=permu))
}




# ---- EML1 for M3pl model -----------------------------------------------------
M3pl_EML1 <- function(
    y,                # N*J mat, observed item responses.
    A_init,           # K*J mat, initial value of A (discrimination parameter).
    b_init,           # J*1 vec, initial value of b (difficulty parameter).
    c_init,           # J*1 vec, initial value of c (guessing parameter).
    Sigma_init,       # K*K mat, initial value of Sigma.
    eta,              # tuning param
    fixed,            # fixed item for identification.
    n_nodes = 5,      # number of quadrature points per dimension
    
    MaxIter.EM = 100,
    MaxIter.CD = 50,
    Tol.EM     = 1e-4,
    Tol.CD     = 1e-6,
    c_max = 0.25,     # max value of c.
    reestimate = TRUE,# re-estimate the parameter by EM
    gamma
){
  
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent trait
  
  # ---- GH nodes and weights ----
  GH <- createNIGrid(dim=K, type="GHN", level=n_nodes)
  x  <- GH$nodes
  wx <- GH$weights
  
  # ---- initial beta ----
  beta_init <- matrix(0, K+1, J)
  beta_init[ 1,] <- b_init
  beta_init[-1,] <- A_init
  
  # ---- lambda mat ----
  lambda_mat <- matrix(-1, K+1, J)
  lambda_mat[1,] <- 1
  # lambda_mat[2:(K+1),fixed] <- diag(1,K)
  lambda_mat[2:(K+1),fixed] <- 1*(A_init[,fixed]!=0) # modified in 2024.04.07
  
  # ---- call m3plm_eml1 ----
  eml1_time <- proc.time()
  eml1_output <- m3plm_eml1(
    y  = y,
    x  = x,
    wx = wx,
    beta  = beta_init,
    c     = c_init,
    sigma = Sigma_init,
    lambda_mat = lambda_mat,
    eta = eta,
    
    maxiter_em = MaxIter.EM,
    maxiter_cd = MaxIter.CD,
    tol_em     = Tol.EM,
    tol_cd     = Tol.CD,
    
    c_max = c_max
  )
  
  eml1_time <- proc.time() - eml1_time
  eml1_time <- as.numeric(eml1_time[3])
  
  eml1_beta  <- eml1_output$beta
  eml1_c     <- eml1_output$c
  eml1_sigma <- eml1_output$sigma
  eml1_qfcn_seq <- eml1_output$qfcn_seq
  eml1_qtt_seq  <- eml1_output$qtt_seq
  eml1_iter     <- eml1_output$iter
  
  eml1_obs_ebic <- calcu_obs_ebic(y = y, x = x, wx = wx, beta = eml1_beta,
                                  c = eml1_c, sigma = eml1_sigma, gamma = gamma)
  
  eml1_cpu_time <- summary(cpu.time)
  
  output <- list()
  output$eta        <- eta
  output$eml1_beta  <- eml1_beta
  output$eml1_c     <- eml1_c
  output$eml1_sigma <- eml1_sigma
  output$eml1_qfcn_seq <- eml1_qfcn_seq
  output$eml1_qtt_seq  <- eml1_qtt_seq
  output$eml1_iter     <- eml1_iter
  output$eml1_obs_ebic <- eml1_obs_ebic
  output$eml1_cpu_time <- eml1_cpu_time
  
  
  # ---- re-estimate the params by EM algorithm ----
  if(reestimate){
    
    lambda_mat <- (eml1_beta!=0)*1
    
    em_time <- proc.time()
    em_output <- m3plm_eml1(
      y  = y,
      x  = x,
      wx = wx,
      beta  = eml1_beta,
      c     = eml1_c,
      sigma = eml1_sigma,
      lambda_mat = lambda_mat,
      eta = 0,
      
      maxiter_em = MaxIter.EM,
      maxiter_cd = MaxIter.CD,
      tol_em     = Tol.EM,
      tol_cd     = Tol.CD,
      
      c_max = c_max
    )
    
    em_time <- proc.time() - em_time
    em_time <- as.numeric(em_time[3])
    
    em_time <- proc.time() - eml1_time
    em_time <- as.numeric(eml1_time[3])
    
    em_beta  <- em_output$beta
    em_c     <- em_output$c
    em_sigma <- em_output$sigma
    em_qfcn_seq <- em_output$qfcn_seq
    em_qtt_seq  <- em_output$qtt_seq
    em_iter     <- em_output$iter
    
    em_obs_ebic <- calcu_obs_ebic(y = y, x = x, wx = wx, beta = em_beta,
                                  c = em_c, sigma = em_sigma, gamma = gamma)
    
    em_cpu_time <- summary(cpu.time)
    
    output$em_beta  <- em_beta
    output$em_c     <- em_c
    output$em_sigma <- em_sigma
    output$em_qfcn_seq <- em_qfcn_seq
    output$em_qtt_seq  <- em_qtt_seq
    output$em_iter     <- em_iter
    output$em_obs_ebic <- em_obs_ebic
    output$em_cpu_time <- em_cpu_time
    
  }
  
  return(output)
  
}


M3pl_EML1_path <- function(
    y,                # N*J mat, observed item responses.
    A_init,           # K*J mat, initial value of A (discrimination parameter).
    b_init,           # J*1 vec, initial value of b (difficulty parameter).
    c_init,           # J*1 vec, initial value of c (guessing parameter).
    Sigma_init,       # K*K mat, initial value of Sigma.
    eta_list,         # tuning param
    fixed,            # fixed item for identification.
    n_nodes = 5,      # number of quadrature points per dimension
    
    MaxIter.EM = 100,
    MaxIter.CD = 50,
    Tol.EM     = 1e-4,
    Tol.CD     = 1e-6,
    c_max = 0.25,     # max value of c.
    reestimate = TRUE,# re-estimate the parameter by EM
    gamma = 0
){
  
  path <- list()
  
  time_total <- proc.time()
  for(i in 1:length(eta_list)){
    
    eta <- eta_list[i]
    cat(sprintf("eta: %02d %.3f\r", i, eta))

    output <- M3pl_EML1(
      y = y,                  # N*J mat, observed item responses.
      A_init = A_init,        # K*J mat, initial value of A.
      b_init = b_init,        # K*1 vec, initial value of b.
      c_init = c_init,
      Sigma_init = Sigma_init,# K*K mat, initial value of Sigma.
      eta = eta,
      fixed = fixed,          # fixed item for identifiability.
      n_nodes = n_nodes,      # number of quadrature points per dimension
      
      MaxIter.EM = MaxIter.EM,
      MaxIter.CD = MaxIter.CD,
      Tol.EM = Tol.EM,
      Tol.CD = Tol.CD,
      c_max = c_max,
      reestimate = reestimate,# re-estimate the parameter by EM
      gamma = gamma
    )
    
    eta_i <- sprintf("eta_%d",i)
    path[[eta_i]] <- output
  }
  time_total <- proc.time() - time_total
  time_total <- as.numeric(time_total[3])
  
  ebic_vec <- rep(Inf, length(path))
  if(reestimate){
    for(i in 1:length(path)){
      cat(sprintf("eta: %.3f\r", eta))
      eta_i <- sprintf("eta_%d",i)
      ebic_vec[i] <- path[[eta_i]]$em_obs_ebic
      
      opt <- which.min(ebic_vec)
      
      opt_eta   <- path[[sprintf("eta_%d",opt)]]$eta
      opt_A     <- path[[sprintf("eta_%d",opt)]]$em_beta[-1,]
      opt_b     <- path[[sprintf("eta_%d",opt)]]$em_beta[1,]
      opt_c     <- path[[sprintf("eta_%d",opt)]]$em_c
      opt_Sigma <- path[[sprintf("eta_%d",opt)]]$em_sigma
    }
  }
  else{
    for(i in 1:length(path)){
      eta_i <- sprintf("eta_%d",i)
      ebic_vec[i] <- path[[eta_i]]$eml1_obs_ebic
      
      opt <- which.min(ebic_vec)
      
      opt_eta   <- path[[sprintf("eta_%d",opt)]]$eta
      opt_A     <- path[[sprintf("eta_%d",opt)]]$eml1_beta[-1,]
      opt_b     <- path[[sprintf("eta_%d",opt)]]$eml1_beta[1,]
      opt_c     <- path[[sprintf("eta_%d",opt)]]$eml1_c
      opt_Sigma <- path[[sprintf("eta_%d",opt)]]$eml1_sigma
    }
  }

  return(list(
    opt        = opt,
    opt_eta    = opt_eta,
    opt_A      = opt_A,
    opt_b      = opt_b,
    opt_c      = opt_c,
    opt_Sigma  = opt_Sigma,
    time_total = time_total,
    path       = path
  ))
}


# ---- A simple test -----------------------------------------------------------
if(sys.nframe() == 0L){
  
  library(magrittr)
  library(MASS)
  
  # ---- true model ----
  a <- seq(2.4,0.6,by=-0.6)
  A_t <- matrix(0,3,12)
  A_t[1,1:4] <- a
  A_t[2,5:8] <- a
  A_t[3,9:12] <- a
  fixed <- c(1,5,9)  # for identification
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 1000       # no. of subjects
  
  # ---- true parameter setting ----
  b_t     <- rnorm(J,0,1)
  c_t     <- runif(J,0.01,0.1)
  Sigma_t <- matrix(0.1,K,K); diag(Sigma_t) <- 1
  
  # ---- generate random sample ----
  set.seed(1)
  x <- mvrnorm(n=N, mu=rep(0,K), Sigma=Sigma_t) # latent traits
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    `*` (matrix(data=1-c_t,nrow=N,ncol=J,byrow=T)) %>% 
    `+` (matrix(data=c_t,nrow=N,ncol=J,byrow=T)) %>% 
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  # ---- EML1 initial parameter ----
  A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE)
  A_init[,fixed] <- diag(1,K)
  b_init <- rep(0,J)
  c_init <- rep(0.05,J)
  Sigma_init <- diag(1,K)
  
  eta_list = (10:1)/100*N
  n_nodes = 5
  MaxIter.EM = 100
  MaxIter.CD = 50
  Tol.EM = 1e-4
  Tol.CD = 1e-6
  c_max = 0.25
  reestimate = TRUE
  gamma = 0
  # ---- call the EML1 algorithm ----
  output <- M3pl_EML1_path(
    y = y,                  # N*J mat, observed item responses.
    A_init = A_init,        # K*J mat, initial value of A.
    b_init = b_init,        # K*1 vec, initial value of b.
    c_init = c_init,
    Sigma_init = Sigma_init,# K*K mat, initial value of Sigma.
    eta_list = eta_list,
    fixed = fixed,          # fixed item for identifiability.
    n_nodes = n_nodes,      # number of quadrature points per dimension
    
    MaxIter.EM = 100,
    MaxIter.CD = 50,
    Tol.EM = 1e-4,
    Tol.CD = 1e-6,
    c_max = 0.25,
    reestimate = TRUE,      # re-estimate the parameter by EM
    gamma = 0
  )
  
  CR <- calcu_CR(output$opt_A, A_t, fixed, col_swap=F)$CR

  cat("\n")
  cat("A_opt:\n");  print(output$opt_A);     cat("\n")
  cat("b_opt:\n");  print(output$opt_b);     cat("\n")
  cat("c_opt:\n");  print(output$opt_c);     cat("\n")
  cat("S_opt:\n");  print(output$opt_Sigma); cat("\n")

  cat("time.total:", output$time_total, "\n")
  cat("CR:", CR, "\n")
  
}
