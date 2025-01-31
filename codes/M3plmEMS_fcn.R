#' -----------------------------------------------------------------------------
#' title:  "Latent variable selection for M3pl model by EMS algorithm"
#' author: "Laixu3107"
#' date:   "2023.11.05"
#' -----------------------------------------------------------------------------

if(sys.nframe() == 0L){rm(list=ls()); gc()}

# ---- required packages -------------------------------------------------------
library(mvQuad)     # GH quadrature points and weights
Rcpp::sourceCpp("./M3plmEMS_algorithm.cpp")


# ---- some useful function ----------------------------------------------------
list_all_candidate_submod <- function(K){
  # List all candidate sub-models.
  # Laixu3107, Aug 9, 2022
  # K:      number of latent traits.
  # Output: (K+1) * (2^K) mat, each col denotes a sub-model, include intercept.
  
  models <- matrix(data=0, nrow=K, ncol=2^K)
  s <- 0
  for(k in 0:K){
    com_k <- combn(K, k)
    for(l in 1:ncol(com_k)){
      s <- s + 1
      models[com_k[,l],s] <- 1
    }
  }
  return(rbind(1,models))
}


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


# ---- EMS for M3pl model ------------------------------------------------------
M3pl_EMS <- function(
    y,                # N*J mat, observed item responses.
    A_init,           # K*J mat, initial value of A (discrimination parameter).
    b_init,           # J*1 vec, initial value of b (difficulty parameter).
    c_init,           # J*1 vec, initial value of c (guessing parameter).
    Sigma_init,       # K*K mat, initial value of Sigma.
    beta_list = NULL, # J  list, initial value of beta for all item j and sub-model s.
    fixed,            # fixed item for identification.
    n_nodes = 5,      # number of quadrature points per dimension
    gamma = 0,        # param for EBIC, 0 <= gamma <= 1. if 0, BIC.  
    
    is.sigmaknown  = 0,
    MaxIter.EMS    = 50,
    MaxIter.Newton = 50,
    Tol.para   = 1e-3,
    Tol.qfcn   = 1e-4,
    Tol.newton = 1e-4,
    c_max = 0.25      # max value of c.
){
  
  J <- ncol(A_init)  # number of items
  K <- nrow(A_init)  # number of latent trait
  
  # ---- GH nodes and weights ----
  GH <- createNIGrid(dim=K, type="GHN", level=n_nodes)
  x  <- GH$nodes
  wx <- GH$weights
  
  # ---- initial model & beta ----
  beta_init <- matrix(0, K+1, J)
  beta_init[ 1,] <- b_init
  beta_init[-1,] <- A_init
  
  mod_init <- matrix(0, K+1, J)
  mod_init[ 1,] <- 1
  mod_init[-1,] <- (A_init!=0)*1
  
  # ---- all sub-models ----
  sub_mods <- list_all_candidate_submod(K)
  
  sub_mods_list <- vector(mode="list", length=J)
  for(j in 1:J){
    if(j %in% fixed){
      sub_mods_list[[j]] <- matrix(mod_init[,j], nrow=K+1, ncol=1)
    }
    else{
      sub_mods_list[[j]] <- sub_mods
    }
  }
  
  # ---- beta list ----
  if(is.null(beta_list)){
    beta_list <- vector(mode="list", length=J)
    for(j in 1:J){
      if(j %in% fixed){
        beta_list[[j]] <- matrix(beta_init[,j], nrow=K+1, ncol=1)
      }
      else{
        beta_list[[j]] <- sub_mods*matrix(beta_init[,j], K+1, ncol(sub_mods))
      }
    }
  }
  
  # ---- call m3plm_ems ----
  time_ems <- proc.time()
  output <- m3plm_ems(
    y  = y,
    x  = x,
    wx = wx,
    beta  = beta_init,
    c     = c_init,
    sigma = Sigma_init,
    mod   = mod_init,
    sub_mods_list = sub_mods_list,
    beta_list     = beta_list,
    gamma = gamma,
    
    is_sigmaknown  = is.sigmaknown,
    maxiter_ems    = MaxIter.EMS,
    maxiter_newton = MaxIter.Newton,
    tol_newton     = Tol.newton,
    tol_param      = Tol.para,
    tol_qfcn       = Tol.qfcn,
    
    c_max = c_max
  )
  
  time_ems <- proc.time() - time_ems
  time_ems <- as.numeric(time_ems[3])
  
  obs_ebic <- calcu_obs_ebic(y = y, x = x, wx = wx, beta = output$beta_opt,
                             c = output$c_opt, sigma = output$sigma_opt, gamma = gamma)
  
  return(list(A_opt     = output$beta_opt[-1,],
              b_opt     = output$beta_opt[1 ,],
              c_opt     = output$c_opt,
              Sigma_opt = output$sigma_opt,
              qfcn_seq  = output$qfcn_seq,
              qtt_seq   = output$qtt_seq,
              iter_ems  = output$iter,
              time_ems  = time_ems,
              obs_ebic   = obs_ebic,
              all_cpu_time = summary(cpu.time)
              
  ))
  
}


# ---- A simple test -----------------------------------------------------------
if(sys.nframe() == 0L){
  
  library(magrittr)
  library(MASS)
  
  # ---- true model ----
  a <- seq(2.4,0.6,by=-0.2)
  A_t <- matrix(0,3,30)
  A_t[1, 1:10] <- a
  A_t[2,11:20] <- a
  A_t[3,21:30] <- a
  fixed <- c(1,11,21)  # for identification
  
  J <- ncol(A_t)  # no. of items
  K <- nrow(A_t)  # no. of latent traits
  N <- 1000       # no. of subjects
  
  # ---- true parameter setting ----
  b_t     <- rnorm(J,0,1)
  c_t     <- runif(J,0.01,0.1)
  Sigma_t <- matrix(0.1,K,K); diag(Sigma_t) <- 1
  
  # ---- generate random sample ----
  set.seed(627)
  x <- mvrnorm(n=N, mu=rep(0,K), Sigma=Sigma_t) # latent traits
  y <- x %>%
    `%*%` (A_t) %>%
    `+` (matrix(data=b_t,nrow=N,ncol=J,byrow=T)) %>%
    plogis(q=.) %>%
    `*` (matrix(data=1-c_t,nrow=N,ncol=J,byrow=T)) %>% 
    `+` (matrix(data=c_t,nrow=N,ncol=J,byrow=T)) %>% 
    rbinom(n=N*J, size=1, prob=.) %>%
    matrix(data=., nrow=N, ncol=J, byrow=F)
  
  # ---- EMS initial parameter ----
  A_init <- matrix(data=1/J, nrow=K, ncol=J, byrow=TRUE)
  A_init[,fixed] <- diag(1,K)
  b_init <- rep(0,J)
  c_init <- rep(0.05,J)
  Sigma_init <- diag(1,K)
  
  n_nodes = 5
  gamma   = 0
  MaxIter.EMS    = 50
  MaxIter.Newton = 50
  Tol.para   = 1e-3
  Tol.qfcn   = 1e-4
  Tol.newton = 1e-4
  c_max = 0.25
  
  # ---- call the EMS algorithm ----
  output <- M3pl_EMS(
    y = y,                  # N*J mat, observed item responses.
    A_init = A_init,        # K*J mat, initial value of A.
    b_init = b_init,        # K*1 vec, initial value of b.
    c_init = c_init,
    Sigma_init = Sigma_init,# K*K mat, initial value of Sigma.
    beta_list = NULL,
    fixed = fixed,          # fixed item for identifiability.
    n_nodes = n_nodes,      # number of quadrature points per dimension
    gamma = gamma,
    
    MaxIter.EMS    = MaxIter.EMS,
    MaxIter.Newton = MaxIter.Newton,
    Tol.para   = Tol.para,
    Tol.qfcn   = Tol.qfcn,
    Tol.newton = Tol.newton,
    c_max = 0.25
  )
  
  CR <- calcu_CR(output$A_opt, A_t, fixed, col_swap=F)$CR
  
  cat("\n")
  cat("A_opt:\n");  print(output$A_opt);     cat("\n")
  cat("b_opt:\n");  print(output$b_opt);     cat("\n")
  cat("c_opt:\n");  print(output$c_opt);     cat("\n")
  cat("S_opt:\n");  print(output$Sigma_opt); cat("\n")
  
  cat("qttplus_seq:\n");  print(output$qfcn_seq); cat("\n")
  cat("qtt_seq:\n");      print(output$qtt_seq);  cat("\n")
  
  cat("iter:",        output$iter_ems,    "\n")
  cat("time.total:",  output$time_ems,    "\n")
  cat("obs_ebic:",    output$obs_ebic,    "\n")
  cat("CR:", CR, "\n")
  
  print(output$all_cpu_time)
  
}
