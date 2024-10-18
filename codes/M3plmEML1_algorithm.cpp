#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#include <RcppArmadillo.h>
#include <RcppClock.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppClock)]]
using namespace Rcpp;


// The E-step of EMS for M3pl model.
// [[Rcpp::export]]
void m3plm_estep(const arma::mat &y,
                 const arma::mat &x,    // fixed points matrix contain a column of 1.
                 const arma::vec &wx,   // weights or prior of x.
                 const arma::mat &beta, // including a & b.
                 const arma::rowvec &c,
                 
                 arma::vec &w0_xg,       // f_g.
                 arma::mat &wj_xg_v1,    // r_gj. (also as denominator of d hat and 1-wj_xg_v1 is denominator of c hat)
                 arma::rowvec &c_hat     // estimated c.
)
{
  int N = y.n_rows; // number of subjects
  int J = y.n_cols; // number of items
  int G = x.n_rows; // number of grid points
  
  arma::mat log_Fj_xg(G,J);
  arma::mat log_Hj_xg(G,J); // Hj = 1 - Fj
  arma::mat log_Pj_xg(G,J);
  arma::mat log_Qj_xg(G,J); // Qj = 1 - Pj
  
  arma::mat Rj_xg(G,J);     // Rj = dj*Fj/Pj
  arma::vec tp_xg(G);
  
  Rj_xg = x * beta;
  Rj_xg = 1 - 1/(1+exp(Rj_xg));               // obtain Fj_xg
  log_Pj_xg = Rj_xg;
  log_Pj_xg.each_row() %= (1 - c);
  log_Pj_xg.each_row() += c;                  // obtain Pj_xg
  log_Qj_xg = 1 - log_Pj_xg;                  // obtain Qj_xg
  
  Rj_xg /= log_Pj_xg;                         // obtain Rj_xg
  
  log_Pj_xg = log(log_Pj_xg);                 // obtain log_Pj_xg
  log_Qj_xg = log(log_Qj_xg);                 // obtain log_Qj_xg
  
  w0_xg.zeros();
  wj_xg_v1.zeros();
  
  int i,j;
  for(i=0;i<N;i++){
    
    tp_xg.zeros();
    
    for(j=0;j<J;j++){
      if(y(i,j)==1){
        tp_xg += log_Pj_xg.col(j);
      }
      else{
        tp_xg += log_Qj_xg.col(j);
      }
    }
    
    tp_xg = exp(tp_xg) % wx;
    tp_xg /= sum(tp_xg);
    
    w0_xg += tp_xg;
    
    for(j=0;j<J;j++){
      if(y(i,j)==1){
        wj_xg_v1.col(j) += tp_xg;  //wj_xg_y1
      }
    }
  }
  
  wj_xg_v1 %= Rj_xg;
  c_hat = 1 - (N - arma::sum(y,0))/(N - arma::sum(wj_xg_v1, 0));
  
}


// Obtain S*.
arma::mat calcu_sigma_hat(const arma::mat &x, const arma::vec &w0, const int N){
  arma::mat sigma = x.t() * (x.each_col() % w0) / N;
  return(sigma);
}


double obj_func_cpp(arma::mat sigma, arma::mat sigma_hat){
  arma::mat sigma_inv = arma::inv(sigma);
  return arma::accu( sigma_inv % sigma_hat ) + log(arma::det(sigma));
}


arma::mat calcu_sigma_cmle_cpp(arma::mat scov, arma::mat scor, double tol){
  
  arma::mat sigma_hat = scov;
  arma::mat sigma0 = scor;
  arma::mat sigma1 = sigma0;
  arma::mat tmp = sigma0;
  double eps = 1;
  double step = 1;
  while(eps > tol){
    step = 1;
    tmp = arma::inv(sigma0);
    sigma1 = sigma0 - step * ( - tmp * sigma_hat * tmp + tmp );
    sigma1.diag().ones();
    sigma1 = arma::symmatu(sigma1);   // add 2021.04.25
    while(obj_func_cpp(sigma0, sigma_hat) < obj_func_cpp(sigma1, sigma_hat) ||
          min(arma::eig_sym(sigma1)) < 0){
      step *= 0.5;
      sigma1 = sigma0 - step * ( - tmp * sigma_hat * tmp + tmp );
      sigma1.diag().ones();
      sigma1 = arma::symmatu(sigma1);   // add 2021.04.25
    }
    eps = obj_func_cpp(sigma0, sigma_hat) - obj_func_cpp(sigma1, sigma_hat);
    // Rprintf("eps= %f\n", eps);
    sigma0 = sigma1;
  }
  return sigma0;
}


// Calculation of Qj(cj|psi(t)).
double calcu_qc(const double    c_est,
                const double    c_hat,
                const arma::vec wj_xg_v0){
  if(c_est == 0){
    return(0.0);
  }
  else{
    return( -2*(c_hat*log(c_est) + (1-c_hat)*log(1-c_est)) * sum(wj_xg_v0) );
  }
}


// Calculation of Qj(aj,bj|psi(t)).
double calcu_qab(const arma::vec &beta, // param for item j.
                 const arma::mat &x,    // contain a column of 1.
                 const arma::vec &w0,
                 const arma::vec &wj,
                 const double lambda,
                 const int K){
  
  double obj_fcn = - sum((x*beta)%wj - log(1+exp(x*beta))%w0) + lambda * arma::accu(beta.subvec(1,K));
  return(obj_fcn);
}


// Calculation of the approximated Q-function ( i.e., Q(psi(t+1)|psi(t)) or Q(psi(t)|psi(t)) )
// given current parameters.
double calcu_qfcn(const arma::mat      &beta,
                  const arma::rowvec   &c,
                  const arma::mat      &sigma,
                  const arma::rowvec   &c_hat,
                  const arma::mat      &sigma_hat,
                  const arma::mat      &xs,
                  const arma::vec      &w0,
                  const arma::mat      &w,
                  const double lambda,
                  const int N,
                  const int K,
                  const int J
){
  int j;
  double qfcn = 0.0;
  qfcn += N * obj_func_cpp(sigma, sigma_hat);
  for(j=0;j<J;j++){
    qfcn += calcu_qab(beta.col(j), xs, w0, w.col(j), lambda, K);
    qfcn += calcu_qc(c(j), c_hat(j), w0-w.col(j));
  }
  return(qfcn);
}



double soft_thres(double delta, double lambda){
  if(abs(delta) <= lambda){
    return(0);
  }
  else if(delta > 0){
    return(delta - lambda);
  }
  else{
    return(delta + lambda);
  }
}





void coordinate_ascent(
    arma::vec &beta_j,
    const arma::Mat<int> &lambda_j,
    const arma::mat &xs,
    const arma::vec &f_g,
    const arma::vec &r_gj,
    const double eta,
    int maxit  = 50,
    double tol = 1e-6
){
  
  int G = xs.n_rows;
  int K = xs.n_cols - 1;
  
  arma::vec Fj_xg(G);
  double deriv1;
  double deriv2;
  double eps = 1.0;
  double old_qab = calcu_qab(beta_j, xs, f_g, r_gj, eta, K); 
  double new_qab;

  int k;
  for(k=0;k<(K+1);k++){
    if(lambda_j(k)==0){
      beta_j(k) = 0;
    }
  }
  
  int it = 0;
  
  while(eps > tol && it < maxit){
    
    it +=1;
    
    Fj_xg  = 1 - 1 / (1 + exp(xs * beta_j));
    deriv1 = sum((r_gj - Fj_xg % f_g));
    deriv2 = - sum(Fj_xg % (1-Fj_xg) % f_g);
    beta_j(0) -= deriv1/deriv2;

    for(k=1;k<(K+1);k++){
      
      if(lambda_j(k) == 0){
        beta_j(k) = 0;
      }
      else{
        Fj_xg  = 1 - 1 / (1 + exp(xs * beta_j));
        deriv1 = sum((r_gj - Fj_xg % f_g) % xs.col(k));
        deriv2 = - sum(Fj_xg % (1-Fj_xg) % f_g % xs.col(k) % xs.col(k));
        
        if(lambda_j(k) == 1){
          beta_j(k) = - soft_thres(deriv1 - deriv2 * beta_j(k), 0)/deriv2;
        }
        else{ // lambda_j(k) == -1
          beta_j(k) = - soft_thres(deriv1 - deriv2 * beta_j(k), eta)/deriv2;
        }
      }
    }
    
    new_qab = calcu_qab(beta_j, xs, f_g, r_gj, eta, K); 
    eps = pow(new_qab - old_qab, 2) / pow(new_qab + new_qab, 2);
    old_qab = new_qab;
  }
  
}


arma::mat update_beta(arma::mat beta,
                 const arma::Mat<int> lambda_mat,
                 const double eta,
                 const arma::mat &xs,
                 const arma::vec &f_g,
                 const arma::mat &r_gj,
                 const int maxit_cd  = 50,
                 const double tol_cd = 1e-6){
  
  
  int J = beta.n_cols;
  int K = beta.n_rows - 1;
  
  arma::vec      beta_j(K+1);
  arma::Col<int> lambda_j(K+1);
  arma::mat      new_beta = beta;

  int j;
  for(j=0;j<J;j++){
    
    beta_j = beta.col(j);
    lambda_j = lambda_mat.col(j);

    coordinate_ascent(beta_j, lambda_j, xs, f_g, r_gj.col(j), eta, maxit_cd, tol_cd);
    new_beta.col(j) = beta_j;
  }

  return(new_beta);
  
}


// [[Rcpp::export]]
Rcpp::List m3plm_eml1(
    const arma::mat &y,
    const arma::mat &x,
    const arma::vec &wx,
    
    arma::mat beta,
    arma::rowvec c,
    arma::mat sigma,

    const arma::Mat<int> lambda_mat, // 0 -> not include, 1 -> include but not penalty, -1 -> include and penalty
    const double eta,
    
    const int maxiter_em = 100,
    const int maxiter_cd = 50,
    const double tol_em = 1e-4,
    const double tol_cd = 1e-6,
    
    const double c_max = 0.25
){
  
  Rcpp::Clock clock; // record CPU time.
  
  int N = y.n_rows; // number of subjects
  int J = y.n_cols; // number of items
  int G = x.n_rows; // number of grid points
  int K = x.n_cols; // number of traits
  
  arma::mat    new_beta  = beta;
  arma::rowvec new_c     = c;
  arma::mat    new_sigma = sigma;
  double       new_qfcn  = 0.0, qfcn_tt = 0.0; // Q(psi(t+1)|psi(t)) and Q(psi(t)|psi(t))
  
  // -- working matrix (or vector) for E-step --
  arma::mat xs = arma::ones(G, K+1);
  arma::vec w0(G);
  arma::mat w(G,J);
  arma::rowvec c_hat(J);
  
  // -- working matrix (or vector) for MS-step --
  int j;
  arma::mat sigma_hat(K,K);
  
  
  // -----------------------
  // ---- EM iteration ----
  // -----------------------
  
  int iter_em = 0;
  double err_qfcn = 1.0;
  arma::vec qfcn_seq(maxiter_em); // record Q(psi(t+1)|psi(t)) for all iterations
  arma::vec qtt_seq(maxiter_em);  // record Q(psi(t)|psi(t))   for all iterations
  
  while(iter_em < maxiter_em){
    
    iter_em += 1; // Rprintf("it:%03d.\r", iter_em);
    
    // ---- E-step: update xs, w0 & w ----
    clock.tick("e-step");
    sigma = 0.5 * sigma + 0.5 * sigma.t();
    xs.cols(1,K) = x * chol(sigma);
    m3plm_estep(y, xs, wx, beta, c, w0, w, c_hat);
    clock.tock("e-step");

    // ---- M-step: ----
    // ---- 1. update sigma ----
    clock.tick("update-sig");
    sigma_hat = calcu_sigma_hat(xs.cols(1,K), w0, N);
    new_sigma = calcu_sigma_cmle_cpp(sigma_hat, sigma, 1e-4);
    new_sigma = arma::symmatu(new_sigma);
    clock.tock("update-sig");

    // ---- 2. Update c for each item ----
    clock.tick("update-c");
    for(j=0;j<J;j++){
      new_c(j) = std::min(c_hat(j), c_max);
    }
    clock.tock("update-c");
    
    // ---- 3. Update A and b ----
    clock.tick("update-beta");
    new_beta = update_beta(beta, lambda_mat, eta, xs, w0, w, maxiter_cd, tol_cd);
    clock.tock("update-beta");

    // ---- stop criterion ----
    clock.tick("calcu_qfcn");
    qfcn_tt  = calcu_qfcn(beta, c, sigma, c_hat, sigma_hat,
                          xs, w0, w, eta, N, K, J);
    new_qfcn = calcu_qfcn(new_beta, new_c, new_sigma, c_hat, sigma_hat,
                          xs, w0, w, eta, N, K, J);
    qtt_seq(iter_em-1)  = qfcn_tt;
    qfcn_seq(iter_em-1) = new_qfcn;
    clock.tock("calcu_qfcn");
    
    err_qfcn = abs((new_qfcn - qfcn_tt)/qfcn_tt);
    
    if( err_qfcn < tol_em ){
      break;
    }
    
    // ---- replace params ----
    beta  = new_beta;
    c     = new_c;
    sigma = new_sigma;

  } // end while, i.e., the em iteration.
  
  clock.stop("cpu.time");
  
  // -- return output --
  List output = List::create(Rcpp::Named("beta")     = new_beta,
                             Rcpp::Named("c")        = new_c,
                             Rcpp::Named("sigma")    = new_sigma,
                             Rcpp::Named("qfcn_seq") = qfcn_seq.subvec(0,iter_em-1), // sequence of q_ttplus
                             Rcpp::Named("qtt_seq")  = qtt_seq.subvec(0,iter_em-1),  // sequence of q_tt
                             Rcpp::Named("iter")     = iter_em,
                             
                             // test term
                             Rcpp::Named("xs") = xs,
                             Rcpp::Named("w0") = w0,
                             Rcpp::Named("w")  = w
                             );
  return output;
  
}


// Calculation of observed ebic(bic) for M3PL model.
// [[Rcpp::export]]
double calcu_obs_ebic(
    const arma::mat &y,
    const arma::mat &x,
    const arma::mat &wx,
    arma::mat &beta,
    arma::rowvec &c,
    arma::mat &sigma,
    double gamma = 0.0
){
  
  int N = y.n_rows;
  int J = y.n_cols;
  int G = x.n_rows;
  int K = x.n_cols;
  int df = arma::accu(beta!=0) - J;
  
  arma::mat xs = arma::ones(G, K+1);
  xs.cols(1,K) = x * chol(sigma);
  
  arma::mat log_Pj_xg(G,J);
  arma::mat log_Qj_xg(G,J);
  arma::vec tp(G);
  
  log_Pj_xg = xs * beta;
  log_Pj_xg = 1 - 1/(1+exp(log_Pj_xg));
  log_Pj_xg = log_Pj_xg.each_row() % (1-c);
  log_Pj_xg = log_Pj_xg.each_row() + c;
  log_Qj_xg = 1 - log_Pj_xg;
  
  log_Pj_xg = log(log_Pj_xg);
  log_Qj_xg = log(log_Qj_xg);
  
  arma::vec py(G);
  double obs_ebic = 0.0;
  
  int i, j;
  for(i=0;i<N;i++){
    
    py.zeros();
    
    for(j=0;j<J;j++){
      if(y(i,j)==1){
        py += log_Pj_xg.col(j);
      }
      else{
        py += log_Qj_xg.col(j);
      }
    }
    
    py = exp(py) % wx;
    obs_ebic += log(sum(py));
  }
  
  obs_ebic = -2 * obs_ebic + df*log(N) + 2*gamma*df*log(K);
  
  return(obs_ebic);
}
