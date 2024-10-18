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


// Updating beta (i.e., (aj, bj)) by Newton method for each sub-model.
// [[Rcpp::export]]
arma::vec update_beta(arma::vec beta,
                      const arma::mat &x,
                      const arma::vec &w0,
                      const arma::vec &wj,
                      const int maxit,
                      const double tol = 1e-4
){

  int G = x.n_rows;
  int K = x.n_cols;

  arma::vec Fj(G);
  arma::vec u(G);
  arma::vec d1(K);
  arma::mat d2(K,K);

  int it = 0;
  double eps = 1.0;

  while(it < maxit && eps > tol){

    it += 1;
    Fj = 1 - 1/(1 + exp(x*beta));
    d1 = x.t() * (wj - w0 % Fj);

    eps = sqrt(sum(d1 % d1));

    u  = Fj % (1 - Fj) % w0;
    d2 = - x.t() * (x.each_col() % u);

    beta = beta - inv(d2)*d1;
  }

  return(beta);
}


// Calculation of Q(aj,bj|psi(t)), i.e., the bic or extended bic.
double calcu_exp_ebic(const arma::vec &beta,
                      const arma::mat &x,
                      const arma::vec &w0,
                      const arma::vec &wj,
                      const int df,
                      const int N,
                      const int K,
                      const double gamma){

  double obj_fcn = -2*sum((x*beta)%wj - log(1+exp(x*beta))%w0) + df*log(N) + 2*gamma*df*log(K);
  return(obj_fcn);
}


// Optimal sub-model selection for each item j.
// [[Rcpp::export]]
void submod_selection(arma::Col<int> &mod,           // output, optimal sub-model for item j.
                      arma::vec &beta,               // output, estimated beta under optimal sub-model for item j.
                      double &qfcn,                  // output, minimal q function value for item j.
                      arma::mat &param_mat,          // output, estimated beta for all sub-models for item j.
                      arma::Mat<int> sub_mod_mat,    // sub-model space, each col a sub-model.
                      const arma::mat xs,
                      const arma::vec w0,
                      const arma::vec wj,
                      const int N,
                      const double gamma,            // tuning parameter for ebic, default = 0 for bic.
                      const int maxiter_newton
){

  int S = sub_mod_mat.n_cols;     // no. of sub-models for item j
  int K = sub_mod_mat.n_rows - 1; // no. of latent traits

  arma::vec exbic_vec(S);
  arma::vec param_s(K+1);


  int s, df;
  for(s=0; s<S; s++){

    arma::uvec ind = find(sub_mod_mat.col(s) != 0);
    df = ind.n_elem;

    param_s = param_mat.col(s);
    param_s(ind) = update_beta(param_s(ind), xs.cols(ind), w0, wj, maxiter_newton);
    param_mat.col(s) = param_s;
    exbic_vec(s) = calcu_exp_ebic(param_s, xs, w0, wj, df, N, K, gamma);

  }

  int opt = arma::index_min(exbic_vec);
  mod  = sub_mod_mat.col(opt);
  beta = param_mat.col(opt);
  qfcn = exbic_vec(opt);

}


// Calculation of the approximated Q-function ( i.e., Q(psi(t+1)|psi(t)) or Q(psi(t)|psi(t)) )
// given current model and its parameters.
double calcu_qfcn(const arma::Mat<int> &mod,
                  const arma::mat      &beta,
                  const arma::rowvec   &c,
                  const arma::mat      &sigma,
                  const arma::rowvec   &c_hat,
                  const arma::mat      &sigma_hat,
                  const arma::mat      &xs,
                  const arma::vec      &w0,
                  const arma::mat      &w,
                  const int N,
                  const int J,
                  const int K,
                  const double gamma
){
  int j;
  double qfcn = 0.0;
  qfcn += N * obj_func_cpp(sigma, sigma_hat);
  for(j=0;j<J;j++){
    qfcn += calcu_exp_ebic(beta.col(j), xs, w0, w.col(j), sum(mod.col(j)), N, K, gamma);
    qfcn += calcu_qc(c(j), c_hat(j), w0-w.col(j));
  }
  return(qfcn);
}


// Latent variable selection for M3pl model by EMS (GEMS) algorithm.
// [[Rcpp::export]]
List m3plm_ems(
    const arma::mat &y,
    const arma::mat &x,
    const arma::mat &wx,

    arma::mat  beta,
    arma::rowvec c,
    arma::mat  sigma,
    arma::Mat<int> mod,

    const arma::field<arma::Mat<int>> &sub_mods_list,
    arma::field<arma::mat>             beta_list,

    const double gamma = 0.0, // tuning parameter for EBIC, gamma = 0 to BIC.
    const int is_sigmaknown  = 0,
    const int maxiter_ems    = 50,
    const int maxiter_newton = 50,
    const double tol_newton = 1e-4,
    const double tol_param  = 1e-3,
    const double tol_qfcn   = 1e-4,

    const double c_max = 0.25
)
{

  Rcpp::Clock clock; // record CPU time.

  int N = y.n_rows;
  int J = y.n_cols;
  int G = x.n_rows;
  int K = x.n_cols;

  arma::mat      new_beta  = beta;
  arma::rowvec   new_c     = c;
  arma::mat      new_sigma = sigma;
  arma::Mat<int> new_mod   = mod;
  double         new_qfcn  = 0.0, qfcn_tt = 0.0; // Q(psi(t+1)|psi(t)) and Q(psi(t)|psi(t))

  // -- working matrix (or vector) for E-step --
  arma::mat xs = arma::ones(G, K+1);
  arma::vec w0(G);
  arma::mat w(G,J);
  arma::rowvec c_hat(J);

  // -- working matrix (or vector) for MS-step --
  int j;
  arma::mat sigma_hat(K,K);
  arma::Col<int> mod_j(K+1);
  arma::vec      beta_j(K+1);
  double         qfcn_j;
  arma::vec      qfcn_vec(J);

  // -----------------------
  // ---- EMS iteration ----
  // -----------------------

  int iter_ems = 0;
  double err_qfcn  = 1.0;
  arma::vec qfcn_seq(maxiter_ems); // record Q(psi(t+1)|psi(t)) for all iterations
  arma::vec qtt_seq(maxiter_ems);  // record Q(psi(t)|psi(t))   for all iterations

  while(iter_ems < maxiter_ems){

    iter_ems += 1; Rprintf("it:%d, e.\r", iter_ems);

    // ---- E-step: update xs, w0 & w ----
    clock.tick("e-step");
    sigma = 0.5 * sigma + 0.5 * sigma.t();
    xs.cols(1,K) = x * chol(sigma);
    m3plm_estep(y, xs, wx, beta, c, w0, w, c_hat);
    clock.tock("e-step");

    // ---- MS-step: ----
    // ---- 1. update sigma -----
    if(is_sigmaknown == 0){
      clock.tick("ms-update-sig");
      sigma_hat = calcu_sigma_hat(xs.cols(1,K), w0, N);
      new_sigma = calcu_sigma_cmle_cpp(sigma_hat, sigma, 1e-4);
      new_sigma = arma::symmatu(new_sigma);
      clock.tock("ms-update-sig");
    }
    
    // ---- 2. update c for each item ----
    for(j=0;j<J;j++){
      new_c(j) = std::min(c_hat(j), c_max);
    }

    // ---- 3. update beta for each item ----
    clock.tick("ms-update-beta");
    for(j=0;j<J;j++){
      Rprintf("it:%d, j:%d.\r", iter_ems, j+1);
      arma::mat beta_mat = beta_list(j);
      submod_selection(mod_j, beta_j, qfcn_j, beta_mat, sub_mods_list(j), xs, w0, w.col(j), N, gamma, maxiter_newton);
      new_mod.col(j)  = mod_j;
      new_beta.col(j) = beta_j;
      beta_list(j)    = beta_mat;
    }
    clock.tock("ms-update-beta");

    // ---- display ----
    // Rprintf("iter_ems:%d\n", iter_ems);
    // new_beta.print("new_beta:");


    // ---- stop criterion ----
    if((new_mod - mod).is_zero()){

      // compute Q(psi(t)|psi(t)) and Q(psi(t+1)|psi(t)).
      clock.tick("calcu_qfcn");
      qfcn_tt  = calcu_qfcn(mod, beta, c, sigma, c_hat, sigma_hat,
                            xs, w0, w, N, J, K, gamma);
      new_qfcn = calcu_qfcn(new_mod, new_beta, new_c, new_sigma, c_hat, sigma_hat,
                            xs, w0, w, N, J, K, gamma);
      qtt_seq(iter_ems-1)  = qfcn_tt;
      qfcn_seq(iter_ems-1) = new_qfcn;
      clock.tock("calcu_qfcn");

      err_qfcn = abs((new_qfcn - qfcn_tt)/qfcn_tt);

      if( err_qfcn < tol_qfcn ){
        break;
      }
    } // end stop criterion

    beta  = new_beta;
    c     = new_c;
    sigma = new_sigma;
    mod   = new_mod;

  } // end while, i.e., the ems iteration.

  Rcpp::List beta_Rcpp_list = Rcpp::wrap( Rcpp::RcppArmadillo::FieldImporter< arma::Mat<double> >( beta_list ) );

  clock.stop("cpu.time");

  // results
  List output = List::create(Rcpp::Named("beta_opt")  = new_beta,
                             Rcpp::Named("c_opt")     = new_c,
                             Rcpp::Named("sigma_opt") = new_sigma,
                             Rcpp::Named("qfcn_seq")  = qfcn_seq.subvec(0,iter_ems-1), // sequence of q_ttplus
                             Rcpp::Named("qtt_seq")   = qtt_seq.subvec(0,iter_ems-1),  // sequence of q_tt
                             Rcpp::Named("iter")      = iter_ems,

                             // test term
                             Rcpp::Named("xs") = xs,
                             Rcpp::Named("w0") = w0,
                             Rcpp::Named("w")  = w,
                             Rcpp::Named("beta_list") = beta_Rcpp_list
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



