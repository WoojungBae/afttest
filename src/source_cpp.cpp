#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @useDynLib afttest, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+"

double target_score2_mis(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec resid = log(Time) + Covari*beta;
  
  arma::uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  arma::mat tempmat_np = zeros(n,p); arma::vec tempvec_n = zeros(n); arma::vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    tempvec_n = sqrt(sum(tempmat_np%tempmat_np,1));
    tempvec_n.replace(0,1);
    
    NumericVector tempnumvec_n = wrap(sqrt(n)*(resid-resid(it))/tempvec_n);
    tempnumvec_n = pnorm(tempnumvec_n);
    
    F_vec += sum(tempmat_np.each_col()%(as<vec>(tempnumvec_n)),0).t()*Delta(it);
  }
  F_vec = F_vec/n - targetvector;
  
  double SumOfSqure = conv_to<double>::from(sum(pow(F_vec,2)));
  
  return SumOfSqure;
}

double target_score2_mns(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec resid = log(Time) + Covari*beta;
  
  arma::uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  arma::mat tempmat_np = zeros(n,p); arma::vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    F_vec += sum(tempmat_np.each_col()%conv_to<vec>::from((resid>=resid(it))),0).t()*Delta(it);
  }
  F_vec = F_vec/n - targetvector;
  
  double SumOfSqure = conv_to<double>::from(sum(pow(F_vec,2)));
  
  return SumOfSqure;
}

arma::vec target_score_mis(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec resid = log(Time) + Covari*beta;
  
  arma::uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  arma::mat tempmat_np = zeros(n,p); arma::vec tempvec_n = zeros(n); arma::vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    tempvec_n = sqrt(sum(tempmat_np%tempmat_np,1));
    tempvec_n.replace(0,1);
    
    NumericVector tempnumvec_n = wrap(sqrt(n)*(resid-resid(it))/tempvec_n);
    tempnumvec_n = pnorm(tempnumvec_n);
    
    F_vec += sum(tempmat_np.each_col()%(as<vec>(tempnumvec_n)),0).t()*Delta(it);
  }
  F_vec = F_vec/n - targetvector;
  
  return F_vec;
}

arma::vec target_score_mns(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec resid = log(Time) + Covari*beta;
  
  arma::uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  arma::mat tempmat_np = zeros(n,p); arma::vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    F_vec += sum(tempmat_np.each_col()%conv_to<vec>::from((resid>=resid(it))),0).t()*Delta(it);
  }
  F_vec = F_vec/n - targetvector;
  
  return F_vec;
}

List dfsane_mis(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){

  arma::vec b_old = beta;
  arma::vec F_old = target_score_mis(b_old,Time,Delta,Covari,targetvector);

  double a_k = (1/sqrt(sum(F_old%F_old))); if (a_k>1){a_k=1;}
  
  arma::vec b_new; arma::vec F_new; double tol_f;
  double optim_tol=1e-7; double tolerance=optim_tol+1; int ittt=0; int maxit=200;
  while(tolerance>optim_tol){
    
    b_new = b_old - a_k*F_old;
    
    F_new = target_score_mis(b_new,Time,Delta,Covari,targetvector);
    
    arma::vec s_k = b_new - b_old;
    arma::vec y_k = F_new - F_old;
    
    double tol_y = sum(y_k%y_k);
    tol_f = sum(F_new%F_new);
    
    if (tol_y>0) {a_k = (sum(s_k%y_k))/tol_y;} else {a_k = a_k;}
    
    tolerance = tol_f;
    
    b_old = b_new;
    F_old = F_new;
    
    if(b_new.has_nan()){tolerance = 0;}
    if(b_new.has_inf()){tolerance = 0;}
    if (ittt>maxit){tolerance = 0;}
    ittt += 1;
  }

  // -----------------------------------------------------------

  return List::create(tol_f,b_new);
}

// List dfsane_mis(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){
// 
//   arma::vec b_old = beta;
//   arma::vec F_old = target_score_mis(b_old,Time,Delta,Covari,targetvector);
//   double sig_k = (1/sqrt(sum(F_old%F_old))); if (sig_k>1){sig_k=1;}
// 
//   arma::vec b_new = b_old - sig_k*F_old;
// 
//   arma::vec F_new = target_score_mis(b_new,Time,Delta,Covari,targetvector);
// 
//   arma::vec s_k = b_new - b_old;
//   arma::vec y_k = F_new - F_old;
// 
//   double tol_0 = sum(F_old%F_old);
//   double tol_s = sum(s_k%s_k);
//   double tol_y = sum(y_k%y_k);
//   double tol_f = sum(F_new%F_new);
// 
//   double tolerance=tol_f+1; double optim_tol=1e-7; double tau_min=0.1; double tau_max=0.5;
//   double sig_min=1e-10; double sig_max=1e+10; double alp_p=1; double alp_m=1; double gam=1e-4;
//   double M=1; arma::vec f_bar=zeros(M); double it=1; double maxit=500;
// 
//   while(tolerance>optim_tol){
// 
//     // STEP 1
//     double eta_k = tol_0/pow(1+it,2);
// 
//     if (tol_y>0) {
//       sig_k = (sum(s_k%y_k))/tol_y;
//     } else {
//       sig_k = sig_k;
//     }
// 
//     if ((sig_min>abs(sig_k)) && (sig_max<abs(sig_k))){
//       if (tol_f<1e-10){
//         sig_k = 1e+5;
//       } else if (tol_f>1){
//         sig_k = 1;
//       } else {
//         sig_k = 1/sqrt(tol_f);
//       }
//     }
// 
//     arma::vec d_k = - sig_k * F_new;
// 
//     // STEP 2
//     int step_tol = 0; int itt = it - M * floor(it/M); f_bar(itt) = tol_f; double a_k =0;
//     while(step_tol == 0){
// 
//       // alpha_plus
//       arma::vec b_new_p = b_new + alp_p * d_k;
//       arma::vec F_new_p = target_score_mis(b_new_p,Time,Delta,Covari,targetvector);
// 
//       double RHS_p = sum(F_new_p%F_new_p);
//       double LHS_p = f_bar.max() + eta_k - gam * pow(alp_p,2) * tol_f;
// 
//       // alpha_minus
//       arma::vec b_new_m = b_new - alp_m * d_k;
//       arma::vec F_new_m = target_score_mis(b_new_m,Time,Delta,Covari,targetvector);
// 
//       double RHS_m = sum(F_new_m%F_new_m);
//       double LHS_m = f_bar.max() + eta_k - gam * pow(alp_m,2) * tol_f;
// 
//       if (RHS_p<=LHS_p){
//         d_k = d_k;
//         a_k = alp_p;
//         b_new = b_old + a_k * d_k;
//         step_tol = 1;
//       } else if (RHS_m<=LHS_m){
//         d_k = - d_k;
//         a_k = alp_m;
//         b_new = b_old + a_k * d_k;
//         step_tol = 1;
//       } else {
// 
//         double alp_p_t = (pow(alp_p,2) * tol_f)/(RHS_p + (2 * alp_p - 1) * tol_f);
// 
//         if (alp_p_t>(tau_max*alp_p)){
//           alp_p = tau_max * alp_p;
//         } else if (alp_p_t<(tau_min*alp_p)){
//           alp_p = tau_min * alp_p;
//         } else {
//           alp_p = alp_p_t;
//         }
// 
//         double alp_m_t = (pow(alp_m,2) * tol_f)/(RHS_m + (2 * alp_m - 1) * tol_f);
// 
//         if (alp_m_t>tau_max*alp_m){
//           alp_m = tau_max * alp_m;
//         } else if (alp_m_t<tau_min*alp_m){
//           alp_m = tau_min * alp_m;
//         } else {
//           alp_m = alp_m_t;
//         }
//       }
//     }
// 
//     // STEP 3
//     F_new = target_score_mis(b_new,Time,Delta,Covari,targetvector);
// 
//     s_k = b_new - b_old;
//     y_k = F_new - F_old;
// 
//     b_old=b_new;
//     F_old=F_new;
// 
//     tol_s = sum(s_k%s_k);
//     tol_y = sum(y_k%y_k);
//     tol_f = sum(F_new%F_new);
// 
//     tolerance = tol_f;
//     if (tol_f>tol_s){tolerance = tol_s;}
//     if (it>maxit){tolerance = 0;}
//     it += 1;
//   }
// 
//   return List::create(tol_f,b_new);
// }

List dfsane_mns(arma::vec beta, arma::vec Time, arma::vec Delta, arma::mat Covari, arma::vec targetvector){
  
  arma::vec b_old = beta;
  arma::vec F_old = target_score_mns(b_old,Time,Delta,Covari,targetvector);
  double sig_k = (1/sqrt(sum(F_old%F_old))); if (sig_k>1){sig_k=1;}
  
  arma::vec b_new = b_old - sig_k*F_old;
  
  arma::vec F_new = target_score_mns(b_new,Time,Delta,Covari,targetvector);
  
  arma::vec s_k = b_new - b_old;
  arma::vec y_k = F_new - F_old;
  
  double tol_0 = sum(F_old%F_old);
  double tol_s = sum(s_k%s_k);
  double tol_y = sum(y_k%y_k);
  double tol_f = sum(F_new%F_new);
  
  double tolerance=tol_f+1; double optim_tol=1e-7; double tau_min=0.1; double tau_max=0.5; 
  double sig_min=1e-10; double sig_max=1e+10; double alp_p=1; double alp_m=1; double gam=1e-4; 
  double M=1; arma::vec f_bar=zeros(M); double it=1; double maxit=500;
  
  while(tolerance>optim_tol){
    
    // STEP 1
    double eta_k = tol_0/pow(1+it,2);
    
    if (tol_y>0) {
      sig_k = (sum(s_k%y_k))/tol_y;
    } else {
      sig_k = sig_k;
    }
    
    if ((sig_min>abs(sig_k)) && (sig_max<abs(sig_k))){
      if (tol_f<1e-10){
        sig_k = 1e+5;
      } else if (tol_f>1){
        sig_k = 1;
      } else {
        sig_k = 1/sqrt(tol_f);
      }
    }
    
    arma::vec d_k = - sig_k * F_new;
    
    // STEP 2
    int step_tol = 0; int itt = it - M * floor(it/M); f_bar(itt) = tol_f; double a_k =0;
    while(step_tol == 0){
      
      // alpha_plus
      arma::vec b_new_p = b_new + alp_p * d_k;
      arma::vec F_new_p = target_score_mns(b_new_p,Time,Delta,Covari,targetvector);
      
      double RHS_p = sum(F_new_p%F_new_p);
      double LHS_p = f_bar.max() + eta_k - gam * pow(alp_p,2) * tol_f;
      
      // alpha_minus
      arma::vec b_new_m = b_new - alp_m * d_k;
      arma::vec F_new_m = target_score_mns(b_new_m,Time,Delta,Covari,targetvector);
      
      double RHS_m = sum(F_new_m%F_new_m);
      double LHS_m = f_bar.max() + eta_k - gam * pow(alp_m,2) * tol_f;
      
      if (RHS_p<=LHS_p){
        d_k = d_k;
        a_k = alp_p;
        b_new = b_old + a_k * d_k;
        step_tol = 1;
      } else if (RHS_m<=LHS_m){
        d_k = - d_k;
        a_k = alp_m;
        b_new = b_old + a_k * d_k;
        step_tol = 1;
      } else {
        
        double alp_p_t = (pow(alp_p,2) * tol_f)/(RHS_p + (2 * alp_p - 1) * tol_f);
        
        if (alp_p_t>(tau_max*alp_p)){
          alp_p = tau_max * alp_p;
        } else if (alp_p_t<(tau_min*alp_p)){
          alp_p = tau_min * alp_p;
        } else {
          alp_p = alp_p_t;
        }
        
        double alp_m_t = (pow(alp_m,2) * tol_f)/(RHS_m + (2 * alp_m - 1) * tol_f);
        
        if (alp_m_t>tau_max*alp_m){
          alp_m = tau_max * alp_m;
        } else if (alp_m_t<tau_min*alp_m){
          alp_m = tau_min * alp_m;
        } else {
          alp_m = alp_m_t;
        }
      }
    }
    
    // STEP 3
    F_new = target_score_mns(b_new,Time,Delta,Covari,targetvector);
    
    s_k = b_new - b_old;
    y_k = F_new - F_old;
    
    b_old=b_new;
    F_old=F_new;    
    
    tol_s = sum(s_k%s_k);
    tol_y = sum(y_k%y_k);
    tol_f = sum(F_new%F_new);
    
    tolerance = tol_f;
    if (tol_f>tol_s){tolerance = tol_s;}
    if (it>maxit){tolerance = 0;}
    it += 1;
  }
  
  return List::create(tol_f,b_new);
}

List omni_mis_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); arma::mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    Mhat_i_t(it) = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    obs_path += (as<arma::vec>(Mhat_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)))/sqrt(n);
  }
  
  List dMhat_i_t(n);
  for(int it=0; it<n; it++){
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,as<arma::vec>(Mhat_i_t(it))));
  }
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mis),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      arma::vec b_s = as<arma::vec>(b_s_opt_results[0]);
      
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::mat term1 = U_pi_phi_t_z/sqrt(n);
    
    arma::mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::mat>(fhat_t_z(it))+cumsum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    
    arma::mat term3 = cumsum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = reshape(stddev(as<arma::mat>(tempmat_n2path),0,1),n,n);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::mat>(app_path(it)))).max();
    app_std_path(it) = (as<arma::mat>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::mat>(app_std_path(it)))).max();
  }
  
  arma::mat obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Omni",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List omni_mns_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); arma::mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    Mhat_i_t(it) = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    obs_path += (as<arma::vec>(Mhat_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)))/sqrt(n);
  }
  
  List dMhat_i_t(n);
  for(int it=0; it<n; it++){
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,as<arma::vec>(Mhat_i_t(it))));
  }
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mns),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      arma::vec b_s = as<arma::vec>(b_s_opt_results[0]);
      
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::mat term1 = U_pi_phi_t_z/sqrt(n);
    
    arma::mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::mat>(fhat_t_z(it))+cumsum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    
    arma::mat term3 = cumsum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = reshape(stddev(as<arma::mat>(tempmat_n2path),0,1),n,n);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::mat>(app_path(it)))).max();
    app_std_path(it) = (as<arma::mat>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::mat>(app_std_path(it)))).max();
  }
  
  arma::mat obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Omni",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List link_mis_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mis),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      arma::vec b_s = as<arma::vec>(b_s_opt_results[0]);
      
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Link",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List link_mns_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mns),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      arma::vec b_s = as<arma::vec>(b_s_opt_results[0]);
      
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Link",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List form_mis_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType, int form){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::vec form_Covari = Covari.col(form-1);
  arma::vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mis),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      arma::vec b_s = as<arma::vec>(b_s_opt_results[0]);
      
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Form",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List form_mns_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType, int form){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::vec form_Covari = Covari.col(form-1);
  arma::vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mns),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      arma::vec b_s = as<arma::vec>(b_s_opt_results[0]);
      
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Form",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List omni_mis_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); arma::mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    Mhat_i_t(it) = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    obs_path += (as<arma::vec>(Mhat_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)))/sqrt(n);
  }
  
  List dMhat_i_t(n);
  for(int it=0; it<n; it++){
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,as<arma::vec>(Mhat_i_t(it))));
  }
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      List b_s_result = dfsane_mis(b, Time, Delta, Covari, U_phi_inf);
      
      b_s = as<arma::vec>(b_s_result[1]);
      
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::mat term1 = U_pi_phi_t_z/sqrt(n);
    
    arma::mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::mat>(fhat_t_z(it))+cumsum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    
    arma::mat term3 = cumsum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = reshape(stddev(as<arma::mat>(tempmat_n2path),0,1),n,n);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::mat>(app_path(it)))).max();
    app_std_path(it) = (as<arma::mat>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::mat>(app_std_path(it)))).max();
  }
  
  arma::mat obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Omni",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List omni_mns_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); arma::mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    Mhat_i_t(it) = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    obs_path += (as<arma::vec>(Mhat_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)))/sqrt(n);
  }
  
  List dMhat_i_t(n);
  for(int it=0; it<n; it++){
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,as<arma::vec>(Mhat_i_t(it))));
  }
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      List b_s_result = dfsane_mns(b, Time, Delta, Covari, U_phi_inf);
      
      b_s = as<arma::vec>(b_s_result[1]);
      
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::mat term1 = U_pi_phi_t_z/sqrt(n);
    
    arma::mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::mat>(fhat_t_z(it))+cumsum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    
    arma::mat term3 = cumsum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = reshape(stddev(as<arma::mat>(tempmat_n2path),0,1),n,n);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::mat>(app_path(it)))).max();
    app_std_path(it) = (as<arma::mat>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::mat>(app_std_path(it)))).max();
  }
  
  arma::mat obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Omni",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List link_mis_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      List b_s_result = dfsane_mis(b, Time, Delta, Covari, U_phi_inf);
      
      b_s = as<arma::vec>(b_s_result[1]);
      
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Link",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List link_mns_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      List b_s_result = dfsane_mns(b, Time, Delta, Covari, U_phi_inf);
      
      b_s = as<arma::vec>(b_s_result[1]);
      
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Link",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List form_mis_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, int form){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::vec form_Covari = Covari.col(form-1);
  arma::vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      List b_s_result = dfsane_mis(b, Time, Delta, Covari, U_phi_inf);
      
      b_s = as<arma::vec>(b_s_result[1]);
      
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Form",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}

List form_mns_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, int form){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  arma::vec zero_vec_1 = zeros(1);
  arma::vec zero_vec_p = zeros(p);
  arma::vec zero_vec_n = zeros(n);
  arma::mat zero_mat_np = zeros(n,p);
  arma::mat zero_mat_nn = zeros(n,n);
  
  arma::vec ones_vec_1 = ones(1);
  
  arma::vec tempvec_p(p);
  arma::vec tempvec_n(n);
  arma::mat tempmat_np(n,p);
  arma::mat tempmat_nn(n,n);
  
  arma::vec resid = log(Time) + Covari*b;
  
  arma::uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); arma::vec S_0_t = zero_vec_n; arma::mat S_1_t = zero_mat_np; arma::mat S_pi_t_z = zero_mat_nn;
  arma::vec form_Covari = Covari.col(form-1);
  arma::vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<arma::vec>(Y_i_t(it));
    S_1_t += as<arma::vec>(Y_i_t(it))*((Covari.row(it)));
    S_pi_t_z += (as<arma::vec>(Y_i_t(it)))*(as<arma::rowvec>(pi_i_z(it)));
  }
  
  arma::vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  
  arma::vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  arma::mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  List Mhat_i_t(n); List dMhat_i_t(n); arma::vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<arma::vec>(N_i_t(it))-(cumsum(as<arma::vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<arma::vec>(pi_i_z(it)));
  }
  obs_path = obs_path/sqrt(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  arma::vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  arma::vec given_data_g = exp(resid);
  
  double bw_g = 1.06 * stddev(given_data_g) * pow(n,-0.2);
  
  arma::vec ghat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += R::dnorm(pred_data(it),given_data_g(itt),bw_g,FALSE);
    }
  }
  ghat_0_t = ghat_0_t/n;
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<arma::rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  arma::vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  arma::vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  arma::vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  arma::vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  
  arma::vec given_data_f = exp(rhat_i);
  
  double bw_f = 1.06 * stddev(given_data_f) * pow(n,-0.2);
  
  arma::vec fhat_0_t = zeros(n);
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += R::dnorm(pred_data(it),given_data_f(itt),bw_f,FALSE);
    }
  }
  fhat_0_t = fhat_0_t/n;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    
    arma::vec Covari_col = Covari.col(itt);
    
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<arma::vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    arma::vec phi_i(n); arma::vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      
      phi_i = rnorm(n);
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<arma::vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<arma::vec>(dMhat_i_t(it))*((Covari.row(it))))*phi_i(it);
      }
      arma::vec U_phi_inf = (sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0)/n).t();
      
      List b_s_result = dfsane_mns(b, Time, Delta, Covari, U_phi_inf);
      
      b_s = as<arma::vec>(b_s_result[1]);
      
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<arma::rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<arma::vec>(dMhat_i_t(it))))*phi_i(it));
    }
    arma::mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    arma::vec resid_s = log(Time) + Covari*b_s;
    
    arma::uvec index_resid_s = sort_index(resid_s);
    
    arma::vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    List N_i_t_s(n); List Y_i_t_s(n); arma::vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      N_i_t_s(it) = (resid_s>=resid_s(it))*Delta_s(it);
      Y_i_t_s(it) = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<arma::vec>(Y_i_t_s(it));
    }
    
    arma::vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    
    arma::vec term1 = U_pi_phi_inf_z/sqrt(n);
    
    arma::vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrt(n);
    for(int it=0; it<p; it++){
      term2 += (as<arma::vec>(fhat_inf_z(it))+(sum((as<arma::mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    
    arma::vec term3 = (sum((S_pi_t_z.each_col())%diff(join_cols(zero_vec_1,Lambdahat_0_t-Lambdahat_0_t_s)))).t()/sqrt(n);
    
    app_path(itt) = term1 - term2 - term3;
    
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  arma::mat se_boot = stddev(as<arma::mat>(tempmat_npath),0,1);
  
  List app_std_path(path); arma::vec absmax_app_path(path); arma::vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    absmax_app_path(it) = (abs(as<arma::vec>(app_path(it)))).max();
    app_std_path(it) = (as<arma::vec>(app_path(it))/se_boot);
    absmax_app_std_path(it) = (abs(as<arma::vec>(app_std_path(it)))).max();
  }
  
  arma::vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  arma::uvec find_unstd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (find_unstd.size()); p_value = p_value/path;
  
  arma::uvec find_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (find_std.size()); p_std_value = p_std_value/path;
  
  return List::create(_["TestType"]="Form",_["path"]=path,_["beta"]=b,_["Time"]=Time,
                      _["Delta"]=Delta,_["Covari"]=Covari,_["Resid"]=resid,_["SE_boot"]=se_boot,
                        _["app_path"]=app_path,_["app_std_path"]=app_std_path,_["p_std_value"]=p_std_value,
                          _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,_["p_value"]=p_value);
}


// omni_mis_optim
List omni_mis_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType);
RcppExport SEXP _afttest_omni_mis_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mis_optim(path, b, Time, Delta, Covari, optimType));
  return rcpp_result_gen;
  END_RCPP
}
// omni_mns_optim
List omni_mns_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType);
RcppExport SEXP _afttest_omni_mns_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mns_optim(path, b, Time, Delta, Covari, optimType));
  return rcpp_result_gen;
  END_RCPP
}
// link_mis_optim
List link_mis_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType);
RcppExport SEXP _afttest_link_mis_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mis_optim(path, b, Time, Delta, Covari, optimType));
  return rcpp_result_gen;
  END_RCPP
}
// link_mns_optim
List link_mns_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType);
RcppExport SEXP _afttest_link_mns_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mns_optim(path, b, Time, Delta, Covari, optimType));
  return rcpp_result_gen;
  END_RCPP
}
// form_mis_optim
List form_mis_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType, int form);
RcppExport SEXP _afttest_form_mis_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP formSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mis_optim(path, b, Time, Delta, Covari, optimType, form));
  return rcpp_result_gen;
  END_RCPP
}
// form_mns_optim
List form_mns_optim(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, String optimType, int form);
RcppExport SEXP _afttest_form_mns_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP formSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mns_optim(path, b, Time, Delta, Covari, optimType, form));
  return rcpp_result_gen;
  END_RCPP
}
// omni_mis_DFSANE
List omni_mis_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari);
RcppExport SEXP _afttest_omni_mis_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mis_DFSANE(path, b, Time, Delta, Covari));
  return rcpp_result_gen;
  END_RCPP
}
// omni_mns_DFSANE
List omni_mns_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari);
RcppExport SEXP _afttest_omni_mns_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mns_DFSANE(path, b, Time, Delta, Covari));
  return rcpp_result_gen;
  END_RCPP
}
// link_mis_DFSANE
List link_mis_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari);
RcppExport SEXP _afttest_link_mis_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mis_DFSANE(path, b, Time, Delta, Covari));
  return rcpp_result_gen;
  END_RCPP
}
// link_mns_DFSANE
List link_mns_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari);
RcppExport SEXP _afttest_link_mns_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mns_DFSANE(path, b, Time, Delta, Covari));
  return rcpp_result_gen;
  END_RCPP
}
// form_mis_DFSANE
List form_mis_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, int form);
RcppExport SEXP _afttest_form_mis_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP formSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mis_DFSANE(path, b, Time, Delta, Covari, form));
  return rcpp_result_gen;
  END_RCPP
}
// form_mns_DFSANE
List form_mns_DFSANE(int path, arma::vec b, arma::vec Time, arma::vec Delta, arma::mat Covari, int form);
RcppExport SEXP _afttest_form_mns_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP formSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< arma::vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< arma::mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mns_DFSANE(path, b, Time, Delta, Covari, form));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_afttest_omni_mis_optim", (DL_FUNC) &_afttest_omni_mis_optim, 6},
  {"_afttest_omni_mns_optim", (DL_FUNC) &_afttest_omni_mns_optim, 6},
  {"_afttest_link_mis_optim", (DL_FUNC) &_afttest_link_mis_optim, 6},
  {"_afttest_link_mns_optim", (DL_FUNC) &_afttest_link_mns_optim, 6},
  {"_afttest_form_mis_optim", (DL_FUNC) &_afttest_form_mis_optim, 7},
  {"_afttest_form_mns_optim", (DL_FUNC) &_afttest_form_mns_optim, 7},
  {"_afttest_omni_mis_DFSANE", (DL_FUNC) &_afttest_omni_mis_DFSANE, 5},
  {"_afttest_omni_mns_DFSANE", (DL_FUNC) &_afttest_omni_mns_DFSANE, 5},
  {"_afttest_link_mis_DFSANE", (DL_FUNC) &_afttest_link_mis_DFSANE, 5},
  {"_afttest_link_mns_DFSANE", (DL_FUNC) &_afttest_link_mns_DFSANE, 5},
  {"_afttest_form_mis_DFSANE", (DL_FUNC) &_afttest_form_mis_DFSANE, 6},
  {"_afttest_form_mns_DFSANE", (DL_FUNC) &_afttest_form_mns_DFSANE, 6},
  {NULL, NULL, 0}
};

RcppExport void R_init_afttest(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}