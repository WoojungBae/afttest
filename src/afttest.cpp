#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @useDynLib afttest, .registration = TRUE
//' @importFrom Rcpp evalCpp
//' @exportPattern "^[[:alpha:]]+"

double target_score2_mis(vec b, vec Time, vec Delta, mat Covari, vec targetvector){
 
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec resid = log(Time) + Covari*b;
  uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  mat tempmat_np = zeros(n,p); vec tempvec_n = zeros(n); vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    tempvec_n = sqrt(sum(tempmat_np%tempmat_np,1));
    tempvec_n.replace(0,1);
    
    tempvec_n = normcdf(sqrtn*(resid-resid(it))/tempvec_n);
    F_vec += sum(tempmat_np.each_col()%tempvec_n,0).t()*Delta(it);
  }
  F_vec -= targetvector;
  F_vec /= n;
  
  double SumOfSqure = norm(F_vec);
  
  return SumOfSqure;
}

double target_score2_mns(vec b, vec Time, vec Delta, mat Covari, vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  mat tempmat_np = zeros(n,p); vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    F_vec += sum(tempmat_np.each_col()%conv_to<vec>::from((resid>=resid(it))),0).t()*Delta(it);
  }
  F_vec -= targetvector;
  F_vec /= n;
  
  double SumOfSqure = norm(F_vec);
  
  return SumOfSqure;
}

vec target_score_mis(vec b, vec Time, vec Delta, mat Covari, vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec resid = log(Time) + Covari*b;
  uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  mat tempmat_np = zeros(n,p); vec tempvec_n = zeros(n); vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    tempvec_n = sqrt(sum(tempmat_np%tempmat_np,1));
    tempvec_n.replace(0,1);
    tempvec_n = normcdf(sqrtn*(resid-resid(it))/tempvec_n);
    F_vec += sum(tempmat_np.each_col()%tempvec_n,0).t()*Delta(it);
  }
  F_vec -= targetvector;
  F_vec /= n;
  
  return F_vec;
}

vec target_score_mns(vec b, vec Time, vec Delta, mat Covari, vec targetvector){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  mat tempmat_np = zeros(n,p); vec F_vec = zeros(p);
  for(int it=0; it<n; it++){
    tempmat_np = Covari.row(it) - Covari.each_row();
    F_vec += sum(tempmat_np.each_col()%conv_to<vec>::from((resid>=resid(it))),0).t()*Delta(it);
  }
  F_vec -= targetvector;
  F_vec /= n;
  
  return F_vec;
}

List dfsane_mis(int n, vec b, vec Time, vec Delta, mat Covari, vec targetvector){
  
  vec b_old = b;
  vec F_old = target_score_mis(b_old,Time,Delta,Covari,targetvector);
  double sig_k = (1/sqrt(sum(F_old%F_old))); if (sig_k>1){sig_k=1;}
  
  vec b_new = b_old - sig_k*F_old;
  
  vec F_new = target_score_mis(b_new,Time,Delta,Covari,targetvector);
  
  vec s_k = b_new - b_old;
  vec y_k = F_new - F_old;
  
  double tol_0 = sum(F_old%F_old);
  double tol_s = sum(s_k%s_k);
  double tol_y = sum(y_k%y_k);
  double tol_f = sum(F_new%F_new);
  
  // Stop sqrt(tol_f)/sqrt(n) <= e_a + e_r * sqrt(tol_0)/sqrt(n)
  // Stop tol_f <= (e_a * sqrt(n) + e_r * tol_0)^{2}
  double e_a = 1e-4; double e_r = 1e-3; 
  double optim_tol = pow(e_a * sqrt(n) + e_r * sqrt(tol_0),2);
  
  double tolerance=tol_f+1; double tau_min=0.1; double tau_max=0.5; 
  double sig_min=0.1; double sig_max=0.5; double alp_p=1; double alp_m=1; double gam=1e-4; 
  double M=1; vec f_bar=zeros(M); double it=1; double maxit=500;
  double eta_k, abssig_k, RHS_p, LHS_p, RHS_m, LHS_m, alp_p_t, alp_m_t;
  vec b_new_p, F_new_p, b_new_m, F_new_m;
  while(tolerance>optim_tol){
    
    // STEP 1
    eta_k = tol_0/pow(1+it,2);
    
    if (tol_y>0) {
      sig_k = (sum(s_k%y_k))/tol_y;
    } 
    // else {
    //   sig_k = sig_k;
    // }
    
    abssig_k = fabs(sig_k);
    if ((sig_min>abssig_k) || (sig_max<abssig_k)){
      if (tol_f<1e-10){
        sig_k = 1e+5;
      } else if (tol_f>1){
        sig_k = 1;
      } else {
        sig_k = 1/sqrt(tol_f);
      }
    }
    
    vec d_k = - sig_k * F_new;
    
    // STEP 2
    int step_tol = 0; int itt = it - M * floor(it/M); f_bar(itt) = tol_f; double a_k =0;
    while(step_tol == 0){
      
      // alpha_plus
      b_new_p = b_new + alp_p * d_k;
      F_new_p = target_score_mis(b_new_p,Time,Delta,Covari,targetvector);
      
      RHS_p = sum(F_new_p%F_new_p);
      LHS_p = f_bar.max() + eta_k - gam * pow(alp_p,2) * tol_f;
      
      // alpha_minus
      b_new_m = b_new - alp_m * d_k;
      F_new_m = target_score_mis(b_new_m,Time,Delta,Covari,targetvector);
      
      RHS_m = sum(F_new_m%F_new_m);
      LHS_m = f_bar.max() + eta_k - gam * pow(alp_m,2) * tol_f;
      
      if (RHS_p<=LHS_p){
        // d_k = d_k;
        a_k = alp_p;
        b_new = b_old + a_k * d_k;
        step_tol = 1;
      } else if (RHS_m<=LHS_m){
        d_k = - d_k;
        a_k = alp_m;
        b_new = b_old + a_k * d_k;
        step_tol = 1;
      } else {
        
        alp_p_t = (pow(alp_p,2) * tol_f)/(RHS_p + (2 * alp_p - 1) * tol_f);
        
        if (alp_p_t>(tau_max*alp_p)){
          alp_p = tau_max * alp_p;
        } else if (alp_p_t<(tau_min*alp_p)){
          alp_p = tau_min * alp_p;
        } else {
          alp_p = alp_p_t;
        }
        
        alp_m_t = (pow(alp_m,2) * tol_f)/(RHS_m + (2 * alp_m - 1) * tol_f);
        
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
    F_new = target_score_mis(b_new,Time,Delta,Covari,targetvector);
    
    s_k = b_new - b_old;
    y_k = F_new - F_old;
    
    b_old = b_new;
    F_old = F_new;    
    
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

List dfsane_mns(int n, vec b, vec Time, vec Delta, mat Covari, vec targetvector){
  
  vec b_old = b;
  vec F_old = target_score_mns(b_old,Time,Delta,Covari,targetvector);
  double sig_k = (1/sqrt(sum(F_old%F_old))); if (sig_k>1){sig_k=1;}
  
  vec b_new = b_old - sig_k*F_old;
  
  vec F_new = target_score_mns(b_new,Time,Delta,Covari,targetvector);
  
  vec s_k = b_new - b_old;
  vec y_k = F_new - F_old;
  
  double tol_0 = sum(F_old%F_old);
  double tol_s = sum(s_k%s_k);
  double tol_y = sum(y_k%y_k);
  double tol_f = sum(F_new%F_new);
  
  // Stop sqrt(tol_f)/sqrt(n) <= e_a + e_r * sqrt(tol_0)/sqrt(n)
  // Stop tol_f <= (e_a * sqrt(n) + e_r * tol_0)^{2}
  double e_a = 1e-4; double e_r = 1e-3; 
  double optim_tol = pow(e_a * sqrt(n) + e_r * sqrt(tol_0),2);
  
  double tolerance=tol_f+1; double tau_min=0.1; double tau_max=0.5; 
  double sig_min=0.1; double sig_max=0.5; double alp_p=1; double alp_m=1; double gam=1e-4; 
  double M=1; vec f_bar=zeros(M); double it=1; double maxit=500;
  double eta_k, abssig_k, RHS_p, LHS_p, RHS_m, LHS_m, alp_p_t, alp_m_t;
  vec b_new_p, F_new_p, b_new_m, F_new_m;
  while(tolerance>optim_tol){
    
    // STEP 1
    eta_k = tol_0/pow(1+it,2);
    
    if (tol_y>0) {
      sig_k = (sum(s_k%y_k))/tol_y;
    } 
    // else {
    //   sig_k = sig_k;
    // }
    
    abssig_k = fabs(sig_k);
    if ((sig_min>abssig_k) || (sig_max<abssig_k)){
      if (tol_f<1e-10){
        sig_k = 1e+5;
      } else if (tol_f>1){
        sig_k = 1;
      } else {
        sig_k = 1/sqrt(tol_f);
      }
    }
    
    vec d_k = - sig_k * F_new;
    
    // STEP 2
    int step_tol = 0; int itt = it - M * floor(it/M); f_bar(itt) = tol_f; double a_k =0;
    while(step_tol == 0){
      
      // alpha_plus
      b_new_p = b_new + alp_p * d_k;
      F_new_p = target_score_mns(b_new_p,Time,Delta,Covari,targetvector);
      
      RHS_p = sum(F_new_p%F_new_p);
      LHS_p = f_bar.max() + eta_k - gam * pow(alp_p,2) * tol_f;
      
      // alpha_minus
      b_new_m = b_new - alp_m * d_k;
      F_new_m = target_score_mns(b_new_m,Time,Delta,Covari,targetvector);
      
      RHS_m = sum(F_new_m%F_new_m);
      LHS_m = f_bar.max() + eta_k - gam * pow(alp_m,2) * tol_f;
      
      if (RHS_p<=LHS_p){
        // d_k = d_k;
        a_k = alp_p;
        b_new = b_old + a_k * d_k;
        step_tol = 1;
      } else if (RHS_m<=LHS_m){
        d_k = - d_k;
        a_k = alp_m;
        b_new = b_old + a_k * d_k;
        step_tol = 1;
      } else {
        
        alp_p_t = (pow(alp_p,2) * tol_f)/(RHS_p + (2 * alp_p - 1) * tol_f);
        
        if (alp_p_t>(tau_max*alp_p)){
          alp_p = tau_max * alp_p;
        } else if (alp_p_t<(tau_min*alp_p)){
          alp_p = tau_min * alp_p;
        } else {
          alp_p = alp_p_t;
        }
        
        alp_m_t = (pow(alp_m,2) * tol_f)/(RHS_m + (2 * alp_m - 1) * tol_f);
        
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
    
    b_old = b_new;
    F_old = F_new;    
    
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

List omni_mis_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  E_pi_t_z.replace(arma::datum::nan,0);
  
  // obs_path; t by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += tempvec_n*(as<rowvec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  Fhat_0_e.replace(arma::datum::nan,0);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      List b_s_result = dfsane_mis(n, b, Time, Delta, Covari, U_phi_inf);
      b_s = as<vec>(b_s_result[1]);
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    mat term1 = U_pi_phi_t_z/sqrtn;
    mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<mat>(fhat_t_z(it))+cumsum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    mat term3 = cumsum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))/sqrtn;
    
    tempmat_nn = term1 - term2 - term3;
    app_path(itt) = tempmat_nn;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec mat_se_boot = stddev(as<mat>(tempmat_n2path),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem,
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min = censoring; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  double censoring = 1-sum(Delta)/n;
  double kappa_min = sqrt(censoring);
  double kappa_max = 1;
  if(kappa_min<0.5){kappa_min = 0.5;}
  
  vec kappa = {kappa_min, kappa_max};
  kappa = quantile(mat_se_boot, kappa);
  mat_se_boot.clamp(kappa(0),kappa(1));
  mat se_boot = reshape(mat_se_boot,n,n);
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempmat_nn = as<mat>(app_path(it));
    absmax_app_path(it) = abs(tempmat_nn).max();
    
    tempmat_nn /= se_boot;
    tempmat_nn.replace(arma::datum::nan,0);
    app_std_path(it) = tempmat_nn;
    absmax_app_std_path(it) = abs(tempmat_nn).max();
  }
  
  mat obs_std_path = obs_path/se_boot;
  obs_std_path.replace(arma::datum::nan,0);
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List omni_mns_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  E_pi_t_z.replace(arma::datum::nan,0);
  
  // obs_path; t by x matrix
  List Mhat_i_t(n); List dMhat_i_t(n); mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += tempvec_n*(as<rowvec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  Fhat_0_e.replace(arma::datum::nan,0);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      List b_s_result = dfsane_mns(n, b, Time, Delta, Covari, U_phi_inf);
      b_s = as<vec>(b_s_result[1]);
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    mat term1 = U_pi_phi_t_z/sqrtn;
    mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<mat>(fhat_t_z(it))+cumsum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    mat term3 = cumsum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))/sqrtn;
    
    tempmat_nn = term1 - term2 - term3;
    app_path(itt) = tempmat_nn;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec mat_se_boot = stddev(as<mat>(tempmat_n2path),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min = censoring; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  double censoring = 1-sum(Delta)/n;
  double kappa_min = sqrt(censoring);
  double kappa_max = 1;
  if(kappa_min<0.5){kappa_min = 0.5;}
  
  vec kappa = {kappa_min, kappa_max};
  kappa = quantile(mat_se_boot, kappa);
  mat_se_boot.clamp(kappa(0),kappa(1));
  mat se_boot = reshape(mat_se_boot,n,n);
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempmat_nn = as<mat>(app_path(it));
    absmax_app_path(it) = abs(tempmat_nn).max();
    
    tempmat_nn /= se_boot;
    tempmat_nn.replace(arma::datum::nan,0);
    app_std_path(it) = tempmat_nn;
    absmax_app_std_path(it) = abs(tempmat_nn).max();
  }
  
  mat obs_std_path = obs_path/se_boot;
  obs_std_path.replace(arma::datum::nan,0);
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List link_mis_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      List b_s_result = dfsane_mis(n, b, Time, Delta, Covari, U_phi_inf);
      b_s = as<vec>(b_s_result[1]);
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List link_mns_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      List b_s_result = dfsane_mns(n, b, Time, Delta, Covari, U_phi_inf);
      b_s = as<vec>(b_s_result[1]);
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List form_mis_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int form, int pathsave){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  vec form_Covari = Covari.col(form-1);
  vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      List b_s_result = dfsane_mis(n, b, Time, Delta, Covari, U_phi_inf);
      b_s = as<vec>(b_s_result[1]);
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List form_mns_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int form, int pathsave){
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  vec form_Covari = Covari.col(form-1);
  vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      List b_s_result = dfsane_mns(n, b, Time, Delta, Covari, U_phi_inf);
      b_s = as<vec>(b_s_result[1]);
      tolerance = as<double>(b_s_result[0]);
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List omni_mis_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n);
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; t by x matrix
  List Mhat_i_t(n); mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    Mhat_i_t(it) = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    obs_path += (as<vec>(Mhat_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path.each_row() -= obs_path.row(0);
  // obs_path.each_col() -= obs_path.col(0);
  
  List dMhat_i_t(n);
  for(int it=0; it<n; it++){
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,as<vec>(Mhat_i_t(it))));
  }
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mis),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      vec b_s = as<vec>(b_s_opt_results[0]);
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    mat term1 = U_pi_phi_t_z/sqrtn;
    mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<mat>(fhat_t_z(it))+cumsum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    mat term3 = cumsum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))/sqrtn;
    
    tempmat_nn = term1 - term2 - term3;
    app_path(itt) = tempmat_nn;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec mat_se_boot = stddev(as<mat>(tempmat_n2path),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min = censoring; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  double censoring = 1-sum(Delta)/n;
  double kappa_min = sqrt(censoring);
  double kappa_max = 1;
  if(kappa_min<0.5){kappa_min = 0.5;}
  
  vec kappa = {kappa_min, kappa_max};
  kappa = quantile(mat_se_boot, kappa);
  mat_se_boot.clamp(kappa(0),kappa(1));
  mat se_boot = reshape(mat_se_boot,n,n);
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempmat_nn = as<mat>(app_path(it));
    absmax_app_path(it) = abs(tempmat_nn).max();
    
    tempmat_nn /= se_boot;
    app_std_path(it) = tempmat_nn;
    absmax_app_std_path(it) = abs(tempmat_nn).max();
  }
  
  mat obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List omni_mns_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; t by x matrix
  List Mhat_i_t(n); mat obs_path = zero_mat_nn;
  for(int it=0; it<n; it++){
    Mhat_i_t(it) = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    obs_path += (as<vec>(Mhat_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path.each_row() -= obs_path.row(0);
  // obs_path.each_col() -= obs_path.col(0);
  
  List dMhat_i_t(n);
  for(int it=0; it<n; it++){
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,as<vec>(Mhat_i_t(it))));
  }
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_t_z(p);
  tempvec_n = fhat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mns),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      vec b_s = as<vec>(b_s_opt_results[0]);
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_t_z = cumsum(tempmat_nn);
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    mat term1 = U_pi_phi_t_z/sqrtn;
    mat term2 = zero_mat_nn;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<mat>(fhat_t_z(it))+cumsum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
    }
    mat term3 = cumsum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))/sqrtn;
    
    tempmat_nn = term1 - term2 - term3;
    app_path(itt) = tempmat_nn;
  }
  
  NumericMatrix tempmat_n2path(pow(n,2),path);
  for(int it=0; it<path; it++){
    tempmat_n2path(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec mat_se_boot = stddev(as<mat>(tempmat_n2path),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min = censoring; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  double censoring = 1-sum(Delta)/n;
  double kappa_min = sqrt(censoring);
  double kappa_max = 1;
  if(kappa_min<0.5){kappa_min = 0.5;}
  
  vec kappa = {kappa_min, kappa_max};
  kappa = quantile(mat_se_boot, kappa);
  mat_se_boot.clamp(kappa(0),kappa(1));
  mat se_boot = reshape(mat_se_boot,n,n);
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempmat_nn = as<mat>(app_path(it));
    absmax_app_path(it) = abs(tempmat_nn).max();
    
    tempmat_nn /= se_boot;
    app_std_path(it) = tempmat_nn;
    absmax_app_std_path(it) = abs(tempmat_nn).max();
  }
  
  mat obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List link_mis_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mis),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      vec b_s = as<vec>(b_s_opt_results[0]);
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List link_mns_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  mat sorted_Covari = sort(Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (prod(Covari.row(it)<=sorted_Covari.row(itt))*1);
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mns),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      vec b_s = as<vec>(b_s_opt_results[0]);
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List form_mis_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int form, int pathsave){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  vec form_Covari = Covari.col(form-1);
  vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mis),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      vec b_s = as<vec>(b_s_opt_results[0]);
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

List form_mns_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int form, int pathsave){
  
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int n = Covari.n_rows;
  int p = Covari.n_cols;
  
  double sqrtn = sqrt(n);
  
  vec zero_vec_1 = zeros(1);
  vec zero_vec_p = zeros(p);
  vec zero_vec_n = zeros(n);
  mat zero_mat_np = zeros(n,p);
  mat zero_mat_nn = zeros(n,n);
  
  vec tempvec_p(p);
  vec tempvec_n(n);
  mat tempmat_np(n,p);
  mat tempmat_nn(n,n);
  
  vec resid = log(Time) + Covari*b;
  
  uvec index_resid = sort_index(resid);
  
  Time = Time(index_resid);
  Delta = Delta(index_resid);
  Covari = Covari.rows(index_resid);
  resid = resid(index_resid);
  
  List pi_i_z(n); List N_i_t(n); List Y_i_t(n); 
  vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
  vec form_Covari = Covari.col(form-1);
  vec sorted_form_Covari = sort(form_Covari);
  tempvec_n = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
    }
    pi_i_z(it) = tempvec_n;
    N_i_t(it) = (resid>=resid(it))*Delta(it);
    Y_i_t(it) = (resid<=resid(it))*1;
    S_0_t += as<vec>(Y_i_t(it));
    S_1_t += as<vec>(Y_i_t(it))*(Covari.row(it));
    S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
  }
  
  vec Lambdahat_0_t = cumsum(Delta/S_0_t);
  vec dLambdahat_0_t = diff(join_cols(zero_vec_1,Lambdahat_0_t));
  
  mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
  
  // obs_path; 1 by x vector
  List Mhat_i_t(n); List dMhat_i_t(n); vec obs_path = zero_vec_n;
  for(int it=0; it<n; it++){
    tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
    Mhat_i_t(it) = tempvec_n;
    dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
    obs_path += (tempvec_n(n-1))*(as<vec>(pi_i_z(it)));
  }
  obs_path /= sqrtn;
  // obs_path -= obs_path(0) * ones(n);
  
  // -----------------------------------------------------------
  // ----------------------Kernel Smoothing---------------------
  // -----------------------------------------------------------
  double bw_base = 1.06 * pow(n,-0.2);
  vec pred_data = exp(resid);
  
  // -----------------------------g0----------------------------
  // vec given_data_g = exp(resid);
  vec given_data_g = pred_data;
  double bw_gn = bw_base * stddev(given_data_g);
  vec ghat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      ghat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn);
    }
  }
  ghat_0_t /= (n*bw_gn);
  
  List ghat_t_z(p);
  tempvec_n = ghat_0_t%Time;
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
    }
    ghat_t_z(itt) = tempmat_nn/n;
  }
  
  // -----------------------------f0----------------------------
  vec Fhat_0_e = 1-cumprod(1-Delta/S_0_t);
  vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
  
  vec Condi_Ehat = zero_vec_n;
  for(int it=0; it<n; it++){
    Condi_Ehat(it) = sum(join_cols(zeros(it+1),ones(n-it-1))%resid%dFhat_0_e)/(1-Fhat_0_e(it));
  }
  Condi_Ehat.replace(datum::nan,0);
  
  vec rhat_i = Delta%resid+(1-Delta)%Condi_Ehat;
  vec given_data_f = exp(rhat_i);
  double bw_fn = bw_base * stddev(given_data_f);
  vec fhat_0_t = zero_vec_n;
  for(int it=0; it<n; it++){
    for(int itt=0; itt<n; itt++){
      fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn) * dFhat_0_e(itt);
    }
  }
  fhat_0_t /= bw_fn;
  
  List fhat_inf_z(p);
  double tempvec_1 = fhat_0_t(n-1)*Time(n-1);
  for(int itt=0; itt<p; itt++){
    vec Covari_col = Covari.col(itt);
    tempvec_n = zero_vec_n;
    for(int it=0; it<n; it++){
      tempvec_n += tempvec_1*((as<vec>(pi_i_z(it)))*(Delta(it)*Covari_col(it)));
    }
    fhat_inf_z(itt) = tempvec_n/n;
  }
  
  // -----------------------------------------------------------
  // ------------------------Sample Path------------------------
  // -----------------------------------------------------------
  // vec one_vec_n = ones(n);
  
  List app_path(path);
  for(int itt=0; itt<path; itt++){
    
    vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
    while(tolerance>tol){
      phi_i = randn(n);
      // phi_i = randg(n) - one_vec_n;
      
      tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
      for(int it=0; it<n; it++){
        tempvec_n += as<vec>(dMhat_i_t(it))*phi_i(it);
        tempmat_np += (as<vec>(dMhat_i_t(it))*(Covari.row(it)))*phi_i(it);
      }
      vec U_phi_inf = sum(((S_0_t%tempmat_np.each_col())-(S_1_t.each_col()%tempvec_n)),0).t();
      
      Rcpp::List b_s_opt_results = optim(Rcpp::_["par"]    = b,
                                         Rcpp::_["fn"]     = Rcpp::InternalFunction(&target_score2_mns),
                                         Rcpp::_["method"] = optimType,
                                         Rcpp::_["Time"] = Time,
                                         Rcpp::_["Delta"] = Delta,
                                         Rcpp::_["Covari"] = Covari,
                                         Rcpp::_["targetvector"] = U_phi_inf);
      
      vec b_s = as<vec>(b_s_opt_results[0]);
      tolerance = as<double>(b_s_opt_results[1]);;
    }
    
    tempmat_nn = zero_mat_nn;
    for(int it=0; it<n; it++){
      tempmat_nn += ((((as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col())%(as<vec>(dMhat_i_t(it))))*phi_i(it));
    }
    mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
    
    vec resid_s = log(Time) + Covari*b_s;
    uvec index_resid_s = sort_index(resid_s);
    
    vec Delta_s = Delta(index_resid_s);
    resid_s = resid_s(index_resid_s);
    
    NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
    for(int it=0; it<n; it++){
      Y_i_t_s = (resid_s<=resid_s(it))*1;
      S_0_t_s += as<vec>(Y_i_t_s);
    }
    vec Lambdahat_0_t_s = cumsum(Delta_s/S_0_t_s);
    vec dLambdahat_0_t_s = diff(join_cols(zero_vec_1,Lambdahat_0_t_s));
    
    vec term1 = U_pi_phi_inf_z/sqrtn;
    vec term2 = zero_vec_n;
    tempvec_p = (b-b_s)*sqrtn;
    for(int it=0; it<p; it++){
      term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
    }
    vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
    
    tempvec_n = term1 - term2 - term3;
    app_path(itt) = tempvec_n;
  }
  
  NumericMatrix tempmat_npath(n,path);
  for(int it=0; it<path; it++){
    tempmat_npath(_,it) = (as<NumericVector>(app_path(it)));
  }
  vec se_boot = stddev(as<mat>(tempmat_npath),0,1);
  // too low values which are 0 or computationally 0 of se_boot makes a problem, 
  // so we adjust them to have kappa = quantile of mat_se_boot
  // e.g., kappa_min =; sqrt(censoring)/2; quantile(mat_se_boot) = {0.2, 1};
  vec kappa = {0.2, 1};
  kappa = quantile(se_boot, kappa);
  se_boot.clamp(kappa(0),kappa(1));
  
  List app_std_path(path); vec absmax_app_path(path); vec absmax_app_std_path(path);
  for(int it=0; it<path; it++){
    tempvec_n = as<vec>(app_path(it));
    absmax_app_path(it) = abs(tempvec_n).max();
    
    tempvec_n /= se_boot;
    app_std_path(it) = tempvec_n;
    absmax_app_std_path(it) = abs(tempvec_n).max();
  }
  
  vec obs_std_path = obs_path/se_boot;
  double absmax_obs_path = (abs(obs_path)).max();
  double absmax_obs_std_path = (abs(obs_std_path)).max();
  
  uvec ind_usd = (find(absmax_app_path>absmax_obs_path));
  double p_value = (ind_usd.size()); p_value = p_value/path;
  
  uvec ind_std = (find(absmax_app_std_path>absmax_obs_std_path));
  double p_std_value = (ind_std.size()); p_std_value = p_std_value/path;
  
  if(path > pathsave){
    pathsave = pathsave - 1;
    app_path = app_path[Range(0,pathsave)];
    app_std_path = app_std_path[Range(0,pathsave)];
  }
  
  return List::create(_["beta"]=b,_["SE_boot"]=se_boot,
                      _["obs_path"]=obs_path,_["obs_std_path"]=obs_std_path,
                      _["app_path"]=app_path,_["app_std_path"]=app_std_path,
                      _["p_value"]=p_value,_["p_std_value"]=p_std_value);
}

// omni_mis_DFSANE
List omni_mis_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave);
RcppExport SEXP _afttest_omni_mis_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mis_DFSANE(path, b, Time, Delta, Covari, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// omni_mns_DFSANE
List omni_mns_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave);
RcppExport SEXP _afttest_omni_mns_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mns_DFSANE(path, b, Time, Delta, Covari, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// link_mis_DFSANE
List link_mis_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave);
RcppExport SEXP _afttest_link_mis_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mis_DFSANE(path, b, Time, Delta, Covari, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// link_mns_DFSANE
List link_mns_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int pathsave);
RcppExport SEXP _afttest_link_mns_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mns_DFSANE(path, b, Time, Delta, Covari, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// form_mis_DFSANE
List form_mis_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int form, int pathsave);
RcppExport SEXP _afttest_form_mis_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP formSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mis_DFSANE(path, b, Time, Delta, Covari, form, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// form_mns_DFSANE
List form_mns_DFSANE(int path, vec b, vec Time, vec Delta, mat Covari, int form, int pathsave);
RcppExport SEXP _afttest_form_mns_DFSANE(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP formSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mns_DFSANE(path, b, Time, Delta, Covari, form, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// omni_mis_optim
List omni_mis_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave);
RcppExport SEXP _afttest_omni_mis_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mis_optim(path, b, Time, Delta, Covari, optimType, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// omni_mns_optim
List omni_mns_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave);
RcppExport SEXP _afttest_omni_mns_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(omni_mns_optim(path, b, Time, Delta, Covari, optimType, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// link_mis_optim
List link_mis_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave);
RcppExport SEXP _afttest_link_mis_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mis_optim(path, b, Time, Delta, Covari, optimType, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// link_mns_optim
List link_mns_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int pathsave);
RcppExport SEXP _afttest_link_mns_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(link_mns_optim(path, b, Time, Delta, Covari, optimType, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// form_mis_optim
List form_mis_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int form, int pathsave);
RcppExport SEXP _afttest_form_mis_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP formSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mis_optim(path, b, Time, Delta, Covari, optimType, form, pathsave));
  return rcpp_result_gen;
  END_RCPP
}
// form_mns_optim
List form_mns_optim(int path, vec b, vec Time, vec Delta, mat Covari, String optimType, int form, int pathsave);
RcppExport SEXP _afttest_form_mns_optim(SEXP pathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP optimTypeSEXP, SEXP formSEXP, SEXP pathsaveSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::RNGScope rcpp_rngScope_gen;
  Rcpp::traits::input_parameter< int >::type path(pathSEXP);
  Rcpp::traits::input_parameter< vec >::type b(bSEXP);
  Rcpp::traits::input_parameter< vec >::type Time(TimeSEXP);
  Rcpp::traits::input_parameter< vec >::type Delta(DeltaSEXP);
  Rcpp::traits::input_parameter< mat >::type Covari(CovariSEXP);
  Rcpp::traits::input_parameter< String >::type optimType(optimTypeSEXP);
  Rcpp::traits::input_parameter< int >::type form(formSEXP);
  Rcpp::traits::input_parameter< int >::type pathsave(pathsaveSEXP);
  rcpp_result_gen = Rcpp::wrap(form_mns_optim(path, b, Time, Delta, Covari, optimType, form, pathsave));
  return rcpp_result_gen;
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_afttest_omni_mis_DFSANE", (DL_FUNC) &_afttest_omni_mis_DFSANE, 6},
  {"_afttest_omni_mns_DFSANE", (DL_FUNC) &_afttest_omni_mns_DFSANE, 6},
  {"_afttest_link_mis_DFSANE", (DL_FUNC) &_afttest_link_mis_DFSANE, 6},
  {"_afttest_link_mns_DFSANE", (DL_FUNC) &_afttest_link_mns_DFSANE, 6},
  {"_afttest_form_mis_DFSANE", (DL_FUNC) &_afttest_form_mis_DFSANE, 7},
  {"_afttest_form_mns_DFSANE", (DL_FUNC) &_afttest_form_mns_DFSANE, 7},
  {"_afttest_omni_mis_optim", (DL_FUNC) &_afttest_omni_mis_optim, 7},
  {"_afttest_omni_mns_optim", (DL_FUNC) &_afttest_omni_mns_optim, 7},
  {"_afttest_link_mis_optim", (DL_FUNC) &_afttest_link_mis_optim, 7},
  {"_afttest_link_mns_optim", (DL_FUNC) &_afttest_link_mns_optim, 7},
  {"_afttest_form_mis_optim", (DL_FUNC) &_afttest_form_mis_optim, 8},
  {"_afttest_form_mns_optim", (DL_FUNC) &_afttest_form_mns_optim, 8},
  {NULL, NULL, 0}
};

RcppExport void R_init_afttest(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
