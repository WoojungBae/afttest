#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
//' @useDynLib afttest, .registration = TRUE
 //' @importFrom Rcpp evalCpp
 //' @exportPattern "^[[:alpha:]]+"
 
#ifdef RCPP_USE_GLOBAL_ROSTREAM
 Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
 Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif
 
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 mat inv_cpp(mat const& MAT) {
   int p_temp = MAT.n_cols;
   int rank_temp = rank(MAT);
   
   mat MATinv;
   if(p_temp > rank_temp){
     MATinv = pinv(MAT);
   } else if(MAT.is_sympd()){
     MATinv = inv_sympd(MAT);
   } else {
     MATinv = inv(MAT);
   }
   
   return MATinv;
 }
 
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 inline double get_rikjl(const arma::rowvec& xdif, const arma::mat& sigma) {
   double val = as_scalar(xdif * sigma * xdif.t());
   return (val > 1e-12) ? sqrt(val) : 0.0;
 }
 
 arma::mat abar_gehan_cpp(arma::vec beta, 
                          arma::vec Y, 
                          arma::mat X, 
                          arma::vec delta, 
                          arma::mat sigma, 
                          arma::vec weights, 
                          arma::vec gehanWeights) {
   int N = X.n_rows;
   int p = X.n_cols;
   double sqrtn = sqrt((double)N);
   
   arma::vec e = Y + X * beta;
   arma::mat abar = zeros<mat>(p, p);
   
   for (int i = 0; i < N; i++) {
     if (delta(i) != 0) {
       for (int j = 0; j < N; j++) {
         
         arma::rowvec xdif = X.row(i) - X.row(j);
         double rikjl = get_rikjl(xdif, sigma);
         
         if (rikjl > 1e-10) {
           // Derivative logic
           double edif = e(j) - e(i); // Direction matches score
           double z = sqrtn * edif / rikjl;
           double h = R::dnorm(z, 0.0, 1.0, 0); 
           
           double coef = gehanWeights(i) * weights(i) * weights(j) * h * sqrtn / rikjl;
           
           // CRITICAL FIX: Subtract outer product (Negative Definite Slope)
           abar -= coef * (xdif.t() * xdif);
         }
       }
     }
   }
   return abar/N;
 }
 
 arma::mat omega_gehan_cpp(arma::vec beta, 
                           arma::vec Y, 
                           arma::mat X, 
                           arma::vec delta, 
                           arma::vec weights) {
   int N = X.n_rows;
   int p = X.n_cols;
   arma::vec e = Y + X * beta;
   arma::mat ksi = zeros<mat>(N, p);
   
   for (int i = 0; i < N; i++) {
     for (int j = 0; j < N; j++) {
       
       // Logic 1: e[i] < e[j]
       if (delta(i) != 0) {
         if (e(i) < e(j)) {
           ksi.row(i) += delta(i) * (X.row(i) - X.row(j)) * weights(j) / (double)N;
         }
       }
       
       // Logic 2: e[i] >= e[j]
       if (delta(j) != 0) {
         if (e(i) >= e(j)) {
           arma::rowvec xdif_sum = zeros<rowvec>(p);
           double emk = 0;
           for (int r = 0; r < N; r++) { 
             if (e(r) >= e(j)) {
               xdif_sum += weights(r) * (X.row(i) - X.row(r));
               emk++;
             }
           }
           if (emk > 0) {
             ksi.row(i) -= xdif_sum / ((double)N * emk);
           }
         }
       }
     }
   }
   return ksi.t() * ksi;
 }
 
 arma::vec score_gehan_is_cpp(arma::vec beta, 
                              arma::vec Y, 
                              arma::mat X, 
                              arma::vec delta, 
                              arma::mat sigma, 
                              arma::vec weights) {
   int N = X.n_rows;
   int p = X.n_cols;
   double sqrtn = sqrt((double)N);
   
   arma::vec e = Y + X * beta;
   arma::vec score = zeros<vec>(p);
   
   for (int i = 0; i < N; i++) {
     if (delta(i) != 0) {
       for (int j = 0; j < N; j++) {
         
         arma::rowvec xdif = X.row(i) - X.row(j);
         double rikjl = get_rikjl(xdif, sigma);
         
         if (rikjl > 1e-10) {
           // Score logic: Approx Indicator(e_i < e_j)
           // When e_i < e_j, (e_j - e_i) is positive -> pnorm ~ 1
           double edif = e(j) - e(i); 
           double z = sqrtn * edif / rikjl;
           double phi = R::pnorm(z, 0.0, 1.0, 1, 0); 
           
           score += weights(i) * weights(j) * xdif.t() * phi;
         }
       }
     }
   }
   return score;
 }
 
 arma::vec score_gehan_ns_cpp(arma::vec beta, 
                              arma::vec Y, 
                              arma::mat X, 
                              arma::vec delta, 
                              arma::vec weights) {
   int N = X.n_rows;
   int p = X.n_cols;
   
   arma::vec e = Y + X * beta;
   arma::vec score = zeros<vec>(p);
   
   for (int i = 0; i < N; i++) {
     if (delta(i) != 0) {
       for (int j = 0; j < N; j++) {
         // Exact Indicator I(e_i < e_j)
         if (e(i) < e(j)) {
           arma::rowvec xdif = X.row(i) - X.row(j);
           score += weights(i) * weights(j) * xdif.t();
         }
       }
     }
   }
   return score;
 }
 
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // Weighted Kaplan-Meier & Imputation
 vec impute_residuals_weighted(vec e, vec delta, vec weights) {
   int n = e.n_elem;
   vec zero_vec_n = zeros(n);
   
   // 1. Sort Data by Residuals
   uvec index_e = sort_index(e);
   vec e_sorted = e(index_e);
   vec d_sorted = delta(index_e);
   vec w_sorted = weights(index_e);
   
   // 2. Identify Unique Times for KM
   // We use a temporary std::vector for dynamic push_back, then convert to arma
   std::vector<double> t_unique_std;
   std::vector<double> w_event_std;
   std::vector<double> w_risk_std;
   
   double current_t = e_sorted(0);
   double sum_w_risk = sum(w_sorted); 
   double block_w_event = 0.0;
   double block_w_loss = 0.0;
   
   for(int it=0; it<n; it++){
     // If new time point found (and not the first iteration)
     if (abs(e_sorted(it) - current_t) > 1e-9) {
       t_unique_std.push_back(current_t);
       w_event_std.push_back(block_w_event);
       w_risk_std.push_back(sum_w_risk);
       
       // Update Risk for next step
       sum_w_risk -= block_w_loss;
       
       // Reset for next block
       current_t = e_sorted(it);
       block_w_event = 0.0;
       block_w_loss = 0.0;
     }
     
     if (d_sorted(it) == 1) block_w_event += w_sorted(it);
     block_w_loss += w_sorted(it);
   }
   // Push last block
   t_unique_std.push_back(current_t);
   w_event_std.push_back(block_w_event);
   w_risk_std.push_back(sum_w_risk);
   
   // Convert to Armadillo vectors
   vec time_unique = conv_to<vec>::from(t_unique_std);
   vec n_event = conv_to<vec>::from(w_event_std);
   vec n_risk = conv_to<vec>::from(w_risk_std);
   
   int n_times = time_unique.n_elem;
   
   // 3. Compute Survival Function S(t)
   vec surv_unique = ones(n_times);
   double S_prev = 1.0;
   
   for(int it=0; it<n_times; it++){
     double haz = (n_risk(it) > 1e-9) ? (n_event(it) / n_risk(it)) : 0.0;
     double S_curr = S_prev * (1.0 - haz);
     surv_unique(it) = S_curr;
     S_prev = S_curr;
   }
   
   // 4. Compute Area (Restricted Mean)
   vec area = zeros(n_times);
   double cum_area = 0.0;
   
   // Calculate dt (difference between times)
   vec dt = zeros(n_times);
   if (n_times > 1) {
     dt.head(n_times-1) = diff(time_unique);
     // Handle tail: assume last dt is same as previous min non-zero dt
     double min_dt = dt(0); 
     for(int k=0; k<n_times-1; k++) if(dt(k) > 0 && dt(k) < min_dt) min_dt = dt(k);
     dt(n_times-1) = min_dt;
   }
   
   // Backward loop for cumulative area
   for(int it=n_times-1; it>=0; it--){
     cum_area += surv_unique(it) * dt(it);
     area(it) = cum_area;
   }
   
   // 5. Impute
   vec e_star = e; // Use original unsorted e
   
   for(int it=0; it<n; it++){
     if (delta(it) == 0) {
       double val = e(it);
       
       // Binary search for time index
       // std::lower_bound works on pointers to arma vector data
       const double* ptr_start = time_unique.memptr();
       const double* ptr_end = ptr_start + n_times;
       const double* ptr_found = std::lower_bound(ptr_start, ptr_end, val);
       int idx = ptr_found - ptr_start;
       
       if (idx >= n_times) idx = n_times - 1;
       
       double S_val = surv_unique(idx);
       double Area_val = area(idx);
       
       if (S_val > 1e-9) {
         e_star(it) += Area_val / S_val;
       }
     }
   }
   
   return e_star;
 }
 
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 vec target_score_ls(vec b, vec time, vec delta, mat covariates, 
                     vec targetvector, vec weights, int n, int p, double sqrtn) {
   
   vec resid = log(time) + covariates * b;
   vec resid_star = impute_residuals_weighted(resid, delta, ones(n));
   
   rowvec X_bar = mean(covariates, 0);
   mat Xc = covariates;
   Xc.each_row() -= X_bar;
   
   vec F_val = ((Xc.t() * resid_star) - targetvector)/n;
   
   return F_val;
 }
 
 vec target_score_is(vec b, vec time, vec delta, mat covariates, 
                      vec targetvector, vec weights, int n, int p, double sqrtn){
   
   vec resid = log(time) + covariates*b;
   uvec index_resid = sort_index(resid);
   
   delta = delta(index_resid);
   covariates = covariates.rows(index_resid);
   resid = resid(index_resid);
   
   mat tempmat_np = zeros(n,p); vec tempvec_n = zeros(n); vec F_vec = zeros(p);
   for(int it=0; it<n; it++){
     if (delta(it)==1){
       tempmat_np = covariates.row(it) - covariates.each_row();
       tempvec_n = sqrt(sum(tempmat_np%tempmat_np,1));
       tempvec_n.replace(0,1);
       
       tempvec_n = normcdf(sqrtn*(resid-resid(it))/tempvec_n);
       F_vec += sum(tempmat_np.each_col()%tempvec_n).t();
     }
   }
   F_vec = (F_vec - targetvector)/n;
   
   return F_vec;
 }
 
 vec target_score_ns(vec b, vec time, vec delta, mat covariates, 
                      vec targetvector, vec weights, int n, int p, double sqrtn){
   
   vec resid = log(time) + covariates*b;
   
   uvec index_resid = sort_index(resid);
   
   delta = delta(index_resid);
   covariates = covariates.rows(index_resid);
   resid = resid(index_resid);
   
   mat tempmat_np = zeros(n,p); vec F_vec = zeros(p);
   for(int it=0; it<n; it++){
     if (delta(it)==1){
       tempmat_np = covariates.row(it) - covariates.each_row();
       F_vec += sum(tempmat_np.each_col()%conv_to<vec>::from((resid>=resid(it)))).t();
     }
   }
   F_vec = (F_vec - targetvector)/n;
   
   return F_vec;
 }
 
 List dfsane(vec b, vec time, vec delta, mat covariates, 
             vec targetvector, vec weights, int n, int p, double sqrtn, 
             std::string eqType){
   
   typedef vec (*score_func_ptr)(vec, vec, vec, mat, vec, vec, int, int, double);
   score_func_ptr get_score;
   
   if (eqType == "ls") {
     get_score = &target_score_ls;
   } else if (eqType == "is") {
     get_score = &target_score_is;
   } else { // "ns"
     get_score = &target_score_ns;
   }
   
   vec b_old = b;
   vec F_old = target_score_ns(b_old,time,delta,covariates,targetvector,weights,n,p,sqrtn);
   double sig_k = (1/sqrt(dot(F_old, F_old))); if (sig_k>1){sig_k=1;}
   
   vec b_new = b_old - sig_k*F_old;
   
   vec F_new = get_score(b_new,time,delta,covariates,targetvector,weights,n,p,sqrtn);
   
   vec s_k = b_new - b_old;
   vec y_k = F_new - F_old;
   
   double tol_0 = dot(F_old, F_old);
   double tol_s = dot(s_k, s_k);
   double tol_y = dot(y_k, y_k);
   double tol_f = dot(F_new, F_new);
   
   // Stop sqrt(tol_f)/sqrt(n) <= e_a + e_r * sqrt(tol_0)/sqrt(n)
   // Stop tol_f <= (e_a * sqrt(n) + e_r * tol_0)^{2}
   // double e_a = 1e-6; double e_r = 1e-5;
   double e_a = 1e-5; double e_r = 1e-4;
   double optim_tol = pow(e_a * sqrt(n) + e_r * sqrt(tol_0),2);
   
   double tolerance=tol_f+1; double tau_min=0.1; double tau_max=0.5; 
   double sig_min=0.1; double sig_max=0.5; double alp_p=1; double alp_m=1; double gam=1e-4; 
   double M=1; vec f_bar=zeros(M); double it=1; double maxit=1e3;
   double eta_k, abssig_k, RHS_p, LHS_p, RHS_m, LHS_m, alp_p_t, alp_m_t;
   vec b_new_p, F_new_p, b_new_m, F_new_m;
   while(tolerance>optim_tol){
     
     // STEP 1
     eta_k = tol_0/pow(1+it,2);
     
     if (tol_y>0) {
       sig_k = dot(s_k, y_k)/tol_y;
     } 
     
     abssig_k = std::abs(sig_k);
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
     int step_tol = 0; int itt = it - M * floor(it/M); f_bar(itt) = tol_f; double a_k = 0;
     while(step_tol == 0){
       
       // alpha_plus
       b_new_p = b_new + alp_p * d_k;
       F_new_p = get_score(b_new_p,time,delta,covariates,targetvector,weights,n,p,sqrtn);
       
       RHS_p = dot(F_new_p, F_new_p);
       LHS_p = f_bar.max() + eta_k - gam * pow(alp_p,2) * tol_f;
       
       // alpha_minus
       b_new_m = b_new - alp_m * d_k;
       F_new_m = get_score(b_new_m,time,delta,covariates,targetvector,weights,n,p,sqrtn);
       
       RHS_m = dot(F_new_m, F_new_m);
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
     F_new = get_score(b_new,time,delta,covariates,targetvector,weights,n,p,sqrtn);
     
     s_k = b_new - b_old;
     y_k = F_new - F_old;
     
     b_old = b_new;
     F_old = F_new;    
     
     tol_s = dot(s_k, s_k);
     tol_y = dot(y_k, y_k);
     tol_f = dot(F_new, F_new);
     
     tolerance = tol_f;
     if (tol_f>tol_s){tolerance = tol_s;}
     if (it>maxit){tolerance = 0;}
     it += 1;
   }
   
   return List::create(tol_f,b_new);
 }
 
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 List omni_cpp(int npath, vec b, vec time, vec delta, mat covariates, 
                  int npathsave, std::string eqType, 
                  bool fast_boot, mat invOmega){
   
   int n = covariates.n_rows;
   int p = covariates.n_cols;
   
   double sqrtn = sqrt(n);
   
   vec zero_vec_1 = zeros(1);
   vec zero_vec_p = zeros(p);
   vec zero_vec_n = zeros(n);
   mat zero_mat_np = zeros(n,p);
   mat zero_mat_nn = zeros(n,n);
   
   vec one_vec_n = ones(n);
   
   vec tempvec_p(p);
   vec tempvec_n(n);
   mat tempmat_np(n,p);
   mat tempmat_nn(n,n);
   
   vec resid = log(time) + covariates*b;
   uvec index_resid = sort_index(resid);
   
   time = time(index_resid);
   delta = delta(index_resid);
   covariates = covariates.rows(index_resid);
   resid = resid(index_resid);
   
   // -----------------------------------------------------------
   List pi_i_z(n); List N_i_t(n); List Y_i_t(n);
   vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
   mat sorted_Covari = sort(covariates);
   tempvec_n = zero_vec_n;
   for(int it=0; it<n; it++){
     for(int itt=0; itt<n; itt++){
       tempvec_n(itt) = (prod(covariates.row(it)<=sorted_Covari.row(itt))*1);
     }
     pi_i_z(it) = tempvec_n;
     if (delta(it)==1){
       N_i_t(it) = (resid>=resid(it));
     } else {
       N_i_t(it) = zero_vec_n;
     }
     Y_i_t(it) = (resid<=resid(it))*1;
     S_0_t += as<vec>(Y_i_t(it));
     S_1_t += as<vec>(Y_i_t(it))*(covariates.row(it));
     S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
   }
   
   vec dLambdahat_0_t = delta/S_0_t;
   dLambdahat_0_t.replace(datum::nan,0);
   
   mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
   E_pi_t_z.replace(datum::nan,0);
   
   // obs_npath; t by x vector
   List Mhat_i_t(n); List dMhat_i_t(n); mat obs_npath = zero_mat_nn;
   for(int it=0; it<n; it++){
     tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
     Mhat_i_t(it) = tempvec_n;
     dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
     obs_npath += tempvec_n*(as<rowvec>(pi_i_z(it)));
   }
   obs_npath /= sqrtn;
   // rowvec obs_npath_row0 = obs_npath.row(0);
   // obs_npath.each_row() -= obs_npath_row0;
   // obs_npath.each_col() -= obs_npath.col(0);
   
   // -----------------------------------------------------------
   // ----------------------Kernel Smoothing---------------------
   // -----------------------------------------------------------
   double bw_base = pow((n*3/4),-0.2);
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
   ghat_0_t /= n;
   
   List ghat_t_z(p);
   tempvec_n = ghat_0_t % time;
   for(int itt=0; itt<p; itt++){
     vec Covari_col = covariates.col(itt);
     tempmat_nn = zero_mat_nn;
     for(int it=0; it<n; it++){
       tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
     }
     ghat_t_z(itt) = tempmat_nn/n;
   }
   
   // -----------------------------f0----------------------------
   vec Shat_0_e = cumprod(one_vec_n - dLambdahat_0_t);
   vec Fhat_0_e = one_vec_n - Shat_0_e;
   vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
   vec fhat_0_t = zero_vec_n;
   for(int it=0; it<n; it++){
     for(int itt=0; itt<n; itt++){
       fhat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn) * dFhat_0_e(itt);
     }
   }
   
   // vec Condi_Ehat = zero_vec_n;
   // for(int it=0; it<n; it++){
   //   Condi_Ehat(it) = sum(resid.tail(n-it-1)%dFhat_0_e.tail(n-it-1))/Shat_0_e(it);
   // }
   // Condi_Ehat.replace(datum::nan,0);
   // 
   // vec rhat_i = delta % resid + (one_vec_n - delta) % Condi_Ehat;
   // vec given_data_f = exp(rhat_i);
   // double bw_fn = bw_base * stddev(given_data_f);
   // vec fhat_0_t = zero_vec_n;
   // for(int it=0; it<n; it++){
   //   for(int itt=0; itt<n; itt++){
   //     fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn);
   //   }
   // }
   // fhat_0_t /= n;
   
   List fhat_t_z(p);
   tempvec_n = fhat_0_t % time;
   for(int itt=0; itt<p; itt++){
     vec Covari_col = covariates.col(itt);
     tempmat_nn = zero_mat_nn;
     for(int it=0; it<n; it++){
       if (delta(it)==1){
         tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
       }
     }
     fhat_t_z(itt) = tempmat_nn/n;
   }
   
   List app_npath(npath);
   if (fast_boot) {
     // -----------------------------------------------------------
     // --------------------------dkappa_t-------------------------
     // -----------------------------------------------------------
     vec lambda_hat_0_t = fhat_0_t / Shat_0_e;
     lambda_hat_0_t.replace(datum::nan, 0);
     lambda_hat_0_t.replace(datum::inf, 0);

     vec lambda_hat_0_tTIMESt = lambda_hat_0_t % pred_data;
     vec dlambda_hat_0_tTIMESt = diff(join_cols(zero_vec_1, lambda_hat_0_tTIMESt));
     mat dkappa_t = - S_1_t;
     dkappa_t.each_col() %= dlambda_hat_0_tTIMESt/S_0_t;

     List term2(p);
     for(int it=0; it<p; it++){
       term2(it) = (as<mat>(fhat_t_z(it)) +
         cumsum((as<mat>(ghat_t_z(it)).each_col()) % dLambdahat_0_t) +
         cumsum((S_pi_t_z.each_col()) % (dkappa_t.col(it)))/n);
     }
     
     // -----------------------------------------------------------
     // ------------------------Sample Path------------------------
     // -----------------------------------------------------------
     for(int itt=0; itt<npath; itt++){
       vec phi_i = randg(n) - one_vec_n; // randn(n);
       
       tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
       for(int i=0; i<n; i++){
         vec dM_phi = as<vec>(dMhat_i_t(i)) * phi_i(i);
         tempvec_n += dM_phi;
         tempmat_np += dM_phi * covariates.row(i);
       }
       tempvec_p = invOmega * ((tempmat_np.t() * S_0_t) - (S_1_t.t() * tempvec_n));
       
       mat term3 = zero_mat_nn;
       for (int it=0; it<p; it++) {
         tempmat_nn += as<mat>(term2(it)) * tempvec_p(it)/sqrtn;
       }
       
       tempmat_nn = zero_mat_nn;
       for(int it=0; it<n; it++){
         tempmat_nn += ((as<vec>(dMhat_i_t(it)))*phi_i(it))%(as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col();
       }
       mat U_pi_phi_t_z = cumsum(tempmat_nn);
       
       app_npath(itt) = (U_pi_phi_t_z - term3)/sqrtn;
     }
   } else {
     // -----------------------------------------------------------
     // --------------------------U_inf_ls-------------------------
     // -----------------------------------------------------------
     mat U_inf_ls;
     if (eqType == "ls") {
       mat E_t = S_1_t.each_col()/S_0_t;
       vec resid_star = impute_residuals_weighted(resid, delta, one_vec_n);
       rowvec X_bar = mean(covariates, 0);
       mat Xc = covariates;
       Xc.each_row() -= X_bar;
       
       U_inf_ls = Xc;
       U_inf_ls.each_col() %= resid_star/n;
       
       vec lambda_hat_0_t = fhat_0_t / Shat_0_e;
       lambda_hat_0_t.replace(datum::nan, 0);
       lambda_hat_0_t.replace(datum::inf, 0);
       
       vec lambda_hat_0_tTIMESt = lambda_hat_0_t % pred_data;
       vec dlambda_hat_0_tTIMESt = diff(join_cols(zero_vec_1, lambda_hat_0_tTIMESt));
       mat dkappa_t = - E_t;
       dkappa_t.each_col() %= dlambda_hat_0_tTIMESt;
     }
     
     // -----------------------------------------------------------
     // ------------------------Sample Path------------------------
     // -----------------------------------------------------------
     for(int itt=0; itt<npath; itt++){
       
       vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
       while(tolerance>tol){
         phi_i = randg(n) - one_vec_n; // randn(n);
         
         vec U_phi_inf(n); 
         if (eqType == "ls") {
           U_phi_inf = U_inf_ls.t() * phi_i;
         } else {
           tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
           for(int i=0; i<n; i++){
             vec dM_phi = as<vec>(dMhat_i_t(i)) * phi_i(i);
             tempvec_n += dM_phi;
             tempmat_np += dM_phi * covariates.row(i);
           }
           U_phi_inf = (tempmat_np.t() * S_0_t) - (S_1_t.t() * tempvec_n);
         }
         
         List b_s_result = dfsane(b, time, delta, covariates, U_phi_inf, phi_i, n, p, sqrtn, eqType);
         tolerance = as<double>(b_s_result[0]);
         b_s = as<vec>(b_s_result[1]);
       }
       
       tempmat_nn = zero_mat_nn;
       for(int it=0; it<n; it++){
         tempmat_nn += ((as<vec>(dMhat_i_t(it)))*phi_i(it))%(as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col();
       }
       mat U_pi_phi_t_z = cumsum(tempmat_nn);
       
       vec resid_s = log(time) + covariates*b_s;
       uvec index_resid_s = sort_index(resid_s);
       
       vec Delta_s = delta(index_resid_s);
       resid_s = resid_s(index_resid_s);
       
       NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
       for(int it=0; it<n; it++){
         Y_i_t_s = (resid_s<=resid_s(it))*1;
         S_0_t_s += as<vec>(Y_i_t_s);
       }
       
       vec dLambdahat_0_t_s = Delta_s/S_0_t_s;
       dLambdahat_0_t_s.replace(datum::nan,0);
       
       mat term1 = U_pi_phi_t_z/sqrtn;
       mat term2 = zero_mat_nn;
       tempvec_p = (b-b_s)*sqrtn;
       for(int it=0; it<p; it++){
         term2 += (as<mat>(fhat_t_z(it))+cumsum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t))*(tempvec_p(it));
       }
       mat term3 = cumsum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))/sqrtn;
       
       app_npath(itt) = term1 - term2 - term3;
     }
   }
   
   NumericMatrix tempmat_n2npath(pow(n,2),npath);
   for(int it=0; it<npath; it++){
     tempmat_n2npath(_,it) = (as<NumericVector>(app_npath(it)));
   }
   vec mat_se_boot = stddev(as<mat>(tempmat_n2npath),0,1);
   // too low values which are 0 or computationally 0 of se_boot makes a problem,
   // so we adjust them to have kappa = quantile of mat_se_boot
   // e.g., kappa_min = censoring; quantile(mat_se_boot) = {censoring, 1};
   // double censoring = 1-sum(delta)/n;
   // double kappa_min = censoring;
   // double kappa_max = 1;
   // if (kappa_min<0.1){kappa_min = 0.1;}
   // vec kappa = {kappa_min, kappa_max};
   
   // double min_mat_se_boot = mat_se_boot(find(mat_se_boot > 0)).min();
   // mat_se_boot.replace(0, min_mat_se_boot); 
   // mat se_boot = reshape(mat_se_boot,n,n);
   
   double censoring = 1-sum(delta)/n;
   vec kappa = {sqrt(censoring), 1};
   kappa = quantile(mat_se_boot, kappa);
   mat_se_boot.clamp(kappa(0),kappa(1));
   // if (kappa(0)<0.1){kappa(0) = 0.1;}
   mat se_boot = reshape(mat_se_boot,n,n);
   
   List app_std_npath(npath); vec absmax_app_npath(npath); vec absmax_app_std_npath(npath);
   for(int it=0; it<npath; it++){
     tempmat_nn = as<mat>(app_npath(it));
     absmax_app_npath(it) = abs(tempmat_nn).max();
     
     tempmat_nn /= se_boot;
     app_std_npath(it) = tempmat_nn;
     absmax_app_std_npath(it) = abs(tempmat_nn).max();
   }
   
   mat obs_std_npath = obs_npath/se_boot;
   double absmax_obs_npath = (abs(obs_npath)).max();
   double absmax_obs_std_npath = (abs(obs_std_npath)).max();
   
   uvec ind_unstd = (find(absmax_app_npath>absmax_obs_npath));
   double p_value = (ind_unstd.size()); p_value = p_value/npath;
   
   uvec ind_std = (find(absmax_app_std_npath>absmax_obs_std_npath));
   double p_std_value = (ind_std.size()); p_std_value = p_std_value/npath;
   
   if (npathsave<1){
     return List::create(_["p_std_value"]=p_std_value,_["p_value"]=p_value);
   } else if (npathsave > npath) {
     return List::create(_["SE_boot"]=se_boot,
                         _["obs_npath"]=obs_npath,_["obs_std_npath"]=obs_std_npath,
                         _["apprx_npath"]=app_npath,_["apprx_std_npath"]=app_std_npath,
                         _["p_value"]=p_value,_["p_std_value"]=p_std_value);
   } else {
     npathsave = npathsave - 1;
     app_npath = app_npath[Range(0,npathsave)];
     app_std_npath = app_std_npath[Range(0,npathsave)];
     return List::create(_["SE_boot"]=se_boot,
                         _["obs_npath"]=obs_npath,_["obs_std_npath"]=obs_std_npath,
                         _["apprx_npath"]=app_npath,_["apprx_std_npath"]=app_std_npath,
                         _["p_value"]=p_value,_["p_std_value"]=p_std_value);
   }
 }
 
 List link_cpp(int npath, vec b, vec time, vec delta, mat covariates, 
                  int npathsave, std::string eqType, 
                  bool fast_boot, mat invOmega){
   
   int n = covariates.n_rows;
   int p = covariates.n_cols;
   
   double sqrtn = sqrt(n);
   
   vec zero_vec_1 = zeros(1);
   vec zero_vec_p = zeros(p);
   vec zero_vec_n = zeros(n);
   mat zero_mat_np = zeros(n,p);
   mat zero_mat_nn = zeros(n,n);
   
   vec one_vec_n = ones(n);
   
   vec tempvec_p(p);
   vec tempvec_n(n);
   mat tempmat_np(n,p);
   mat tempmat_nn(n,n);
   
   vec resid = log(time) + covariates*b;
   uvec index_resid = sort_index(resid);
   
   time = time(index_resid);
   delta = delta(index_resid);
   covariates = covariates.rows(index_resid);
   resid = resid(index_resid);
   
   // -----------------------------------------------------------
   List pi_i_z(n); List N_i_t(n); List Y_i_t(n);
   vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
   mat sorted_Covari = sort(covariates);
   tempvec_n = zero_vec_n;
   for(int it=0; it<n; it++){
     for(int itt=0; itt<n; itt++){
       tempvec_n(itt) = (prod(covariates.row(it)<=sorted_Covari.row(itt))*1);
     }
     pi_i_z(it) = tempvec_n;
     if (delta(it)==1){
       N_i_t(it) = (resid>=resid(it));
     } else {
       N_i_t(it) = zero_vec_n;
     }
     Y_i_t(it) = (resid<=resid(it))*1;
     S_0_t += as<vec>(Y_i_t(it));
     S_1_t += as<vec>(Y_i_t(it))*(covariates.row(it));
     S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
   }
   
   vec dLambdahat_0_t = delta/S_0_t;
   dLambdahat_0_t.replace(datum::nan,0);
   
   mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
   E_pi_t_z.replace(datum::nan,0);
   
   // obs_npath; 1 by x vector
   List Mhat_i_t(n); List dMhat_i_t(n); vec obs_npath = zero_vec_n;
   for(int it=0; it<n; it++){
     tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
     Mhat_i_t(it) = tempvec_n;
     dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
     obs_npath += (tempvec_n(n-1))*as<vec>(pi_i_z(it));
   }
   obs_npath /= sqrtn;
   // vec obs_npath_0 = obs_npath(0) * one_vec_n;
   // obs_npath -= obs_npath_0;
   
   // -----------------------------------------------------------
   // ----------------------Kernel Smoothing---------------------
   // -----------------------------------------------------------
   double bw_base = pow((n*3/4),-0.2);
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
   ghat_0_t /= n;
   
   List ghat_t_z(p);
   tempvec_n = ghat_0_t % time;
   for(int itt=0; itt<p; itt++){
     vec Covari_col = covariates.col(itt);
     tempmat_nn = zero_mat_nn;
     for(int it=0; it<n; it++){
       tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
     }
     ghat_t_z(itt) = tempmat_nn/n;
   }
   
   // -----------------------------f0----------------------------
   vec Shat_0_e = cumprod(one_vec_n - dLambdahat_0_t);
   vec Fhat_0_e = one_vec_n - Shat_0_e;
   vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
   vec fhat_0_t = zero_vec_n;
   for(int it=0; it<n; it++){
     for(int itt=0; itt<n; itt++){
       fhat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn) * dFhat_0_e(itt);
     }
   }
   
   // vec Condi_Ehat = zero_vec_n;
   // for(int it=0; it<n; it++){
   //   Condi_Ehat(it) = sum(resid.tail(n-it-1)%dFhat_0_e.tail(n-it-1))/Shat_0_e(it);
   // }
   // Condi_Ehat.replace(datum::nan,0);
   // 
   // vec rhat_i = delta % resid + (one_vec_n - delta) % Condi_Ehat;
   // vec given_data_f = exp(rhat_i);
   // double bw_fn = bw_base * stddev(given_data_f);
   // vec fhat_0_t = zero_vec_n;
   // for(int it=0; it<n; it++){
   //   for(int itt=0; itt<n; itt++){
   //     fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn);
   //   }
   // }
   // fhat_0_t /= n;
   
   List fhat_inf_z(p);
   double tempvec_1 = fhat_0_t(n-1) * time(n-1);
   for(int itt=0; itt<p; itt++){
     vec Covari_col = covariates.col(itt);
     tempvec_n = zero_vec_n;
     for(int it=0; it<n; it++){
       if (delta(it)==1){
         tempvec_n += tempvec_1*(as<vec>(pi_i_z(it))*Covari_col(it));
       }
     }
     fhat_inf_z(itt) = tempvec_n/n;
   }
   
   List app_npath(npath);
   if (fast_boot) {
     // -----------------------------------------------------------
     // --------------------------dkappa_t-------------------------
     // -----------------------------------------------------------
     vec lambda_hat_0_t = fhat_0_t / Shat_0_e;
     lambda_hat_0_t.replace(datum::nan, 0);
     lambda_hat_0_t.replace(datum::inf, 0);
     
     vec lambda_hat_0_tTIMESt = lambda_hat_0_t % pred_data;
     vec dlambda_hat_0_tTIMESt = diff(join_cols(zero_vec_1, lambda_hat_0_tTIMESt));
     mat dkappa_t = - S_1_t;
     dkappa_t.each_col() %= dlambda_hat_0_tTIMESt/S_0_t;
     
     mat term2 = zero_mat_np;
     for(int it=0; it<p; it++){
       term2.col(it) = (as<vec>(fhat_inf_z(it)) +
         sum((as<mat>(ghat_t_z(it)).each_col()) % dLambdahat_0_t).t() +
         sum((S_pi_t_z.each_col()) % (dkappa_t.col(it))).t()/n);
     }
     term2 = term2 * invOmega;
     
     // -----------------------------------------------------------
     // ------------------------Sample Path------------------------
     // -----------------------------------------------------------
     for(int itt=0; itt<npath; itt++){
       vec phi_i = randg(n) - one_vec_n; // randn(n);
       
       tempvec_n = zero_vec_n; 
       tempmat_np = zero_mat_np;
       for(int i=0; i<n; i++){
         vec dM_phi = as<vec>(dMhat_i_t(i)) * phi_i(i);
         tempvec_n += dM_phi;
         tempmat_np += dM_phi * covariates.row(i);
       }
       mat term3 = term2 * ((tempmat_np.t() * S_0_t) - (S_1_t.t() * tempvec_n))/sqrtn;
       
       tempmat_nn = zero_mat_nn;
       for(int it=0; it<n; it++){
         tempmat_nn += ((as<vec>(dMhat_i_t(it)))*phi_i(it))%(as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col();
       }
       mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
       
       app_npath(itt) = (U_pi_phi_inf_z - term3)/sqrtn;
     }
   } else {
     // -----------------------------------------------------------
     // --------------------------U_inf_ls-------------------------
     // -----------------------------------------------------------
     mat U_inf_ls;
     if (eqType == "ls") {
       mat E_t = S_1_t.each_col()/S_0_t;
       vec resid_star = impute_residuals_weighted(resid, delta, one_vec_n);
       rowvec X_bar = mean(covariates, 0);
       mat Xc = covariates;
       Xc.each_row() -= X_bar;
       
       U_inf_ls = Xc;
       U_inf_ls.each_col() %= resid_star/n;
       
       vec lambda_hat_0_t = fhat_0_t / Shat_0_e;
       lambda_hat_0_t.replace(datum::nan, 0);
       lambda_hat_0_t.replace(datum::inf, 0);
       
       vec lambda_hat_0_tTIMESt = lambda_hat_0_t % pred_data;
       vec dlambda_hat_0_tTIMESt = diff(join_cols(zero_vec_1, lambda_hat_0_tTIMESt));
       mat dkappa_t = - E_t;
       dkappa_t.each_col() %= dlambda_hat_0_tTIMESt;
     }
     
     // -----------------------------------------------------------
     // ------------------------Sample Path------------------------
     // -----------------------------------------------------------
     for(int itt=0; itt<npath; itt++){
       
       vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
       while(tolerance>tol){
         phi_i = randg(n) - one_vec_n; // randn(n);
         
         vec U_phi_inf(n); 
         if (eqType == "ls") {
           U_phi_inf = U_inf_ls.t() * phi_i;
         } else {
           tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
           for(int i=0; i<n; i++){
             vec dM_phi = as<vec>(dMhat_i_t(i)) * phi_i(i);
             tempvec_n += dM_phi;
             tempmat_np += dM_phi * covariates.row(i);
           }
           U_phi_inf = (tempmat_np.t() * S_0_t) - (S_1_t.t() * tempvec_n);
         }
         
         List b_s_result = dfsane(b, time, delta, covariates, U_phi_inf, phi_i, n, p, sqrtn, eqType);
         tolerance = as<double>(b_s_result[0]);
         b_s = as<vec>(b_s_result[1]);
       }
       
       tempmat_nn = zero_mat_nn;
       for(int it=0; it<n; it++){
         tempmat_nn += ((as<vec>(dMhat_i_t(it)))*phi_i(it))%(as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col();
       }
       mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
       
       vec resid_s = log(time) + covariates*b_s;
       uvec index_resid_s = sort_index(resid_s);
       
       vec Delta_s = delta(index_resid_s);
       resid_s = resid_s(index_resid_s);
       
       NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
       for(int it=0; it<n; it++){
         Y_i_t_s = (resid_s<=resid_s(it))*1;
         S_0_t_s += as<vec>(Y_i_t_s);
       }
       
       vec dLambdahat_0_t_s = Delta_s/S_0_t_s;
       dLambdahat_0_t_s.replace(datum::nan,0);
       
       vec term1 = U_pi_phi_inf_z/sqrtn;
       vec term2 = zero_vec_n;
       tempvec_p = (b-b_s)*sqrtn;
       for(int it=0; it<p; it++){
         term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
       }
       vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
       
       app_npath(itt) = term1 - term2 - term3;
     }
   }
   
   NumericMatrix tempmat_nnpath(n,npath);
   for(int it=0; it<npath; it++){
     tempmat_nnpath(_,it) = (as<NumericVector>(app_npath(it)));
   }
   vec se_boot = stddev(as<mat>(tempmat_nnpath),0,1);
   // too low values which are 0 or computationally 0 of se_boot makes a problem,
   // so we adjust them to have kappa = quantile of mat_se_boot
   // e.g., kappa_min = censoring; quantile(mat_se_boot) = {censoring, 1};
   // double censoring = 1-sum(delta)/n;
   // double kappa_min = censoring;
   // double kappa_max = 1;
   // if (kappa_min<0.1){kappa_min = 0.1;}
   // vec kappa = {kappa_min, kappa_max};
   
   // double min_se_boot = se_boot(find(se_boot > 0)).min();
   // se_boot.replace(0, min_se_boot); 
   
   double censoring = 1-sum(delta)/n;
   vec kappa = {censoring, 1};
   kappa = quantile(se_boot, kappa);
   // if (kappa(0)<0.1){kappa(0) = 0.1;}
   se_boot.clamp(kappa(0),kappa(1));
   
   List app_std_npath(npath); vec absmax_app_npath(npath); vec absmax_app_std_npath(npath);
   for(int it=0; it<npath; it++){
     tempvec_n = as<vec>(app_npath(it));
     absmax_app_npath(it) = abs(tempvec_n).max();
     
     tempvec_n /= se_boot;
     app_std_npath(it) = tempvec_n;
     absmax_app_std_npath(it) = abs(tempvec_n).max();
   }
   
   vec obs_std_npath = obs_npath/se_boot;
   double absmax_obs_npath = (abs(obs_npath)).max();
   double absmax_obs_std_npath = (abs(obs_std_npath)).max();
   
   uvec ind_unstd = (find(absmax_app_npath>absmax_obs_npath));
   double p_value = (ind_unstd.size()); p_value = p_value/npath;
   
   uvec ind_std = (find(absmax_app_std_npath>absmax_obs_std_npath));
   double p_std_value = (ind_std.size()); p_std_value = p_std_value/npath;
   
   if (npathsave<1){
     return List::create(_["p_std_value"]=p_std_value,_["p_value"]=p_value);
   } else if (npathsave > npath) {
     return List::create(_["SE_boot"]=se_boot,
                         _["obs_npath"]=obs_npath,_["obs_std_npath"]=obs_std_npath,
                         _["apprx_npath"]=app_npath,_["apprx_std_npath"]=app_std_npath,
                         _["p_value"]=p_value,_["p_std_value"]=p_std_value);
   } else {
     npathsave = npathsave - 1;
     app_npath = app_npath[Range(0,npathsave)];
     app_std_npath = app_std_npath[Range(0,npathsave)];
     return List::create(_["SE_boot"]=se_boot,
                         _["obs_npath"]=obs_npath,_["obs_std_npath"]=obs_std_npath,
                         _["apprx_npath"]=app_npath,_["apprx_std_npath"]=app_std_npath,
                         _["p_value"]=p_value,_["p_std_value"]=p_std_value);
   }
 }
 
 List form_cpp(int npath, vec b, vec time, vec delta, mat covariates, 
                  int covTested, int npathsave, std::string eqType, 
                  bool fast_boot, mat invOmega){
   
   int n = covariates.n_rows;
   int p = covariates.n_cols;
   
   double sqrtn = sqrt(n);
   
   vec zero_vec_1 = zeros(1);
   vec zero_vec_p = zeros(p);
   vec zero_vec_n = zeros(n);
   mat zero_mat_np = zeros(n,p);
   mat zero_mat_nn = zeros(n,n);
   
   vec one_vec_n = ones(n);
   
   vec tempvec_p(p);
   vec tempvec_n(n);
   mat tempmat_np(n,p);
   mat tempmat_nn(n,n);
   
   vec resid = log(time) + covariates*b;
   uvec index_resid = sort_index(resid);
   
   time = time(index_resid);
   delta = delta(index_resid);
   covariates = covariates.rows(index_resid);
   resid = resid(index_resid);
   
   // -----------------------------------------------------------
   List pi_i_z(n); List N_i_t(n); List Y_i_t(n);
   vec S_0_t = zero_vec_n; mat S_1_t = zero_mat_np; mat S_pi_t_z = zero_mat_nn;
   vec form_Covari = covariates.col(covTested-1);
   vec sorted_form_Covari = sort(form_Covari);
   tempvec_n = zero_vec_n;
   for(int it=0; it<n; it++){
     for(int itt=0; itt<n; itt++){
       tempvec_n(itt) = (form_Covari(it)<=sorted_form_Covari(itt))*1;
     }
     pi_i_z(it) = tempvec_n;
     if (delta(it)==1){
       N_i_t(it) = (resid>=resid(it));
     } else {
       N_i_t(it) = zero_vec_n;
     }
     Y_i_t(it) = (resid<=resid(it))*1;
     S_0_t += as<vec>(Y_i_t(it));
     S_1_t += as<vec>(Y_i_t(it))*(covariates.row(it));
     S_pi_t_z += (as<vec>(Y_i_t(it)))*(as<rowvec>(pi_i_z(it)));
   }
   
   vec dLambdahat_0_t = delta/S_0_t;
   dLambdahat_0_t.replace(datum::nan,0);
   
   mat E_pi_t_z = S_pi_t_z.each_col()/S_0_t;
   E_pi_t_z.replace(datum::nan,0);
   
   // obs_npath; 1 by x vector
   List Mhat_i_t(n); List dMhat_i_t(n); vec obs_npath = zero_vec_n;
   for(int it=0; it<n; it++){
     tempvec_n = as<vec>(N_i_t(it))-(cumsum(as<vec>(Y_i_t(it))%(dLambdahat_0_t)));
     Mhat_i_t(it) = tempvec_n;
     dMhat_i_t(it) = diff(join_cols(zero_vec_1,tempvec_n));
     obs_npath += (tempvec_n(n-1))*as<vec>(pi_i_z(it));
   }
   obs_npath /= sqrtn;
   // vec obs_npath_0 = obs_npath(0) * one_vec_n;
   // obs_npath -= obs_npath_0;
   
   // -----------------------------------------------------------
   // ----------------------Kernel Smoothing---------------------
   // -----------------------------------------------------------
   double bw_base = pow((n*3/4),-0.2);
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
   ghat_0_t /= n;
   
   List ghat_t_z(p);
   tempvec_n = ghat_0_t % time;
   for(int itt=0; itt<p; itt++){
     vec Covari_col = covariates.col(itt);
     tempmat_nn = zero_mat_nn;
     for(int it=0; it<n; it++){
       tempmat_nn += tempvec_n*((as<rowvec>(pi_i_z(it)))*Covari_col(it));
     }
     ghat_t_z(itt) = tempmat_nn/n;
   }
   
   // -----------------------------f0----------------------------
   vec Shat_0_e = cumprod(one_vec_n - dLambdahat_0_t);
   vec Fhat_0_e = one_vec_n - Shat_0_e;
   vec dFhat_0_e = diff(join_cols(zero_vec_1,Fhat_0_e));
   vec fhat_0_t = zero_vec_n;
   for(int it=0; it<n; it++){
     for(int itt=0; itt<n; itt++){
       fhat_0_t(it) += normpdf(pred_data(it),given_data_g(itt),bw_gn) * dFhat_0_e(itt);
     }
   }
   // vec Condi_Ehat = zero_vec_n;
   // for(int it=0; it<n; it++){
   //   Condi_Ehat(it) = sum(resid.tail(n-it-1)%dFhat_0_e.tail(n-it-1))/Shat_0_e(it);
   // }
   // Condi_Ehat.replace(datum::nan,0);
   // 
   // vec rhat_i = delta % resid + (one_vec_n - delta) % Condi_Ehat;
   // vec given_data_f = exp(rhat_i);
   // double bw_fn = bw_base * stddev(given_data_f);
   // vec fhat_0_t = zero_vec_n;
   // for(int it=0; it<n; it++){
   //   for(int itt=0; itt<n; itt++){
   //     fhat_0_t(it) += normpdf(pred_data(it),given_data_f(itt),bw_fn);
   //   }
   // }
   // fhat_0_t /= n;
   
   List fhat_inf_z(p);
   double tempvec_1 = fhat_0_t(n-1) * time(n-1);
   for(int itt=0; itt<p; itt++){
     vec Covari_col = covariates.col(itt);
     tempvec_n = zero_vec_n;
     for(int it=0; it<n; it++){
       if (delta(it)==1){
         tempvec_n += tempvec_1*(as<vec>(pi_i_z(it))*Covari_col(it));
       }
     }
     fhat_inf_z(itt) = tempvec_n/n;
   }
   
   List app_npath(npath);
   if (fast_boot) {
     // -----------------------------------------------------------
     // --------------------------dkappa_t-------------------------
     // -----------------------------------------------------------
     vec lambda_hat_0_t = fhat_0_t / Shat_0_e;
     lambda_hat_0_t.replace(datum::nan, 0);
     lambda_hat_0_t.replace(datum::inf, 0);
     
     vec lambda_hat_0_tTIMESt = lambda_hat_0_t % pred_data;
     vec dlambda_hat_0_tTIMESt = diff(join_cols(zero_vec_1, lambda_hat_0_tTIMESt));
     mat dkappa_t = - S_1_t;
     dkappa_t.each_col() %= dlambda_hat_0_tTIMESt/S_0_t;
     
     mat term2 = zero_mat_np;
     for(int it=0; it<p; it++){
       term2.col(it) = (as<vec>(fhat_inf_z(it)) +
         sum((as<mat>(ghat_t_z(it)).each_col()) % dLambdahat_0_t).t() +
         sum((S_pi_t_z.each_col()) % (dkappa_t.col(it))).t()/n);
     }
     term2 = term2 * invOmega;
     
     // -----------------------------------------------------------
     // ------------------------Sample Path------------------------
     // -----------------------------------------------------------
     for(int itt=0; itt<npath; itt++){
       vec phi_i = randg(n) - one_vec_n; // randn(n);
       
       tempvec_n = zero_vec_n; 
       tempmat_np = zero_mat_np;
       for(int i=0; i<n; i++){
         vec dM_phi = as<vec>(dMhat_i_t(i)) * phi_i(i);
         tempvec_n += dM_phi;
         tempmat_np += dM_phi * covariates.row(i);
       }
       mat term3 = term2 * ((tempmat_np.t() * S_0_t) - (S_1_t.t() * tempvec_n))/sqrtn;
       
       tempmat_nn = zero_mat_nn;
       for(int it=0; it<n; it++){
         tempmat_nn += ((as<vec>(dMhat_i_t(it)))*phi_i(it))%(as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col();
       }
       mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
       
       app_npath(itt) = (U_pi_phi_inf_z - term3)/sqrtn;
     }
   } else {
     // -----------------------------------------------------------
     // --------------------------U_inf_ls-------------------------
     // -----------------------------------------------------------
     mat U_inf_ls;
     if (eqType == "ls") {
       mat E_t = S_1_t.each_col()/S_0_t;
       vec resid_star = impute_residuals_weighted(resid, delta, one_vec_n);
       rowvec X_bar = mean(covariates, 0);
       mat Xc = covariates;
       Xc.each_row() -= X_bar;
       
       U_inf_ls = Xc;
       U_inf_ls.each_col() %= resid_star/n;
       
       vec lambda_hat_0_t = fhat_0_t / Shat_0_e;
       lambda_hat_0_t.replace(datum::nan, 0);
       lambda_hat_0_t.replace(datum::inf, 0);
       
       vec lambda_hat_0_tTIMESt = lambda_hat_0_t % pred_data;
       vec dlambda_hat_0_tTIMESt = diff(join_cols(zero_vec_1, lambda_hat_0_tTIMESt));
       mat dkappa_t = - E_t;
       dkappa_t.each_col() %= dlambda_hat_0_tTIMESt;
     }
     
     // -----------------------------------------------------------
     // ------------------------Sample Path------------------------
     // -----------------------------------------------------------
     for(int itt=0; itt<npath; itt++){
       
       vec phi_i(n); vec b_s(p); double tol = pow(p,2); double tolerance = tol+1;
       while(tolerance>tol){
         phi_i = randg(n) - one_vec_n; // randn(n);
         
         vec U_phi_inf(n); 
         if (eqType == "ls") {
           U_phi_inf = U_inf_ls.t() * phi_i;
         } else {
           tempvec_n = zero_vec_n; tempmat_np = zero_mat_np;
           for(int i=0; i<n; i++){
             vec dM_phi = as<vec>(dMhat_i_t(i)) * phi_i(i);
             tempvec_n += dM_phi;
             tempmat_np += dM_phi * covariates.row(i);
           }
           U_phi_inf = (tempmat_np.t() * S_0_t) - (S_1_t.t() * tempvec_n);
         }
         
         List b_s_result = dfsane(b, time, delta, covariates, U_phi_inf, phi_i, n, p, sqrtn, eqType);
         tolerance = as<double>(b_s_result[0]);
         b_s = as<vec>(b_s_result[1]);
       }
       
       tempmat_nn = zero_mat_nn;
       for(int it=0; it<n; it++){
         tempmat_nn += ((as<vec>(dMhat_i_t(it)))*phi_i(it))%(as<rowvec>(pi_i_z(it))-E_pi_t_z.each_row()).each_col();
       }
       mat U_pi_phi_inf_z = (sum(tempmat_nn)).t();
       
       vec resid_s = log(time) + covariates*b_s;
       uvec index_resid_s = sort_index(resid_s);
       
       vec Delta_s = delta(index_resid_s);
       resid_s = resid_s(index_resid_s);
       
       NumericVector Y_i_t_s(n); vec S_0_t_s = zero_vec_n;
       for(int it=0; it<n; it++){
         Y_i_t_s = (resid_s<=resid_s(it))*1;
         S_0_t_s += as<vec>(Y_i_t_s);
       }
       
       vec dLambdahat_0_t_s = Delta_s/S_0_t_s;
       dLambdahat_0_t_s.replace(datum::nan,0);
       
       vec term1 = U_pi_phi_inf_z/sqrtn;
       vec term2 = zero_vec_n;
       tempvec_p = (b-b_s)*sqrtn;
       for(int it=0; it<p; it++){
         term2 += (as<vec>(fhat_inf_z(it))+(sum((as<mat>(ghat_t_z(it)).each_col())%dLambdahat_0_t)).t())*(tempvec_p(it));
       }
       vec term3 = (sum((S_pi_t_z.each_col())%(dLambdahat_0_t - dLambdahat_0_t_s))).t()/sqrtn;
       
       app_npath(itt) = term1 - term2 - term3;
     }
   }
   
   NumericMatrix tempmat_nnpath(n,npath);
   for(int it=0; it<npath; it++){
     tempmat_nnpath(_,it) = (as<NumericVector>(app_npath(it)));
   }
   vec se_boot = stddev(as<mat>(tempmat_nnpath),0,1);
   // too low values which are 0 or computationally 0 of se_boot makes a problem,
   // so we adjust them to have kappa = quantile of mat_se_boot
   // e.g., kappa_min = censoring; quantile(mat_se_boot) = {censoring, 1};
   // double censoring = 1-sum(delta)/n;
   // double kappa_min = censoring;
   // double kappa_max = 1;
   // if (kappa_min<0.1){kappa_min = 0.1;}
   // vec kappa = {kappa_min, kappa_max};
   
   // double min_se_boot = se_boot(find(se_boot > 0)).min();
   // se_boot.replace(0, min_se_boot); 
   
   double censoring = 1-sum(delta)/n;
   vec kappa = {censoring, 1};
   kappa = quantile(se_boot, kappa);
   // if (kappa(0)<0.1){kappa(0) = 0.1;}
   se_boot.clamp(kappa(0),kappa(1));
   
   List app_std_npath(npath); vec absmax_app_npath(npath); vec absmax_app_std_npath(npath);
   for(int it=0; it<npath; it++){
     tempvec_n = as<vec>(app_npath(it));
     absmax_app_npath(it) = abs(tempvec_n).max();
     
     tempvec_n /= se_boot;
     app_std_npath(it) = tempvec_n;
     absmax_app_std_npath(it) = abs(tempvec_n).max();
   }
   
   vec obs_std_npath = obs_npath/se_boot;
   double absmax_obs_npath = (abs(obs_npath)).max();
   double absmax_obs_std_npath = (abs(obs_std_npath)).max();
   
   uvec ind_unstd = (find(absmax_app_npath>absmax_obs_npath));
   double p_value = (ind_unstd.size()); p_value = p_value/npath;
   
   uvec ind_std = (find(absmax_app_std_npath>absmax_obs_std_npath));
   double p_std_value = (ind_std.size()); p_std_value = p_std_value/npath;
   
   if (npathsave<1){
     return List::create(_["p_std_value"]=p_std_value,_["p_value"]=p_value);
   } else if (npathsave > npath) {
     return List::create(_["SE_boot"]=se_boot,
                         _["obs_npath"]=obs_npath,_["obs_std_npath"]=obs_std_npath,
                         _["apprx_npath"]=app_npath,_["apprx_std_npath"]=app_std_npath,
                         _["p_value"]=p_value,_["p_std_value"]=p_std_value);
   } else {
     npathsave = npathsave - 1;
     app_npath = app_npath[Range(0,npathsave)];
     app_std_npath = app_std_npath[Range(0,npathsave)];
     return List::create(_["SE_boot"]=se_boot,
                         _["obs_npath"]=obs_npath,_["obs_std_npath"]=obs_std_npath,
                         _["apprx_npath"]=app_npath,_["apprx_std_npath"]=app_std_npath,
                         _["p_value"]=p_value,_["p_std_value"]=p_std_value);
   }
 }
 
 // inv_cpp
 arma::mat inv_cpp(arma::mat const& MAT);
 RcppExport SEXP _afttest_inv_cpp(SEXP MATSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< arma::mat >::type MAT(MATSEXP);
   rcpp_result_gen = Rcpp::wrap(inv_cpp(MAT));
   return rcpp_result_gen;
   END_RCPP
 }
 // abar_gehan_gehan_cpp
 arma::mat abar_gehan_cpp(arma::vec beta, arma::vec Y, arma::mat X, arma::vec delta, arma::mat sigma, arma::vec weights, arma::vec gehanWeights);
 RcppExport SEXP _afttest_abar_gehan_cpp(SEXP betaSEXP, SEXP YSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP sigmaSEXP, SEXP weightsSEXP, SEXP gehanWeightsSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
   Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
   Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type gehanWeights(gehanWeightsSEXP);
   rcpp_result_gen = Rcpp::wrap(abar_gehan_cpp(beta, Y, X, delta, sigma, weights, gehanWeights));
   return rcpp_result_gen;
   END_RCPP
 }
 // omega_cpp
 arma::mat omega_gehan_cpp(arma::vec beta, arma::vec Y, arma::mat X, arma::vec delta, arma::vec weights);
 RcppExport SEXP _afttest_omega_gehan_cpp(SEXP betaSEXP, SEXP YSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP weightsSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
   Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
   rcpp_result_gen = Rcpp::wrap(omega_gehan_cpp(beta, Y, X, delta, weights));
   return rcpp_result_gen;
   END_RCPP
 }
 // score_gehan_is_cpp
 arma::vec score_gehan_is_cpp(arma::vec beta, arma::vec Y, arma::mat X, arma::vec delta, arma::mat sigma, arma::vec weights);
 RcppExport SEXP _afttest_score_gehan_is_cpp(SEXP betaSEXP, SEXP YSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP sigmaSEXP, SEXP weightsSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
   Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
   Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
   rcpp_result_gen = Rcpp::wrap(score_gehan_is_cpp(beta, Y, X, delta, sigma, weights));
   return rcpp_result_gen;
   END_RCPP
 }
 // score_gehan_ns_cpp
 arma::vec score_gehan_ns_cpp(arma::vec beta, arma::vec Y, arma::mat X, arma::vec delta, arma::vec weights);
 RcppExport SEXP _afttest_score_gehan_ns_cpp(SEXP betaSEXP, SEXP YSEXP, SEXP XSEXP, SEXP deltaSEXP, SEXP weightsSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
   Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type delta(deltaSEXP);
   Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
   rcpp_result_gen = Rcpp::wrap(score_gehan_ns_cpp(beta, Y, X, delta, weights));
   return rcpp_result_gen;
   END_RCPP
 }
 // omni_cpp
 List omni_cpp(int npath, vec b, vec time, vec delta, mat covariates, int npathsave, std::string eqType, bool fast_boot, mat invOmega);
 RcppExport SEXP _afttest_omni_cpp(SEXP npathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP npathsaveSEXP, SEXP eqTypeSEXP, SEXP fast_bootSEXP, SEXP invOmegaSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< int >::type npath(npathSEXP);
   Rcpp::traits::input_parameter< vec >::type b(bSEXP);
   Rcpp::traits::input_parameter< vec >::type time(TimeSEXP);
   Rcpp::traits::input_parameter< vec >::type delta(DeltaSEXP);
   Rcpp::traits::input_parameter< mat >::type covariates(CovariSEXP);
   Rcpp::traits::input_parameter< int >::type npathsave(npathsaveSEXP);
   Rcpp::traits::input_parameter< std::string >::type eqType(eqTypeSEXP);
   Rcpp::traits::input_parameter< bool >::type fast_boot(fast_bootSEXP);
   Rcpp::traits::input_parameter< mat >::type invOmega(invOmegaSEXP);
   rcpp_result_gen = Rcpp::wrap(omni_cpp(npath, b, time, delta, covariates, npathsave, eqType, fast_boot, invOmega));
   return rcpp_result_gen;
   END_RCPP
 }
 // link_cpp
 List link_cpp(int npath, vec b, vec time, vec delta, mat covariates, int npathsave, std::string eqType, bool fast_boot, mat invOmega);
 RcppExport SEXP _afttest_link_cpp(SEXP npathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP npathsaveSEXP, SEXP eqTypeSEXP, SEXP fast_bootSEXP, SEXP invOmegaSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< int >::type npath(npathSEXP);
   Rcpp::traits::input_parameter< vec >::type b(bSEXP);
   Rcpp::traits::input_parameter< vec >::type time(TimeSEXP);
   Rcpp::traits::input_parameter< vec >::type delta(DeltaSEXP);
   Rcpp::traits::input_parameter< mat >::type covariates(CovariSEXP);
   Rcpp::traits::input_parameter< int >::type npathsave(npathsaveSEXP);
   Rcpp::traits::input_parameter< std::string >::type eqType(eqTypeSEXP);
   Rcpp::traits::input_parameter< bool >::type fast_boot(fast_bootSEXP);
   Rcpp::traits::input_parameter< mat >::type invOmega(invOmegaSEXP);
   rcpp_result_gen = Rcpp::wrap(link_cpp(npath, b, time, delta, covariates, npathsave, eqType, fast_boot, invOmega));
   return rcpp_result_gen;
   END_RCPP
 }
 // form_cpp
 List form_cpp(int npath, vec b, vec time, vec delta, mat covariates, int covTested, int npathsave, std::string eqType, bool fast_boot, mat invOmega);
 RcppExport SEXP _afttest_form_cpp(SEXP npathSEXP, SEXP bSEXP, SEXP TimeSEXP, SEXP DeltaSEXP, SEXP CovariSEXP, SEXP cov_testedSEXP, SEXP npathsaveSEXP, SEXP eqTypeSEXP, SEXP fast_bootSEXP, SEXP invOmegaSEXP) {
   BEGIN_RCPP
   Rcpp::RObject rcpp_result_gen;
   Rcpp::RNGScope rcpp_rngScope_gen;
   Rcpp::traits::input_parameter< int >::type npath(npathSEXP);
   Rcpp::traits::input_parameter< vec >::type b(bSEXP);
   Rcpp::traits::input_parameter< vec >::type time(TimeSEXP);
   Rcpp::traits::input_parameter< vec >::type delta(DeltaSEXP);
   Rcpp::traits::input_parameter< mat >::type covariates(CovariSEXP);
   Rcpp::traits::input_parameter< int >::type covTested(cov_testedSEXP);
   Rcpp::traits::input_parameter< int >::type npathsave(npathsaveSEXP);
   Rcpp::traits::input_parameter< std::string >::type eqType(eqTypeSEXP);
   Rcpp::traits::input_parameter< bool >::type fast_boot(fast_bootSEXP);
   Rcpp::traits::input_parameter< mat >::type invOmega(invOmegaSEXP);
   rcpp_result_gen = Rcpp::wrap(form_cpp(npath, b, time, delta, covariates, covTested, npathsave, eqType, fast_boot, invOmega));
   return rcpp_result_gen;
   END_RCPP
 }
 
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 // -----------------------------------------------------------
 static const R_CallMethodDef CallEntries[] = {
   // New helper functions
   {"_afttest_inv_cpp", (DL_FUNC) &_afttest_inv_cpp, 1},
   {"_afttest_abar_gehan_cpp", (DL_FUNC) &_afttest_abar_gehan_cpp, 7},
   {"_afttest_score_gehan_is_cpp", (DL_FUNC) &_afttest_score_gehan_is_cpp, 6},
   {"_afttest_score_gehan_ns_cpp", (DL_FUNC) &_afttest_score_gehan_ns_cpp, 5},
   {"_afttest_omega_gehan_cpp", (DL_FUNC) &_afttest_omega_gehan_cpp, 5},
   
   {"_afttest_omni_cpp", (DL_FUNC) &_afttest_omni_cpp, 9},
   {"_afttest_link_cpp", (DL_FUNC) &_afttest_link_cpp, 9},
   {"_afttest_form_cpp", (DL_FUNC) &_afttest_form_cpp, 10},
   {NULL, NULL, 0}
 };
 
 RcppExport void R_init_afttest(DllInfo *dll) {
   R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
   R_useDynamicSymbols(dll, FALSE);
 }
 