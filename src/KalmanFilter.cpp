#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


//' @description Derives the normal density
//' @param x diagonal mat with the quantile values
//' @param var diagonal mat with the variances
//' @returns row vector with the densities

 vec dnormVec_fctn(vec x, mat Var) {
   int n = x.size();
   vec dens(n);
   for(int i = 0; i < n; i++) {
     dens[i] = R::dnorm4(x[i], 0, sqrt(Var(i,i)), FALSE);
   }
   return dens;
 }


//' @description Function that applies parameter constraints
//' @param spec character indicating the model specification
//' @param par vector structural parameters
//' @return constrained parameter vector
// [[Rcpp::export]]
vec ParConstrain(vec paramVec, bool dgp1) {
  vec constrPar = paramVec;

  if (dgp1 == true){
// Constrain AR parameters of the cycle to result in a stable process
double phi_1 = 2 * paramVec(0) / (1 + abs(paramVec(0)));
double phi_2 = -(1 - abs(phi_1)) * paramVec(1) / (1 + abs(paramVec(1))) - abs(phi_1);
constrPar(0) = phi_1;
constrPar(1) = phi_2;
// Constrain system St.Devs to be positive
constrPar(span(2, 5)) = exp(paramVec(span(2, 5)));
  } else {
    constrPar(0) = exp(paramVec(0));
  }
  return constrPar;
}

//' @description Function constructs the system matrices for the Kalman filter
//' @param paramVec vector structural parameters
//' @param constrainParam boolean. If true, the parameter constraining function is applied.
//' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
//' @return list with the system matrices
 // [[Rcpp::export]]
 List SystemMat_fctn(vec paramVec, bool constrainParam, bool dgp1) {
   vec constrPar;
   if (constrainParam == true){
     constrPar = ParConstrain(paramVec, dgp1); 
   } else {
     constrPar = paramVec;
   }
   int transObs;
   if (dgp1 == true){
     transObs = 4;
   } else {
     transObs = 2;
   }
   
   vec SdVec;
   double SdEta;
   mat A(transObs, transObs);
   mat B;
   mat C(1, transObs);
   mat D(1, 1, fill::zeros); 
   mat Q;
   
   if (dgp1 == true){
        double phi_1 = constrPar[0];
        double phi_2 = constrPar[1];
        SdEta = constrPar[2];
       double SdU = constrPar[3];
      double SdE = constrPar[4];
      double SdEpsilon = constrPar[5];
// Transition eq
A << 1 << 1 << 0 << 0 << endr
  << 0 << 1 << 0 << 0 << endr
  << 0 << 0 << phi_1 << phi_2 << endr
  << 0 << 0 << 1 << 0 << endr;
 B = join_cols(eye(3, 3), zeros<rowvec>(3));
 SdVec = {pow(SdEta, 2), pow(SdU, 2), pow(SdE, 2)};
Q = diagmat(SdVec);
// Measurement eq
 C << 1 << 0 << 1 << 0 << endr;
 D.col(0).row(0) = pow(SdEpsilon, 2);
   } else {
     SdEta = constrPar[0];
     // Transition eq
     A << 1 << 1 << endr
       << 0 << 1 << endr;
     B << 1 << endr
       << 0 << endr;
     SdVec = {pow(SdEta, 2)};
     Q = diagmat(SdVec);
     // Measurement eq
     C << 1 << 0 << endr;
   }
        List outputList = List::create(
                        Named("A") = A,
                        _["B"] = B,
                        _["C"] = C,
                        _["D"] = D,
                        _["Q"] = Q);                

        return outputList;
}


//' @description Function that executes the Kalman recursions. Notation as in Angelini et al. (2022).
//' @param paramVec vector of structural parameters
//' @param data matrix with the data. One column per observed variable
//' @param outLogLik boolean. If true, function outputs a likelihood value. Else the filter output
//' @param constrainParam boolean. If true, the parameter constraining function is applied.
//' @param dgp1 boolean. If True, the simulation and estimation concern DGP 1. Else DGP 2
//' @return constrained parameter vector
// [[Rcpp::export]]
List KalmanRecursions(vec paramVec, mat data, bool outLogLik, bool constrainParam, bool dgp1)
{
        // Produce system matrices
      List systemList = SystemMat_fctn(paramVec, constrainParam, dgp1);
        
        // Transition eq
        mat A = systemList["A"];
        mat transA = trans(A);
        mat B = systemList["B"];
        mat Q = systemList["Q"];
        mat Q_expand = B * Q * trans(B);
        // Measurement eq
        mat C = systemList["C"];
        mat transC = trans(C);
        mat D = systemList["D"];

        // Dimensions
        int transDim = A.n_cols;
        int obsDim = data.n_cols;
        // Number of time periods
        int nPeriods = data.n_rows;

        // Initialize the output
        vec logLik_vec(nPeriods, fill::zeros);
        mat Z_tt_mat(nPeriods, transDim, fill::zeros);
        cube P_tt_array(transDim, transDim, nPeriods, fill::zeros);
        cube K_t_array(transDim, obsDim, nPeriods, fill::zeros);
        mat e_hat_mat(nPeriods, obsDim, fill::zeros);
        cube Sigma_sqrt_t(obsDim, obsDim, nPeriods, fill::zeros);

        // Initialize the filter routine (diffusely)
        // State vector with zeros
        vec Z_t1(transDim, fill::zeros);
        // State vector var-cov matrix with very high variance
        vec P_t1_vec(transDim, fill::value(1000));
        mat P_t1 = diagmat(P_t1_vec);
        vec epsilon_t;
        mat Sigma_t;
        mat Sigma_inv_t;
        vec e_hat_t;
        mat K_t;
        vec Z_tt;
        mat P_tt;

        // Run the recursions (until nPeriods-1 bc state vector enters measurement eq with lag (eq. 2))
        for (int i = 0; i < nPeriods; i++)
        {

                //-------------------//
                // Get innovations
                //-------------------//

                // One step ahead prediction error
                epsilon_t = data.row(i).t() - C * Z_t1;
                // One step ahead prediction error
                Sigma_t = C * P_t1 * transC + D;
                Sigma_inv_t = inv(Sigma_t);
                // standardized prediction errors eq. (13)
                // Bootstrap algorithm step 1
                e_hat_t = real(sqrtmat(Sigma_inv_t)) * epsilon_t;
                Sigma_sqrt_t.slice(i) = real(sqrtmat(Sigma_t));
                
                //-------------------//
                // Updating step
                //-------------------//

                // Kalman gain
                K_t = P_t1 * transC * Sigma_inv_t;
                // State vector
                Z_tt = Z_t1 + K_t * epsilon_t;
                // State vector var-cov matrix
                P_tt = P_t1 - K_t * C * P_t1;
                // Either store the likelihood or the filter output
                if (outLogLik == true)
                {
                        // Calculation and storage of the log likelihood
                       logLik_vec[i] = -log(as_scalar(dnormVec_fctn(epsilon_t, Sigma_t)));
                }
                else
                {
                        K_t_array.slice(i) = K_t;
                        Z_tt_mat.row(i) = Z_tt.t();
                        P_tt_array.slice(i) = P_tt;
                        e_hat_mat.row(i) = e_hat_t.t();
                }

                //-------------------//
                // Prediction step
                //-------------------//

                // State vector
                Z_t1 = A * Z_tt;
                // State vector var-cov matrix
                P_t1 = A * P_tt * transA + Q_expand;
        }
        // Set the output
        if (outLogLik == TRUE)
        {
                double logLik = sum(logLik_vec);
                List logLikList = List::create(logLik);
                return logLikList;
        }
        else
        {
                // Note that first innovation is dropped (eq 13)
                vec e_hat_0(obsDim);
                e_hat_0.fill(datum::nan);
                e_hat_mat.row(0) = e_hat_0.t();
                // Construct the output list
                List outputList = List::create(
                    Named("Z_tt") = Z_tt_mat,
                    _["P_tt"] = P_tt_array,
                    _["K_t"] = K_t_array,
                    _["e_hat"] = e_hat_mat,
                    _["DGP1"] = dgp1,
                    _["C_theta"] = C,
                    _["Sigma_sqrt"] = Sigma_sqrt_t);
                return outputList;
         }
}