#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' @description Function that executes the Kalman recursions. Notation as in Angelini et al. (2022).
//' @param data matrix with the data. One column per observed variable
//' @param paramVec vector of structural parameters
//' @param systemList List with system matrices
//' @param outLogLik boolean. If true, function outputs a likelihood value. Else the filter output
//' @return constrained parameter vector
// [[Rcpp::export]]
List KalmanRecursions(mat data, vec paramVec, List systemList, bool outLogLik)
{
        // Transition eq
        mat A = systemList["A"];
        mat B = systemList["B"];
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
        mat logLik_mat(nPeriods, obsDim, fill::zeros);
        mat Z_tt_mat(nPeriods, transDim, fill::zeros);
        cube P_tt_array(transDim, transDim, nPeriods, fill::zeros);
        cube K_t_array(transDim, transDim, nPeriods, fill::zeros);
        mat e_hat_mat(nPeriods, transDim, fill::zeros);

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
        for (int i = 0; i < (nPeriods - 1); i++)
        {
                //-------------------//
                // Get innovations
                //-------------------//

                // One step ahead prediction error
                epsilon_t = data.row(i) - C * Z_t1;
                // One step ahead prediction error
                Sigma_t = C * P_t1 * transC + D;
                Sigma_inv_t = inv(Sigma_t);
                // standardized prediction errors eq. (13)
                // Bootstrap algorithm step 1
                e_hat_t = pow(Sigma_t, -.5) * epsilon_t;

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
                        logLik_mat.row(i) = .5 * log(2 * datum::pi) - .5 * log(abs(det(Sigma_t))) + (-.5 * epsilon_t * Sigma_inv_t * epsilon_t);
                }
                else
                {
                        K_t_array.slice(i) = K_t;
                        Z_tt_mat.row(i) = Z_tt;
                        P_tt_array.slice(i) = P_tt;
                        e_hat_mat.row(i) = e_hat_t;
                }

                //-------------------//
                // Prediction step
                //-------------------//

                // State vector
                Z_t1 = C * Z_tt;
                // State vector var-cov matrix
                P_t1 = C * P_tt * transC + B;
        }
        // Set the output
        if (outLogLik == TRUE)
        {
                mat logLik = -sum(logLik_mat);
                List logLikList = List::create(logLik);
                return logLikList;
        }
        else
        {
                // Note that first innovation is dropped (eq 13)
                vec e_hat_0(obsDim);
                e_hat_0.fill(datum::nan);
                e_hat_t.row(0) = e_hat_0.t();
                // Construct the output list
                List outputList = List::create(
                    Named("Z_tt") = Z_tt_mat,
                    _["P_tt"] = P_tt_array,
                    _["K_t"] = K_t_array,
                    _["e_hat"] = e_hat_mat);
                return outputList;
        }
}