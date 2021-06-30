#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Read in data
  // DATA_VECTOR(t); // ORDERED survival times (Don't need this; input the ranks directly)
  DATA_IVECTOR(cens); // censoring indicators, with 0 denotes right-censoring
  DATA_IVECTOR(ranks); // rank of each observation, correcting for ties by Breslow method
  int n = ranks.size(); // Sample size
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  DATA_SPARSE_MATRIX(D); // Differencing matrix to compute delta;
  DATA_SPARSE_MATRIX(design); // eta = design * W
  DATA_SCALAR(nu); // Matern shape
  DATA_SCALAR(rho_u);
  DATA_SCALAR(rho_alpha);
  DATA_SCALAR(sigma_u);
  DATA_SCALAR(sigma_alpha);
  DATA_SCALAR(betaprec);
  DATA_MATRIX(DS); // Distance matrix for calculating matern
  
  int d = P.cols(); // Number of B-Spline coefficients
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U), eta = design * U
  PARAMETER(theta); // theta = -2log(sigma)
  PARAMETER(logkappa); // Transformed matern params
  PARAMETER(logtau);

  double pi = 3.141592653589793115998;
  
  // rho and sigma
  Type kappa = exp(logkappa);
  Type tau = exp(logtau);
  Type rho = sqrt(8.0*nu) / kappa;
  Type sigma = tau / ( pow(kappa,nu) * sqrt( exp(lgamma(nu + 1.0)) * (4.0*pi) / exp(lgamma(nu))));
  
  
  // Split the param into B-Spline coefficients and polynomial coefficients
  int Wdim = W.size();
  vector<Type> U(d);
  for (int i=0;i<d;i++) U(i) = W(i+n);
  REPORT(U);
  int betadim = 3; // for the leuk example
  vector<Type> beta(betadim);
  for (int i=0;i<betadim;i++) beta(i) = W(i+d+n);

  // Transformations
  vector<Type> eta = design * W;
  REPORT(eta); // Check this works
  vector<Type> delta_red = D * eta;
  REPORT(delta_red);
  vector<Type> delta(n);
  delta(0) = 0;
  for (int i=1;i<n;i++) delta(i) = delta_red(i-1);
  REPORT(delta);
  
  // Log likelihood
  Type ll = 0;
  for (int i=0;i<n;i++){
    int nn = n-ranks(i)+1;
    vector<Type> delta_vec_i(nn); //rank starts at 1!!!
    for(int j=0;j<nn;j++) {
      delta_vec_i(j) = delta(n - nn + j);
      }
    vector<Type> diffvec(nn);
    for (int j=0;j<nn;j++) {
      diffvec(j) = exp(delta(i) - delta_vec_i(j));
      }
    ll += -cens(i) * log(diffvec.sum());
  }
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;
  // Cross product
  vector<Type> PU = P*U;
  Type UPU = (U * PU).sum();
  lpW += -0.5 * exp(theta) * UPU; // U part
  // also the fixed effect part!
  Type bb = (beta * beta).sum();
  REPORT(beta);
  REPORT(bb);
  lpW += -0.5 * betaprec * bb; // Beta part

  // Log determinant
  Type logdet1 = d * theta + logPdet;
  lpW += 0.5 * logdet1; // P part
  Type logdet2 = betadim * log(betaprec);
  lpW += 0.5 * logdet2; // beta part
  REPORT(logdet1);
  REPORT(lpW);
  REPORT(PU);
  REPORT(UPU);
  
  // Log prior for theta
  Type lpT = 0;
  Type phi = -log(alpha) / u;
  lpT += log(0.5 * phi) - phi*exp(-0.5*theta) - 0.5*theta;
  REPORT(lpT);
  
  Type lp = 0;
  // Prior for sigma,rho. with dimension 2 fixed, formula simplifies
  // the logkappa + logtau at the end is the jacobian
  Type lambda1 = -1.0 * (rho_u / sqrt(8.0*nu)) * log(rho_alpha);
  Type lambda2 = ( -1.0 * pow(kappa,-1.0 * nu) * sqrt( exp(lgamma(nu))  / ( exp(lgamma(nu + 1.0)) * (4.0*pi) ) ) ) * log(sigma_alpha) / sigma_u;
  Type lpt = log(lambda1) + log(lambda2) - lambda1 * kappa - lambda2 * tau + logkappa + logtau;
  lp += lpt;
  
  // prior for the matern family
  matrix<Type> C(DS);
  for(int i=0; i<C.rows(); i++)
    for(int j=0; j<C.cols(); j++)
      C(i,j) = pow(sigma,2.0) * matern(DS(i,j), rho / sqrt(8.0 * nu), nu);
  vector<Type> W1(n);
  for (int i=0;i<n;i++) W1(i) = W(i);
  Type nll1 = density::MVNORM_t<Type>(C)(W1);
  REPORT(nll1);
  REPORT(C);
  lp -= nll1; // Negative prior
  REPORT(lp);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT + lp);
  
  return logpost;
}
