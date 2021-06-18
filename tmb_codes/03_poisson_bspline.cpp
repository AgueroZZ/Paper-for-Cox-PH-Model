#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Read in data
  DATA_VECTOR(y); // Response
  int n = y.size(); // Sample size
  DATA_SPARSE_MATRIX(BX); // Design matrix- cbind(BX,X)
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  int d = P.cols(); // Number of B-Spline coefficients
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u); // pc prior, u param
  DATA_SCALAR(alpha); // pc prior, alpha param
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)

  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta), eta = BX * U + X * beta
  PARAMETER(theta); // theta = -2log(sigma)
  
  // Split the param into B-Spline coefficients and polynomial coefficients
  int Wdim = W.size();
  int betadim = Wdim - d;
  vector<Type> U(d);
  vector<Type> beta(betadim);
  for (int i=0;i<d;i++) U(i) = W(i);
  for (int i=0;i<betadim;i++) beta(i) = W(i+d);
  REPORT(U);
  REPORT(beta);
  
  // Transformations
  vector<Type> eta = BX * W;
  REPORT(eta); // Check this works
  Type sigma = exp(-0.5*theta);
  REPORT(sigma);
  
  // Log likelihood
  Type ll = 0;
  for (int i=0;i<n;i++) ll += y(i) * eta(i) - exp(eta(i));
  REPORT(ll);
  
  // Log prior on W
  Type lpW = 0;
  // Cross product
  vector<Type> PU = P*U;
  Type UPU = (U * PU).sum();
  lpW += -0.5 * exp(theta) * UPU; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  
  // Log determinant
  Type logdet1 = d * theta + logPdet;
  lpW += 0.5 * logdet1; // P part
  Type logdet2 = betadim * log(betaprec);
  lpW += 0.5 * logdet2; // beta part
  REPORT(logdet1);
  REPORT(logdet2);
  REPORT(lpW);
  REPORT(PU);
  REPORT(UPU);
  
  // Log prior for theta
  Type lpT = 0;
  Type phi = -log(alpha) / u;
  lpT += log(0.5 * phi) - phi*exp(-0.5*theta) - 0.5*theta;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (ll + lpW + lpT);

  return logpost;
}
