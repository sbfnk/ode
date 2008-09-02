#ifndef DIMINISH_HH
#define DIMINISH_HH

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>
#include <math.h>

//------------------------------------------------------------

struct DiminishParams
{
  DiminishParams() : nvars(0) {};
  
  // No. of equations
  unsigned int nvars;
  
  double alpha, beta, gamma, delta, lambda;
  double omega, rho;

  unsigned int N, M;
  
  // overloading operator<<
  friend std::ostream& operator <<
    (std::ostream& os, const DiminishParams& x);
      
  // overloading operator=
  void operator= (const DiminishParams& x);
      
}; // DiminishParams

//------------------------------------------------------------

void DiminishParams::operator= (const DiminishParams& x)
{
  nvars = x.nvars;

  alpha = x.alpha;
  beta  = x.beta;
  gamma = x.gamma;
  delta = x.delta;
  lambda = x.lambda;
  M = x.M;
  N = x.N;
  omega = x.omega;
  rho = x.rho;
  
} // operator=

//------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, const DiminishParams& x)
{
  os << std::endl
     << "Model parameters:" << std::endl
     << "-----------------" << std::endl
     << "alpha  = " << x.alpha << std::endl
     << "beta   = " << x.beta << std::endl
     << "gamma  = " << x.gamma << std::endl
     << "delta  = " << x.delta << std::endl
     << "lambda = " << x.lambda << std::endl
     << "omega  = " << x.omega << std::endl
     << "rho    = " << x.rho << std::endl
     << "M      = " << x.M << std::endl
     << "N      = " << x.N << std::endl
     << "nvars  = " << x.nvars << std::endl
     << std::endl;
  
  return os;
   
} // operator<<

//------------------------------------------------------------

struct Diminish
{

  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {         
    DiminishParams p = *(static_cast<DiminishParams*>(params));
    
    // local readable short variables
    double alpha=p.alpha, beta=p.beta, gamma=p.gamma, delta=p.delta;
    double lambda=p.lambda, omega=p.omega, rho=p.rho;
    unsigned int M=p.M, N=p.N;

    double Nsum=0, Ssum=0, Isum=0, Rsum = 0;

    for (unsigned int i = 0; i <= M; ++i) {
      Ssum += y[3*i];
      Isum += y[3*i+1];
      Rsum += y[3*i+2];
      std::cout << "S[" << i << "]=" << y[3*i] << "  ";
      std::cout << "I[" << i << "]=" << y[3*i+1] << "  ";
      std::cout << "R[" << i << "]=" << y[3*i+2] << std::endl;
    }

    std::cout << "S=" << Ssum << "  ";
    std::cout << "I=" << Isum << "  ";
    std::cout << "R=" << Rsum << std::endl;

    double IfullSum = Isum;

    for (unsigned int i = 0; i <= M; ++i) {
      unsigned int Sindex = 3*i, Iindex = Sindex + 1, Rindex = Sindex + 2;

      double S_im1 = 0, I_im1 = 0, R_im1 = 0;

      if (i > 0) {
        S_im1 = y[Sindex-3]; I_im1 = y[Iindex-3]; R_im1 = y[Rindex-3];
      }
      
      double N_im1 = S_im1 + I_im1 + R_im1;
      double S_i = y[Sindex], I_i = y[Iindex], R_i = y[Rindex];


      rhs[Sindex] = -(1-pow(rho, i))*beta*S_i/N*IfullSum +
        delta*R_i -
        alpha*S_i/N*Nsum +
        alpha*N_im1/N*Ssum + 
        lambda*S_im1;
      rhs[Iindex] = +(1-pow(rho, i))*beta*S_i/N*IfullSum -
        gamma*I_i -
        alpha*I_i/N*Nsum +
        alpha*N_im1/N*Isum +
	lambda*I_im1;
      rhs[Rindex] = - delta*R_i +
        gamma*I_i -
        alpha*R_i/N*Nsum +
        alpha*N_im1/N*Rsum +
	lambda*R_im1;

      if (i < M) {
        rhs[Sindex] -= lambda * S_i;
        rhs[Iindex] -= lambda * I_i;
        rhs[Rindex] -= lambda * R_i;
      }

      if (i == 0) {
        rhs[Iindex] += omega * (IfullSum - I_i);
      } else {
        rhs[Iindex] -= omega * I_i;
      }

      std::cout << "S'[" << i << "]=" << rhs[Sindex] << " ";
      std::cout << "I'[" << i << "]=" << rhs[Iindex] << " ";
      std::cout << "R'[" << i << "]=" << rhs[Rindex] << std::endl;

      Nsum += S_i + I_i + R_i;
      Ssum -= S_i;
      Isum -= I_i;
      Rsum -= R_i;
    }

    return GSL_SUCCESS;         
  }
      
}; // Diminish

//------------------------------------------------------------
         
#endif /* DIMINISH_HH */
