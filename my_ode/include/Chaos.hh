#ifndef CHAOS_HH
#define CHAOS_HH

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>
#include <math.h>

//------------------------------------------------------------

struct ChaosParams
{
  ChaosParams() : nvars(0) {};
  
  // No. of equations
  unsigned int nvars;
  
  double alpha, beta, gamma, delta;

  unsigned int N, n;
  
  // overloading operator<<
  friend std::ostream& operator <<
    (std::ostream& os, const ChaosParams& x);
      
  // overloading operator=
  void operator= (const ChaosParams& x);
      
}; // ChaosParams

//------------------------------------------------------------

void ChaosParams::operator= (const ChaosParams& x)
{
  nvars = x.nvars;

  alpha  = x.alpha;
  beta  = x.beta;
  gamma = x.gamma;
  delta = x.delta;
  n = x.n;
  N = x.N;
  
} // operator=

//------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, const ChaosParams& x)
{
  os << std::endl
     << "Model parameters:" << std::endl
     << "-----------------" << std::endl
     << "alpha   = " << x.alpha << std::endl
     << "beta   = " << x.beta << std::endl
     << "gamma  = " << x.gamma << std::endl
     << "delta  = " << x.delta << std::endl
     << "n      = " << x.n << std::endl
     << "N      = " << x.N << std::endl
     << "nvars  = " << x.nvars << std::endl
     << std::endl;
  
  return os;
   
} // operator<<

//------------------------------------------------------------

struct Chaos
{

  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {         
    ChaosParams p = *(static_cast<ChaosParams*>(params));
    
    // local readable short variables
    double alpha=p.alpha, beta=p.beta, gamma=p.gamma, delta=p.delta;
    unsigned int n=p.n, N=p.N;

    unsigned int Sindex = 0, Iindex = 1, Rindex = 2;

    double S = y[Sindex], I = y[Iindex], R = y[Rindex];


    rhs[Sindex] = -S/N*(beta+alpha*pow(I, n)) + delta*R;
    rhs[Iindex] = +S/N*(beta+alpha*pow(I, n)) - gamma*I;
    rhs[Rindex] = + gamma*I - delta*R;

    return GSL_SUCCESS;         
  }
      
}; // Chaos

//------------------------------------------------------------
         
#endif /* CHAOS_HH */
