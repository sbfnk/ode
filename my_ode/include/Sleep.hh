#ifndef SLEEP_HH
#define SLEEP_HH

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>

//------------------------------------------------------------

struct SleepParams
{
  SleepParams() : nvars(3) {};
  
  // No. of equations
  unsigned int nvars;
  
  double beta, gamma, delta;
//   double sigma, tau;
  
  // overloading operator<<
  friend std::ostream& operator <<
    (std::ostream& os, const SleepParams& x);
      
  // overloading operator=
  void operator= (const SleepParams& x);
      
}; // SleepParams

//------------------------------------------------------------

void SleepParams::operator= (const SleepParams& x)
{
  nvars = x.nvars;
  
  beta  = x.beta;
  gamma = x.gamma;
  delta = x.delta;
//   sigma = x.sigma;
//   tau = x.tau;
  
} // operator=

//------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, const SleepParams& x)
{
  os << std::endl
     << "Model parameters:" << std::endl
     << "-----------------" << std::endl
     << "nvars = " << x.nvars << std::endl
     << "gamma = " << x.gamma << std::endl
     << "beta  = " << x.beta << std::endl
     << "delta = " << x.delta << std::endl
//      << "sigma = " << x.sigma << std::endl
//      << "tau   = " << x.tau << std::endl
     << std::endl;
  
  return os;
   
} // operator<<

//------------------------------------------------------------

struct Sleep
{

//   static double screening_rate(double sigma, double tau, double t)
//   {
//     if (t >= tau) {
//       return sigma * (t - tau);
//     } else {
//       return 0;
//     }
//   }
  
  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {         
    SleepParams p = *(static_cast<SleepParams*>(params));
    
    // local readable short variables
    double beta=p.beta, gamma=p.gamma, delta=p.delta;
//     double tau=p.tau, sigma=p.sigma;

//     rhs[0] = mu*y[2] - beta*y[0]*(y[1]+y[2]) + screening_rate(sigma, tau, t) * y[1];
//     rhs[1] = beta*y[0]*(y[1]+y[2]) - gamma*y[1] - screening_rate(sigma, tau, t) * y[1];
//     rhs[2] = gamma*y[1] - mu*y[2];
//    rhs[3] = screening_rate(sigma, tau, t)*y[1];

    rhs[0] = delta*y[2] - beta*y[0]*y[1];
    rhs[1] = beta*y[0]*y[1] - gamma*y[1];
    rhs[2] = gamma*y[1] - delta*y[2];

    return GSL_SUCCESS;         
  }
      
}; // Sleep

//------------------------------------------------------------
         
#endif /* SLEEP_HH */
