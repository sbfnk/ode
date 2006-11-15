/*******************************************************************/
//
// Class MeanField
// ---------------
//
// This is the MeanField class interface.
//
// See van_der_pol.hh for documentation.
//
/******************************************************************/

#ifndef MFA_H
#define MFA_H

/******************************************************************/

#include "ode.hh"

using namespace std;

/******************************************************************/

struct ModelParams
{
      double beta_d, gamma_d, delta_d;
      double beta_i, gamma_i, delta_i;
      double beta_m1, beta_m2, alpha, nu, lambda;
      double Qd, Qi, N;
      size_t njac;
};

class MeanField : public Ode
{
   public:
      MeanField();
      ~MeanField();
      
      // user supplied functions 
      static int derivs(double, const double *, double *, void *);
      static int jac(double, const double *, double *, double *, void *);
      
      // utility function 
      void PrtModelPrms();
      
}; // class MeanField

/******************************************************************/

#endif // MFA_H
