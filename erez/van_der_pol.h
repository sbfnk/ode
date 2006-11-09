/********************************************************************/

#ifndef VAN_DER_POL_H
#define VAN_DER_POL_H

/********************************************************************/

#include "ode.h"

using namespace std;

/********************************************************************/

struct ModelParams
{
      double beta_d, gamma_d, delta_d;
      double beta_i, gamma_i, delta_i;
      double beta_m, alpha, nu, lambda;
      double mu1, mu2;
      size_t njac;
};

class VanDerPol : public Ode
{
  public:
   VanDerPol();
   ~VanDerPol();

   /* user supplied functions */
   static int derivs(double, const double *, double *, void *);
   static int jac(double, const double *, double *, double *, void *);

   /* utility function */
   void set_model_params(ModelParams *model_params);
   ModelParams *get_model_params() { return model_params; };
   void plugin_model_params();      
   void prt_model_prms();
   
  private:
   ModelParams *model_params;

}; // class VanDerPol

/********************************************************************/

#endif // VAN_DER_POL_H
