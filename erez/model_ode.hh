/*******************************************************************/
//
// Class MeanField
// ---------------
//
// This is the MeanFiled class interface.
//
// See van_der_pol.hh for documentation.
//
/******************************************************************/

#ifndef MODEL_ODE_H
#define MODEL_ODE_H

/******************************************************************/

#include "ode.hh"

using namespace std;

/******************************************************************/

struct ModelParams
{
      // rates
      double beta_d, gamma_d, delta_d;
      double beta_i, gamma_i, delta_i;
      double beta_m1, beta_m2, alpha, nu, lambda;

      // Graph properties
      double Qd, Qi, N;

      // clustering
      double Cidi, Cidd, Cdii, Cdid;
      double Cddi, Cddd, Ciii, Ciid;

      // Jacobian size (if needed)
      size_t njac;      
};

class ModelOde : public Ode
{
   public:
      ModelOde(); 
      ~ModelOde();
      
      // mutators
      void SetCddi(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Cddi=C; }
      void SetCddd(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Cddd=C; }
      void SetCdii(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Cdii=C; }
      void SetCdid(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Cdid=C; }
      void SetCidi(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Cidi=C; }
      void SetCidd(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Cidd=C; }
      void SetCiii(const double C)
      { static_cast<ModelParams *>(GetModelParams())->Ciii=C; }
      void SetCiid(const double C) 
      { static_cast<ModelParams *>(GetModelParams())->Ciid=C; }
      void SetQd(const double Qd)
      { static_cast<ModelParams *>(GetModelParams())->Qd=Qd; }      
      void SetQi(const double Qi)
      { static_cast<ModelParams *>(GetModelParams())->Qi=Qi; }

      // user supplied functions      
      static int MFderivs(double, const double *, double *, void *);
      static int PAderivs(double, const double *, double *, void *);
      
      // helper function for PAderivs
      static double PA(double kk, double ij, double jk, double j);
      static double CC(double C_i, double C_d, double ki_i, double ki_d,
                       double ki, double N_Qd, double N_Qi);

      //static int MFjac(double, const double *, double *, double *, void *);
      
      // utility function
      void PluginMFderivs() { PluginFuncs(&MFderivs, NULL); }
      void PluginPAderivs() { PluginFuncs(&PAderivs, NULL); }
      
      void PrtModelPrms() const;
      void PrtGraphPrms() const;      
            
}; // class ModelOde

/******************************************************************/

#endif // MODEL_ODE_H
