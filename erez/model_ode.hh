/*******************************************************************/
//
// Class MeanField
// ---------------
//
// This is the VanDerPol class interface. This class is AN EXAMPLE for
// using the Ode class. This class is a derived class of Ode, that
// provides the Ode class with the specific Model equations and
// parameters. The user should define the I/O functionality to read
// the relevant parameters (e.g. ode_io_utils.hh).
// The class contains two function, required by gsl/odeiv, for
// calculating the rhs and the Jacobian. Pointers to these functions
// are passed to the Ode object in the class constructor, using
// Ode::OdePluginFuncs. These functions MUST be declared static,
// so the pointers to them will NOT be of type Ode::*p2func, which
// is not accepted by the GSL functions. These function can call
// other (static) help functions.
//
// I/O: See Ode.hh documentation.
//
// Parameters:
// -----------
// The class's parameters are Model related. The specific model
// parameters should be plugged-in from derived classes using the
// Ode member function Ode::SetModelParams. This is done at the end
// the ReadModelParams function, after all parameters were read.
// The parameters are stored in a structure defined as
//
// struct ModelParams
// {
//     ....
// };
//
// Functions:
// ----------
// Constructor - The constructor initialize the pointers to the
//               derivs and jac functions, that claculate the RHS
//               and Jacobian of the system. It uses OdePluginFuncs.
//
// destructor - The destructor free the memory allocated for the Model
//              parameters structure.
//
// derivs - Compute the RHS of the system. Declared as static and may
//          call help functions. The signature of this function is
//          dictated by the GSL functions.
//
// jac - Compute the Jacobian of the system. Declared as static and may
//       call help functions. The signature of this function is dictated
//       by the GSL functions.
//
// PrtModelPrms - Print the Ode and Model parameters. It first calls
//                Ode::PrtOdePrms and then prints out the Model
//                parameters to the std::out.
//
// Special notices:
// ----------------
// - If a Jacobian is not supplied, the pointer to the function that
// calculates it should be set to NULL (see the constructor).
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
