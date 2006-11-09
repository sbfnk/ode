//*******************************************************************
//
// Class Ode
// ---------
//
// This is the Ode class interface. This class solves ODE systems
// using GSL - GNU Scientific Library (www.gnu.org/software/gsl).
// This is a baes class, i.e. in order to solve a specific system
// the user must provide the class with pointers to two functions,
// one for calculating the RHS and one for the Jacobian. The user
// must also provide the relevant parameters of the model. See as
// an example the derived class VanDerPol in van_der_pol.h .
//
// Parameters:
// -----------
// The class's parameters are ODE related. The specific model
// parameters should be plugged-in from derived classes using an
// Ode member function.
//
// some noteworthy parameters are listed below:
//
// nvars - # of variables.
// step_algo - stepping algorithm type (see GSL docs).
// abs_tol, rel_tol - absolute/relative tolerances (see GSL docs).
// dt - initial step size.
// *rhs - store derivatives, i.e. du/dt=rhs(u,t).
// *_params - pointer to model parameters.
// *_p2derivs - pointer to function derivs, which computes the rhs.
// *_p2jac - pointer to function jac, which computes the Jacobian.
//
// Functions:
// ----------
// Constructor, destructor, mutators and accessors - straight forward.
// OdePluginFuncs - Once the functions for computing the derivatives and 
//    the Jacobian, this function enables the user to set the Ode class
//    pointers to those functions. It is used in the derived class
//    constructor. 
// SetParams - Once the structure of parameters has been defined by the
//    user in the derived class and has been initialized, a pointer to
//    this structure has to be passed to initialize *_params. This is done
//    using SetParams. The user defined derived class member function
//    PluginModelParams() passes the required pointer to *_params after 
//    casting it to void, i.e. static_cast<void *>(...).
//    using
// InitRhsFromFile - allocate rhs as rhs_ic ...
// OdeSolve - runs InitRhsFromFile, open and close dat file, allocate all
//    necessary objects for gsl/odevi and perform the main loop ...
// 
//  
//
// Special notices:
// ----------------
// - rhs is actually allocated in private member function
//   Ode::InitRhsFromFile(), as rhs_ic. The array rhs_ic is then initialized
//   from ic_file, and a pointer to this array is copied to rhs by
//   Ode::SetRhs(rhs_ic). The array is set free at Ode::~Ode().
// - step_algo can be: rk2, rk4, rkf45, rkck, rk8pd, rk2imp, rk4imp, bsimp,
//                     gear1, gear2.
// - The signature of the functions derivs and jac, and therefor the form of
//   the pointers *_p2derivs ans *_p2jac, is dictated by GSL (see GSL docs).
//
/********************************************************************/

#ifndef ODE_H
#define ODE_H

/********************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <exception>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

using namespace std;

/********************************************************************/

const size_t MAX_STR_LEN=64; 

/********************************************************************/

class Ode
{
  public:
   
   /* constructors and destructors */
   Ode();
   ~Ode();
   
   /* mutators */
   void SetNvars(const size_t nvars);
   void SetNsave(const int nsave);
   void SetStepAlgo(const char *step_algo);
   void SetTol(const double abs_tol, const double rel_tol);
   void SetTmax(const double tmax);
   void SetDt(const double dt);
   void SetParams(void *params);
   void SetoFileName(const char *ofile_name);
   void SeticFileName(const char *ic_file_name);
   void SetRhs(double *rhs_ic);
   
   /* accessors */
   size_t GetNvars() const;
   int GetNsave() const;
   char *GetStepAlgo();
   double GetAtol() const;
   double GetRtol() const;
   double GetTmax() const;
   double GetDt() const;
   double *GetRhs() const;
   char *GetoFileName();
   char *GeticFileName();
   
   /* gsl/odeiv solve stuff */
   void OdeSolve();
   void OdePluginFuncs(int (*p2derivs)(double,const double *,double *,void *),
                       int (*p2jac)(double,const double *,double *,double *,void *));

   /* i/o stuff */
   void PrtOdePrms();
   
  private:
   
   /* gsl/odeiv solve stuff */
   const gsl_odeiv_step_type *SetStepType();
   void OpenoFile(ofstream *ofile);
   void CloseoFile(ofstream *ofile);
   void InitRhsFromFile();
   
   /* ode parameters */
   size_t nvars;                // size of the ode system
   int nsave;                   // save solution every nsave time steps
   char step_algo[MAX_STR_LEN]; // name of stepping algorithm type
   double abs_tol, rel_tol;     // absolute/relative tolerances
   double tmax, dt;             // ...
   double *rhs;                 // rhs[navrs] for the variables
   void *_params;               // pointer to model parameters. will be set
                                // by the user to some specific model
                                // parameters using the member function
                                // Ode::SetParams(void *params);
   
   /* output file */
   char ofile_name[MAX_STR_LEN];   // output file name
   char ic_file_name[MAX_STR_LEN]; // ic input file name
   
   /* pointers to rhs and jacobian declaration */
   int (*_p2derivs)(double, const double *, double *, void *);
   int (*_p2jac)(double, const double *, double *, double *, void *);
   
}; // class Ode

#endif // ODE_H
