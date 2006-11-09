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
   void SetStep();
   void OpenoFile(ofstream *ofile);
   void CloseoFile(ofstream *ofile);
   void InitRhsFromFile();
   
   /* ode parameters */
   size_t nvars;                // size of the ode system
   int nsave;                   // save solution every nsave time steps
   char step_algo[MAX_STR_LEN]; // type of stepping algorithm
   double abs_tol, rel_tol;     // absolute/relative tolerances
   double tmax, dt;             //
   double *rhs;                 // rhs[navrs] for the variables
   void *_params;               // model parameters
   
   /* output file */
   char ofile_name[MAX_STR_LEN];   // output file name
   char ic_file_name[MAX_STR_LEN]; // ic input file name
   
   /* gsl/odeiv parameter */
   const gsl_odeiv_step_type *step_type;
   
   /* pointers to rhs and jacobian declaration */
   int (*_p2derivs)(double, const double *, double *, void *);
   int (*_p2jac)(double, const double *, double *, double *, void *);
   
}; // class Ode

#endif // ODE_H
