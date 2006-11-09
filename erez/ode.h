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
   void set_nvars(const size_t nvars);
   void set_nsave(const int nsave);
   void set_step_algo(const char *step_algo);
   void set_tol(const double abs_tol, const double rel_tol);
   void set_tmax(const double tmax);
   void set_dt(const double dt);
   void set_params(void *params);
   void set_ofile_name(const char *ofile_name);
   void set_ic_file_name(const char *ic_file_name);
   void set_rhs(double *rhs_ic);
   
   /* accessors */
   size_t get_nvars() const;
   int get_nsave() const;
   char *get_step_algo();
   double get_atol() const;
   double get_rtol() const;
   double get_tmax() const;
   double get_dt() const;
   double *get_rhs() const;
   char *get_ofile_name();
   char *get_ic_file_name();
   
   /* gsl/odeiv solve stuff */
   void ode_solve();
   void ode_plugin_funcs(int (*p2derivs)(double,const double *,double *,void *),
                         int (*p2jac)(double,const double *,double *,double *,void *));

   /* i/o stuff */
   void prt_ode_prms();
   
  private:
   
   /* gsl/odeiv solve stuff */
   void set_step();
   void open_ofile(ofstream *ofile);
   void close_ofile(ofstream *ofile);
   void init_rhs_from_file();
   
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
