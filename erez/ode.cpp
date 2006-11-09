/********************************************************************/

#include "ode.h"

using namespace std;

/********************************************************************/

/* constructors and destructors */

/********************************************************************/

Ode::Ode() : nvars(2), nsave(1), abs_tol(1e-6), rel_tol(1e-6), tmax(1.5),
             dt(1e-6) 
{
   strcpy(Ode::step_algo,"rkf45");
   strcpy(Ode::ofile_name,"no_ofile_name");
   strcpy(Ode::ic_file_name,"no_ic_file_name");
   Ode::rhs=NULL;
   Ode::_params=NULL;
   Ode::_p2derivs=NULL;
   Ode::_p2jac=NULL;
}

/********************************************************************/

Ode::~Ode()
{
   delete [] rhs;
}

/********************************************************************/

/* mutators */

/********************************************************************/

void Ode::SetNvars(const size_t nvars)
{
   Ode::nvars=nvars;
}

/********************************************************************/

void Ode::SetNsave(const int nsave)
{
   Ode::nsave=nsave;
}

/********************************************************************/

void Ode::SetStepAlgo(const char *step_algo)
{
   strcpy(Ode::step_algo,step_algo);
}

/********************************************************************/

void Ode::SetTol(const double abs_tol, const double rel_tol)
{
   Ode::abs_tol=abs_tol;
   Ode::rel_tol=rel_tol;
}

/********************************************************************/

void Ode::SetTmax(const double tmax)
{
   Ode::tmax=tmax;
}

/********************************************************************/

void Ode::SetDt(const double dt)
{
   Ode::dt=dt;
}

/********************************************************************/

void Ode::SetParams(void *params)
{
   Ode::_params=params;
}

/*******************************************************************/

void Ode::SetoFileName(const char *ofile_name)
{
   strcpy(Ode::ofile_name,ofile_name);
}

/********************************************************************/

void Ode::SeticFileName(const char *ic_file_name)
{
   strcpy(Ode::ic_file_name,ic_file_name);
}

/********************************************************************/

void Ode::SetRhs(double *rhs_ic)
{
   Ode::rhs=rhs_ic;
}

/*******************************************************************/

/* accessors */

/********************************************************************/

inline size_t Ode::Getnvars() const { return nvars; }

/********************************************************************/

inline int Ode::Getnsave() const { return nsave; }

/********************************************************************/

inline char *Ode::Getstep_algo() { return step_algo; }

/********************************************************************/

inline double Ode::Getatol() const { return abs_tol; }

/********************************************************************/

inline double Ode::Getrtol() const { return rel_tol; }

/********************************************************************/

inline double Ode::Gettmax() const { return tmax; }

/********************************************************************/

inline double Ode::Getdt() const { return dt; }

/********************************************************************/

inline double *Ode::Getrhs() const { return rhs; }

/********************************************************************/

inline char *Ode::Getofile_name() {return ofile_name; }

/********************************************************************/

inline char *Ode::Getic_file_name() {return ic_file_name; }

/********************************************************************/

/* gsl/odeiv solve stuff */

/********************************************************************/

void Ode::open_ofile(ofstream *ofile)
{
   cout << "... opening ode output file: " << ofile_name << endl;
   
   try
   {
      (*ofile).open(ofile_name, ios::out);
   }
   
   catch (exception &e)
   {
      cout << "... unable to open Ode output file " 
           << ofile_name << endl;
      cout << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }

   cout << "... done\n\n";   
}

/********************************************************************/

void Ode::close_ofile(ofstream *ofile)
{
   cout << "... closing ode output file: " << ofile_name << endl;
   
   try
   {
      (*ofile).close();
   }

   catch (exception &e)
   {
      cout << "... unable to close Ode output file " 
           << ofile_name << endl;
      cout << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   cout << "... done\n\n";
}

/********************************************************************/

void Ode::init_rhs_from_file()
{
   ifstream ic_file;
   double *rhs_ic;

   /* allocating rhs */
   try
   {      
      rhs_ic=new double[nvars]; // will be deleted in destructor
   }
   
   catch (exception &e)
   {
      cout << "... unable to alloc rhs\n"; 
      cout << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   /* opening ic_file */
   cout << "... opening ic_file: " << ic_file_name << endl;
   
   try
   {
      ic_file.open(ic_file_name, ios::in);
   }
   
   catch (exception &e)
   {
      cout << "... unable to open ic_file " 
           << ic_file_name << endl;
      cout << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   cout << "... done\n\n";
   
   /* init rhs from ic_file */
   for(unsigned int i=0; i<nvars; i++)
      ic_file >> rhs_ic[i];
   
   Setrhs(rhs_ic);
   
   /* closing ic_file */
   cout << "... closing ode output file: " << ofile_name << endl;
   
   try
   {
      ic_file.close();
   }
   
   catch (exception &e) {
      cout << "... unable to close ic_file " 
           << ic_file_name << endl;
      cout << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   cout << "... done\n\n";
}

/********************************************************************/

/* gsl/odeiv solve stuff */

/********************************************************************/

void Ode::ode_solve()
{
   double t;
   int status;
   ofstream ofile;
   
   Setstep();
   init_rhs_from_file();
   open_ofile(&ofile);
   
   gsl_odeiv_step *step = gsl_odeiv_step_alloc(step_type,nvars);
   gsl_odeiv_control *control = gsl_odeiv_control_y_new(abs_tol,rel_tol);
   gsl_odeiv_evolve *evolve  = gsl_odeiv_evolve_alloc(nvars);
   
   gsl_odeiv_system sys = {_p2derivs, _p2jac, nvars, _params};
   
   while (t < tmax)
   {
      status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tmax,
                                       &dt, rhs);      
      
      if(status != GSL_SUCCESS)
         exit(1);
      
      ofile << t << '\t' << rhs[0] << '\t' << rhs[1] << '\t' << dt  <<endl;
   }

   close_ofile(&ofile);
   
   gsl_odeiv_evolve_free (evolve);
   gsl_odeiv_control_free (control);
   gsl_odeiv_step_free (step);
}

/********************************************************************/

void Ode::Setstep()
{
   cout << "... setting stepping algorithm to "
        << step_algo << endl;
   
   if(!strcmp(step_algo,"rk2"))
   {   
      step_type = gsl_odeiv_step_rk2;
   }
   else if(!strcmp(step_algo,"rk4"))
   {
      step_type = gsl_odeiv_step_rk2;
   }
   else if(!strcmp(step_algo,"rkf45"))
   {
      step_type = gsl_odeiv_step_rkf45;
   }
   else if(!strcmp(step_algo,"rkck"))
   {
      step_type = gsl_odeiv_step_rkck;
   }
   else if(!strcmp(step_algo,"rk8pd"))
   {
      step_type = gsl_odeiv_step_rk8pd;
   }
   else if(!strcmp(step_algo,"rk2imp"))
   {
      step_type = gsl_odeiv_step_rk2imp;
   }
   else if(!strcmp(step_algo,"rk4imp"))
   {
      step_type = gsl_odeiv_step_rk4imp;
   }
   else if(!strcmp(step_algo,"bsimp"))
   {
      step_type = gsl_odeiv_step_bsimp;
   }
   else if(!strcmp(step_algo,"gear1"))
   {
      step_type = gsl_odeiv_step_gear1;
   }
   else if(!strcmp(step_algo,"gear2"))
   {
      step_type = gsl_odeiv_step_gear2;
   }
   else
   {
      cout << "... wrong step_algo value, set to default rkf45\n\n";
   }
   
   cout << "... done\n\n";
}

/********************************************************************/

void Ode::ode_plugin_funcs(int (*p2derivs)(double,const double *,double *,void *),
                      int (*p2jac)(double,const double *,double *,double *,void *))
{
   Ode::_p2derivs=p2derivs;
   Ode::_p2jac=p2jac;
}

/********************************************************************/

void Ode::prt_ode_prms() 
{
   cout << "Ode parameters:" << endl
        << "---------------" << endl
        << "nvars        = " << Getnvars() << endl
        << "nsave        = " << Getnsave() << endl
        << "step type    = " << Getstep_algo() << endl
        << "abs_tol      = " << Getatol() << endl
        << "rel_tol      = " << Getrtol() << endl
        << "tmax         = " << Gettmax() << endl
        << "dt           = " << Getdt() << endl
        << "ofile_name   = " << Getofile_name() << endl
        << "ic_file_name = " << Getic_file_name() << endl
        << endl;
}

/********************************************************************/
