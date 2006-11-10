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

void Ode::SetModelParams(void *params)
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

size_t Ode::GetNvars() const { return nvars; }

/********************************************************************/

int Ode::GetNsave() const { return nsave; }

/********************************************************************/

char *Ode::GetStepAlgo() { return step_algo; }

/********************************************************************/

double Ode::GetAtol() const { return abs_tol; }

/********************************************************************/

double Ode::GetRtol() const { return rel_tol; }

/********************************************************************/

double Ode::GetTmax() const { return tmax; }

/********************************************************************/

double Ode::GetDt() const { return dt; }

/********************************************************************/

void *Ode::GetModelParams() const { return _params; }

/********************************************************************/

double *Ode::GetRhs() const { return rhs; }

/********************************************************************/

char *Ode::GetoFileName() {return ofile_name; }

/********************************************************************/

char *Ode::GeticFileName() {return ic_file_name; }

/********************************************************************/

/* gsl/odeiv solve stuff */

/********************************************************************/

void Ode::OpenoFile(ofstream *ofile)
{
   cout << "... opening ode output file: " << ofile_name;
   
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

   cout << " ... done\n";   
}

/********************************************************************/

void Ode::CloseoFile(ofstream *ofile)
{
   cout << "... closing ode output file: " << ofile_name;
   
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
   
   cout << " ... done\n";
}

/********************************************************************/

void Ode::InitRhsFromFile()
{
   ifstream ic_file;
   double *rhs_ic;

   
   /* allocating rhs */
   cout << "... allocating rhs memory";
   
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

   cout << " ... done\n";
   
   /* opening ic_file */
   cout << "... opening ode ic_file: " << ic_file_name;
   
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
   
   cout << " ... done\n";
   
   /* init rhs from ic_file */
   cout << "... reading content of ic_file";
   
   for(unsigned int i=0; i<nvars; i++)
      ic_file >> rhs_ic[i];
   
   SetRhs(rhs_ic);

   cout << " ... done\n";
   
   /* closing ic_file */
   cout << "... closing ode ic_file: " << ofile_name;
   
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
   
   cout << " ... done\n";
}

/********************************************************************/

/* gsl/odeiv solve stuff */

/********************************************************************/

void Ode::OdeSolve()
{
   double t;
   int status;
   ofstream ofile;

   cout << "Solving the ode system\n"
        << "----------------------\n";

   /* setting step type */
   const gsl_odeiv_step_type *step_type;
   step_type = SetStepType();

   /* allocating rhs and initialize it from ic_file */
   InitRhsFromFile();

   /* opening output file */
   OpenoFile(&ofile);

   /* allocating stepping function */
   gsl_odeiv_step *step = gsl_odeiv_step_alloc(step_type,nvars);

   /* allocating control function */
   gsl_odeiv_control *control = gsl_odeiv_control_y_new(abs_tol,rel_tol);

   /* allocating evolution function */
   gsl_odeiv_evolve *evolve  = gsl_odeiv_evolve_alloc(nvars);

   /* defining the system */
   gsl_odeiv_system sys = {_p2derivs, _p2jac, nvars, _params};

   /*** main loop ***/
   cout << "... in main loop";

   /* write t=0 rhs */
   t=0;
   ofile << t << '\t' << rhs[0] << '\t' << rhs[1] << '\t' << dt  <<endl;
   
   while (t < tmax)
   {
      status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tmax,
                                       &dt, rhs);      
      
      if(status != GSL_SUCCESS)
         exit(1);

      // USE NSAVE !!!!!!!
      ofile << t << '\t' << rhs[0] << '\t' << rhs[1] << '\t' << dt  <<endl;
   }
   
   cout << " .................. done\n";
   
   CloseoFile(&ofile);

   /* free allocated memory */
   gsl_odeiv_evolve_free (evolve);
   gsl_odeiv_control_free (control);
   gsl_odeiv_step_free (step);
}

/********************************************************************/

const gsl_odeiv_step_type *Ode::SetStepType()
{
   const gsl_odeiv_step_type *tmp_step_type;

   cout << "... setting stepping algorithm to " << step_algo;
   
   if(!strcmp(step_algo,"rk2"))
   {   
      tmp_step_type = gsl_odeiv_step_rk2;
   }
   else if(!strcmp(step_algo,"rk4"))
   {
      tmp_step_type = gsl_odeiv_step_rk2;
   }
   else if(!strcmp(step_algo,"rkf45"))
   {
      tmp_step_type = gsl_odeiv_step_rkf45;
   }
   else if(!strcmp(step_algo,"rkck"))
   {
      tmp_step_type = gsl_odeiv_step_rkck;
   }
   else if(!strcmp(step_algo,"rk8pd"))
   {
      tmp_step_type = gsl_odeiv_step_rk8pd;
   }
   else if(!strcmp(step_algo,"rk2imp"))
   {
      tmp_step_type = gsl_odeiv_step_rk2imp;
   }
   else if(!strcmp(step_algo,"rk4imp"))
   {
      tmp_step_type = gsl_odeiv_step_rk4imp;
   }
   else if(!strcmp(step_algo,"bsimp"))
   {
      tmp_step_type = gsl_odeiv_step_bsimp;
   }
   else if(!strcmp(step_algo,"gear1"))
   {
      tmp_step_type = gsl_odeiv_step_gear1;
   }
   else if(!strcmp(step_algo,"gear2"))
   {
      tmp_step_type = gsl_odeiv_step_gear2;
   }
   else
   {
      cout << "... wrong step_algo value, set to default rkf45\n";
   }
   
   cout << " ... done\n";

   return tmp_step_type;
}

/********************************************************************/

void Ode::OdePluginFuncs(int (*p2derivs)(double,const double *,double *,void *),
                         int (*p2jac)(double,const double *,double *,double *,void *))
{
   Ode::_p2derivs=p2derivs;
   Ode::_p2jac=p2jac;
}

/********************************************************************/

void Ode::PrtOdePrms() 
{
   cout << endl
        << "Ode parameters:" << endl
        << "---------------" << endl
        << "nvars        = " << GetNvars() << endl
        << "nsave        = " << GetNsave() << endl
        << "step type    = " << GetStepAlgo() << endl
        << "abs_tol      = " << GetAtol() << endl
        << "rel_tol      = " << GetRtol() << endl
        << "tmax         = " << GetTmax() << endl
        << "dt           = " << GetDt() << endl
        << "ofile_name   = " << GetoFileName() << endl
        << "ic_file_name = " << GeticFileName() << endl
        << endl;
}

/********************************************************************/
