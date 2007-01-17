/******************************************************************/

#include "ode.hh"

using namespace std;

/******************************************************************/

// constructors and destructors 

/******************************************************************/
Ode::Ode(bool v) : nvars(2), nsave(1), abs_tol(1e-6), rel_tol(1e-6),
                   tmax(1.5), dt(1e-6), verbose(v)
{
   strcpy(Ode::step_algo,"rkf45");
   strcpy(Ode::ofile_name,"no_ofile_name");
   strcpy(Ode::ic_file_name,"no_ic_file_name");
   Ode::rhs=NULL;
   Ode::params=NULL;
   Ode::p2derivs=NULL;
   Ode::p2jac=NULL;
}

/******************************************************************/

//Ode::~Ode()
//{
//   delete [] rhs;
//}

/******************************************************************/

// gsl/odeiv solve stuff 

/******************************************************************/

void Ode::OpenoFile(ofstream& ofile)
{
  if (verbose) cout << "... opening ode output file: " << ofile_name;
   
   try
   {
      ofile.open(ofile_name, ios::out);
   }
   
   catch (exception &e)
   {
      cerr << "... unable to open Ode output file " 
           << ofile_name << endl;
      cerr << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }

   if (verbose) cout << " ... done\n";   
}

/******************************************************************/

void Ode::CloseoFile(ofstream& ofile)
{
   if (verbose) cout << "... closing ode output file: " << ofile_name;
   
   try
   {
      ofile.close();
   }

   catch (exception &e)
   {
      cerr << "... unable to close Ode output file " 
           << ofile_name << endl;
      cerr << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   if (verbose) cout << " ... done\n";
}

/******************************************************************/

void Ode::InitRhsFromFile()
{
   ifstream ic_file;
   double *rhs_ic;

   
   // allocating rhs 
   if (verbose) cout << "... allocating rhs memory";
   
   try
   {      
      rhs_ic=new double[nvars]; // will be deleted in destructor
   }
   
   catch (exception &e)
   {
      cerr << "... unable to alloc rhs\n"; 
      cerr << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }

   if (verbose) cout << " ... done\n";
   
   // opening ic_file 
   if (verbose) cout << "... opening ode ic_file: " << ic_file_name;
   
   try
   {
      ic_file.open(ic_file_name, ios::in);
   }
   
   catch (exception &e)
   {
      cerr << "... unable to open ic_file " 
           << ic_file_name << endl;
      cerr << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   if (verbose) cout << " ... done\n";
   
   // init rhs from ic_file 
   if (verbose) cout << "... reading content of ic_file";
   
   for(unsigned int i=0; i<nvars; i++)
      ic_file >> rhs_ic[i];
   
   SetRhs(rhs_ic);

   if (verbose) cout << " ... done\n";
   
   // closing ic_file 
   if (verbose) cout << "... closing ode ic_file: " << ofile_name;
   
   try
   {
      ic_file.close();
   }
   
   catch (exception &e)
   {
      cerr << "... unable to close ic_file " 
           << ic_file_name << endl;
      cerr << "... Standard exception: " << e.what() << endl;      
      exit(1); // uses cstdlib.h
   }
   
   if (verbose) cout << " ... done\n";
}

/******************************************************************/

// gsl/odeiv solve stuff 

/******************************************************************/

void Ode::Solve()
{
   double t;
   int status;
   ofstream ofile;

   if (verbose) cout << "Solving the ode system\n"
                     << "----------------------\n";

   // setting step type 
   const gsl_odeiv_step_type *step_type;
   step_type = SetStepType();

   // allocating rhs and initialize it from ic_file 
   InitRhsFromFile();

   // opening output file 
   OpenoFile(ofile);

   // allocating stepping function 
   gsl_odeiv_step *step = gsl_odeiv_step_alloc(step_type,nvars);

   // allocating control function 
   gsl_odeiv_control *control = gsl_odeiv_control_y_new(abs_tol,rel_tol);

   // allocating evolution function 
   gsl_odeiv_evolve *evolve  = gsl_odeiv_evolve_alloc(nvars);

   // defining the system 
   gsl_odeiv_system sys = {p2derivs, p2jac, nvars, params};

   // write t=0 rhs 
   t=0;
   WriteRHS(ofile, t);

   unsigned int o_count=1;
   
   ///// main loop /////
   if (verbose) cout << "... in main loop";   
   
   while (t < tmax)
   {
      // stepping solution 
      status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tmax,
                                       &dt, rhs);      
      
      if(status != GSL_SUCCESS)
         exit(1);
      
      // writing RHS to ofile 
      if(!(o_count++ % nsave))
         WriteRHS(ofile, t);
   }
   
   if (verbose) cout << ".................. done\n";
   
   CloseoFile(ofile);
   
   // free allocated memory 
   gsl_odeiv_evolve_free(evolve);
   gsl_odeiv_control_free(control);
   gsl_odeiv_step_free(step);

   delete [] rhs;
   rhs = NULL;
}

/******************************************************************/

const gsl_odeiv_step_type *Ode::SetStepType()
{
   const gsl_odeiv_step_type *tmp_step_type;

   if (verbose) cout << "... setting stepping algorithm to " << step_algo;
   
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
     if (verbose) cout << "... wrong step_algo value, set to default rkf45\n";
   }
   
   if (verbose) cout << " ... done\n";

   return tmp_step_type;
}

/******************************************************************/

void Ode::PluginFuncs(int (*p2derivs)(double,const double *,double *,void *),
                      int (*p2jac)(double,const double *,double *,double *,void *))
{
   Ode::p2derivs=p2derivs;
   Ode::p2jac=p2jac;
}

/******************************************************************/

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

/******************************************************************/

void Ode::WriteRHS(ofstream& ofile, double t)
{
   ofile << t;
   for(unsigned int i=0; i<nvars; i++)
      ofile << '\t' << rhs[i];
   ofile << endl;
}

/******************************************************************/
