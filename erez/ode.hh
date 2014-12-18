/*******************************************************************/
//
// Class Ode
// ---------
//
// This is the Ode class interface. This class solves ODE systems
// using GSL - GNU Scientific Library (www.gnu.org/software/gsl).
// This is a baes class, i.e. in order to solve a specific system
// the user must provide a derived class with pointers to two
// functions, one for calculating the RHS and one for the Jacobian.
// The user must also provide the relevant parameters of the model.
// See as an example the derived class VanDerPol in van_der_pol.h .
//
// I/O:
// ----
// Input
//
// *.prm - Parameters input file. The class reads the required
//    parameters
//    from a file with .prm suffix, e.g. model.prm. Both Ode related
//    and Model related parameters are to be supplied using this file.
//    The specific IO implementation is user dependent, using mutators,
//    e.g. ode_io_utils.h. The IO functions look for the string "Model
//    parameters" and "Ode parameters" for reading the model and ode
//    parameters, respectively.
//    The reading is done by order, not by name !
//
// *.init - Initial condition file. This file containes the values of
//    the system variables at t=0. Each value is in a separate line.
//    This file is read by order ! The function that reads this file
//    is Ode::InitRhsFromFile().
//
// Output
//
// *.dat - Solutions output file. The solutions for the system is
//    written to this file every nsave time steps. The format is
//    copetible with Gnuplot, i.e.
//
//    time y[0] y[1] y[2] ....
//
// std::out - The object will write messages to the standard output
//    reporting general operations like opening and closing files, reading
//    and writing from them, and parameters values. It is recommended
//    to output these messages into a *.log file, e.g.
//
//    ./a.out > file.log 
//
// Parameters:
// -----------
// The class's parameters are ODE related. The specific model
// parameters should be plugged-in from derived classes using the
// Ode member function Ode::SetModelParams.
//
// some noteworthy parameters are listed below:
//
// nvars - # of variables.
// step_algo - stepping algorithm type (see GSL docs).
// abs_tol, rel_tol - absolute/relative tolerances (see GSL docs).
// dt - initial step size.
// *rhs - store derivatives, i.e. du/dt=rhs(u,t).
// *params - pointer to model parameters structure.
// *p2derivs - pointer to function derivs, which computes the rhs.
// *p2jac - pointer to function jac, which computes the Jacobian.
//
// Functions:
// ----------
// Constructor, destructor, mutators and accessors - straight forward.
//
// OdePluginFuncs - Once the functions for computing the derivatives and 
//    the Jacobian are defined by the user, this function enables the user
//    to set the Ode class pointers to those functions. It is used in the
//    derived class constructor.
//
// SetModelParams - Once the structure of parameters has been defined by the
//    user in the derived class and has been initialized, a pointer to
//    this structure has to be passed to initialize *_params. This is done
//    using SetModelParams. For example, once the user has defined and
//    initialized the parameters structure, this member function is used
//    to initialize *_params, after casting it to void. See ReadModelParams
//    in ode_io_utils.h. The structure is deleted in the derived class
//    destructor.
//
// InitRhsFromFile - The memory allocation for *rhs is actually done in this
//    function, as *rhs_ic. The array rhs_ic is then initialized from a file
//    (ic_file), and the pointer to this array is assigned to rhs using the
//    member function Ode::SetRhs. The array is deleted by the object
//    destructor. 
//
// OdeSolve - This is the main function that actually SOLVE the ODE system.
//    It Sets the stepping algorithm type by invoking Ode::SetStepType, and
//    allocate the stepping, controlling, evolution and system objects (see
//    GSL docs). The output file is opened, results from the main loop are
//    written to it, then it is closed and the allocated memory is returned.
//
// Special notices:
// ----------------
// - rhs is actually allocated in private member function Ode::InitRhsFromFile(),
//   as rhs_ic. The array rhs_ic is then initialized from ic_file, and a pointer
//   to this array is copied to rhs by Ode::SetRhs(rhs_ic). The array is set free
//   at Ode::~Ode().
// - step_algo can be: rk2, rk4, rkf45, rkck, rk8pd, rk2imp, rk4imp, bsimp,
//                     gear1, gear2.
// - The signature of the functions derivs and jac, and therefor the form of
//   the pointers *p2derivs ans *p2jac, is dictated by GSL (see GSL docs).
//
/******************************************************************/

#ifndef ODE_H
#define ODE_H

/******************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <exception>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

using namespace std;

/******************************************************************/

const size_t MAX_STR_LEN=128; 

/******************************************************************/

class Ode
{
   public:
      
      // constructors and destructors 
      Ode(bool v = false);
      virtual ~Ode() {;}
  
      
      // mutators 
      void SetNvars(const size_t nvars) { Ode::nvars=nvars; }
      void SetNinit(const size_t ninit) { Ode::ninit=ninit; }
      void SetNsave(const unsigned int nsave) { Ode::nsave=nsave; }
      void SetStepAlgo(const char *step_algo) { strcpy(Ode::step_algo,step_algo); }
      void SetAbsTol(const double abs_tol) { Ode::abs_tol=abs_tol; }
      void SetRelTol(const double rel_tol) { Ode::rel_tol=rel_tol; } 
      void SetTmax(const double tmax) { Ode::tmax=tmax; }
      void SetDt(const double dt) { Ode::dt=dt; }
      void SetModelParams(void *params) { Ode::params=params; }
      void SetFileId(const char *file_id) { strcpy(Ode::file_id,file_id); }
      void SetoFileName(const char *ofile_name) { strcpy(Ode::ofile_name,ofile_name); }
      void SeticFileName(const char *ic_file_name) { strcpy(Ode::ic_file_name,ic_file_name);} 
      void SetRhs(double *rhs_ic) { Ode::rhs=rhs_ic; }
      
      // accessors 
      size_t GetNvars() const { return nvars; }
      size_t GetNinit() const { return ninit; }
      unsigned int GetNsave() const { return nsave; }
      const char *GetStepAlgo() const { return step_algo; }
      double GetAtol() const { return abs_tol; }
      double GetRtol() const { return rel_tol; }
      double GetTmax() const { return tmax; }
      double GetDt() const { return dt; }
      void *GetModelParams() const { return params; }
      double *GetRhs() const { return rhs; }
      char *GetFileId()  {return file_id; }
      char *GetoFileName()  {return ofile_name; }
      char *GeticFileName() {return ic_file_name; }
      
      // gsl/odeiv solve stuff 
      void Solve();
      void PluginFuncs(int (*p2derivs)(double,const double *,double *,void *),
                       int (*p2jac)(double,const double *,double *,double *,void *));
      
      // i/o stuff 
      void PrtOdePrms();

   protected:

      double *rhs;                 // rhs[navrs] for the variables
      
   private:
      
      // gsl/odeiv solve stuff 
      const gsl_odeiv_step_type *SetStepType();
      void OpenoFile(ofstream& ofile);
      void CloseoFile(ofstream& ofile);
      void InitRhsFromFile();
      void WriteRHS(ofstream& ofile, double t);

      //this can be overridden by derived classes if needed 
      virtual void InitParameters() {;}

  
      // ode parameters 
      size_t nvars;                // size of the ode system
      size_t ninit;                // number of initial conditions
      unsigned int nsave;          // save solution every nsave time steps
      char step_algo[MAX_STR_LEN]; // name of stepping algorithm type
      double abs_tol, rel_tol;     // absolute/relative tolerances
      double tmax, dt;             // ...
      void *params;                // pointer to model parameters. will be set
                                   // by the user to some specific model
                                   // parameters using the member function
                                   // Ode::SetParams(void *params);
      
      // output file
      char file_id[MAX_STR_LEN];      // file identifier 
      char ofile_name[MAX_STR_LEN];   // output file name
      char ic_file_name[MAX_STR_LEN]; // ic input file name

      bool verbose;                // print verbose output
      
      // pointers to rhs and jacobian declaration 
      int (*p2derivs)(double, const double *, double *, void *);
      int (*p2jac)(double, const double *, double *, double *, void *);
      
}; // class Ode

#endif // ODE_H
