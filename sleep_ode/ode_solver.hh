#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

//------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <exception>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "file_io.hh"
#include "convergence.hh"

//------------------------------------------------------------

const size_t MAX_STR_LEN=128; 

namespace fio = file_io_utils;
namespace cv = convergence;

//------------------------------------------------------------

namespace ode
{
  //------------------------------------------------------------
  // declarations (needed for template friend function
  // see C++ FAQ lite [35.16]:
  // http://www.parashift.com/c++-faq-lite/templates.html
   
  template <class ModelParams, class ModelDerivs> class OdeSolver;
  template <class ModelParams, class ModelDerivs> std::ostream& operator<<
    (std::ostream& os, const OdeSolver<ModelParams, ModelDerivs>& x);
   
  //------------------------------------------------------------
   
  template <class ModelParams, class ModelDerivs> class OdeSolver
  {
  public:
         
    // constructors and destructors 
    OdeSolver();
    ~OdeSolver();
         
    // mutators 
    void set_nsave(const size_t nsave) { OdeSolver::nsave = nsave; }
    void set_step_algo(const char *step_algo)
    { strcpy(OdeSolver::step_algo, step_algo); }
    void set_abs_tol(const double abs_tol) { OdeSolver::abs_tol = abs_tol; }
    void set_rel_tol(const double rel_tol) { OdeSolver::rel_tol = rel_tol; } 
    void set_t0(const double t0) { OdeSolver::t0 = t0; }
    void set_tmax(const double tmax) { OdeSolver::tmax = tmax; }
    void set_dt(const double dt) { OdeSolver::dt = dt; }
    void set_convergence_check(const bool b) { OdeSolver::check_if_converged = b; }
    void set_file_id(const char *file_id)
    { std::strcpy(OdeSolver::file_id, file_id); }
    void set_output_file_name(const char *output_file_name)
    { std::strcpy(OdeSolver::output_file_name, output_file_name); }
    void set_ic_file_name(const char *ic_file_name)
    { std::strcpy(OdeSolver::ic_file_name, ic_file_name); } 
    void set_rhs(double *rhs_ic) { OdeSolver::rhs = rhs_ic; }
    void set_verbose(const bool verbose)
    { OdeSolver::verbose = verbose; }
         
    // accessors 
    size_t get_nsave() const { return nsave; }
    const char* get_step_algo() const { return step_algo; }
    double get_abs_tol() const { return abs_tol; }
    double get_rel_tol() const { return rel_tol; }
    double get_t0() const { return t0; }
    double get_tmax() const { return tmax; }
    double get_dt() const { return dt; }
    bool get_convergence_check() const { return check_if_converged; }
    double* get_rhs() const { return rhs; }
    const char* get_file_id() const {return file_id; }
    const char* get_output_file_name() const {return output_file_name; }
    const char* get_ic_file_name() const {return ic_file_name; }
    bool get_verbose() const { return verbose; }
    ModelParams* get_model_params() const { return model_params; }
    ModelDerivs get_model_derivs() const { return model_derivs; }
         
    // gsl/odeiv solve stuff 
    void solve();
         
    // overloding operators
    friend std::ostream& operator<< <ModelParams, ModelDerivs>
    (std::ostream& os, const OdeSolver<ModelParams, ModelDerivs>& x);

    template <class OtherDerivs>
    void operator= (const OdeSolver<ModelParams, OtherDerivs>& x);
         
  private:
    double *rhs;                 // rhs[navrs] for the variables
         
    // gsl/odeiv solve stuff 
    const gsl_odeiv_step_type* set_step_type();
    void init_rhs();
    void init_Q();
    void write_rhs(std::ofstream& output_file, double t);
    void write_last_line();
         
    // ode parameters 
    size_t nsave;          // save solution every nsave time steps
    char step_algo[MAX_STR_LEN]; // name of stepping algorithm type
    double abs_tol, rel_tol;     // absolute/relative tolerances
    double t0, tmax, dt;             // ...
         
    // output file
    char file_id[MAX_STR_LEN];      // file identifier 
         
    // CREATE FROM FILE_ID --- IN MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    char output_file_name[MAX_STR_LEN];   // output file name
    char ic_file_name[MAX_STR_LEN]; // ic input file name
         
    bool verbose;                // print verbose output
    bool check_if_converged;
         
    // pointers to rhs and jacobian declaration 
    //int (*p2derivs)(double, const double *, double *, void *);
    //int (*p2jac)(double, const double *, double *, double *, void *);
         
    // Model object
    ModelParams* model_params;
    ModelDerivs model_derivs;
         
  }; /* class OdeSolver */

  //------------------------------------------------------------
   
  // constructors and destructors 
   
  template <class ModelParams, class ModelDerivs>
  OdeSolver<ModelParams, ModelDerivs>::OdeSolver() : nsave(1), abs_tol(1e-6),
                                                     rel_tol(1e-6), t0(0),
                                                     tmax(10), dt(1e-6),
                                                     verbose(false),
                                                     check_if_converged(false)
  {
    std::strcpy(step_algo,"rkf45");
    std::strcpy(file_id,"no_file_id");
    std::strcpy(output_file_name,"no_ofile_name");
    std::strcpy(ic_file_name,"no_ic_file_name");
    rhs=NULL;
    model_params = new ModelParams;
      
  } // OdeSolver 
   
  //------------------------------------------------------------
   
  template <class ModelParams, class ModelDerivs>
  OdeSolver<ModelParams, ModelDerivs>::~OdeSolver()
  {
    delete model_params;
    delete [] rhs;
      
  } // ~OdeSolver 
   
  //------------------------------------------------------------
   
  template <class ModelParams, class ModelDerivs>
  void OdeSolver<ModelParams, ModelDerivs>::init_rhs()
  {
    // opening ic_file 
    fio::File<std::ifstream> ic(ic_file_name, verbose);
      
    double *rhs_ic;
    unsigned int nv = model_params->nvars;
      
    // allocating rhs 
    if (verbose) std::cout << "... allocating rhs memory";      
    try {      
      rhs_ic = new double[nv]; // will be deleted in destructor
    }      
    catch (std::exception &e) {
      std::cerr << "... unable to alloc rhs\n"
                << "... Standard exception: " << e.what() << std::endl;      
      std::exit(1); 
    }      
    if (verbose) std::cout << " ... done\n";
      
    // init rhs from ic_file 
    if (verbose) std::cout << "... reading content of ic-file";
      
    for(unsigned int i = 0; i < nv; i++)
      ic.fs >> rhs_ic[i];

    if (nv == 69) // dib model
      for (unsigned int i = 48; i < 69; i++) {
        rhs_ic[i - 21] -= rhs_ic[i];
        rhs_ic[i - 42] -= rhs_ic[i];
      }

    set_rhs(rhs_ic);
      
    if (verbose) std::cout << " ... done\n";
      
    // closing ic_file 
    ic.close(ic_file_name);
      
  } // init_rhs
   
  //------------------------------------------------------------
   
  template <class ModelParams, class ModelDerivs>
  void OdeSolver<ModelParams, ModelDerivs>::solve()
  {
    double t;
    int status;
    int nv = model_params->nvars;

    // convergence check object
    cv::ConvergenceCheck conv(nv);

    // print message    
    if (verbose) std::cout << std::endl 
                           << "\n----------------------\n"
                           << " FILE ID = " << file_id 
                           << "\n----------------------\n"
                           << "\nSolving the ode system\n"
                           << "----------------------\n";
    
    // opening output file 
    fio::File<std::ofstream> data(output_file_name, verbose);
      
    // setting step type 
    const gsl_odeiv_step_type *step_type;
    step_type = set_step_type();
      
    // allocating rhs and initialize it from ic_file 
    init_rhs();
    
    // printing
    if (verbose)
      std::cout << *this ;
      
    // allocating stepping function 
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(step_type, nv);
      
    // allocating control function 
    gsl_odeiv_control *control = gsl_odeiv_control_y_new(abs_tol,rel_tol);
      
    // allocating evolution function 
    gsl_odeiv_evolve *evolve  = gsl_odeiv_evolve_alloc(nv);

    // defining the system 
    gsl_odeiv_system sys =
      {&(model_derivs.rhs_eval), NULL, nv, static_cast<void*>(model_params) };
      
    // write t=t0 rhs
    t=t0;
    write_rhs(data.fs, t);
      
    size_t o_count = 1;

    //--- main loop ---//
    if (verbose) std::cout << "... doing main loop !!!\n";   
    
    while (t < tmax)
      {
        // stepping solution 
        status = gsl_odeiv_evolve_apply (evolve, control, step, &sys, &t, tmax,
                                         &dt, rhs);      
        
        if(status != GSL_SUCCESS)
          exit(1);
        
        // writing RHS to ofile 
        if(!(o_count % nsave))
          write_rhs(data.fs, t);

        // check convergence
        if (check_if_converged)
          if (!(o_count % conv.get_samples_interval()))
            if(conv.check(rhs)) {
              std::cout << "... Ode_Solver converged at t = "
                        << t << " !!!" << std::endl;
              break;
            }
        
        // advance counter
        ++o_count;        
        
      }
    
    if (verbose) std::cout << "... writing last line ";
    write_last_line();
    if (verbose) std::cout << "... done\n";
    
    if (verbose) std::cout << "... done ode\n";
    
    // closing output file
    data.close(output_file_name);
    
    // free allocated memory 
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    
  } // solve 
  
  //------------------------------------------------------------
   
  template <class ModelParams, class ModelDerivs>
  const gsl_odeiv_step_type* OdeSolver<ModelParams, ModelDerivs>::set_step_type()
  {
    if (verbose) std::cout << "... setting stepping algorithm to " << step_algo
                           << std::endl;
      
    if(!strcmp(step_algo,"rk2")) {
      return gsl_odeiv_step_rk2; }
    else if(!strcmp(step_algo,"rk4")) {         
      return gsl_odeiv_step_rk2; }
    else if(!strcmp(step_algo,"rkf45")) {
      return gsl_odeiv_step_rkf45; }
    else if(!strcmp(step_algo,"rkck")) {
      return gsl_odeiv_step_rkck; }
    else if(!strcmp(step_algo,"rk8pd")) {
      return gsl_odeiv_step_rk8pd; }
    else if(!strcmp(step_algo,"rk2imp")) {
      return gsl_odeiv_step_rk2imp; }
    else if(!strcmp(step_algo,"rk4imp")) {
      return gsl_odeiv_step_rk4imp; }
    else if(!strcmp(step_algo,"bsimp")) {
      return gsl_odeiv_step_bsimp; }
    else if(!strcmp(step_algo,"gear1")) {
      return gsl_odeiv_step_gear1; }
    else if(!strcmp(step_algo,"gear2")) {
      return gsl_odeiv_step_gear2; }
    else {
      if (verbose)
        std::cout << "... wrong step_algo value, set to default rkf45\n"; 

      return gsl_odeiv_step_rkf45;
    }
              
  } // get_step_type

  //------------------------------------------------------------
   
  template <class ModelParams, class ModelDerivs>
  void OdeSolver<ModelParams, ModelDerivs>::write_rhs(std::ofstream& ofile, double t)
  {
    int nv = model_params->nvars;
    
    ofile << t; 
    for(int i = 0; i < nv; i++)
      ofile << '\t' << rhs[i];
    ofile << std::endl;
      
  } // write_rhs 
   
  //------------------------------------------------------------

  template <class ModelParams, class ModelDerivs>
  void OdeSolver<ModelParams, ModelDerivs>::write_last_line()
  {
    int nv = model_params->nvars;
    std::string lastLineFileName(file_id);

    lastLineFileName += ".final";

    // open file_id.final
    std::ofstream lastLine(lastLineFileName.c_str(), std::ios::out);
    
    // writing to stringstream
    for(int i = 0; i < nv; i++)
      lastLine << rhs[i] << '\t';
    lastLine << std::endl;

    lastLine.close();
    
  } // write_last_line
  
  //------------------------------------------------------------
   
  // overloding ode::operator<<
  template <class ModelParams, class ModelDerivs>
  std::ostream& operator<< (std::ostream& output,
                            const OdeSolver<ModelParams, ModelDerivs>& x)   
  {
    output << std::endl
           << "Ode solver parameters:" << std::endl
           << "----------------------" << std::endl
           << "step type .......... " << x.get_step_algo() << std::endl
           << "abs tol ............ " << x.get_abs_tol() << std::endl
           << "rel tol ............ " << x.get_rel_tol() << std::endl
           << "t0.. ............... " << x.get_t0() << std::endl
           << "tmax ............... " << x.get_tmax() << std::endl
           << "dt ................. " << x.get_dt() << std::endl
           << "nsave .............. " << x.get_nsave() << std::endl
           << "file id ............ " << x.get_file_id() << std::endl
           << "output file ........ " << x.get_output_file_name() << std::endl
           << "ic file ............ " << x.get_ic_file_name() << std::endl
           << "verbose mode ....... " << x.get_verbose() << std::endl
           << "convergence check .. " << x.get_convergence_check() << std::endl
           << *(x.model_params);
            
    return output;
      
  } /* operator<< */
   
  //------------------------------------------------------------
   
  // overloding ode::operator=
  template <class ModelParams, class ModelDerivs>
  template <class OtherDerivs>
  void OdeSolver<ModelParams, ModelDerivs>::operator=
  (const OdeSolver<ModelParams, OtherDerivs>& x)
  {

    // ode parameters
    nsave   = x.get_nsave();
    abs_tol = x.get_abs_tol();
    rel_tol = x.get_rel_tol();
    t0      = x.get_t0();
    tmax    = x.get_tmax();
    dt      = x.get_dt();
    verbose = x.get_verbose();
    check_if_converged = x.get_convergence_check();
      
    strcpy(step_algo, x.get_step_algo());
    strcpy(file_id, x.get_file_id());
    strcpy(ic_file_name, x.get_ic_file_name());
      
    // model parameters
      
    *model_params = *(x.get_model_params());
      
  } /* operator= */
   
  //------------------------------------------------------------
   
} /* namespace ode */

//------------------------------------------------------------

#endif // ODE_SOLVER_H

