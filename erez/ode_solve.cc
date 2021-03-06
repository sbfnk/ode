/******************************************************************/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <boost/program_options.hpp>

#include "ode.hh"
#include "model_ode.hh"
#include "ode_io_utils.hh"

//#define DEBUG

using namespace std;
namespace po = boost::program_options;

/******************************************************************/

int main(int argc, char *argv[])
{
   // variables
   const size_t mf_nvars = 6;
   const size_t pa_nvars = 48;

   // initial conditions
   const size_t mf_ninit = 48;
   const size_t pa_ninit = 49;

   // suffixes
   const char mf_dat[] = ".mf.dat";
   const char pa_dat[] = ".pa.dat";
   
// use common init file for mean field and pair approximation
   const char mf_init[] = ".init";
   const char pa_init[] = ".init";
//    const char mf_init[] = ".mf.init";
//    const char pa_init[] = ".pa.init";

   const char gp[] = ".gp";
   const char params[] = ".prm.txt";

   // for file name manipulations
   char file_id[MAX_STR_LEN];
   char fname[MAX_STR_LEN];
   
   bool verbose = false;
   
   // checking command line arguments
   po::options_description command_line_options
      ("Usage: ode_solve -p params_file -m model_file [options]... \n\nAllowed options");
   
   command_line_options.add_options()
      ("help,h", "produce help message")
      ("verbose,v", "produce verbose output")
      ("params-file,p",po::value<std::string>(),
       "file containing ode parameters")
      ("model-file,m",po::value<std::string>(),
       "file containing model parameters")
      ;
   
   po::options_description ode_options
      ("ODE parameters");
   
   ode_options.add_options()
      ("file_id", po::value<std::string>(),
       "file id for output files")
      ("tmax", po::value<double>(),
       "stopping time")
      ("dt", po::value<double>(),
       "size of time step")
      ("nsave", po::value<unsigned int>(),
       "save solution every nsave time steps")
      ("step_algo", po::value<std::string>(),
       "name of stepping algorithm type")
      ("abs_tol", po::value<double>(),
       "absolute tolerance")
      ("rel_tol", po::value<double>(),
       "relative tolerance")
      ;
   
   po::options_description model_options
      ("Model parameters");
   
   model_options.add_options()
      ("beta--", po::value<double>(),
       "disease transmission rate uninformed->uninformed")
      ("beta+-", po::value<double>(),
       "disease transmission rate informed->uninformed")
      ("beta-+", po::value<double>(),
       "disease transmission rate uninformed->informed")
      ("beta++", po::value<double>(),
       "disease transmission rate informed->informed")
      ("gamma-", po::value<double>(),
       "recovery rate of uninformed")
      ("gamma+", po::value<double>(),
       "recovery rate of informed")
      ("delta-", po::value<double>(),
       "loss of immunity rate of uninformed")
      ("delta+", po::value<double>(),
       "loss of immunity rate of informed")
      ("alpha", po::value<double>(),
       "information transmission rate")
      ("nu", po::value<double>(),
       "information generation rate")
      ("omega", po::value<double>(),
       "local information generation rate")
      ("lambda", po::value<double>(),
       "loss of information rate")
      ("vertices,N", po::value<double>(),
       "total number of individuals")
      ("njac", po::value<double>(),
       "size of Jacobian (if needed)")
      ;
   
   po::options_description all_options;
   all_options.add(command_line_options).add(ode_options).add(model_options);
   
   po::variables_map vm;

   // parse command line options and allow unregistered options
   po::parsed_options parsed = po::command_line_parser(argc, argv).
     options(all_options).allow_unregistered().run();

   // remove all unregistered options
   for (unsigned i = 0; i < parsed.options.size(); i++) {
     if (parsed.options[i].unregistered) {
       parsed.options.erase(parsed.options.begin()+i);
       --i;
     }
   }
   
   po::store(parsed, vm);
   po::notify(vm);
   
   if (vm.count("help"))
   {
      std::cout << all_options << std::endl;
      return 1;
   }
   
   if (vm.count("verbose"))
   {
     verbose = true;
   }
   
   if (vm.count("params-file"))
   {
      std::ifstream ifs(vm["params-file"].as<std::string>().c_str());
      try
      {
         po::store(po::parse_config_file(ifs, all_options), vm);
      }
      catch (std::exception& e)
      {
         std::cout << "Error parsing params file: " << e.what() << std::endl;
         return 1;
      }
   }
   
   if (vm.count("model-file"))
   {
      std::ifstream ifs(vm["model-file"].as<std::string>().c_str());
      try
      {
         po::store(po::parse_config_file(ifs, all_options), vm);
      }
      catch (std::exception& e)
      {
         std::cout << "Error parsing model file: " << e.what() << std::endl;
         return 1;
      }
   }
   
   po::notify(vm);
   
   //MeanField mfa;
   ModelOde Model(verbose);

   // reading ode parameters
   if (ReadOdeParams(vm, Model) != 0)
   {
      po::options_description options;
      options.add(command_line_options).add(ode_options);
      std::cout << command_line_options << ode_options << std::endl;
      return 1;
   }
   
   // reading model parameters 
   if (ReadModelParams(vm, Model) != 0)
   {
      std::cout << command_line_options << model_options << std::endl;
      return 1;
   }
   
   if (verbose) {
     cout << endl
          << "Starting ODE solver\n"
          << "-------------------\n";
   }
   
   // SOLVING MFA //
   
   // init nvars 
   Model.SetNvars(mf_nvars);
   Model.SetNinit(mf_ninit);

   // init file names
   strcpy(file_id, Model.GetFileId());
   
   strcpy(fname, file_id);
   strcat(fname, mf_dat);
   Model.SetoFileName(fname);

   strcpy(fname, file_id);
   strcat(fname, mf_init);
   Model.SeticFileName(fname); 
   
   // plugin MF derivs   
   Model.PluginMFderivs();
   
   // solve the system    
   Model.Solve();
   
   if (verbose) {
     // printing ode + model parameters 
     Model.PrtOdePrms();
     Model.PrtModelPrms();
     Model.PrtGraphPrms();
   }
   
   if (verbose) std::cout << std::endl;
   
   // SOLVING PA //
   
   // init nvars 
   Model.SetNvars(pa_nvars);
   Model.SetNinit(pa_ninit);
   
   // init file names
   strcpy(fname, file_id);
   strcat(fname, pa_dat);
   Model.SetoFileName(fname);
   
   strcpy(fname, file_id);
   strcat(fname, pa_init);
   Model.SeticFileName(fname); 
   
   // plugin MF derivs   
   Model.PluginPAderivs();
   
   // solve the system    
   Model.Solve();

   // printing ode + model parameters 
//    if (verbose) {
//      Model.PrtOdePrms();
//      Model.PrtModelPrms();
//      Model.PrtGraphPrms();
//    }
   
   // write parameters for gnuplot
   strcpy(fname, file_id);
   strcat(fname, gp);
   WriteGnuPlot(fname, Model, verbose);

   // write model parameters to file
   strcpy(fname, file_id);
   strcat(fname, params);
   WriteModelParams(fname, Model, verbose);
   
   if (verbose) std::cout << std::endl;
   
   return 0;
}

   
   
