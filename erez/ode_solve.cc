/******************************************************************/

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <boost/program_options.hpp>

#include "ode.hh"
#include "model_ode.hh"
#include "ode_io_utils.hh"

using namespace std;
namespace po = boost::program_options;

/******************************************************************/

int main(int argc, char *argv[])
{
   // constants
   const size_t mf_nvars = 6;
   const size_t pa_nvars = 48;

   // suffixes
   const char mf_dat[] = ".mf.dat";
   const char pa_dat[] = ".pa.dat";
   
   const char mf_init[] = ".mf.init";
   const char pa_init[] = ".pa.init";

   // for file name manipulations
   char file_id[MAX_STR_LEN];
   char fname[MAX_STR_LEN];
   
   //MeanField mfa;
   ModelOde Model;
   
   // checking command line arguments
   po::options_description command_line_options
     ("Usage: ode_solve -p params_file [options]... \n\nAllowed options");

   command_line_options.add_options()
     ("help,h",
      "produce help message")
     ("params-file,p",po::value<std::string>(),
      "file containing model and ode parameters")
     ;

   po::options_description ode_options
     ("ODE parameters");

   ode_options.add_options()
     ("file_id", po::value<std::string>(),
      "file-id for output files")
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
     ("N", po::value<double>(),
      "total number of individuals")
     ("njac", po::value<double>(),
      "size of Jacobian (if needed)")
     ;

   po::options_description all_options;
   all_options.add(command_line_options).add(ode_options).add(model_options);
  
   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, all_options), vm);
   po::notify(vm);

   if (vm.count("help")) {
     std::cout << all_options << std::endl;
     return 1;
   }
  
   if (vm.count("params-file")) {
     std::ifstream ifs(vm["params-file"].as<std::string>().c_str());
     try {
       po::store(po::parse_config_file(ifs, all_options), vm);
     }
     catch (std::exception& e) {
       std::cout << "Error parsing config file: " << e.what() << std::endl;
       return 1;
     }
   }

   po::notify(vm);

   // reading ode parameters
   if (ReadOdeParams(vm, Model) != 0) {
     po::options_description options;
     options.add(command_line_options).add(ode_options);
     std::cout << options << std::endl;
     return 1;
   }
   
   // reading model parameters 
   if (ReadModelParams(vm, Model) != 0) {
     po::options_description options;
     options.add(command_line_options).add(model_options);
     std::cout << options << std::endl;
     return 1;
   }
    
   cout << endl
        << "Starting ODE solver\n"
        << "-------------------\n";
   
   // init stuff from Graph object //
   
   // init Qd, Qd, C...
   Model.SetQd(3);
   Model.SetQi(3);
   
   Model.SetCddi(0.0);
   Model.SetCddd(0.0);
   Model.SetCdii(0.0);
   Model.SetCdid(0.0);
   Model.SetCidi(0.0);
   Model.SetCidd(0.0);
   Model.SetCiii(0.0);
   Model.SetCiid(0.0);
   
   // SOLVING MFA //
   
   // init nvars 
   Model.SetNvars(mf_nvars);

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
   
   // printing ode + model parameters 
   Model.PrtOdePrms();
   Model.PrtModelPrms();
   Model.PrtGraphPrms();
   
   // solve the system    
   Model.Solve();

    // SOLVING PA //
   
   // init nvars 
   Model.SetNvars(pa_nvars);

   // init file names
   strcpy(fname, file_id);
   strcat(fname, pa_dat);
   Model.SetoFileName(fname);

   strcpy(fname, file_id);
   strcat(fname, pa_init);
   Model.SeticFileName(fname); 
   
   // plugin MF derivs   
   Model.PluginPAderivs();
   
   // printing ode + model parameters 
   //Model.PrtOdePrms();
   //Model.PrtModelPrms();
   //Model.PrtGraphPrms();
   
   // solve the system    
   Model.Solve();
   
   return 0;
}

   
   
