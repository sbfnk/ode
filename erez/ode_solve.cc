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

   po::variables_map vm;
   po::store(po::parse_command_line(argc, argv, command_line_options), vm);
   po::notify(vm);

   if (vm.count("help")) {
     std::cout << command_line_options << std::endl;
     return 1;
   }

   std::string paramsFile = "";
   
   if (vm.count("params-file")) {
     paramsFile = vm["params-file"].as<std::string>();
   }

   // reading ode parameters 
   ReadOdeParams(argc, argv, paramsFile.c_str(), Model);
   
   // reading model parameters 
   ReadModelParams(argc, argv, paramsFile.c_str(), Model);
   
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

   
   
