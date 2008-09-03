#ifndef ODE_IO_UTILS_H
#define ODE_IO_UTILS_H

//------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <string>
#include <exception>
#include <boost/program_options.hpp>

#include "ode_solver.hh"
#include "file_io.hh"

//------------------------------------------------------------

namespace po = boost::program_options;
namespace fio = file_io_utils;

//------------------------------------------------------------

int read_comm_line_args(int argc, char* argv[], po::variables_map& vm)
{
   
  // checking command line arguments
  po::options_description command_line_options
    ("Usage: ode_solve -p ode_params_file -m model_params_file [options]...\n\nAllowed options");
   
  // general command line options
  command_line_options.add_options()
    ("help,h", "produce help message")
    ("verbose,v", "produce verbose output")
    ("ode-params-file,p",po::value<std::string>(),
     "file containing ode parameters")
    ("model-params-file,m",po::value<std::string>(),
     "file containing model parameters")
    ("check-convergence",
     "check for convergence and stop before Tmax")

    ;
   
  // ode command line options
  po::options_description ode_options
    ("ODE parameters");
   
  ode_options.add_options()
    ("file-id", po::value<std::string>(),
     "file id for output files")
    ("ic-file", po::value<std::string>(),
     "initial conditions file")
    ("tmax", po::value<double>(),
     "stopping time")
    ("dt", po::value<double>(),
     "size of time step")
    ("nsave", po::value<unsigned int>(),
     "save solution every nsave time steps")
    ("step-algo", po::value<std::string>(),
     "name of stepping algorithm type")
    ("abs-tol", po::value<double>(),
     "absolute tolerance")
    ("rel-tol", po::value<double>(),
     "relative tolerance")
    ;

  // model command line options
  po::options_description model_options
    ("Model parameters");
   
  model_options.add_options()
    ("alpha", po::value<double>(),
     "information transmission rate")
    ("beta", po::value<double>(),
     "disease transmission rate")
    ("gamma", po::value<double>(),
     "rate of transition stage 1->2")
    ("delta", po::value<double>(),
     "loss of immunity rate")
    ("lambda", po::value<double>(),
     "rate of forgetting")
    ("omega", po::value<double>(),
     "rate of information generation")
    ("rho", po::value<double>(),
     "fractional loss of information at transmission of forgetting")
    ("N", po::value<unsigned int>(),
     "size of the population")
    ("M", po::value<unsigned int>(),
     "number of generations to consider")
    ;
  
  // all options
  po::options_description all_options;
  all_options.add(command_line_options).add(ode_options).add(model_options);
  
  //po::variables_map vm;
  
  // parse command line options and allow unregistered options
  po::parsed_options parsed = po::command_line_parser(argc, argv).
    options(all_options).allow_unregistered().run();
  
  // remove all unregistered options
  unsigned int i;
  for (i = 0; i < parsed.options.size(); i++)
    {
      if (parsed.options[i].unregistered)
        {
          parsed.options.erase(parsed.options.begin()+i);
          --i;
        }
    }
  
  // store and notify
  po::store(parsed, vm);
  po::notify(vm);
  
  // print help messages
  if (vm.count("help"))
    {
      std::cout << all_options << std::endl;
      return 1;
    }
  
  // parsing ode params file
  if (vm.count("ode-params-file"))
    {
      std::ifstream ifs(vm["ode-params-file"].as<std::string>().c_str());
      try
        {
          po::store(po::parse_config_file(ifs, all_options), vm);
        }
      catch (std::exception& e)
        {
          std::cout << "Error parsing ode params file: " << e.what() << std::endl;
          return 1;
        }
      ifs.close();
    }
  
  // parsing model params file
  if (vm.count("model-params-file"))
    {
      std::ifstream ifs(vm["model-params-file"].as<std::string>().c_str());
      try
        {
          po::store(po::parse_config_file(ifs, all_options), vm);
        }
      catch (std::exception& e)
        {
          std::cout << "Error parsing model params file: " << e.what() << std::endl;
          return 1;
        }
      ifs.close();
    }
  
  // notify
  po::notify(vm);
  
  return 0;
}

//------------------------------------------------------------

template <class ModelParams, class ModelDerivs>
int init_ode_params_from_comm_line(po::variables_map& vm,
                                   ode::OdeSolver<ModelParams, ModelDerivs>& x)
{
  if (vm.count("file-id")) {
    // ode file id name 
    x.set_file_id(vm["file-id"].as<std::string>().c_str());      
  } else  {
    std::cerr << "ERROR: no file-id given" << std::endl;
    return 1;
  }
  if (vm.count("ic-file")) {
    // ode ic file name 
    x.set_ic_file_name(vm["ic-file"].as<std::string>().c_str());      
  } else  {
    std::cerr << "ERROR: no ic-file given" << std::endl;
    return 1;
  }   
  if (vm.count("tmax")) {
    // tmax 
    x.set_tmax(vm["tmax"].as<double>());
  } else {
    std::cerr << "ERROR: no tmax given" << std::endl;
    return 1;
  }
  if (vm.count("dt")) {
    // dt 
    x.set_dt(vm["dt"].as<double>());
  } else {
    std::cerr << "ERROR: no dt given" << std::endl;
    return 1;
  }
  if (vm.count("nsave")) {
    // nsave 
    x.set_nsave(vm["nsave"].as<unsigned int>());
  } else {
    std::cerr << "ERROR: no nsave given" << std::endl;
    return 1;
  }
  if (vm.count("step-algo")) {
    // step_algo 
    x.set_step_algo(vm["step-algo"].as<std::string>().c_str());
  } else {
    std::cerr << "ERROR: no step-algo given" << std::endl;
    return 1;
  }
  if (vm.count("abs-tol")) {
    // abs_tol 
    x.set_abs_tol(vm["abs-tol"].as<double>());
  } else {
    std::cerr << "ERROR: no abs-tol given" << std::endl;
    return 1;
  }
  if (vm.count("rel-tol")) {
    // rel_tol 
    x.set_rel_tol(vm["rel-tol"].as<double>());
  } else {
    std::cerr << "ERROR: no rel-tol given" << std::endl;
    return 1;
  }
  if (vm.count("verbose")) {
    // verbose 
    x.set_verbose(true);
  }
  if (vm.count("check-convergence")) {
    // convergence
    x.set_convergence_check(true);
  }
  
  return 0;
}

//------------------------------------------------------------

template <class ModelParams, class ModelDerivs>
int init_model_params_from_comm_line(po::variables_map& vm,
                                     ode::OdeSolver<ModelParams, ModelDerivs>& x)
{
  ModelParams* model_params = x.get_model_params(); 
  
  if (vm.count("alpha")) {
    model_params->alpha=vm["alpha"].as<double>();
  } else {
    std::cerr << "WARNING: no alpha given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->alpha=0;
  }
  if (vm.count("beta")) {
    model_params->beta=vm["beta"].as<double>();
  } else {
    std::cerr << "WARNING: no beta given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->beta=0;
  }
  if (vm.count("gamma")) {
    model_params->gamma=vm["gamma"].as<double>();
  } else {
    std::cerr << "WARNING: no gamma given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->gamma=0;
  }
  if (vm.count("delta")) {
    model_params->delta=vm["delta"].as<double>();
  } else {
    std::cerr << "WARNING: no delta given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->delta=0;
  }
  if (vm.count("lambda")) {
    model_params->lambda=vm["lambda"].as<double>();
  } else {
    std::cerr << "WARNING: no lambda given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->lambda=0;
  }
  if (vm.count("omega")) {
    model_params->omega=vm["omega"].as<double>();
  } else {
    std::cerr << "WARNING: no omega given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->omega=0;
  }
  if (vm.count("rho")) {
    model_params->rho=vm["rho"].as<double>();
  } else {
    std::cerr << "WARNING: no rho given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->rho=0;
  }
  if (vm.count("M")) {
    model_params->M=vm["M"].as<unsigned int>();
  } else {
    std::cerr << "WARNING: no M given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->M=0;
  }
  if (vm.count("N")) {
    model_params->N=vm["N"].as<unsigned int>();
  } else {
    std::cerr << "WARNING: no N given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->N=0;
  }
  model_params->nvars=(model_params->M+1) * 3;

  return 0;
}

#endif // ODE_IO_UTILS_H
