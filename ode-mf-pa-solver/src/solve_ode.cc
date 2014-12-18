#include <iostream>

#include <boost/program_options.hpp>

#include "ode_solver.hh"
//#include "InfoSIRSmf.hh"
#include "InfoSIRSpa.hh"
#include "ode_io_utils.hh"

//------------------------------------------------------------

int main(int argc, char* argv[])
{
  
  // generate main options
  po::options_description* main_options = generate_main_options();
  
  // generate ode options
  po::options_description* ode_options = generate_ode_options();

  // generate specific model options
  po::options_description* model_options = InfoSIRSpa::generate_model_options();
  
  // gathering options
  po::options_description all_options, all_ode_options;
  all_ode_options.add(*main_options).add(*ode_options);
  all_options.add(*main_options).add(*ode_options).add(*model_options);
  
  //------------------------------------------------------------
  // parsing command line and files
  
  // variables map for command line args
  po::variables_map vm;
  
  // parse command line arguments
  parse_comm_line_args(argc, argv, all_options, vm);
  
  // parse ode-params file
  parse_ode_args_file(all_ode_options, vm);
  
  // parse model-params file
  parse_model_args_file(*model_options, vm);
  
  //------------------------------------------------------------
  // help messages
  
  // print help messages
  if (vm.count("help"))
    {
      std::cout << *main_options << std::endl;
      return 1;
    }
  if (vm.count("longhelp"))
    {
      std::cout << all_options << std::endl;
      return 1;
    }
  
  //------------------------------------------------------------
  
  // define Model variable
  ode::OdeSolver<InfoSIRSpa::Params, InfoSIRSpa::Eqs> ode_system;
  
  // initializing ode parameters
  init_ode_params(vm, ode_system);
  
  // initializing model parameters
  InfoSIRSpa::init_model_params(vm, ode_system);
  
  // print ode_system
  std::cout << ode_system << std::endl;
  
  // init rhs
  ode_system.init_rhs();
  
  // FOR MF - RESCALE RATES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  // solve
  ode_system.solve();
  
  //------------------------------------------------------------
  // write parameters to Gnuplot script
  
  if (!vm.count("no-gp-file")) {
    write_gp_script(vm["file-id"].as<std::string>(),
                    vm["Qd"].as<double>(),
                    vm["Qi"].as<double>(),
                    vm["tmax"].as<double>(),
                    vm["vertices"].as<double>(),
                    vm.count("verbose"));
  }
  
  
  return 0;
}

//------------------------------------------------------------
