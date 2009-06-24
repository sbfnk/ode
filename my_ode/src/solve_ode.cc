#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "ode_solver.hh"
#include "InfoSIRSmf.hh"
#include "InfoSIRSpa.hh"
#include "ode_io_utils.hh"

namespace mf = InfoSIRSmf;
namespace pa = InfoSIRSpa;

// for overloading generate_..._options functions
using namespace mf;
using namespace pa;

//------------------------------------------------------------

template <class T>
int solve_ode_system(int argc, char* argv[]);

//------------------------------------------------------------

int main(int argc, char* argv[])
{
  int status;
  
  // parse command line looking just for --ode-model=mf/pa  
  std::string ode_type = find_ode_type(argc, argv);
  
  if ((ode_type != "mf") && (ode_type != "pa")) {
    std::cerr << "ERROR: ode-type is neither mf nor pa" << std::endl;
    return 1;
  }
  
  // define Model variable
  if (ode_type == "mf") {
    status = solve_ode_system<ode::OdeSolver<mf::Params, mf::Eqs> >(argc, argv);
  } else {    
    status = solve_ode_system<ode::OdeSolver<pa::Params, pa::Eqs> >(argc, argv);
  }
  
  return status;
  
}

//------------------------------------------------------------

template <class T>
int solve_ode_system(int argc, char* argv[])
{
  // ode system variable
  T ode_sys;
  
  // generate main options
  po::options_description* main_options = generate_main_options();
  
  // generate ode options
  po::options_description* ode_options = generate_ode_options();
  
  // generate specific model options
  po::options_description* model_options = generate_model_options(ode_sys);
  
  // generate graph options
  po::options_description* graph_options = generate_graph_options(ode_sys);
  
  // generate graph options
  po::options_description* hidden_options = generate_hidden_options(ode_sys);
  
  // gathering options
  po::options_description all_options, all_ode_options, all_graph_options;
  
  all_ode_options.add(*main_options).add(*ode_options);
  all_options.add(*main_options).add(*ode_options).add(*model_options).add(*graph_options).add(*hidden_options);
  all_graph_options.add(*graph_options).add(*hidden_options);
  
  po::options_description help_options, long_help_options;
  help_options.add(*main_options);
  long_help_options.add(*main_options).add(*ode_options).add(*model_options).add(*graph_options);
  
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
  
  // parse graph-params file
  parse_graph_args_file(all_graph_options, vm);
  
  //------------------------------------------------------------
  // help messages
  
  // print help messages
  if (vm.count("help"))
    {
      std::cout << help_options << std::endl;
      return 1;
    }
  if (vm.count("longhelp"))
    {
      std::cout << long_help_options << std::endl;
      return 1;
    }

  //------------------------------------------------------------
  
  // initializing ode parameters
  init_ode_params(vm, ode_sys);
  
  // initializing model parameters
  init_model_params(vm, ode_sys);
  
  // initializing graph parameters
  init_graph_params(vm, ode_sys);
  
  // print ode_system
  std::cout << "print ode system" << std::endl;
  std::cout << ode_sys << std::endl;
  
  // init rhs
  ode_sys.init_rhs();
  
  // solve
  ode_sys.solve();
  
  //------------------------------------------------------------
  // write parameters to Gnuplot script
  
  if (vm.count("gp-file")) {
    write_gp_script(vm["file-id"].as<std::string>(),
                    vm["Qd"].as<double>(),
                    vm["Qi"].as<double>(),
                    vm["tmax"].as<double>(),
                    vm["vertices"].as<double>(),
                    vm.count("verbose"));
  }
  
  // cleaning
  delete main_options;
  delete ode_options;
  delete model_options;
  delete graph_options;    

  return 0;
}

//------------------------------------------------------------
