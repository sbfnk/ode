#include <iostream>

#include <boost/program_options.hpp>

#include "ode_solver.hh"
#include "Sleep.hh"
#include "ode_io_utils.hh"

//------------------------------------------------------------

int main(int argc, char* argv[])
{
  // some constants
  const std::string suffix = ".dat" ;
  
  // ode solver
  ode::OdeSolver<SleepParams, Sleep> solver;
  
  // variables map for command line args
  po::variables_map vm;
  
  // reading command line args
  if (read_comm_line_args(argc, argv, vm, solver)) // something is wrong - return 1
    return 1;
  
  // initializing ode parameters
  if (init_ode_params_from_comm_line(vm, solver)) // something is wrong - return 1
    return 1;
  
  // initializing model parameters
  init_model_params_from_comm_line(vm, solver);    

  // initialize output and nvars files
  std::string str;

  str = vm["file-id"].as<std::string>() + suffix;
  solver.set_output_file_name(str.c_str());
  
  solver.solve();
  
  return 0;
}
   
//------------------------------------------------------------
