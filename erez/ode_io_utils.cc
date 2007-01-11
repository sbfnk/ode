/******************************************************************/  

#include "ode_io_utils.hh"
#include "ode.hh"
#include "model_ode.hh"

#include <boost/program_options.hpp>

using namespace std;

namespace po = boost::program_options;

/******************************************************************/  

int ReadOdeParams(int argc, char* argv[], const char* params_file, Ode& obj)
{

  po::options_description ode_options
    ("ODE parameters");

  ode_options.add_options()
    ("file-id", po::value<std::string>(),
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

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, ode_options), vm);
  std::ifstream ifs(params_file);
  try {
    po::store(po::parse_config_file(ifs, ode_options), vm);
  }
  catch (std::exception& e) {
    std::cout << "Error parsing config file: " << e.what() << std::endl;
    return 1;
  }
  
  po::notify(vm);

  if (vm.count("file-id")) {
    // ode file id name 
    obj.SetFileId(vm["file-id"].as<char*>());
  } else {
    std::cerr << "ERROR: no file-id given" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }
  
  if (vm.count("tmax")) {
    // tmax 
    obj.SetTmax(vm["tmax"].as<double>());
  } else {
    std::cerr << "ERROR: no tmax give" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }
  if (vm.count("dt")) {
    // dt 
    obj.SetDt(vm["dt"].as<double>());
  } else {
    std::cerr << "ERROR: no dt given" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }
  if (vm.count("nsave")) {
    // nsave 
    obj.SetNsave(vm["nsave"].as<unsigned int>());
  } else {
    std::cerr << "ERROR: no nsave given" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }
  if (vm.count("step_algo")) {
    // step_algo 
    obj.SetStepAlgo(vm["step_algo"].as<char*>());
  } else {
    std::cerr << "ERROR: no step_algo given" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }
  if (vm.count("abs_tol")) {
    // abs_tol 
    obj.SetAbsTol(vm["abs_tol"].as<double>());
  } else {
    std::cerr << "ERROR: no abs_tol given" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }
  if (vm.count("rel_tol")) {
    // rel_tol 
    obj.SetRelTol(vm["rel_tol"].as<double>());
  } else {
    std::cerr << "ERROR: no rel_tol given" << std::endl;
    std::cerr << std::endl;
    std::cerr << ode_options << std::endl;
    return 1;
  }

  return 0;
}

/******************************************************************/  

int ReadModelParams(int argc, char* argv[], const char* params_file,
                     ModelOde& obj)
{
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
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, model_options), vm);
  std::ifstream ifs(params_file);
  try {
    po::store(po::parse_config_file(ifs, model_options), vm);
  }
  catch (std::exception& e) {
    std::cout << "Error parsing config file: " << e.what() << std::endl;
    return 1;
  }
  
  po::notify(vm);

  ModelParams *model_params = new ModelParams;
  
//    // beta_d    
//    model_params->beta_d=read_dbl_val(ifile);
  
//    // gamma_d 
//    model_params->gamma_d=read_dbl_val(ifile);
  
//    // delta_d 
//    model_params->delta_d=read_dbl_val(ifile);
      
//    // beta_i 
//    model_params->beta_i=read_dbl_val(ifile);
   
//    // gamma_i 
//    model_params->gamma_i=read_dbl_val(ifile);
   
//    // delta_i 
//    model_params->delta_i=read_dbl_val(ifile);
   
//    // beta_m1 
//    model_params->beta_m1=read_dbl_val(ifile);
   
//    // beta_m2 
//    model_params->beta_m2=read_dbl_val(ifile);
   
//    // alpha 
//    model_params->alpha=read_dbl_val(ifile);
   
//    // nu 
//    model_params->nu=read_dbl_val(ifile);
   
//    // lambda 
//    model_params->lambda=read_dbl_val(ifile);
   
//    // Qd 
//    //model_params->Qd=read_dbl_val(ifile);
   
//    // Qi 
//    //model_params->Qi=read_dbl_val(ifile);
   
//    // N
//    model_params->N=read_dbl_val(ifile);
   
//    // njac 
//    model_params->njac=obj.GetNvars();

  return 0;
  
}

