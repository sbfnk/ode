/******************************************************************/  

#include "ode_io_utils.hh"
#include "ode.hh"
#include "model_ode.hh"

using namespace std;

/******************************************************************/  

int ReadOdeParams(po::variables_map& vm, Ode& obj)
{
  if (vm.count("file_id")) {
    // ode file id name 
    obj.SetFileId(vm["file_id"].as<std::string>().c_str());
  } else {
    std::cerr << "ERROR: no file_id given" << std::endl;
    return 1;
  }
  
  if (vm.count("tmax")) {
    // tmax 
    obj.SetTmax(vm["tmax"].as<double>());
  } else {
    std::cerr << "ERROR: no tmax give" << std::endl;
    return 1;
  }
  if (vm.count("dt")) {
    // dt 
    obj.SetDt(vm["dt"].as<double>());
  } else {
    std::cerr << "ERROR: no dt given" << std::endl;
    return 1;
  }
  if (vm.count("nsave")) {
    // nsave 
    obj.SetNsave(vm["nsave"].as<unsigned int>());
  } else {
    std::cerr << "ERROR: no nsave given" << std::endl;
    return 1;
  }
  if (vm.count("step_algo")) {
    // step_algo 
    obj.SetStepAlgo(vm["step_algo"].as<std::string>().c_str());
  } else {
    std::cerr << "ERROR: no step_algo given" << std::endl;
    return 1;
  }
  if (vm.count("abs_tol")) {
    // abs_tol 
    obj.SetAbsTol(vm["abs_tol"].as<double>());
  } else {
    std::cerr << "ERROR: no abs_tol given" << std::endl;
    return 1;
  }
  if (vm.count("rel_tol")) {
    // rel_tol 
    obj.SetRelTol(vm["rel_tol"].as<double>());
  } else {
    std::cerr << "ERROR: no rel_tol given" << std::endl;
    return 1;
  }

  return 0;
}

/******************************************************************/  

int ReadModelParams(po::variables_map& vm, ModelOde& obj)
{
  ModelParams *model_params = new ModelParams;
  
  if (vm.count("beta--")) {
    model_params->beta[0][0]=vm["beta--"].as<double>();
  } else {
    std::cerr << "WARNING: no beta-- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->beta[0][0]=0;
  }
  if (vm.count("beta+-")) {
    model_params->beta[1][0]=vm["beta+-"].as<double>();
  } else {
    std::cerr << "WARNING: no beta+- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->beta[1][0]=0;
  }
  if (vm.count("beta-+")) {
    model_params->beta[0][1]=vm["beta-+"].as<double>();
  } else {
    std::cerr << "WARNING: no beta-+ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->beta[0][1]=0;
  }
  if (vm.count("beta++")) {
    model_params->beta[1][1]=vm["beta++"].as<double>();
  } else {
    std::cerr << "WARNING: no beta++ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->beta[1][1]=0;
  }
  if (vm.count("gamma-")) {
    model_params->gamma[0]=vm["gamma-"].as<double>();
  } else {
    std::cerr << "WARNING: no gamma- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->gamma[0]=0;
  }
  if (vm.count("gamma+")) {
    model_params->gamma[1]=vm["gamma+"].as<double>();
  } else {
    std::cerr << "WARNING: no gamma+ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->gamma[1]=0;
  }
  if (vm.count("delta-")) {
    model_params->delta[0]=vm["delta-"].as<double>();
  } else {
    std::cerr << "WARNING: no delta- given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->delta[0]=0;
  }
  if (vm.count("delta+")) {
    model_params->delta[1]=vm["delta+"].as<double>();
  } else {
    std::cerr << "WARNING: no delta+ given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->delta[1]=0;
  }
  if (vm.count("alpha")) {
    model_params->alpha=vm["alpha"].as<double>();
  } else {
    std::cerr << "WARNING: no alpha given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->alpha=0;
  }
  if (vm.count("nu")) {
    model_params->nu=vm["nu"].as<double>();
  } else {
    std::cerr << "WARNING: no nu given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->nu=0;
  }
  if (vm.count("omega")) {
    model_params->omega=vm["omega"].as<double>();
  } else {
    std::cerr << "WARNING: no omega given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->omega=0;
  }
  if (vm.count("lambda")) {
    model_params->lambda=vm["lambda"].as<double>();
  } else {
    std::cerr << "WARNING: no lambda given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    model_params->lambda=0;
  }
  if (vm.count("N")) {
    model_params->N=vm["N"].as<double>();
  } else {
    std::cerr << "ERROR: no N given" << std::endl;
    return 1;
  }
  if (vm.count("njac")) {
    model_params->njac=vm["njac"].as<size_t>();
  }

  // passing a pointer to model_params to Ode class member _params 
  obj.SetModelParams(static_cast<void *>(model_params));
  
  return 0;
}

