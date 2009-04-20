#ifndef SLEEP_HH
#define SLEEP_HH

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>
#include <boost/algorithm/string.hpp>

//------------------------------------------------------------

namespace po = boost::program_options;

//------------------------------------------------------------

struct SleepParams
{
  SleepParams() : nvars(3) {};
  
  // No. of equations
  unsigned int nvars;
  
  double beta, gamma, delta;
  double rho, epsilon, sigma;

  unsigned int N;

  std::vector<std::pair<double, double> > sigma_vector;
  unsigned int next_sigma;
  
  // overloading operator<<
  friend std::ostream& operator <<
    (std::ostream& os, const SleepParams& x);
      
  // overloading operator=
  void operator= (const SleepParams& x);

  po::options_description get_command_line_params();
  void init_model_params(po::variables_map& vm);

  bool read_sigma_file(std::string fileName);
      
}; // SleepParams

//------------------------------------------------------------

void SleepParams::operator= (const SleepParams& x)
{
  nvars = x.nvars;
  
  beta  = x.beta;
  gamma = x.gamma;
  delta = x.delta;
//   sigma = x.sigma;
//   tau = x.tau;
  
} // operator=

//------------------------------------------------------------

po::options_description SleepParams::get_command_line_params()
{
  po::options_description model_options
    ("Model parameters");

  model_options.add_options()
    ("beta", po::value<double>(),
     "disease transmission rate")
    ("gamma", po::value<double>(),
     "death rate of infected")
    ("rho", po::value<double>(),
     "rate of recognition in the absence of surveillance")
    ("epsilon", po::value<double>()->default_value(1.),
     "efficiency of screening")
    ("sigma", po::value<double>(),
     "rate of screening")
    ("delta", po::value<double>(),
     "healing rate for recognized infected")
    ("N", po::value<unsigned int>(),
     "size of the population")
    ("sigma-file", po::value<std::string>(),
     "file containing screening rates for different years")
    ;

  return model_options;
}

//------------------------------------------------------------

void SleepParams::init_model_params(po::variables_map& vm)
{
  if (vm.count("beta")) {
    beta=vm["beta"].as<double>();
  } else {
    std::cerr << "WARNING: no beta given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    beta=0;
  }
  if (vm.count("gamma")) {
    gamma=vm["gamma"].as<double>();
  } else {
    std::cerr << "WARNING: no gamma given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    gamma=0;
  }
  if (vm.count("rho")) {
    rho=vm["rho"].as<double>();
  } else {
    std::cerr << "WARNING: no rho given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    rho=0;
  }
  if (vm.count("epsilon")) {
    epsilon=vm["epsilon"].as<double>();
  } else {
    std::cerr << "WARNING: no epsilon given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    epsilon=0;
  }
  if (vm.count("sigma")) {
    sigma=vm["sigma"].as<double>();
    if (vm.count("sigma-file")) {
      std::cerr << "WARNING: sigma and sigma-file given" << std::endl;
      std::cerr << "using only sigma of " << sigma << std::endl;
    }
  } else {
    if (vm.count("sigma-file")) {
      std::string fileName = vm["sigma-file"].as<std::string>();
      if (read_sigma_file(fileName) == false) {
        std::cerr << "WARNING: problem reading" << fileName << std::endl;
        std::cerr << "setting sigma to 0" << std::endl;
      }
    } else {
      std::cerr << "WARNING: no sigma given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      sigma=0;
    }
  }
  if (vm.count("delta")) {
    delta=vm["delta"].as<double>();
  } else {
    std::cerr << "WARNING: no delta given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    delta=0;
  }
  if (vm.count("N")) {
    N=vm["N"].as<unsigned int>();
  } else {
    std::cerr << "WARNING: no N given" << std::endl;
    std::cerr << "setting to 0" << std::endl;
    N=0;
  }
  nvars=3;
  
}

//------------------------------------------------------------

bool SleepParams::read_sigma_file(std::string fileName)
{
  std::ifstream file;
  try {
    file.open(fileName.c_str(), std::ios::in);
  }
  catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    return false;
  }

  if (file.is_open()) {
    while(!file.eof()) {
      std::string line; 
      //read line
      getline(file, line);
      std::vector<std::string> strs;
      boost::split(strs, line, boost::is_any_of("\t "));
      if (strs.size() > 1) {
        double time, sigma;
        std::istringstream iss;
        iss.str(strs[0]);
        iss >> time;
        iss.str(strs[1]);
        iss >> sigma;
        sigma_vector.push_back(std::make_pair(time, sigma));
      }
    }
    next_sigma=0;
    sigma=0;
    
    file.close();
  }

  return true;
  
}

//------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, const SleepParams& x)
{
  os << std::endl
     << "Model parameters:" << std::endl
     << "-----------------" << std::endl
     << "nvars   = " << x.nvars << std::endl
     << "gamma   = " << x.gamma << std::endl
     << "beta    = " << x.beta << std::endl
     << "delta   = " << x.delta << std::endl
     << "sigma   = " << x.sigma << std::endl
     << "rho     = " << x.rho << std::endl
     << "epsilon = " << x.epsilon << std::endl
     << "N       = " << x.N << std::endl
     << std::endl;
  
  return os;
   
} // operator<<

//------------------------------------------------------------

struct Sleep
{

  static double screening_rate(double sigma, double tau, double t)
  {
    if (t >= tau) {
      return sigma * (t - tau);
    } else {
      return 0;
    }
  }
  
  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {         
    SleepParams p = *(static_cast<SleepParams*>(params));
    
    // local readable short variables
    double beta=p.beta, gamma=p.gamma, delta=p.delta;
    double rho=p.rho, sigma=p.sigma, epsilon=p.epsilon;
    unsigned int N=p.N;
    std::vector<std::pair<double, double> > sigma_vector=p.sigma_vector;
    unsigned int next_sigma=p.next_sigma;

    if (sigma_vector.size() > 0) {
      if (sigma_vector[next_sigma].first < t) {
        sigma = sigma_vector[next_sigma].second;
        ++next_sigma;
      }
    }

    rhs[0] = -beta*y[0]*y[1]/N + (rho+sigma+gamma)*y[1];
    rhs[1] = beta*y[0]*y[1]/N - (rho+sigma+gamma)*y[1];
    rhs[2] = (rho+sigma)*y[1] - delta*y[2];

    return GSL_SUCCESS;         
  }
      
}; // Sleep

//------------------------------------------------------------
         
#endif /* SLEEP_HH */
