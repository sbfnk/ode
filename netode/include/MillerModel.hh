#ifndef MILLER_MODEL_HH
#define MILLER_MODEL_HH

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>
#include <math.h>

#include "net_ode.hh"

#include <boost/program_options.hpp>

//------------------------------------------------------------

namespace po = boost::program_options;

//------------------------------------------------------------

namespace MillerModel
{

  //------------------------------------------------------------
  // Params structure
  
  struct Params
  {
    Params() : nvars(2) {};
    
    // No. of equations
    unsigned int nvars;
    
    // rates
    double beta, gamma;

    std::vector<double> degreeDist;

    double (*suscept)(double, std::vector<double> const&);
    double (*suscept_prime)(double, std::vector<double> const&);
    
    // overloading operator<<
    friend std::ostream& operator <<
      (std::ostream& os, const Params& x);
    
  }; // Params
  
  //------------------------------------------------------------
  // overloading operator<<
  
  std::ostream& operator<< (std::ostream& os, const Params& x)
  {
    os << std::endl
       << "Model parameters:\n"
       << "=================\n"
       << "beta     = " << x.beta << std::endl
       << "gamma    = " << x.gamma << std::endl;
    
    return os;
    
  } // operator<<
  
  static double psi (double theta, std::vector<double> const& degDist)
  {
    double res = .0;
    for (unsigned int i = 0; i < degDist.size(); ++i) {
      res += degDist[i] * pow(theta, i);
    }
    return res;
  }

  static double nbpsi (double theta, std::vector<double> const& degDist)
  {
    double res = .0;
    for (unsigned int i = 1; i < degDist.size(); ++i) {
      res += i * degDist[i] * pow(theta, i);
    }
    if (theta < 1) { res /= nbpsi(1, degDist); }
    return res;
  }

  static double psi_prime (double theta, std::vector<double> const& degDist)
  {
    double res = .0;
    for (unsigned int i = 1; i < degDist.size(); ++i) {
      res += i * degDist[i] * pow(theta, i-1);
    }
    return res;
  }

  static double nbpsi_prime (double theta, std::vector<double> const& degDist)
  {
    double res = .0;
    for (unsigned int i = 1; i < degDist.size(); ++i) {
      res += i * i * degDist[i] * pow(theta, i-1);
    }
    if (theta < 1) { res /= nbpsi(1, degDist); }
    return res;
  }

  static double psi_2prime (double theta, std::vector<double> const& degDist)
  {
    double res = .0;
    for (unsigned int i = 1; i < degDist.size(); ++i) {
      res += i * (i-1) * degDist[i] * pow(theta, i-2);
    }
    return res;
  }

  //------------------------------------------------------------
  // Equations structure
  
  struct Eqs
  {

    static void write_rhs(double rhs[], void* params, std::ofstream& os)
    {
      Params p = *(static_cast<Params*>(params));
      std::vector<double> const& degDist = p.degreeDist;
      double beta = p.beta, gamma = p.gamma;

      double R = rhs[1];
      double S = p.suscept(rhs[0], degDist);
      double I = 1 - R - S;
      
      os << rhs[0] << '\t';
      os << S << '\t';
      os << I << '\t';
      os << R << '\t';
      os << (- gamma*I -
             p.suscept_prime(rhs[0], degDist) *
             (-beta * rhs[0] +
              gamma * (1 - rhs[0]) +
              beta * psi_prime(rhs[0], degDist) / psi_prime(1, degDist)));
    }
    
    static void init(double rhs[], void* params)
    {
      Params p = *(static_cast<Params*>(params));
      std::vector<double> const& degDist = p.degreeDist;
      double beta = p.beta, gamma = p.gamma;

      rhs[1] = (1 - p.suscept(rhs[0], degDist))/
        (-beta/gamma * (1-psi_2prime(1, degDist) / psi_prime(1, degDist)));
    }
    
    static void finish(double rhs[], void* params, std::ofstream& os)
    {
      // Params p = *(static_cast<Params*>(params));
      // std::vector<double> const& degDist = p.degreeDist;
      // double beta = p.beta, gamma = p.gamma;

      // os << "Final theta: " << (1-beta/(beta+gamma) +
      //                           beta/(beta+gamma) *
      //                           psi_prime(rhs[0], degDist) /
      //                           psi_prime(1, degDist)) << std::endl; 
    }
    
    // rhs function
    static int rhs_eval (double t, const double y[], double rhs[], void* params)
    {         
      Params p = *(static_cast<Params*>(params));
      
      // local readable short variables
      double beta = p.beta, gamma = p.gamma;
      std::vector<double> const& degDist = p.degreeDist;

      // Equations

      rhs[0] = -beta * y[0] + gamma * (1 - y[0]) +
        beta * psi_prime(y[0], degDist) / psi_prime(1, degDist);

      rhs[1] = gamma * (1 - y[1] - p.suscept(y[0], degDist));
      
      return GSL_SUCCESS;         
    }
    
  }; // Eqs

  //------------------------------------------------------------
  // generating model options
  
  po::options_description* generate_model_options(ode::OdeSolver<Params, Eqs>& dummy)
  {
    // dummy is for overloading
    
    po::options_description* opt =
      new po::options_description("Model parameters");
    
    opt->add_options()
      ("beta,b", po::value<double>(),
       "per contact transmission rate")
      ("gamma,g", po::value<double>(),
       "recovery rate")
      ("transmission,T", po::value<double>(),
       "transmission probability")
      ("degree-file,d", po::value<std::string>(),
       "degree file")
      ("susfunc,s", po::value<int>(),
       "susceptible function (0: #susc, 1:#nb susc")
      ;      
    
    return opt;
  }
  
  //------------------------------------------------------------  
  //  initial model parameters from vm
  
  int init_model_params(po::variables_map& vm,
                        ode::OdeSolver<Params, Eqs>& x)
  {
    Params* model_params = x.get_model_params(); 
    
    if (vm.count("beta")) {
      model_params->beta=vm["beta"].as<double>();
    } else {
      if (!vm.count("transmission")) {
        std::cerr << "ERROR: no beta given" << std::endl;
        return 1;
      }
    }
    if (vm.count("gamma")) {
      model_params->gamma=vm["gamma"].as<double>();
    } else {
      if (!vm.count("transmission")) {
        std::cerr << "ERROR: no gamma given" << std::endl;
        return 1;
      }
    }

    if (vm.count("transmission")) {
      double T = vm["transmission"].as<double>();
      model_params->beta= T/(1-T);
      model_params->gamma=1;
    }
    
    //------------------------------------------------
    if (vm.count("degree-file")) {
      if (!read_single_degree_dist(model_params->degreeDist, vm["degree-file"].
                                   as<std::string>())) {
        std::cerr << "ERROR: problem reading " << vm["degree-file"].
          as<std::string>() << std::endl;
        return 1;
      }
    } else {
      std::cerr << "ERROR: no degree-file" << std::endl;
      return 1;
    }
    
    if (vm.count("susfunc")) {
      int susfunc=vm["susfunc"].as<int>();
      switch (susfunc) {
       case 0:
        model_params->suscept = &psi;
        model_params->suscept_prime = &psi_prime;
        break;
       case 1:
        model_params->suscept = &nbpsi;
        model_params->suscept_prime = &nbpsi_prime;
        break;
       default:
        std::cerr << "WARNING: no susfunc " << susfunc << ", setting to 0"
                  << std::endl;
        model_params->suscept = psi;
        model_params->suscept_prime = psi_prime;
      }
    } else {
      std::cerr << "WARNING: no susfunc given, setting to 0" << std::endl;
      model_params->suscept = psi;
      model_params->suscept_prime = psi_prime;
    }

    return 0;
  }
  
  
} // namespace MillerModel

//------------------------------------------------------------

#endif // MILLER_MODEL_HH
