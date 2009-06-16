#ifndef INFO_SIRS_MF_H
#define INFO_SIRS_MF_H

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>

#include <boost/program_options.hpp>

//------------------------------------------------------------

namespace po = boost::program_options;

//------------------------------------------------------------

namespace InfoSIRSmf
{

  //------------------------------------------------------------
  // Params structure
  
  struct Params
  {
    Params() : nvars(6) {};
    
    // No. of equations
    unsigned int nvars;
    
    // rates
    double gamma[2], delta[2];
    double tau[2][2], chi, mu, lambda;
    double omega;
    double Qd, Qi;
    
    // Graph properties
    double N;
    
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
       << "nvars  = " << x.nvars << std::endl
       << "beta--  = " << x.tau[0][0] << std::endl
       << "beta-+  = " << x.tau[0][1] << std::endl
       << "beta+-  = " << x.tau[1][0] << std::endl
       << "beta++  = " << x.tau[1][1] << std::endl
       << "gamma- = " << x.gamma[0] << std::endl
       << "gamma+ = " << x.gamma[1] << std::endl
       << "chi    = " << x.chi << std::endl
       << "lambda = " << x.lambda << std::endl
       << "mu     = " << x.mu << std::endl
       << "omega  = " << x.omega << std::endl
       << "-----------------" << std::endl
       << "N      = " << x.N << std::endl
       << "Qd     = " << x.Qd << std::endl
       << "Qi     = " << x.Qi << std::endl      
       << "-----------------" << std::endl;
    
    return os;
    
  } // operator<<
  
  //------------------------------------------------------------
  // Equations structure
  
  struct Eqs
  {
    // rhs function
    static int rhs_eval (double t, const double y[], double rhs[], void* params)
    {         
      Params p = *(static_cast<Params*>(params));
      
      // local readable short variables
      double tmm=p.tau[0][0], gm=p.gamma[0], dm=p.delta[0];
      double tpp=p.tau[1][1], gp=p.gamma[1], dp=p.delta[1];
      double tmp=p.tau[0][1], tpm=p.tau[1][0], chi=p.chi, mu=p.mu;
      double lm=p.lambda, om=p.omega;
      double N=p.N, invN=1./N;
      
      // Equations

      rhs[0] = -invN * ( tmm*y[1] + tmp*y[4] ) * y[0] + dm * y[2]
        - chi * invN * ( y[3] + y[4] + y[5] ) * y[0] + lm * y[3]
        - mu * invN * ( y[1] + y[4] ) * y[0];
      
      rhs[1] =  invN * ( tmm*y[1] + tmp*y[4] ) * y[0] - gm * y[1]
        - chi * invN * ( y[3] + y[4] + y[5] ) * y[1] + lm * y[4]
        - mu * invN * ( y[1] + y[4] ) * y[1] - om * y[1];
      
      rhs[2] =  gm * y[1] - dm * y[2] + lm * y[5]
        - chi * invN * ( y[3] + y[4] + y[5] ) * y[2]
        - mu * invN * ( y[1] + y[4] ) * y[2];
      
      rhs[3] = -invN * ( tpp*y[4] + tpm*y[1] ) * y[3] + dp * y[5]
        + chi * invN * ( y[3] + y[4] + y[5] ) * y[0] - lm * y[3]
        + mu * invN * ( y[1] + y[4] ) * y[0];
      
      rhs[4] =  invN * ( tpp*y[4] + tpm*y[1] ) * y[3] - gp * y[4]
        + chi * invN * ( y[3] + y[4] + y[5] ) * y[1] - lm * y[4]
        + mu * invN * ( y[1] + y[4] ) * y[1] + om * y[1];
      
      rhs[5] =  gp * y[4] - dp * y[5] - lm * y[5]
        + chi * invN * ( y[3] + y[4] + y[5] ) * y[2]
        + mu * invN * ( y[1] + y[4] ) * y[2];
       
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
      ("tau--", po::value<double>(),
       "disease transmission rate uninformed->uninformed")
      ("tau+-", po::value<double>(),
       "disease transmission rate informed->uninformed")
      ("tau-+", po::value<double>(),
       "disease transmission rate uninformed->informed")
      ("tau++", po::value<double>(),
       "disease transmission rate informed->informed")
      //------------------------------------------------
      ("gamma-", po::value<double>(),
       "recovery rate of uninformed")
      ("gamma+", po::value<double>(),
       "recovery rate of informed")
      //------------------------------------------------
      ("delta-", po::value<double>(),
       "loss of immunity rate of uninformed")
      ("delta+", po::value<double>(),
       "loss of immunity rate of informed")
      //------------------------------------------------
      ("chi", po::value<double>(),
       "information transmission rate")
      ("mu", po::value<double>(),
       "information generation rate")
      ("omega", po::value<double>(),
       "local information generation rate")
      ("lambda", po::value<double>(),
       "loss of information rate")
      //------------------------------------------------
      ("sigma", po::value<double>(),
       "ratio between uninformed/uninformed susceptibilities")
      ("vertices,N", po::value<double>(),
       "total number of individuals")
      ("nvars", po::value<int>(),
       "total number of equations")
      //------------------------------------------------
      ("beta--", po::value<double>(),       
       "disease transmission rate uninformed->uninformed")
      ("beta+-", po::value<double>(),
       "disease transmission rate informed->uninformed")
      ("beta-+", po::value<double>(),
       "disease transmission rate uninformed->informed")
      ("beta++", po::value<double>(),
       "disease transmission rate informed->informed")
      ;
    
    return opt;
    
  }
  
  //------------------------------------------------------------
  // generating graph options
  
  po::options_description* generate_graph_options(ode::OdeSolver<Params, Eqs>& dummy)
  {
    // dummy is for overloading    
    
    po::options_description* opt =
      new po::options_description("Graph parameters");
    
    opt->add_options()
      ("Qd", po::value<double>(),
       "average d-degree")
      ("Qi", po::value<double>(),
       "average i-degree")
      ;
    
    return opt;
  }

  //------------------------------------------------------------
  // generating hidden options
  
  po::options_description* generate_hidden_options(ode::OdeSolver<Params, Eqs>& dummy)
  {
    // dummy is for overloading
    
    po::options_description* opt =
      new po::options_description("Hidden parameters");
    
    opt->add_options()
      ("vertices,N", po::value<double>(),
       "total number of individuals")
      ("qid", po::value<double>(),
       "conditional probability for i-edge given d-edge")
      ("qdi", po::value<double>(),
       "conditional probability for d-edge given i-edge")
      //------------------------------------------------
      ("Cddd", po::value<double>(),
       "local clustering coefficient")
      ("Cddi", po::value<double>(),
       "local clustering coefficient")
      ("Cdid", po::value<double>(),
       "local clustering coefficient")
      ("Cdii", po::value<double>(),
       "local clustering coefficient")
      ("Ciid", po::value<double>(),
       "local clustering coefficient")
      ("Ciii", po::value<double>(),
       "local clustering coefficient")
      ;
    
    return opt;
  }
  
  //------------------------------------------------------------  
  //  initial model parameters from vm
  
  int init_model_params(po::variables_map& vm,
                        ode::OdeSolver<Params, Eqs>& x)
  {
    Params* model_params = x.get_model_params(); 
    
    if (vm.count("tau--")) {
      model_params->tau[0][0]=vm["tau--"].as<double>();
    } else {
      std::cerr << "WARNING: no tau-- given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->tau[0][0]=0;
    }
    if (vm.count("tau+-")) {
      model_params->tau[1][0]=vm["tau+-"].as<double>();
    } else {
      std::cerr << "WARNING: no tau+- given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->tau[1][0]=0;
    }
    if (vm.count("tau-+")) {
      model_params->tau[0][1]=vm["tau-+"].as<double>();
    } else {
      std::cerr << "WARNING: no tau-+ given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->tau[0][1]=0;
    }
    if (vm.count("tau++")) {
      model_params->tau[1][1]=vm["tau++"].as<double>();
    } else {
      std::cerr << "WARNING: no tau++ given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->tau[1][1]=0;
    }
    //------------------------------------------------
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
    //------------------------------------------------
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
    //------------------------------------------------
    if (vm.count("chi")) {
      model_params->chi=vm["chi"].as<double>();
    } else {
      std::cerr << "WARNING: no chi given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->chi=0;
    }
    if (vm.count("mu")) {
      model_params->mu=vm["mu"].as<double>();
    } else {
      std::cerr << "WARNING: no mu given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->mu=0;
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
    //------------------------------------------------
    if (vm.count("vertices")) {
      model_params->N=vm["vertices"].as<double>();
    } else {
      std::cerr << "ERROR: no N given" << std::endl;
      return 1;
    }
    if (vm.count("nvars")) {
      model_params->nvars=vm["nvars"].as<int>();
    }
    //------------------------------------------------
    if (vm.count("beta--")) {
      std::cout << "WARNING: beta-- is given" << std::endl;
      std::cout << "setting to tau--" << std::endl;
      
      model_params->tau[0][0]=vm["beta--"].as<double>();
    }    
    if (vm.count("beta-+")) {
      std::cout << "WARNING: beta-+ is given" << std::endl;
      std::cout << "setting to tau-+" << std::endl;
      
      model_params->tau[0][1]=vm["beta-+"].as<double>();
    }    
    if (vm.count("beta+-")) {
      std::cout << "WARNING: beta+- is given" << std::endl;
      std::cout << "setting to tau+-" << std::endl;
      
      model_params->tau[1][0]=vm["beta+-"].as<double>();
    }    
    if (vm.count("beta++")) {
      std::cout << "WARNING: beta++ is given" << std::endl;
      std::cout << "setting to tau++" << std::endl;
      
      model_params->tau[1][1]=vm["beta++"].as<double>();
    }
    //------------------------------------------------
    if (vm.count("sigma")) {
      std::cout << "WARNING: sigma is given" << std::endl;
      std::cout << "rescaling tau+-, tau++" << std::endl;
      
      model_params->tau[1][0] = model_params->tau[0][0] * vm["sigma"].as<double>();
      model_params->tau[1][1] = model_params->tau[0][1] * vm["sigma"].as<double>();
    }
    
    return 0;
  }
  
  //------------------------------------------------------------  
  //  initial graph parameters from vm
  
  int init_graph_params(po::variables_map& vm,
                        ode::OdeSolver<Params, Eqs>& x)
  {
    std::cout << "init_graph_params" << std::endl;
    Params* model_params = x.get_model_params(); 
    
    if (vm.count("Qd")) {
      model_params->Qd=vm["Qd"].as<double>();
    } else {
      std::cerr << "WARNING: no Qd given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->Qd=0;
    }
    model_params->tau[0][0] *= model_params->Qd;
    model_params->tau[0][1] *= model_params->Qd;
    model_params->tau[1][0] *= model_params->Qd;
    model_params->tau[1][1] *= model_params->Qd;
    
    if (vm.count("Qi")) {
      model_params->Qi=vm["Qi"].as<double>();
    } else {
      std::cerr << "WARNING: no Qi given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->Qi=0;
    }
    model_params->chi *= model_params->Qi;
    model_params->mu  *= model_params->Qi;

    return 0;
  }

  //------------------------------------------------------------  
  
} // namespace InfoSIRSmf

//------------------------------------------------------------

#endif // INFO_SIRS_MF_H
