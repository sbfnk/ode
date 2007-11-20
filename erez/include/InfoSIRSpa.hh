#ifndef INFO_SIRS_PA_H
#define INFO_SIRS_PA_H

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>

#include "pa_macros.hh"

#include <boost/program_options.hpp>

//------------------------------------------------------------

namespace po = boost::program_options;

//------------------------------------------------------------

namespace InfoSIRSpa
{

  //------------------------------------------------------------
  // Params structure
  
  struct Params
  {
    Params() : nvars(48) {};
    
    // No. of equations
    unsigned int nvars;
    
    // rates
    double gamma[2], delta[2];
    double tau[2][2], chi, mu, lambda;
    double omega;
    
    // Graph properties
    double Qd, Qi, N;
    double qdi, qid;


    // FIX CLUSTERING COEFFICIENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    
    // clustering
    double C[3][3][3];
    
    // claculate Qd, Qi and Qb from i.c.
    //void init_Q(double*, bool verbose);
    
    
    // overloading operator<<
    friend std::ostream& operator <<
      (std::ostream& os, const Params& x);
    
  }; // Params
  
  //------------------------------------------------------------
  // generating model options
  
  po::options_description* generate_model_options()
  {
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
      ("gamma-", po::value<double>(),
       "recovery rate of uninformed")
      ("gamma+", po::value<double>(),
       "recovery rate of informed")
      ("delta-", po::value<double>(),
       "loss of immunity rate of uninformed")
      ("delta+", po::value<double>(),
       "loss of immunity rate of informed")
      ("chi", po::value<double>(),
       "information transmission rate")
      ("mu", po::value<double>(),
       "information generation rate")
      ("omega", po::value<double>(),
       "local information generation rate")
      ("lambda", po::value<double>(),
       "loss of information rate")
      ("sigma", po::value<double>(),
       "ratio between uninformed/uninformed susceptibilities")
      ("vertices,N", po::value<double>(),
       "total number of individuals")
      ("nvars", po::value<int>(),
       "total number of equations")
      ("clustering", po::value<double>(),
       "clustering coefficient value for all coefficients")
      ("Qd", po::value<double>(),
       "average d-degree")
      ("Qi", po::value<double>(),
       "average i-degree")
      ("qid", po::value<double>(),
       "conditional probability for i-edge given d-edge")
      ("qdi", po::value<double>(),
       "conditional probability for d-edge given i-edge")
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
  
  template <class Params, class Eqs>
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
    if (vm.count("Qd")) {
      model_params->Qd=vm["Qd"].as<double>();
    } else {
      std::cerr << "WARNING: no Qd given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->Qd=0;
    }
    if (vm.count("Qi")) {
      model_params->Qi=vm["Qi"].as<double>();
    } else {
      std::cerr << "WARNING: no Qi given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->Qi=0;
    }
    if (vm.count("qid")) {
      model_params->qid=vm["qid"].as<double>();
    } else {
      std::cerr << "WARNING: no qid given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->qid=0;
    }
    if (vm.count("qdi")) {
      model_params->qdi=vm["qdi"].as<double>();
    } else {
      std::cerr << "WARNING: no qdi given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      model_params->qdi=0;
    }
    if (vm.count("sigma")) {
      model_params->lambda=vm["sigma"].as<double>();
    }
    if (vm.count("vertices")) {
      model_params->N=vm["vertices"].as<double>();
    } else {
      std::cerr << "ERROR: no N given" << std::endl;
      return 1;
    }
    if (vm.count("nvars")) {
      model_params->nvars=vm["nvars"].as<int>();
    }
    if (vm.count("clustering")) {
      int i,j,k;
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          for (k = 0; k < 3; k++)
            model_params->C[i][j][k] = vm["clustering"].as<double>();
    } else {
      std::cerr << "WARNING: no cluster-coefficient given" << std::endl;
      std::cerr << "setting to 0" << std::endl;
      int i,j,k;
      for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
          for (k = 0; k < 3; k++)
            model_params->C[i][j][k] = 0;
    }
    
    
    if (vm.count("beta--")) 
      model_params->tau[0][0]=vm["beta--"].as<double>();
    
    if (vm.count("beta+-")) 
      model_params->tau[1][0]=vm["beta+-"].as<double>();
    
    if (vm.count("beta-+")) 
      model_params->tau[0][1]=vm["beta-+"].as<double>();
    
    if (vm.count("beta++")) 
      model_params->tau[1][1]=vm["beta++"].as<double>();
    
    return 0;
  }
  
  //------------------------------------------------------------
  

  


  
  //------------------------------------------------------------
  
  // void InfoSISParams::init_Q(double* rhs, bool verbose)
  // {   
  //   // calculate Qd and Qi in from i.c.
  //   if (nvars > 6) {
    
  //     // calculate total number of disease pairs
  //     if (verbose) std::cout << "... calculating Qd from i.c.";
    
  //     double dPairs =
  //       SS_d_t + SI_d_t + SR_d_t + Ss_d_t + Si_d_t + Sr_d_t +
  //       II_d_t + IR_d_t + Is_d_t + Ii_d_t + Ir_d_t +
  //       RR_d_t + Rs_d_t + Ri_d_t + Rr_d_t +
  //       ss_d_t + si_d_t + sr_d_t +
  //       ii_d_t + ir_d_t +
  //       rr_d_t;
    
  //     Qd = (2 * dPairs / N);
    
  //     if (verbose) std::cout << " ... done "
  //                            << "(Qd=" << Qd << ")" << std::endl;
    
  //     // calculate total number of info pairs
  //     if (verbose) std::cout << "... calculating Qi from i.c.";
    
  //     double iPairs =
  //       SS_i_t + SI_i_t + SR_i_t + Ss_i_t + Si_i_t + Sr_i_t +
  //       II_i_t + IR_i_t + Is_i_t + Ii_i_t + Ir_i_t +
  //       RR_i_t + Rs_i_t + Ri_i_t + Rr_i_t +
  //       ss_i_t + si_i_t + sr_i_t +
  //       ii_i_t + ir_i_t +
  //       rr_i_t;
    
  //     Qi = (2 * iPairs / N);      
  //     if (verbose) std::cout << " ... done "
  //                            << "(Qi=" << Qi << ")" << std::endl;
    
  //   } else {
  //     if (verbose) std::cout << "... initialize Qd and Qi to default values";
  //     if (verbose) std::cout << " ... done\n";
  //   }
  
  //   // dib model
  //   if (nvars > 48) {
      
  //     // calculate total number of parallel pairs
  //     if (verbose) std::cout << "... calculating Qb from i.c.";
      
  //     double bPairs =
  //       SS_b_t + SI_b_t + SR_b_t + Ss_b_t + Si_b_t + Sr_b_t +
  //       II_b_t + IR_b_t + Is_b_t + Ii_b_t + Ir_b_t +
  //       RR_b_t + Rs_b_t + Ri_b_t + Rr_b_t +
  //       ss_b_t + si_b_t + sr_b_t +
  //       ii_b_t + ir_b_t +
  //       rr_b_t;
      
  //     Qb = (2 * bPairs / N);      
  //     if (verbose) std::cout << " ... done "
  //                            << "(Qb=" << Qb << ")" << std::endl;
      
  //   } else {
  //     if (verbose) std::cout << "... initialize Qb to default values";
  //     if (verbose) std::cout << " ... done\n";
  //   }

  // } // init_Q

  //------------------------------------------------------------

  //------------------------------------------------------------
  // overloading operator<<
  
  std::ostream& operator<< (std::ostream& os, const Params& x)
  {
    os << std::endl
       << "Model parameters:\n"
       << "=================\n"
       << "nvars  = " << x.nvars << std::endl
       << "tau--  = " << x.tau[0][0] << std::endl
       << "tau-+  = " << x.tau[0][1] << std::endl
       << "tau+-  = " << x.tau[1][0] << std::endl
       << "tau++  = " << x.tau[1][1] << std::endl
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
       << "qdi    = " << x.qdi << std::endl
       << "qid    = " << x.qid << std::endl
       << "-----------------" << std::endl
       << "C_ddi  = " << x.C[1][1][0] << std::endl
       << "C_ddd  = " << x.C[1][1][1] << std::endl
       << "C_dii  = " << x.C[1][0][0] << std::endl
       << "C_did  = " << x.C[1][0][1] << std::endl
       << "C_idi  = " << x.C[0][1][0] << std::endl
       << "C_idd  = " << x.C[0][1][1] << std::endl
       << "C_iii  = " << x.C[0][0][0] << std::endl
       << "C_iid  = " << x.C[0][0][1] << std::endl
       << "-----------------" << std::endl;
    
    return os;
    
  } // operator<<
  
  //------------------------------------------------------------
  // Equations structure
  
  struct Eqs
  {
    // pair approx      
    static double pa(double kk, double ij, double jk, double j)
    {
      const double eps = 1e-4;
      
      if(j < eps)
        return 0.0;
      else
        return kk*((ij*jk)/j);
    }
    
    // clustering correction
    static double cc(double C_i, double C_d, double ki_i, double ki_d,
                     double ki, double N_Qd, double N_Qi)
    {
      const double eps = 1e-4;
      
      if(ki < eps)
        return ((1.0 - C_i - C_d));
      else
        return ((1.0 - C_i - C_d) + C_i*N_Qi*ki_i/ki + C_d*N_Qd*ki_d/ki);
    }           
    
    // rhs function
    static int rhs_eval (double t, const double y[], double rhs[], void* params)
    {         
      Params p = *(static_cast<Params*>(params));
      
      // local readable short variables
      double tmm=p.tau[0][0], gm=p.gamma[0], dm=p.delta[0];
      double tpp=p.tau[1][1], gp=p.gamma[1], dp=p.delta[1];
      double tmp=p.tau[0][1], tpm=p.tau[1][0], chi=p.chi, mu=p.mu;
      double lm=p.lambda, om=p.omega;
      double N=p.N, Qd=p.Qd, Qi=p.Qi;
      double N_Qd=N/Qd, N_Qi=N/Qi;

      // Clustering corrections
      double idi=p.C[0][1][0], idd=p.C[0][1][1];
      double dii=p.C[1][0][0], did=p.C[1][0][1];
      double ddi=p.C[1][1][0], ddd=p.C[1][1][1];
      double iii=p.C[0][0][0], iid=p.C[0][0][1];
      
      double qdi=p.qdi, qid=p.qid;
      double kdd=(Qd-1.0)/Qd, kii=(Qi-1.0)/Qi;
      double kdi=(Qd-qdi)/Qd, kid=(Qi-qid)/Qi;
      
      // Equations

      S_t = - tmm*SI_d
        - tmp*Si_d
        + dm*R
        + lm*s
        - chi*Ss_i
        - chi*Si_i
        - chi*Sr_i
        - mu*SI_i
        - mu*Si_i;

      I_t = + tmm*SI_d
        + tmp*Si_d
        - gm*I
        + lm*i
        - om*I
        - chi*Is_i
        - chi*Ii_i
        - chi*Ir_i
        - mu*II_i
        - mu*Ii_i;
 
      R_t = + gm*I
        - dm*R
        + lm*r
        - chi*Rs_i
        - chi*Ri_i
        - chi*Rr_i
        - mu*IR_i
        - mu*Ri_i;
 
      s_t = - tpp*si_d
        - tpm*Is_d
        + dp*r
        - lm*s
        + chi*Ss_i
        + chi*Si_i
        + chi*Sr_i
        + mu*SI_i
        + mu*Si_i;
 
      i_t = + tpp*si_d
        + tpm*Is_d
        - gp*i
        - lm*i
        + om*I
        + chi*Is_i
        + chi*Ii_i
        + chi*Ir_i
        + mu*II_i
        + mu*Ii_i;
 
      r_t = + gp*i
        - dp*r
        - lm*r
        + chi*Rs_i
        + chi*Ri_i
        + chi*Rr_i
        + mu*IR_i
        + mu*Ri_i;
 
      SS_d_t = 2*(
                  - tmm*pa(kdd, SS_d, SI_d, S)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
                  - tmp*pa(kdd, SS_d, Si_d, S)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
                  + dm*SR_d
                  + lm*Ss_d
                  - chi*pa(kdi, SS_d, Ss_i, S)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
                  - chi*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
                  - chi*pa(kdi, SS_d, Sr_i, S)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
                  - mu*pa(kdi, SS_d, SI_i, S)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
                  - mu*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
                  );

      SS_i_t = 2*(
                  - tmm*pa(kid, SS_i, SI_d, S)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
                  - tmp*pa(kid, SS_i, Si_d, S)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
                  + dm*SR_i
                  + lm*Ss_i
                  - chi*pa(kii, SS_i, Ss_i, S)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
                  - chi*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
                  - chi*pa(kii, SS_i, Sr_i, S)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
                  - mu*pa(kii, SS_i, SI_i, S)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
                  - mu*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
                  );

      SI_d_t = - tmm*SI_d
        + tmm*pa(kdd, SS_d, SI_d, S)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - tmm*pa(kdd, IS_d, SI_d, S)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
        + tmp*pa(kdd, SS_d, Si_d, S)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - tmp*pa(kdd, iS_d, SI_d, S)*cc(ddi, ddd, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - gm*SI_d
        + dm*IR_d
        + lm*Is_d
        + lm*Si_d
        - om*SI_d
        - chi*pa(kid, sS_i, SI_d, S)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
        - chi*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - chi*pa(kid, rS_i, SI_d, S)*cc(idi, idd, rI_i, rI_d, r*I, N_Qd, N_Qi)
        - chi*pa(kdi, SI_d, Is_i, I)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        - chi*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - chi*pa(kdi, SI_d, Ir_i, I)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kid, IS_i, SI_d, S)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kdi, SI_d, II_i, I)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - mu*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - qid*mu*SI_d;

      SI_i_t = + tmm*pa(kid, SS_i, SI_d, S)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - tmm*pa(kdi, IS_d, SI_i, S)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
        + tmp*pa(kid, SS_i, Si_d, S)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - tmp*pa(kdi, iS_d, SI_i, S)*cc(dii, did, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - gm*SI_i
        + dm*IR_i
        + lm*Is_i
        + lm*Si_i
        - om*SI_i
        - chi*pa(kii, sS_i, SI_i, S)*cc(iii, iid, sI_i, sI_d, s*I, N_Qd, N_Qi)
        - chi*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - chi*pa(kii, rS_i, SI_i, S)*cc(iii, iid, rI_i, rI_d, r*I, N_Qd, N_Qi)
        - chi*pa(kii, SI_i, Is_i, I)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        - chi*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - chi*pa(kii, SI_i, Ir_i, I)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*SI_i
        - mu*pa(kii, IS_i, SI_i, S)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kii, SI_i, II_i, I)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - mu*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - qdi*tmm*SI_i;

      SR_d_t = - tmm*pa(kdd, IS_d, SR_d, S)*cc(ddi, ddd, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - tmp*pa(kdd, iS_d, SR_d, S)*cc(ddi, ddd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + gm*SI_d
        - dm*SR_d
        + dm*RR_d
        + lm*Rs_d
        + lm*Sr_d
        - chi*pa(kid, sS_i, SR_d, S)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
        - chi*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - chi*pa(kid, rS_i, SR_d, S)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kdi, SR_d, Rs_i, R)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        - chi*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - chi*pa(kdi, SR_d, Rr_i, R)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kid, IS_i, SR_d, S)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kdi, SR_d, RI_i, R)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        ;

      SR_i_t = - tmm*pa(kdi, IS_d, SR_i, S)*cc(dii, did, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - tmp*pa(kdi, iS_d, SR_i, S)*cc(dii, did, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + gm*SI_i
        - dm*SR_i
        + dm*RR_i
        + lm*Rs_i
        + lm*Sr_i
        - chi*pa(kii, sS_i, SR_i, S)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
        - chi*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - chi*pa(kii, rS_i, SR_i, S)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kii, SR_i, Rs_i, R)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        - chi*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - chi*pa(kii, SR_i, Rr_i, R)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kii, IS_i, SR_i, S)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kii, SR_i, RI_i, R)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        ;
            
      Ss_d_t = - tmm*pa(kdd, IS_d, Ss_d, S)*cc(ddi, ddd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - tpp*pa(kdd, Ss_d, si_d, s)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - tmp*pa(kdd, iS_d, Ss_d, S)*cc(ddi, ddd, is_i, is_d, i*s, N_Qd, N_Qi)
        - tpm*pa(kdd, Ss_d, sI_d, s)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
        + dm*Rs_d
        + dp*Sr_d
        - lm*Ss_d
        + lm*ss_d
        + chi*pa(kdi, SS_d, Ss_i, S)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        - chi*pa(kid, sS_i, Ss_d, S)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
        + chi*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - chi*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        + chi*pa(kdi, SS_d, Sr_i, S)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - chi*pa(kid, rS_i, Ss_d, S)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kdi, SS_d, SI_i, S)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kid, IS_i, Ss_d, S)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - mu*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        - qid*chi*Ss_d;
               
      Ss_i_t = - tmm*pa(kdi, IS_d, Ss_i, S)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - tpp*pa(kid, Ss_i, si_d, s)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - tmp*pa(kdi, iS_d, Ss_i, S)*cc(dii, did, is_i, is_d, i*s, N_Qd, N_Qi)
        - tpm*pa(kid, Ss_i, sI_d, s)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
        + dm*Rs_i
        + dp*Sr_i
        - lm*Ss_i
        + lm*ss_i
        - chi*Ss_i
        + chi*pa(kii, SS_i, Ss_i, S)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        - chi*pa(kii, sS_i, Ss_i, S)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
        + chi*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - chi*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        + chi*pa(kii, SS_i, Sr_i, S)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - chi*pa(kii, rS_i, Ss_i, S)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kii, SS_i, SI_i, S)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kii, IS_i, Ss_i, S)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - mu*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        ;
            
      Si_d_t = - tmm*pa(kdd, IS_d, Si_d, S)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tpp*pa(kdd, Ss_d, si_d, s)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - tmp*Si_d
        - tmp*pa(kdd, iS_d, Si_d, S)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + tpm*pa(kdd, Ss_d, sI_d, s)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - gp*Si_d
        + dm*Ri_d
        + lm*si_d
        - lm*Si_d
        + om*SI_d
        - chi*pa(kid, sS_i, Si_d, S)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
        - chi*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        - chi*pa(kid, rS_i, Si_d, S)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + chi*pa(kdi, SI_d, Is_i, I)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        + chi*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        + chi*pa(kdi, SI_d, Ir_i, I)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kid, IS_i, Si_d, S)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kdi, SI_d, II_i, I)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + mu*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - qid*chi*Si_d
        - qid*mu*Si_d
        ;
            
      Si_i_t = - tmm*pa(kdi, IS_d, Si_i, S)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tpp*pa(kid, Ss_i, si_d, s)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - tmp*pa(kdi, iS_d, Si_i, S)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + tpm*pa(kid, Ss_i, sI_d, s)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - gp*Si_i
        + dm*Ri_i
        + lm*si_i
        - lm*Si_i
        + om*SI_i
        - chi*pa(kii, sS_i, Si_i, S)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
        - chi*Si_i
        - chi*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        - chi*pa(kii, rS_i, Si_i, S)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + chi*pa(kii, SI_i, Is_i, I)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        + chi*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        + chi*pa(kii, SI_i, Ir_i, I)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kii, IS_i, Si_i, S)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kii, SI_i, II_i, I)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*Si_i
        - mu*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + mu*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - qdi*tmp*Si_i;

      Sr_d_t = - tmm*pa(kdd, IS_d, Sr_d, S)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - tmp*pa(kdd, iS_d, Sr_d, S)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + gp*Si_d
        + dm*Rr_d
        - dp*Sr_d
        + lm*sr_d
        - lm*Sr_d
        - chi*pa(kid, sS_i, Sr_d, S)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
        - chi*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - chi*pa(kid, rS_i, Sr_d, S)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kdi, SR_d, Rs_i, R)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        + chi*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        + chi*pa(kdi, SR_d, Rr_i, R)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kid, IS_i, Sr_d, S)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kdi, SR_d, RI_i, R)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
        - qid*chi*Sr_d;

      Sr_i_t = - tmm*pa(kdi, IS_d, Sr_i, S)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - tmp*pa(kdi, iS_d, Sr_i, S)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + gp*Si_i
        + dm*Rr_i
        - dp*Sr_i
        + lm*sr_i
        - lm*Sr_i
        - chi*pa(kii, sS_i, Sr_i, S)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
        - chi*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - chi*Sr_i
        - chi*pa(kii, rS_i, Sr_i, S)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kii, SR_i, Rs_i, R)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
        + chi*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        + chi*pa(kii, SR_i, Rr_i, R)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
        - mu*pa(kii, IS_i, Sr_i, S)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kii, SR_i, RI_i, R)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
        - mu*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
        ;
            
      II_d_t = 2*(
                  + tmm*SI_d
                  + tmm*pa(kdd, IS_d, SI_d, S)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
                  + tmp*pa(kdd, iS_d, SI_d, S)*cc(ddi, ddd, iI_i, iI_d, i*I, N_Qd, N_Qi)
                  - gm*II_d
                  + lm*Ii_d
                  - om*II_d
                  - chi*pa(kdi, II_d, Is_i, I)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
                  - chi*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  - chi*pa(kdi, II_d, Ir_i, I)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                  - mu*pa(kdi, II_d, II_i, I)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
                  - mu*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  - qid*mu*II_d
                  );

      II_i_t = 2*(
                  + tmm*pa(kdi, IS_d, SI_i, S)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
                  + tmp*pa(kdi, iS_d, SI_i, S)*cc(dii, did, iI_i, iI_d, i*I, N_Qd, N_Qi)
                  - gm*II_i
                  + lm*Ii_i
                  - om*II_i
                  - chi*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
                  - chi*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  - chi*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                  - mu*II_i
                  - mu*pa(kii, II_i, II_i, I)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
                  - mu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  + qdi*tmm*SI_i
                  );

      IR_d_t = + tmm*pa(kdd, IS_d, SR_d, S)*cc(ddi, ddd, IR_i, IR_d, I*R, N_Qd, N_Qi)
        + tmp*pa(kdd, iS_d, SR_d, S)*cc(ddi, ddd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + gm*II_d
        - gm*IR_d
        - dm*IR_d
        + lm*Ri_d
        + lm*Ir_d
        - om*IR_d
        - chi*pa(kid, sI_i, IR_d, I)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
        - chi*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - chi*pa(kid, rI_i, IR_d, I)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kdi, IR_d, Rs_i, R)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - chi*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - chi*pa(kdi, IR_d, Rr_i, R)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - mu*pa(kid, II_i, IR_d, I)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kdi, IR_d, RI_i, R)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - qid*mu*IR_d;

      IR_i_t = + tmm*pa(kdi, IS_d, SR_i, S)*cc(dii, did, IR_i, IR_d, I*R, N_Qd, N_Qi)
        + tmp*pa(kdi, iS_d, SR_i, S)*cc(dii, did, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + gm*II_i
        - gm*IR_i
        - dm*IR_i
        + lm*Ri_i
        + lm*Ir_i
        - om*IR_i
        - chi*pa(kii, sI_i, IR_i, I)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
        - chi*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - chi*pa(kii, rI_i, IR_i, I)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - chi*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - chi*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - mu*pa(kii, II_i, IR_i, I)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*IR_i
        - mu*pa(kii, IR_i, RI_i, R)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        ;
            
      Is_d_t = + tmm*pa(kdd, IS_d, Ss_d, S)*cc(ddi, ddd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - tpp*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tmp*pa(kdd, iS_d, Ss_d, S)*cc(ddi, ddd, is_i, is_d, i*s, N_Qd, N_Qi)
        - tpm*Is_d
        - tpm*pa(kdd, Is_d, sI_d, s)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
        - gm*Is_d
        + dp*Ir_d
        - lm*Is_d
        + lm*si_d
        - om*Is_d
        + chi*pa(kid, sS_i, SI_d, S)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
        + chi*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
        + chi*pa(kid, rS_i, SI_d, S)*cc(idi, idd, rI_i, rI_d, r*I, N_Qd, N_Qi)
        - chi*pa(kid, sI_i, Is_d, I)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
        - chi*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        - chi*pa(kid, rI_i, Is_d, I)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kid, IS_i, SI_d, S)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kid, II_i, Is_d, I)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - mu*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        - qid*chi*Is_d
        + qid*mu*SI_d;

      Is_i_t = + tmm*pa(kdi, IS_d, Ss_i, S)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - tpp*pa(kid, Is_i, si_d, s)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tmp*pa(kdi, iS_d, Ss_i, S)*cc(dii, did, is_i, is_d, i*s, N_Qd, N_Qi)
        - tpm*pa(kid, Is_i, sI_d, s)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
        - gm*Is_i
        + dp*Ir_i
        - lm*Is_i
        + lm*si_i
        - om*Is_i
        + chi*pa(kii, sS_i, SI_i, S)*cc(iii, iid, sI_i, sI_d, s*I, N_Qd, N_Qi)
        + chi*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
        + chi*pa(kii, rS_i, SI_i, S)*cc(iii, iid, rI_i, rI_d, r*I, N_Qd, N_Qi)
        - chi*Is_i
        - chi*pa(kii, sI_i, Is_i, I)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
        - chi*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        - chi*pa(kii, rI_i, Is_i, I)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*SI_i
        + mu*pa(kii, IS_i, SI_i, S)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
        - mu*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        - qdi*tpm*Is_i;

      Ii_d_t = + tmm*pa(kdd, IS_d, Si_d, S)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tpp*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tmp*Si_d
        + tmp*pa(kdd, iS_d, Si_d, S)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + tpm*Is_d
        + tpm*pa(kdd, Is_d, sI_d, s)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
        - gm*Ii_d
        - gp*Ii_d
        - lm*Ii_d
        + lm*ii_d
        + om*II_d
        - om*Ii_d
        + chi*pa(kdi, II_d, Is_i, I)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - chi*pa(kid, sI_i, Ii_d, I)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
        + chi*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - chi*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + chi*pa(kdi, II_d, Ir_i, I)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - chi*pa(kid, rI_i, Ii_d, I)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + mu*pa(kdi, II_d, II_i, I)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kid, II_i, Ii_d, I)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - mu*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        - qid*chi*Ii_d
        + qid*mu*II_d
        - qid*mu*Ii_d;

      Ii_i_t = + tmm*pa(kdi, IS_d, Si_i, S)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tpp*pa(kid, Is_i, si_d, s)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + tmp*pa(kdi, iS_d, Si_i, S)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + tpm*pa(kid, Is_i, sI_d, s)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
        - gm*Ii_i
        - gp*Ii_i
        - lm*Ii_i
        + lm*ii_i
        + om*II_i
        - om*Ii_i
        + chi*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        - chi*pa(kii, sI_i, Ii_i, I)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
        - chi*Ii_i
        + chi*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - chi*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + chi*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - chi*pa(kii, rI_i, Ii_i, I)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + mu*II_i
        + mu*pa(kii, II_i, II_i, I)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - mu*Ii_i
        + mu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - mu*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + qdi*tmp*Si_i
        + qdi*tpm*Is_i;

      Ir_d_t = + tmm*pa(kdd, IS_d, Sr_d, S)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + tmp*pa(kdd, iS_d, Sr_d, S)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - gm*Ir_d
        + gp*Ii_d
        - dp*Ir_d
        + lm*ir_d
        - lm*Ir_d
        - om*Ir_d
        - chi*pa(kid, sI_i, Ir_d, I)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
        - chi*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - chi*pa(kid, rI_i, Ir_d, I)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kdi, IR_d, Rs_i, R)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + chi*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + chi*pa(kdi, IR_d, Rr_i, R)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - mu*pa(kid, II_i, Ir_d, I)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kdi, IR_d, RI_i, R)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - qid*chi*Ir_d
        + qid*mu*IR_d;

      Ir_i_t = + tmm*pa(kdi, IS_d, Sr_i, S)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + tmp*pa(kdi, iS_d, Sr_i, S)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - gm*Ir_i
        + gp*Ii_i
        - dp*Ir_i
        + lm*ir_i
        - lm*Ir_i
        - om*Ir_i
        - chi*pa(kii, sI_i, Ir_i, I)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
        - chi*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - chi*Ir_i
        - chi*pa(kii, rI_i, Ir_i, I)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + chi*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + chi*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        - mu*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*IR_i
        + mu*pa(kii, IR_i, RI_i, R)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
        - mu*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        ;
            
      RR_d_t = 2*(
                  + gm*IR_d
                  - dm*RR_d
                  + lm*Rr_d
                  - chi*pa(kdi, RR_d, Rs_i, R)*cc(dii, did, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
                  - chi*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                  - chi*pa(kdi, RR_d, Rr_i, R)*cc(dii, did, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
                  - mu*pa(kdi, RR_d, RI_i, R)*cc(dii, did, RI_i, RI_d, R*I, N_Qd, N_Qi)
                  - mu*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                  );

      RR_i_t = 2*(
                  + gm*IR_i
                  - dm*RR_i
                  + lm*Rr_i
                  - chi*pa(kii, RR_i, Rs_i, R)*cc(iii, iid, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
                  - chi*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                  - chi*pa(kii, RR_i, Rr_i, R)*cc(iii, iid, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
                  - mu*pa(kii, RR_i, RI_i, R)*cc(iii, iid, RI_i, RI_d, R*I, N_Qd, N_Qi)
                  - mu*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                  );

      Rs_d_t = - tpp*pa(kdd, Rs_d, si_d, s)*cc(ddi, ddd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        - tpm*pa(kdd, Rs_d, sI_d, s)*cc(ddi, ddd, RI_i, RI_d, R*I, N_Qd, N_Qi)
        + gm*Is_d
        - dm*Rs_d
        + dp*Rr_d
        - lm*Rs_d
        + lm*sr_d
        + chi*pa(kid, sS_i, SR_d, S)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
        + chi*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + chi*pa(kid, rS_i, SR_d, S)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kid, sR_i, Rs_d, R)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
        - chi*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        - chi*pa(kid, rR_i, Rs_d, R)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kid, IS_i, SR_d, S)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kid, IR_i, Rs_d, R)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        - qid*chi*Rs_d;

      Rs_i_t = - tpp*pa(kid, Rs_i, si_d, s)*cc(idi, idd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        - tpm*pa(kid, Rs_i, sI_d, s)*cc(idi, idd, RI_i, RI_d, R*I, N_Qd, N_Qi)
        + gm*Is_i
        - dm*Rs_i
        + dp*Rr_i
        - lm*Rs_i
        + lm*sr_i
        + chi*pa(kii, sS_i, SR_i, S)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
        + chi*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + chi*pa(kii, rS_i, SR_i, S)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*Rs_i
        - chi*pa(kii, sR_i, Rs_i, R)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
        - chi*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        - chi*pa(kii, rR_i, Rs_i, R)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kii, IS_i, SR_i, S)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        ;
            
      Ri_d_t = + tpp*pa(kdd, Rs_d, si_d, s)*cc(ddi, ddd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        + tpm*pa(kdd, Rs_d, sI_d, s)*cc(ddi, ddd, RI_i, RI_d, R*I, N_Qd, N_Qi)
        + gm*Ii_d
        - gp*Ri_d
        - dm*Ri_d
        - lm*Ri_d
        + lm*ir_d
        + om*IR_d
        + chi*pa(kid, sI_i, IR_d, I)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
        + chi*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + chi*pa(kid, rI_i, IR_d, I)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kid, sR_i, Ri_d, R)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
        - chi*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        - chi*pa(kid, rR_i, Ri_d, R)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + mu*pa(kid, II_i, IR_d, I)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kid, IR_i, Ri_d, R)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        - qid*chi*Ri_d
        - qid*mu*Ri_d;

      Ri_i_t = + tpp*pa(kid, Rs_i, si_d, s)*cc(idi, idd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        + tpm*pa(kid, Rs_i, sI_d, s)*cc(idi, idd, RI_i, RI_d, R*I, N_Qd, N_Qi)
        + gm*Ii_i
        - gp*Ri_i
        - dm*Ri_i
        - lm*Ri_i
        + lm*ir_i
        + om*IR_i
        + chi*pa(kii, sI_i, IR_i, I)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
        + chi*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        + chi*pa(kii, rI_i, IR_i, I)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
        - chi*pa(kii, sR_i, Ri_i, R)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
        - chi*Ri_i
        - chi*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        - chi*pa(kii, rR_i, Ri_i, R)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + mu*pa(kii, II_i, IR_i, I)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
        - mu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
        - mu*Ri_i
        - mu*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        ;
            
      Rr_d_t = + gm*Ir_d
        + gp*Ri_d
        - dm*Rr_d
        - dp*Rr_d
        - lm*Rr_d
        + lm*rr_d
        + chi*pa(kdi, RR_d, Rs_i, R)*cc(dii, did, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
        - chi*pa(kid, sR_i, Rr_d, R)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
        + chi*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        - chi*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + chi*pa(kdi, RR_d, Rr_i, R)*cc(dii, did, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
        - chi*pa(kid, rR_i, Rr_d, R)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + mu*pa(kdi, RR_d, RI_i, R)*cc(dii, did, RI_i, RI_d, R*I, N_Qd, N_Qi)
        - mu*pa(kid, IR_i, Rr_d, R)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        - mu*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - qid*chi*Rr_d;

      Rr_i_t = + gm*Ir_i
        + gp*Ri_i
        - dm*Rr_i
        - dp*Rr_i
        - lm*Rr_i
        + lm*rr_i
        + chi*pa(kii, RR_i, Rs_i, R)*cc(iii, iid, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
        - chi*pa(kii, sR_i, Rr_i, R)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
        + chi*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        - chi*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - chi*Rr_i
        + chi*pa(kii, RR_i, Rr_i, R)*cc(iii, iid, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
        - chi*pa(kii, rR_i, Rr_i, R)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + mu*pa(kii, RR_i, RI_i, R)*cc(iii, iid, RI_i, RI_d, R*I, N_Qd, N_Qi)
        - mu*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
        - mu*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        ;
            
      ss_d_t = 2*(
                  - tpp*pa(kdd, ss_d, si_d, s)*cc(ddi, ddd, si_i, si_d, s*i, N_Qd, N_Qi)
                  - tpm*pa(kdd, ss_d, sI_d, s)*cc(ddi, ddd, sI_i, sI_d, s*I, N_Qd, N_Qi)
                  + dp*sr_d
                  - lm*ss_d
                  + chi*pa(kid, sS_i, Ss_d, S)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
                  + chi*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
                  + chi*pa(kid, rS_i, Ss_d, S)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
                  + mu*pa(kid, IS_i, Ss_d, S)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
                  + mu*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
                  + qid*chi*Ss_d
                  );

      ss_i_t = 2*(
                  - tpp*pa(kid, ss_i, si_d, s)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
                  - tpm*pa(kid, ss_i, sI_d, s)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
                  + dp*sr_i
                  - lm*ss_i
                  + chi*Ss_i
                  + chi*pa(kii, sS_i, Ss_i, S)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
                  + chi*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
                  + chi*pa(kii, rS_i, Ss_i, S)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
                  + mu*pa(kii, IS_i, Ss_i, S)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
                  + mu*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
                  );

      si_d_t = - tpp*si_d
        + tpp*pa(kdd, ss_d, si_d, s)*cc(ddi, ddd, si_i, si_d, s*i, N_Qd, N_Qi)
        - tpp*pa(kdd, is_d, si_d, s)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + tpm*pa(kdd, ss_d, sI_d, s)*cc(ddi, ddd, sI_i, sI_d, s*I, N_Qd, N_Qi)
        - tpm*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - gp*si_d
        + dp*ir_d
        - lm*si_d
        - lm*si_d
        + om*Is_d
        + chi*pa(kid, sS_i, Si_d, S)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
        + chi*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + chi*pa(kid, rS_i, Si_d, S)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + chi*pa(kid, sI_i, Is_d, I)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
        + chi*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        + chi*pa(kid, rI_i, Is_d, I)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kid, IS_i, Si_d, S)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kid, II_i, Is_d, I)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + mu*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        + qid*chi*Si_d
        + qid*chi*Is_d
        + qid*mu*Si_d;

      si_i_t = + tpp*pa(kid, ss_i, si_d, s)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
        - tpp*pa(kdi, is_d, si_i, s)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + tpm*pa(kid, ss_i, sI_d, s)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
        - tpm*pa(kdi, Is_d, si_i, s)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        - gp*si_i
        + dp*ir_i
        - lm*si_i
        - lm*si_i
        + om*Is_i
        + chi*pa(kii, sS_i, Si_i, S)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
        + chi*Si_i
        + chi*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + chi*pa(kii, rS_i, Si_i, S)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + chi*Is_i
        + chi*pa(kii, sI_i, Is_i, I)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
        + chi*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        + chi*pa(kii, rI_i, Is_i, I)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kii, IS_i, Si_i, S)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*Si_i
        + mu*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + mu*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        - qdi*tpp*si_i;

      sr_d_t = - tpp*pa(kdd, is_d, sr_d, s)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - tpm*pa(kdd, Is_d, sr_d, s)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + gp*si_d
        - dp*sr_d
        + dp*rr_d
        - lm*sr_d
        - lm*sr_d
        + chi*pa(kid, sS_i, Sr_d, S)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
        + chi*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + chi*pa(kid, rS_i, Sr_d, S)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kid, sR_i, Rs_d, R)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
        + chi*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        + chi*pa(kid, rR_i, Rs_d, R)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kid, IS_i, Sr_d, S)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kid, IR_i, Rs_d, R)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
        + qid*chi*Sr_d
        + qid*chi*Rs_d;

      sr_i_t = - tpp*pa(kdi, is_d, sr_i, s)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
        - tpm*pa(kdi, Is_d, sr_i, s)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + gp*si_i
        - dp*sr_i
        + dp*rr_i
        - lm*sr_i
        - lm*sr_i
        + chi*pa(kii, sS_i, Sr_i, S)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
        + chi*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + chi*Sr_i
        + chi*pa(kii, rS_i, Sr_i, S)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*Rs_i
        + chi*pa(kii, sR_i, Rs_i, R)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
        + chi*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        + chi*pa(kii, rR_i, Rs_i, R)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
        + mu*pa(kii, IS_i, Sr_i, S)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
        + mu*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
        ;
            
      ii_d_t = 2*(
                  + tpp*si_d
                  + tpp*pa(kdd, is_d, si_d, s)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
                  + tpm*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  - gp*ii_d
                  - lm*ii_d
                  + om*Ii_d
                  + chi*pa(kid, sI_i, Ii_d, I)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
                  + chi*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
                  + chi*pa(kid, rI_i, Ii_d, I)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
                  + mu*pa(kid, II_i, Ii_d, I)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  + mu*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
                  + qid*chi*Ii_d
                  + qid*mu*Ii_d
                  );

      ii_i_t = 2*(
                  + tpp*pa(kdi, is_d, si_i, s)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
                  + tpm*pa(kdi, Is_d, si_i, s)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  - gp*ii_i
                  - lm*ii_i
                  + om*Ii_i
                  + chi*pa(kii, sI_i, Ii_i, I)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
                  + chi*Ii_i
                  + chi*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
                  + chi*pa(kii, rI_i, Ii_i, I)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
                  + mu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                  + mu*Ii_i
                  + mu*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
                  + qdi*tpp*si_i
                  );

      ir_d_t = + tpp*pa(kdd, is_d, sr_d, s)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + tpm*pa(kdd, Is_d, sr_d, s)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + gp*ii_d
        - gp*ir_d
        - dp*ir_d
        - lm*ir_d
        - lm*ir_d
        + om*Ir_d
        + chi*pa(kid, sI_i, Ir_d, I)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
        + chi*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + chi*pa(kid, rI_i, Ir_d, I)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kid, sR_i, Ri_d, R)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
        + chi*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + chi*pa(kid, rR_i, Ri_d, R)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + mu*pa(kid, II_i, Ir_d, I)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kid, IR_i, Ri_d, R)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + qid*chi*Ir_d
        + qid*chi*Ri_d
        + qid*mu*Ri_d;

      ir_i_t = + tpp*pa(kdi, is_d, sr_i, s)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + tpm*pa(kdi, Is_d, sr_i, s)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + gp*ii_i
        - gp*ir_i
        - dp*ir_i
        - lm*ir_i
        - lm*ir_i
        + om*Ir_i
        + chi*pa(kii, sI_i, Ir_i, I)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
        + chi*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + chi*Ir_i
        + chi*pa(kii, rI_i, Ir_i, I)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
        + chi*pa(kii, sR_i, Ri_i, R)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
        + chi*Ri_i
        + chi*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        + chi*pa(kii, rR_i, Ri_i, R)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
        + mu*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
        + mu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
        + mu*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
        + mu*Ri_i
        + mu*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
        ;
            
      rr_d_t = 2*(
                  + gp*ir_d
                  - dp*rr_d
                  - lm*rr_d
                  + chi*pa(kid, sR_i, Rr_d, R)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
                  + chi*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
                  + chi*pa(kid, rR_i, Rr_d, R)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
                  + mu*pa(kid, IR_i, Rr_d, R)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                  + mu*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
                  + qid*chi*Rr_d
                  );

      rr_i_t = 2*(
                  + gp*ir_i
                  - dp*rr_i
                  - lm*rr_i
                  + chi*pa(kii, sR_i, Rr_i, R)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
                  + chi*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
                  + chi*Rr_i
                  + chi*pa(kii, rR_i, Rr_i, R)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
                  + mu*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                  + mu*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
                  );
      
      return GSL_SUCCESS;         
    }
    
  }; // Eqs

  //------------------------------------------------------------
  
} // namespace InfoSIRS

//------------------------------------------------------------

#endif // INFO_SIRS_PA_H
