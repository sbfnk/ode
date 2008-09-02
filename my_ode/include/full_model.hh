#ifndef FULL_MODEL_H
#define FULL_MODEL_H

//------------------------------------------------------------

#include <iostream>
#include <gsl/gsl_errno.h>

#include "pa_macros.hh"

//------------------------------------------------------------

struct FullModelParams
{
  FullModelParams() : nvars(2) {};
      
  // No. of equations
  unsigned int nvars;
      
  // rates
  double gamma[2], delta[2];
  double beta[2][2], alpha, nu, lambda;
  double omega;
  double sigma;
      
  // Graph properties
  double Qd, Qi, Qb, N;
  double qdi, qid;
      
  // clustering
  double C[3][3][3];

  // claculate Qd, Qi and Qb from i.c.
  void init_Q(double* rhs, bool verbose);
      
  // overloading operator<<
  friend std::ostream& operator <<
    (std::ostream& os, const FullModelParams& x);
      
  // overloading operator=
  void operator= (const FullModelParams& x);
      
}; // FullModelParams

//------------------------------------------------------------

void FullModelParams::init_Q(double* rhs, bool verbose)
{   
  // calculate Qd and Qi in from i.c.
  if (nvars > 6) {
      
    // calculate total number of disease pairs
    if (verbose) std::cout << "... calculating Qd from i.c.";
      
    double dPairs =
      SS_d_t + SI_d_t + SR_d_t + Ss_d_t + Si_d_t + Sr_d_t +
      II_d_t + IR_d_t + Is_d_t + Ii_d_t + Ir_d_t +
      RR_d_t + Rs_d_t + Ri_d_t + Rr_d_t +
      ss_d_t + si_d_t + sr_d_t +
      ii_d_t + ir_d_t +
      rr_d_t;
      
    Qd = (2 * dPairs / N);
      
    if (verbose) std::cout << " ... done "
                           << "(Qd=" << Qd << ")" << std::endl;
      
    // calculate total number of info pairs
    if (verbose) std::cout << "... calculating Qi from i.c.";
      
    double iPairs =
      SS_i_t + SI_i_t + SR_i_t + Ss_i_t + Si_i_t + Sr_i_t +
      II_i_t + IR_i_t + Is_i_t + Ii_i_t + Ir_i_t +
      RR_i_t + Rs_i_t + Ri_i_t + Rr_i_t +
      ss_i_t + si_i_t + sr_i_t +
      ii_i_t + ir_i_t +
      rr_i_t;
      
    Qi = (2 * iPairs / N);      
    if (verbose) std::cout << " ... done "
                           << "(Qi=" << Qi << ")" << std::endl;
      
  } else {
    if (verbose) std::cout << "... initialize Qd and Qi to default values";
    //Qd = 3;
    //Qi = 3;
    if (verbose) std::cout << " ... done\n";
  }

  // dib model
  if (nvars > 48) {
      
    // calculate total number of parallel pairs
    if (verbose) std::cout << "... calculating Qb from i.c.";
      
    double bPairs =
      SS_b_t + SI_b_t + SR_b_t + Ss_b_t + Si_b_t + Sr_b_t +
      II_b_t + IR_b_t + Is_b_t + Ii_b_t + Ir_b_t +
      RR_b_t + Rs_b_t + Ri_b_t + Rr_b_t +
      ss_b_t + si_b_t + sr_b_t +
      ii_b_t + ir_b_t +
      rr_b_t;
      
    Qb = (2 * bPairs / N);      
    if (verbose) std::cout << " ... done "
                           << "(Qb=" << Qb << ")" << std::endl;
      
  } else {
    if (verbose) std::cout << "... initialize Qb to default values";
    //Qb = 3;
    if (verbose) std::cout << " ... done\n";
  }
  
  if (sigma > 0) {
    beta[1][0] = sigma*beta[0][0];
    beta[1][1] = sigma*beta[0][1];
    if (verbose) {
      std::cout << "sigma given, setting beta+-=" << beta[1][0] << " and "
                << "beta++= " << beta[1][1] << std::endl;
    }
  }
   
   
  // ADD INIT_C !!!!!!
   
}

//------------------------------------------------------------

void FullModelParams::operator= (const FullModelParams& x)
{
  int i,j,k;
   
  nvars = x.nvars;
   
  for(i = 0; i < 2; i++)
    {
      gamma[i] = x.gamma[i];
      delta[i] = x.delta[i];
      for(j = 0; j < 2; j++)
        beta[i][j] = x.beta[i][j];
    }
  alpha  = x.alpha;
  nu     = x.nu;
  lambda = x.lambda;
  omega  = x.omega;
   
  N   = x.N;
  Qd  = x.Qd;
  Qi  = x.Qi;
  Qb = x.Qb;
  qdi = x.qdi;
  qid = x.qid;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
        C[i][j][k] = x.C[i][j][k];

} // operator=

//------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, const FullModelParams& x)
{
  os << std::endl
     << "Model parameters:" << std::endl
     << "-----------------" << std::endl
     << "nvars =    " << x.nvars << std::endl
     << "beta--   = " << x.beta[0][0] << std::endl
     << "beta+-   = " << x.beta[1][0] << std::endl
     << "beta-+   = " << x.beta[0][1] << std::endl
     << "beta++   = " << x.beta[1][1] << std::endl
     << "gamma-   = " << x.gamma[0] << std::endl
     << "gamma+   = " << x.gamma[1] << std::endl
     << "delta-   = " << x.delta[0] << std::endl
     << "delta+   = " << x.delta[1] << std::endl
     << "alpha    = " << x.alpha << std::endl
     << "nu       = " << x.nu << std::endl
     << "lambda   = " << x.lambda << std::endl
     << "omega    = " << x.omega << std::endl
     << std::endl
     << "N        = " << x.N << std::endl
     << "Qd       = " << x.Qd << std::endl
     << "Qi       = " << x.Qi << std::endl
     << "Qb      = " << x.Qb << std::endl
     << "qdi      = " << x.qdi << std::endl
     << "qid      = " << x.qid << std::endl
     << std::endl
     << "R_0 d -- = " << (x.beta[0][0])/(x.gamma[0]) << std::endl
     << "R_0 d -+ = " << (x.beta[1][0])/(x.gamma[1]) << std::endl
     << "R_0 d +- = " << (x.beta[0][1])/(x.gamma[0]) << std::endl
     << "R_0 d ++ = " << (x.beta[1][1])/(x.gamma[1]) << std::endl
     << "R_0 i    = " << (x.alpha)/(x.lambda) << std::endl
     << std::endl
     << "C_ddi  = " << x.C[1][1][0] << std::endl
     << "C_ddd  = " << x.C[1][1][1] << std::endl
     << "C_dii  = " << x.C[1][0][0] << std::endl
     << "C_did  = " << x.C[1][0][1] << std::endl
     << "C_idi  = " << x.C[0][1][0] << std::endl
     << "C_idd  = " << x.C[0][1][1] << std::endl
     << "C_iii  = " << x.C[0][0][0] << std::endl
     << "C_iid  = " << x.C[0][0][1] << std::endl
     << std::endl;
   
  return os;
   
} // operator<<

//------------------------------------------------------------

struct FullModelMeanField
{
  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {         
    FullModelParams p = *(static_cast<FullModelParams*>(params));

    // local readable short variables
    double bd=p.beta[0][0], gd=p.gamma[0], dd=p.delta[0];
    double bi=p.beta[1][1], gi=p.gamma[1], di=p.delta[1];
    double bm1=p.beta[0][1], bm2=p.beta[1][0], al=p.alpha, nu=p.nu;
    double lm=p.lambda, om=p.omega;
    double N=p.N, invN=1./N;
    double Qd=p.Qd, Qi=p.Qi;
         
    //bd*=Qd; bi*=Qd; bm1*=Qd; bm2*=Qd;
    //al*=Qi; nu*=Qi;
         
         
    rhs[0] = -invN * ( bd*y[1] + bm1*y[4] ) * y[0] + dd * y[2]
      - al * invN * ( y[3] + y[4] + y[5] ) * y[0] + lm * y[3]
      - nu * invN * ( y[1] + y[4] ) * y[0];
         
    rhs[1] =  invN * ( bd*y[1] + bm1*y[4] ) * y[0] - gd * y[1]
      - al * invN * ( y[3] + y[4] + y[5] ) * y[1] + lm * y[4]
      - nu * invN * ( y[1] + y[4] ) * y[1] - om * y[1];
         
    rhs[2] =  gd * y[1] - dd * y[2] + lm * y[5]
      - al * invN * ( y[3] + y[4] + y[5] ) * y[2]
      - nu * invN * ( y[1] + y[4] ) * y[2];
         
    rhs[3] = -invN * ( bi*y[4] + bm2*y[1] ) * y[3] + di * y[5]
      + al * invN * ( y[3] + y[4] + y[5] ) * y[0] - lm * y[3]
      + nu * invN * ( y[1] + y[4] ) * y[0];
         
    rhs[4] =  invN * ( bi*y[4] + bm2*y[1] ) * y[3] - gi * y[4]
      + al * invN * ( y[3] + y[4] + y[5] ) * y[1] - lm * y[4]
      + nu * invN * ( y[1] + y[4] ) * y[1] + om * y[1];
         
    rhs[5] =  gi * y[4] - di * y[5] - lm * y[5]
      + al * invN * ( y[3] + y[4] + y[5] ) * y[2]
      + nu * invN * ( y[1] + y[4] ) * y[2];
         
    return GSL_SUCCESS;         
  }
      
}; // FullModelMeanField

//------------------------------------------------------------

struct FullModelPairApprox1
{
  // pair approx      
  static double pa(double kk, double ij, double jk, double j)
  {
    const double eps = 1e-8;
         
    if(j < eps)
      return 0.0;
    else
      return kk*((ij*jk)/j);
  }
      
  // clustering correction
  static double  cc(double C_i, double C_d, double ki_i, double ki_d,
                    double ki, double N_Qd, double N_Qi)
  {
    if(ki == 0.0)
      return ((1.0 - C_i - C_d));
    else
      return ((1.0 - C_i - C_d) + C_i*N_Qi*ki_i/ki + C_d*N_Qd*ki_d/ki);
  }           
      
  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {
    FullModelParams p = *(static_cast<FullModelParams*>(params));
         
    // local readable short variables
    double bd=p.beta[0][0], gd=p.gamma[0], dd=p.delta[0];
    double bi=p.beta[1][1], gi=p.gamma[1], di=p.delta[1];
    double b1=p.beta[0][1], b2=p.beta[1][0], al=p.alpha, nu=p.nu;
    double lm=p.lambda, om=p.omega;
    double Qd=p.Qd, Qi=p.Qi, N=p.N;
    double N_Qd=N/Qd, N_Qi=N/Qi;

    //bd/=Qd; bi/=Qd; b1/=Qd; b2/=Qd;
    //al/=Qi; nu/=Qi;
         
    // Clustering corrections
    double idi=p.C[0][1][0], idd=p.C[0][1][1];
    double dii=p.C[1][0][0], did=p.C[1][0][1];
    double ddi=p.C[1][1][0], ddd=p.C[1][1][1];
    double iii=p.C[0][0][0], iid=p.C[0][0][1];
    double kdd=(Qd-1.0)/Qd, kii=(Qi-1.0)/Qi;
    double kdi=1.0, kid=1.0;
         
    double qdi=p.qdi, qid=p.qid;
         

    S_t = - bd*SI_d
      - b1*Si_d
      + dd*R
      + lm*s
      - al*Ss_i
      - al*Si_i
      - al*Sr_i
      - nu*SI_i
      - nu*Si_i;

    I_t = + bd*SI_d
      + b1*Si_d
      - gd*I
      + lm*i
      - om*I
      - al*Is_i
      - al*Ii_i
      - al*Ir_i
      - nu*II_i
      - nu*Ii_i;
 
    R_t = + gd*I
      - dd*R
      + lm*r
      - al*Rs_i
      - al*Ri_i
      - al*Rr_i
      - nu*IR_i
      - nu*Ri_i;
 
    s_t = - bi*si_d
      - b2*Is_d
      + di*r
      - lm*s
      + al*Ss_i
      + al*Si_i
      + al*Sr_i
      + nu*SI_i
      + nu*Si_i;
 
    i_t = + bi*si_d
      + b2*Is_d
      - gi*i
      - lm*i
      + om*I
      + al*Is_i
      + al*Ii_i
      + al*Ir_i
      + nu*II_i
      + nu*Ii_i;
 
    r_t = + gi*i
      - di*r
      - lm*r
      + al*Rs_i
      + al*Ri_i
      + al*Rr_i
      + nu*IR_i
      + nu*Ri_i;
 
    SS_d_t = 2*(
                - bd*pa(kdd, SS_d, SI_d, S)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
                - b1*pa(kdd, SS_d, Si_d, S)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
                + dd*SR_d
                + lm*Ss_d
                - al*pa(kdi, SS_d, Ss_i, S)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
                - al*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
                - al*pa(kdi, SS_d, Sr_i, S)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
                - nu*pa(kdi, SS_d, SI_i, S)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
                - nu*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
                );

    SS_i_t = 2*(
                - bd*pa(kid, SS_i, SI_d, S)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
                - b1*pa(kid, SS_i, Si_d, S)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
                + dd*SR_i
                + lm*Ss_i
                - al*pa(kii, SS_i, Ss_i, S)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
                - al*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
                - al*pa(kii, SS_i, Sr_i, S)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
                - nu*pa(kii, SS_i, SI_i, S)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
                - nu*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
                );

    SI_d_t = - bd*SI_d
      + bd*pa(kdd, SS_d, SI_d, S)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - bd*pa(kdd, IS_d, SI_d, S)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
      + b1*pa(kdd, SS_d, Si_d, S)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - b1*pa(kdd, iS_d, SI_d, S)*cc(ddi, ddd, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - gd*SI_d
      + dd*IR_d
      + lm*Is_d
      + lm*Si_d
      - om*SI_d
      - al*pa(kid, sS_i, SI_d, S)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
      - al*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - al*pa(kid, rS_i, SI_d, S)*cc(idi, idd, rI_i, rI_d, r*I, N_Qd, N_Qi)
      - al*pa(kdi, SI_d, Is_i, I)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      - al*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - al*pa(kdi, SI_d, Ir_i, I)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kid, IS_i, SI_d, S)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kdi, SI_d, II_i, I)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - nu*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - qid*nu*SI_d;

    SI_i_t = + bd*pa(kid, SS_i, SI_d, S)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - bd*pa(kdi, IS_d, SI_i, S)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
      + b1*pa(kid, SS_i, Si_d, S)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - b1*pa(kdi, iS_d, SI_i, S)*cc(dii, did, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - gd*SI_i
      + dd*IR_i
      + lm*Is_i
      + lm*Si_i
      - om*SI_i
      - al*pa(kii, sS_i, SI_i, S)*cc(iii, iid, sI_i, sI_d, s*I, N_Qd, N_Qi)
      - al*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - al*pa(kii, rS_i, SI_i, S)*cc(iii, iid, rI_i, rI_d, r*I, N_Qd, N_Qi)
      - al*pa(kii, SI_i, Is_i, I)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      - al*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - al*pa(kii, SI_i, Ir_i, I)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*SI_i
      - nu*pa(kii, IS_i, SI_i, S)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kii, SI_i, II_i, I)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - nu*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - qdi*bd*SI_i;

    SR_d_t = - bd*pa(kdd, IS_d, SR_d, S)*cc(ddi, ddd, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - b1*pa(kdd, iS_d, SR_d, S)*cc(ddi, ddd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + gd*SI_d
      - dd*SR_d
      + dd*RR_d
      + lm*Rs_d
      + lm*Sr_d
      - al*pa(kid, sS_i, SR_d, S)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
      - al*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - al*pa(kid, rS_i, SR_d, S)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kdi, SR_d, Rs_i, R)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      - al*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - al*pa(kdi, SR_d, Rr_i, R)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kid, IS_i, SR_d, S)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kdi, SR_d, RI_i, R)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      ;

    SR_i_t = - bd*pa(kdi, IS_d, SR_i, S)*cc(dii, did, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - b1*pa(kdi, iS_d, SR_i, S)*cc(dii, did, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + gd*SI_i
      - dd*SR_i
      + dd*RR_i
      + lm*Rs_i
      + lm*Sr_i
      - al*pa(kii, sS_i, SR_i, S)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
      - al*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - al*pa(kii, rS_i, SR_i, S)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kii, SR_i, Rs_i, R)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      - al*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - al*pa(kii, SR_i, Rr_i, R)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kii, IS_i, SR_i, S)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kii, SR_i, RI_i, R)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      ;
            
    Ss_d_t = - bd*pa(kdd, IS_d, Ss_d, S)*cc(ddi, ddd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - bi*pa(kdd, Ss_d, si_d, s)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - b1*pa(kdd, iS_d, Ss_d, S)*cc(ddi, ddd, is_i, is_d, i*s, N_Qd, N_Qi)
      - b2*pa(kdd, Ss_d, sI_d, s)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
      + dd*Rs_d
      + di*Sr_d
      - lm*Ss_d
      + lm*ss_d
      + al*pa(kdi, SS_d, Ss_i, S)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      - al*pa(kid, sS_i, Ss_d, S)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
      + al*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - al*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      + al*pa(kdi, SS_d, Sr_i, S)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - al*pa(kid, rS_i, Ss_d, S)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kdi, SS_d, SI_i, S)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kid, IS_i, Ss_d, S)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kdi, SS_d, Si_i, S)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - nu*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      - qid*al*Ss_d;
               
    Ss_i_t = - bd*pa(kdi, IS_d, Ss_i, S)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - bi*pa(kid, Ss_i, si_d, s)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - b1*pa(kdi, iS_d, Ss_i, S)*cc(dii, did, is_i, is_d, i*s, N_Qd, N_Qi)
      - b2*pa(kid, Ss_i, sI_d, s)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
      + dd*Rs_i
      + di*Sr_i
      - lm*Ss_i
      + lm*ss_i
      - al*Ss_i
      + al*pa(kii, SS_i, Ss_i, S)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      - al*pa(kii, sS_i, Ss_i, S)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
      + al*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - al*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      + al*pa(kii, SS_i, Sr_i, S)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - al*pa(kii, rS_i, Ss_i, S)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kii, SS_i, SI_i, S)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kii, IS_i, Ss_i, S)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kii, SS_i, Si_i, S)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - nu*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      ;
            
    Si_d_t = - bd*pa(kdd, IS_d, Si_d, S)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + bi*pa(kdd, Ss_d, si_d, s)*cc(ddi, ddd, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - b1*Si_d
      - b1*pa(kdd, iS_d, Si_d, S)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + b2*pa(kdd, Ss_d, sI_d, s)*cc(ddi, ddd, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - gi*Si_d
      + dd*Ri_d
      + lm*si_d
      - lm*Si_d
      + om*SI_d
      - al*pa(kid, sS_i, Si_d, S)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
      - al*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      - al*pa(kid, rS_i, Si_d, S)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + al*pa(kdi, SI_d, Is_i, I)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      + al*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      + al*pa(kdi, SI_d, Ir_i, I)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kid, IS_i, Si_d, S)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kdi, SI_d, II_i, I)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + nu*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - qid*al*Si_d
      - qid*nu*Si_d
      ;
            
    Si_i_t = - bd*pa(kdi, IS_d, Si_i, S)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + bi*pa(kid, Ss_i, si_d, s)*cc(idi, idd, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - b1*pa(kdi, iS_d, Si_i, S)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + b2*pa(kid, Ss_i, sI_d, s)*cc(idi, idd, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - gi*Si_i
      + dd*Ri_i
      + lm*si_i
      - lm*Si_i
      + om*SI_i
      - al*pa(kii, sS_i, Si_i, S)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
      - al*Si_i
      - al*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      - al*pa(kii, rS_i, Si_i, S)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + al*pa(kii, SI_i, Is_i, I)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      + al*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      + al*pa(kii, SI_i, Ir_i, I)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kii, IS_i, Si_i, S)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kii, SI_i, II_i, I)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*Si_i
      - nu*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + nu*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - qdi*b1*Si_i;

    Sr_d_t = - bd*pa(kdd, IS_d, Sr_d, S)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - b1*pa(kdd, iS_d, Sr_d, S)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + gi*Si_d
      + dd*Rr_d
      - di*Sr_d
      + lm*sr_d
      - lm*Sr_d
      - al*pa(kid, sS_i, Sr_d, S)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
      - al*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - al*pa(kid, rS_i, Sr_d, S)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kdi, SR_d, Rs_i, R)*cc(dii, did, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      + al*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      + al*pa(kdi, SR_d, Rr_i, R)*cc(dii, did, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kid, IS_i, Sr_d, S)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kdi, SR_d, RI_i, R)*cc(dii, did, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, Si_i, Si_d, S*i, N_Qd, N_Qi)
      - qid*al*Sr_d;

    Sr_i_t = - bd*pa(kdi, IS_d, Sr_i, S)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - b1*pa(kdi, iS_d, Sr_i, S)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + gi*Si_i
      + dd*Rr_i
      - di*Sr_i
      + lm*sr_i
      - lm*Sr_i
      - al*pa(kii, sS_i, Sr_i, S)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
      - al*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - al*Sr_i
      - al*pa(kii, rS_i, Sr_i, S)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kii, SR_i, Rs_i, R)*cc(iii, iid, Ss_i, Ss_d, S*s, N_Qd, N_Qi)
      + al*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      + al*pa(kii, SR_i, Rr_i, R)*cc(iii, iid, Sr_i, Sr_d, S*r, N_Qd, N_Qi)
      - nu*pa(kii, IS_i, Sr_i, S)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kii, SR_i, RI_i, R)*cc(iii, iid, SI_i, SI_d, S*I, N_Qd, N_Qi)
      - nu*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, Si_i, Si_d, S*i, N_Qd, N_Qi)
      ;
            
    II_d_t = 2*(
                + bd*SI_d
                + bd*pa(kdd, IS_d, SI_d, S)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
                + b1*pa(kdd, iS_d, SI_d, S)*cc(ddi, ddd, iI_i, iI_d, i*I, N_Qd, N_Qi)
                - gd*II_d
                + lm*Ii_d
                - om*II_d
                - al*pa(kdi, II_d, Is_i, I)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
                - al*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                - al*pa(kdi, II_d, Ir_i, I)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                - nu*pa(kdi, II_d, II_i, I)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
                - nu*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                - qid*nu*II_d
                );

    II_i_t = 2*(
                + bd*pa(kdi, IS_d, SI_i, S)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
                + b1*pa(kdi, iS_d, SI_i, S)*cc(dii, did, iI_i, iI_d, i*I, N_Qd, N_Qi)
                - gd*II_i
                + lm*Ii_i
                - om*II_i
                - al*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
                - al*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                - al*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                - nu*II_i
                - nu*pa(kii, II_i, II_i, I)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
                - nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                + qdi*bd*SI_i
                );

    IR_d_t = + bd*pa(kdd, IS_d, SR_d, S)*cc(ddi, ddd, IR_i, IR_d, I*R, N_Qd, N_Qi)
      + b1*pa(kdd, iS_d, SR_d, S)*cc(ddi, ddd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + gd*II_d
      - gd*IR_d
      - dd*IR_d
      + lm*Ri_d
      + lm*Ir_d
      - om*IR_d
      - al*pa(kid, sI_i, IR_d, I)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
      - al*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - al*pa(kid, rI_i, IR_d, I)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kdi, IR_d, Rs_i, R)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - al*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - al*pa(kdi, IR_d, Rr_i, R)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - nu*pa(kid, II_i, IR_d, I)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kdi, IR_d, RI_i, R)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - qid*nu*IR_d;

    IR_i_t = + bd*pa(kdi, IS_d, SR_i, S)*cc(dii, did, IR_i, IR_d, I*R, N_Qd, N_Qi)
      + b1*pa(kdi, iS_d, SR_i, S)*cc(dii, did, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + gd*II_i
      - gd*IR_i
      - dd*IR_i
      + lm*Ri_i
      + lm*Ir_i
      - om*IR_i
      - al*pa(kii, sI_i, IR_i, I)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
      - al*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - al*pa(kii, rI_i, IR_i, I)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - al*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - al*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - nu*pa(kii, II_i, IR_i, I)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*IR_i
      - nu*pa(kii, IR_i, RI_i, R)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      ;
            
    Is_d_t = + bd*pa(kdd, IS_d, Ss_d, S)*cc(ddi, ddd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - bi*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + b1*pa(kdd, iS_d, Ss_d, S)*cc(ddi, ddd, is_i, is_d, i*s, N_Qd, N_Qi)
      - b2*Is_d
      - b2*pa(kdd, Is_d, sI_d, s)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
      - gd*Is_d
      + di*Ir_d
      - lm*Is_d
      + lm*si_d
      - om*Is_d
      + al*pa(kid, sS_i, SI_d, S)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
      + al*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
      + al*pa(kid, rS_i, SI_d, S)*cc(idi, idd, rI_i, rI_d, r*I, N_Qd, N_Qi)
      - al*pa(kid, sI_i, Is_d, I)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
      - al*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      - al*pa(kid, rI_i, Is_d, I)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kid, IS_i, SI_d, S)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kid, II_i, Is_d, I)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kid, iS_i, SI_d, S)*cc(idi, idd, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - nu*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      - qid*al*Is_d
      + qid*nu*SI_d;

    Is_i_t = + bd*pa(kdi, IS_d, Ss_i, S)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - bi*pa(kid, Is_i, si_d, s)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + b1*pa(kdi, iS_d, Ss_i, S)*cc(dii, did, is_i, is_d, i*s, N_Qd, N_Qi)
      - b2*pa(kid, Is_i, sI_d, s)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
      - gd*Is_i
      + di*Ir_i
      - lm*Is_i
      + lm*si_i
      - om*Is_i
      + al*pa(kii, sS_i, SI_i, S)*cc(iii, iid, sI_i, sI_d, s*I, N_Qd, N_Qi)
      + al*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
      + al*pa(kii, rS_i, SI_i, S)*cc(iii, iid, rI_i, rI_d, r*I, N_Qd, N_Qi)
      - al*Is_i
      - al*pa(kii, sI_i, Is_i, I)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
      - al*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      - al*pa(kii, rI_i, Is_i, I)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*SI_i
      + nu*pa(kii, IS_i, SI_i, S)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iI_i, iI_d, i*I, N_Qd, N_Qi)
      - nu*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      - qdi*b2*Is_i;

    Ii_d_t = + bd*pa(kdd, IS_d, Si_d, S)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + bi*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + b1*Si_d
      + b1*pa(kdd, iS_d, Si_d, S)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + b2*Is_d
      + b2*pa(kdd, Is_d, sI_d, s)*cc(ddi, ddd, II_i, II_d, I*I, N_Qd, N_Qi)
      - gd*Ii_d
      - gi*Ii_d
      - lm*Ii_d
      + lm*ii_d
      + om*II_d
      - om*Ii_d
      + al*pa(kdi, II_d, Is_i, I)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - al*pa(kid, sI_i, Ii_d, I)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
      + al*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - al*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + al*pa(kdi, II_d, Ir_i, I)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - al*pa(kid, rI_i, Ii_d, I)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + nu*pa(kdi, II_d, II_i, I)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kid, II_i, Ii_d, I)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kdi, II_d, Ii_i, I)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - nu*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      - qid*al*Ii_d
      + qid*nu*II_d
      - qid*nu*Ii_d;

    Ii_i_t = + bd*pa(kdi, IS_d, Si_i, S)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + bi*pa(kid, Is_i, si_d, s)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + b1*pa(kdi, iS_d, Si_i, S)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + b2*pa(kid, Is_i, sI_d, s)*cc(idi, idd, II_i, II_d, I*I, N_Qd, N_Qi)
      - gd*Ii_i
      - gi*Ii_i
      - lm*Ii_i
      + lm*ii_i
      + om*II_i
      - om*Ii_i
      + al*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      - al*pa(kii, sI_i, Ii_i, I)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
      - al*Ii_i
      + al*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - al*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + al*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - al*pa(kii, rI_i, Ii_i, I)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + nu*II_i
      + nu*pa(kii, II_i, II_i, I)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - nu*Ii_i
      + nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - nu*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + qdi*b1*Si_i
      + qdi*b2*Is_i;

    Ir_d_t = + bd*pa(kdd, IS_d, Sr_d, S)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + b1*pa(kdd, iS_d, Sr_d, S)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - gd*Ir_d
      + gi*Ii_d
      - di*Ir_d
      + lm*ir_d
      - lm*Ir_d
      - om*Ir_d
      - al*pa(kid, sI_i, Ir_d, I)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
      - al*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - al*pa(kid, rI_i, Ir_d, I)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kdi, IR_d, Rs_i, R)*cc(dii, did, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + al*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + al*pa(kdi, IR_d, Rr_i, R)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - nu*pa(kid, II_i, Ir_d, I)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kdi, IR_d, RI_i, R)*cc(dii, did, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - qid*al*Ir_d
      + qid*nu*IR_d;

    Ir_i_t = + bd*pa(kdi, IS_d, Sr_i, S)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + b1*pa(kdi, iS_d, Sr_i, S)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - gd*Ir_i
      + gi*Ii_i
      - di*Ir_i
      + lm*ir_i
      - lm*Ir_i
      - om*Ir_i
      - al*pa(kii, sI_i, Ir_i, I)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
      - al*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - al*Ir_i
      - al*pa(kii, rI_i, Ir_i, I)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + al*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + al*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      - nu*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*IR_i
      + nu*pa(kii, IR_i, RI_i, R)*cc(iii, iid, II_i, II_d, I*I, N_Qd, N_Qi)
      - nu*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      ;
            
    RR_d_t = 2*(
                + gd*IR_d
                - dd*RR_d
                + lm*Rr_d
                - al*pa(kdi, RR_d, Rs_i, R)*cc(dii, did, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
                - al*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                - al*pa(kdi, RR_d, Rr_i, R)*cc(dii, did, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
                - nu*pa(kdi, RR_d, RI_i, R)*cc(dii, did, RI_i, RI_d, R*I, N_Qd, N_Qi)
                - nu*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                );

    RR_i_t = 2*(
                + gd*IR_i
                - dd*RR_i
                + lm*Rr_i
                - al*pa(kii, RR_i, Rs_i, R)*cc(iii, iid, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
                - al*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                - al*pa(kii, RR_i, Rr_i, R)*cc(iii, iid, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
                - nu*pa(kii, RR_i, RI_i, R)*cc(iii, iid, RI_i, RI_d, R*I, N_Qd, N_Qi)
                - nu*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
                );

    Rs_d_t = - bi*pa(kdd, Rs_d, si_d, s)*cc(ddi, ddd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      - b2*pa(kdd, Rs_d, sI_d, s)*cc(ddi, ddd, RI_i, RI_d, R*I, N_Qd, N_Qi)
      + gd*Is_d
      - dd*Rs_d
      + di*Rr_d
      - lm*Rs_d
      + lm*sr_d
      + al*pa(kid, sS_i, SR_d, S)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
      + al*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + al*pa(kid, rS_i, SR_d, S)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kid, sR_i, Rs_d, R)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
      - al*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      - al*pa(kid, rR_i, Rs_d, R)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kid, IS_i, SR_d, S)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kid, IR_i, Rs_d, R)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kid, iS_i, SR_d, S)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      - qid*al*Rs_d;

    Rs_i_t = - bi*pa(kid, Rs_i, si_d, s)*cc(idi, idd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      - b2*pa(kid, Rs_i, sI_d, s)*cc(idi, idd, RI_i, RI_d, R*I, N_Qd, N_Qi)
      + gd*Is_i
      - dd*Rs_i
      + di*Rr_i
      - lm*Rs_i
      + lm*sr_i
      + al*pa(kii, sS_i, SR_i, S)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
      + al*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + al*pa(kii, rS_i, SR_i, S)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*Rs_i
      - al*pa(kii, sR_i, Rs_i, R)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
      - al*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      - al*pa(kii, rR_i, Rs_i, R)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kii, IS_i, SR_i, S)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      ;
            
    Ri_d_t = + bi*pa(kdd, Rs_d, si_d, s)*cc(ddi, ddd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      + b2*pa(kdd, Rs_d, sI_d, s)*cc(ddi, ddd, RI_i, RI_d, R*I, N_Qd, N_Qi)
      + gd*Ii_d
      - gi*Ri_d
      - dd*Ri_d
      - lm*Ri_d
      + lm*ir_d
      + om*IR_d
      + al*pa(kid, sI_i, IR_d, I)*cc(idi, idd, sR_i, sR_d, s*R, N_Qd, N_Qi)
      + al*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + al*pa(kid, rI_i, IR_d, I)*cc(idi, idd, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kid, sR_i, Ri_d, R)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
      - al*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      - al*pa(kid, rR_i, Ri_d, R)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + nu*pa(kid, II_i, IR_d, I)*cc(idi, idd, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kid, IR_i, Ri_d, R)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kid, iI_i, IR_d, I)*cc(idi, idd, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      - qid*al*Ri_d
      - qid*nu*Ri_d;

    Ri_i_t = + bi*pa(kid, Rs_i, si_d, s)*cc(idi, idd, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      + b2*pa(kid, Rs_i, sI_d, s)*cc(idi, idd, RI_i, RI_d, R*I, N_Qd, N_Qi)
      + gd*Ii_i
      - gi*Ri_i
      - dd*Ri_i
      - lm*Ri_i
      + lm*ir_i
      + om*IR_i
      + al*pa(kii, sI_i, IR_i, I)*cc(iii, iid, sR_i, sR_d, s*R, N_Qd, N_Qi)
      + al*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      + al*pa(kii, rI_i, IR_i, I)*cc(iii, iid, rR_i, rR_d, r*R, N_Qd, N_Qi)
      - al*pa(kii, sR_i, Ri_i, R)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
      - al*Ri_i
      - al*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      - al*pa(kii, rR_i, Ri_i, R)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + nu*pa(kii, II_i, IR_i, I)*cc(iii, iid, IR_i, IR_d, I*R, N_Qd, N_Qi)
      - nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iR_i, iR_d, i*R, N_Qd, N_Qi)
      - nu*Ri_i
      - nu*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      ;
            
    Rr_d_t = + gd*Ir_d
      + gi*Ri_d
      - dd*Rr_d
      - di*Rr_d
      - lm*Rr_d
      + lm*rr_d
      + al*pa(kdi, RR_d, Rs_i, R)*cc(dii, did, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
      - al*pa(kid, sR_i, Rr_d, R)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
      + al*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      - al*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + al*pa(kdi, RR_d, Rr_i, R)*cc(dii, did, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
      - al*pa(kid, rR_i, Rr_d, R)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + nu*pa(kdi, RR_d, RI_i, R)*cc(dii, did, RI_i, RI_d, R*I, N_Qd, N_Qi)
      - nu*pa(kid, IR_i, Rr_d, R)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      - nu*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - qid*al*Rr_d;

    Rr_i_t = + gd*Ir_i
      + gi*Ri_i
      - dd*Rr_i
      - di*Rr_i
      - lm*Rr_i
      + lm*rr_i
      + al*pa(kii, RR_i, Rs_i, R)*cc(iii, iid, Rs_i, Rs_d, R*s, N_Qd, N_Qi)
      - al*pa(kii, sR_i, Rr_i, R)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
      + al*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      - al*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - al*Rr_i
      + al*pa(kii, RR_i, Rr_i, R)*cc(iii, iid, Rr_i, Rr_d, R*r, N_Qd, N_Qi)
      - al*pa(kii, rR_i, Rr_i, R)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + nu*pa(kii, RR_i, RI_i, R)*cc(iii, iid, RI_i, RI_d, R*I, N_Qd, N_Qi)
      - nu*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, Ri_i, Ri_d, R*i, N_Qd, N_Qi)
      - nu*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      ;
            
    ss_d_t = 2*(
                - bi*pa(kdd, ss_d, si_d, s)*cc(ddi, ddd, si_i, si_d, s*i, N_Qd, N_Qi)
                - b2*pa(kdd, ss_d, sI_d, s)*cc(ddi, ddd, sI_i, sI_d, s*I, N_Qd, N_Qi)
                + di*sr_d
                - lm*ss_d
                + al*pa(kid, sS_i, Ss_d, S)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
                + al*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
                + al*pa(kid, rS_i, Ss_d, S)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
                + nu*pa(kid, IS_i, Ss_d, S)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
                + nu*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
                + qid*al*Ss_d
                );

    ss_i_t = 2*(
                - bi*pa(kid, ss_i, si_d, s)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
                - b2*pa(kid, ss_i, sI_d, s)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
                + di*sr_i
                - lm*ss_i
                + al*Ss_i
                + al*pa(kii, sS_i, Ss_i, S)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
                + al*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
                + al*pa(kii, rS_i, Ss_i, S)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
                + nu*pa(kii, IS_i, Ss_i, S)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
                + nu*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
                );

    si_d_t = - bi*si_d
      + bi*pa(kdd, ss_d, si_d, s)*cc(ddi, ddd, si_i, si_d, s*i, N_Qd, N_Qi)
      - bi*pa(kdd, is_d, si_d, s)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + b2*pa(kdd, ss_d, sI_d, s)*cc(ddi, ddd, sI_i, sI_d, s*I, N_Qd, N_Qi)
      - b2*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - gi*si_d
      + di*ir_d
      - lm*si_d
      - lm*si_d
      + om*Is_d
      + al*pa(kid, sS_i, Si_d, S)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
      + al*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + al*pa(kid, rS_i, Si_d, S)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + al*pa(kid, sI_i, Is_d, I)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
      + al*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      + al*pa(kid, rI_i, Is_d, I)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kid, IS_i, Si_d, S)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kid, II_i, Is_d, I)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kid, iS_i, Si_d, S)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + nu*pa(kid, iI_i, Is_d, I)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      + qid*al*Si_d
      + qid*al*Is_d
      + qid*nu*Si_d;

    si_i_t = + bi*pa(kid, ss_i, si_d, s)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
      - bi*pa(kdi, is_d, si_i, s)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + b2*pa(kid, ss_i, sI_d, s)*cc(idi, idd, sI_i, sI_d, s*I, N_Qd, N_Qi)
      - b2*pa(kdi, Is_d, si_i, s)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      - gi*si_i
      + di*ir_i
      - lm*si_i
      - lm*si_i
      + om*Is_i
      + al*pa(kii, sS_i, Si_i, S)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
      + al*Si_i
      + al*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + al*pa(kii, rS_i, Si_i, S)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + al*Is_i
      + al*pa(kii, sI_i, Is_i, I)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
      + al*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      + al*pa(kii, rI_i, Is_i, I)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kii, IS_i, Si_i, S)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kii, II_i, Is_i, I)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*Si_i
      + nu*pa(kii, iS_i, Si_i, S)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + nu*pa(kii, iI_i, Is_i, I)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      - qdi*bi*si_i;

    sr_d_t = - bi*pa(kdd, is_d, sr_d, s)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - b2*pa(kdd, Is_d, sr_d, s)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + gi*si_d
      - di*sr_d
      + di*rr_d
      - lm*sr_d
      - lm*sr_d
      + al*pa(kid, sS_i, Sr_d, S)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
      + al*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + al*pa(kid, rS_i, Sr_d, S)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kid, sR_i, Rs_d, R)*cc(idi, idd, ss_i, ss_d, s*s, N_Qd, N_Qi)
      + al*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      + al*pa(kid, rR_i, Rs_d, R)*cc(idi, idd, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kid, IS_i, Sr_d, S)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kid, IR_i, Rs_d, R)*cc(idi, idd, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, is_i, is_d, i*s, N_Qd, N_Qi)
      + qid*al*Sr_d
      + qid*al*Rs_d;

    sr_i_t = - bi*pa(kdi, is_d, sr_i, s)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
      - b2*pa(kdi, Is_d, sr_i, s)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + gi*si_i
      - di*sr_i
      + di*rr_i
      - lm*sr_i
      - lm*sr_i
      + al*pa(kii, sS_i, Sr_i, S)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
      + al*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + al*Sr_i
      + al*pa(kii, rS_i, Sr_i, S)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*Rs_i
      + al*pa(kii, sR_i, Rs_i, R)*cc(iii, iid, ss_i, ss_d, s*s, N_Qd, N_Qi)
      + al*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      + al*pa(kii, rR_i, Rs_i, R)*cc(iii, iid, rs_i, rs_d, r*s, N_Qd, N_Qi)
      + nu*pa(kii, IS_i, Sr_i, S)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, Is_i, Is_d, I*s, N_Qd, N_Qi)
      + nu*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, is_i, is_d, i*s, N_Qd, N_Qi)
      ;
            
    ii_d_t = 2*(
                + bi*si_d
                + bi*pa(kdd, is_d, si_d, s)*cc(ddi, ddd, ii_i, ii_d, i*i, N_Qd, N_Qi)
                + b2*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                - gi*ii_d
                - lm*ii_d
                + om*Ii_d
                + al*pa(kid, sI_i, Ii_d, I)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
                + al*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
                + al*pa(kid, rI_i, Ii_d, I)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
                + nu*pa(kid, II_i, Ii_d, I)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                + nu*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
                + qid*al*Ii_d
                + qid*nu*Ii_d
                );

    ii_i_t = 2*(
                + bi*pa(kdi, is_d, si_i, s)*cc(dii, did, ii_i, ii_d, i*i, N_Qd, N_Qi)
                + b2*pa(kdi, Is_d, si_i, s)*cc(dii, did, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                - gi*ii_i
                - lm*ii_i
                + om*Ii_i
                + al*pa(kii, sI_i, Ii_i, I)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
                + al*Ii_i
                + al*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
                + al*pa(kii, rI_i, Ii_i, I)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
                + nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
                + nu*Ii_i
                + nu*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
                + qdi*bi*si_i
                );

    ir_d_t = + bi*pa(kdd, is_d, sr_d, s)*cc(ddi, ddd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + b2*pa(kdd, Is_d, sr_d, s)*cc(ddi, ddd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + gi*ii_d
      - gi*ir_d
      - di*ir_d
      - lm*ir_d
      - lm*ir_d
      + om*Ir_d
      + al*pa(kid, sI_i, Ir_d, I)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
      + al*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + al*pa(kid, rI_i, Ir_d, I)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kid, sR_i, Ri_d, R)*cc(idi, idd, si_i, si_d, s*i, N_Qd, N_Qi)
      + al*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + al*pa(kid, rR_i, Ri_d, R)*cc(idi, idd, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + nu*pa(kid, II_i, Ir_d, I)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kid, IR_i, Ri_d, R)*cc(idi, idd, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + qid*al*Ir_d
      + qid*al*Ri_d
      + qid*nu*Ri_d;

    ir_i_t = + bi*pa(kdi, is_d, sr_i, s)*cc(dii, did, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + b2*pa(kdi, Is_d, sr_i, s)*cc(dii, did, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + gi*ii_i
      - gi*ir_i
      - di*ir_i
      - lm*ir_i
      - lm*ir_i
      + om*Ir_i
      + al*pa(kii, sI_i, Ir_i, I)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
      + al*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + al*Ir_i
      + al*pa(kii, rI_i, Ir_i, I)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
      + al*pa(kii, sR_i, Ri_i, R)*cc(iii, iid, si_i, si_d, s*i, N_Qd, N_Qi)
      + al*Ri_i
      + al*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      + al*pa(kii, rR_i, Ri_i, R)*cc(iii, iid, ri_i, ri_d, r*i, N_Qd, N_Qi)
      + nu*pa(kii, II_i, Ir_i, I)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
      + nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, Ii_i, Ii_d, I*i, N_Qd, N_Qi)
      + nu*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
      + nu*Ri_i
      + nu*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, ii_i, ii_d, i*i, N_Qd, N_Qi)
      ;
            
    rr_d_t = 2*(
                + gi*ir_d
                - di*rr_d
                - lm*rr_d
                + al*pa(kid, sR_i, Rr_d, R)*cc(idi, idd, sr_i, sr_d, s*r, N_Qd, N_Qi)
                + al*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
                + al*pa(kid, rR_i, Rr_d, R)*cc(idi, idd, rr_i, rr_d, r*r, N_Qd, N_Qi)
                + nu*pa(kid, IR_i, Rr_d, R)*cc(idi, idd, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                + nu*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, ir_i, ir_d, i*r, N_Qd, N_Qi)
                + qid*al*Rr_d
                );

    rr_i_t = 2*(
                + gi*ir_i
                - di*rr_i
                - lm*rr_i
                + al*pa(kii, sR_i, Rr_i, R)*cc(iii, iid, sr_i, sr_d, s*r, N_Qd, N_Qi)
                + al*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
                + al*Rr_i
                + al*pa(kii, rR_i, Rr_i, R)*cc(iii, iid, rr_i, rr_d, r*r, N_Qd, N_Qi)
                + nu*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, Ir_i, Ir_d, I*r, N_Qd, N_Qi)
                + nu*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, ir_i, ir_d, i*r, N_Qd, N_Qi)
                );



    return GSL_SUCCESS;         
  }
      
}; // FullModelPairApprox1

//------------------------------------------------------------

struct FullModelPairApprox2
{
  // pair approx      
  static double pa(double kk, double ij, double jk, double j)
  {
    const double eps = 1e-8;
         
    if(j < eps)
      return 0.0;
    else
      return kk*((ij*jk)/j);
  }
      
  // clustering correction
  static double  cc(double C_i, double C_d, double C_b,
                    double ki_i, double ki_d, double ki_b,
                    double ki, double N_Qd, double N_Qi, double N_Qb)
  {
    if(ki == 0.0)
      return ((1.0 - C_i - C_d - C_b));
    else
      return ((1.0 - C_i - C_d - C_b)
              + C_i*N_Qi*ki_i/ki + C_d*N_Qd*ki_d/ki
              + C_b*N_Qb*ki_b/ki);
  }           
      
  // rhs function
  static int rhs_eval (double t, const double y[], double rhs[], void* params)
  {
    FullModelParams p = *(static_cast<FullModelParams*>(params));
         
    // local readable short variables
    double bd=p.beta[0][0], gd=p.gamma[0], dd=p.delta[0];
    double bi=p.beta[1][1], gi=p.gamma[1], di=p.delta[1];
    double b1=p.beta[0][1], b2=p.beta[1][0], al=p.alpha, nu=p.nu;
    double lm=p.lambda, om=p.omega;
    double Qd=p.Qd, Qi=p.Qi, Qb=p.Qb, N=p.N;


    double N_Qd = 0.0, N_Qi = 0.0, N_Qb = 0.0;
    if (Qd > 1.0) N_Qd = N/Qd;
    if (Qi > 1.0) N_Qi = N/Qi;
    if (Qb > 1.0) N_Qb = N/Qb;
         
    // scaling rates -> rates per contact
    double bdbr = bd, bibr = bi, b1b = b1, b2b = b2;
    double alb = al, nub = nu;
         
    //bd/=Qd; bi/=Qd; b1/=Qd; b2/=Qd;
    //al/=Qi; nu/=Qi;
         
    // Clustering corrections         
    double idi=p.C[0][1][0], idd=p.C[0][1][1];
    double dii=p.C[1][0][0], did=p.C[1][0][1];
    double ddi=p.C[1][1][0], ddd=p.C[1][1][1];
    double iii=p.C[0][0][0], iid=p.C[0][0][1];
         
    double ibi=p.C[0][2][0], ibb=p.C[0][2][2];
    double bii=p.C[2][0][0], bib=p.C[2][0][2];
    double bbi=p.C[2][2][0], bbb=p.C[2][2][2];
    double                   iib=p.C[0][0][2];

    double bdb=p.C[2][1][2], bdd=p.C[2][1][1];
    double dbb=p.C[1][2][2], dbd=p.C[1][2][1];
    double ddb=p.C[1][1][2]                  ;
    double                   bbd=p.C[2][2][1];

    double idb=p.C[0][1][2], dib=p.C[1][0][2];
    double ibd=p.C[0][2][1], bid=p.C[2][0][1];
    double bdi=p.C[2][1][0], dbi=p.C[1][2][0];
         
    double kdd = 0.0, kii = 0.0, kbb = 0.0;
    if (Qd > 1.0) kdd = (Qd-1.0)/Qd;
    if (Qi > 1.0) kii = (Qi-1.0)/Qi;
    if (Qb > 1.0) kbb = (Qb-1.0)/Qb;
         
    double kdi = 1.0, kid = 1.0, kdb = 1.0, kbd = 1.0;
    double kbi = 1.0, kib = 1.0;

    // overlaping
    //double qdi=p.qdi, qid=p.qid;



    S_t = - bd*SI_d
      - b1*Si_d
      - bdbr*SI_b
      - b1b*Si_b
      + dd*R
      + lm*s
      - al*Ss_i
      - al*Si_i
      - al*Sr_i
      - nu*SI_i
      - nu*Si_i
      - alb*Ss_b
      - alb*Si_b
      - alb*Sr_b
      - nub*SI_b
      - nub*Si_b;

    I_t = + bd*SI_d
      + b1*Si_d
      + bdbr*SI_b
      + b1b*Si_b
      - gd*I
      + lm*i
      - om*I
      - al*Is_i
      - al*Ii_i
      - al*Ir_i
      - nu*II_i
      - nu*Ii_i
      - alb*Is_b
      - alb*Ii_b
      - alb*Ir_b
      - nub*II_b
      - nub*Ii_b;

    R_t = + gd*I
      - dd*R
      + lm*r
      - al*Rs_i
      - al*Ri_i
      - al*Rr_i
      - nu*IR_i
      - nu*Ri_i
      - alb*Rs_b
      - alb*Ri_b
      - alb*Rr_b
      - nub*IR_b
      - nub*Ri_b;

    s_t = - bi*si_d
      - b2*Is_d
      - bibr*si_b
      - b2b*Is_b
      + di*r
      - lm*s
      + al*Ss_i
      + al*Si_i
      + al*Sr_i
      + nu*SI_i
      + nu*Si_i
      + alb*Ss_b
      + alb*Si_b
      + alb*Sr_b
      + nub*SI_b
      + nub*Si_b;

    i_t = + bi*si_d
      + b2*Is_d
      + bibr*si_b
      + b2b*Is_b
      - gi*i
      - lm*i
      + om*I
      + al*Is_i
      + al*Ii_i
      + al*Ir_i
      + nu*II_i
      + nu*Ii_i
      + alb*Is_b
      + alb*Ii_b
      + alb*Ir_b
      + nub*II_b
      + nub*Ii_b;

    r_t = + gi*i
      - di*r
      - lm*r
      + al*Rs_i
      + al*Ri_i
      + al*Rr_i
      + nu*IR_i
      + nu*Ri_i
      + alb*Rs_b
      + alb*Ri_b
      + alb*Rr_b
      + nub*IR_b
      + nub*Ri_b;

    SS_d_t = 2*(
                - bd*pa(kdd, SS_d, SI_d, S)*cc(ddi, ddd, ddb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - b1*pa(kdd, SS_d, Si_d, S)*cc(ddi, ddd, ddb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - bdbr*pa(kdb, SS_d, SI_b, S)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - b1b*pa(kdb, SS_d, Si_b, S)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                + dd*SR_d
                + lm*Ss_d
                - al*pa(kdi, SS_d, Ss_i, S)*cc(dii, did, dib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kdi, SS_d, Si_i, S)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kdi, SS_d, Sr_i, S)*cc(dii, did, dib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kdi, SS_d, SI_i, S)*cc(dii, did, dib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kdi, SS_d, Si_i, S)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, SS_d, Ss_b, S)*cc(dbi, dbd, dbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, SS_d, Si_b, S)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, SS_d, Sr_b, S)*cc(dbi, dbd, dbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kdb, SS_d, SI_b, S)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kdb, SS_d, Si_b, S)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                );

    SS_i_t = 2*(
                - bd*pa(kid, SS_i, SI_d, S)*cc(idi, idd, idb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - b1*pa(kid, SS_i, Si_d, S)*cc(idi, idd, idb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - bdbr*pa(kib, SS_i, SI_b, S)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - b1b*pa(kib, SS_i, Si_b, S)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                + dd*SR_i
                + lm*Ss_i
                - al*pa(kii, SS_i, Ss_i, S)*cc(iii, iid, iib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kii, SS_i, Si_i, S)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kii, SS_i, Sr_i, S)*cc(iii, iid, iib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kii, SS_i, SI_i, S)*cc(iii, iid, iib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kii, SS_i, Si_i, S)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, SS_i, Ss_b, S)*cc(ibi, ibd, ibb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, SS_i, Si_b, S)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, SS_i, Sr_b, S)*cc(ibi, ibd, ibb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kib, SS_i, SI_b, S)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kib, SS_i, Si_b, S)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                );

    SS_b_t = 2*(
                - bd*pa(kbd, SS_b, SI_d, S)*cc(bdi, bdd, bdb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - b1*pa(kbd, SS_b, Si_d, S)*cc(bdi, bdd, bdb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - bdbr*pa(kbb, SS_b, SI_b, S)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - b1b*pa(kbb, SS_b, Si_b, S)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                + dd*SR_b
                + lm*Ss_b
                - al*pa(kbi, SS_b, Ss_i, S)*cc(bii, bid, bib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kbi, SS_b, Si_i, S)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kbi, SS_b, Sr_i, S)*cc(bii, bid, bib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kbi, SS_b, SI_i, S)*cc(bii, bid, bib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kbi, SS_b, Si_i, S)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, SS_b, Ss_b, S)*cc(bbi, bbd, bbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, SS_b, Si_b, S)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, SS_b, Sr_b, S)*cc(bbi, bbd, bbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kbb, SS_b, SI_b, S)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kbb, SS_b, Si_b, S)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
                );

    SI_d_t = - bd*SI_d
      + bd*pa(kdd, SS_d, SI_d, S)*cc(ddi, ddd, ddb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bd*pa(kdd, IS_d, SI_d, S)*cc(ddi, ddd, ddb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdd, SS_d, Si_d, S)*cc(ddi, ddd, ddb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdd, iS_d, SI_d, S)*cc(ddi, ddd, ddb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kdb, SS_d, SI_b, S)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbd, IS_b, SI_d, S)*cc(bdi, bdd, bdb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kdb, SS_d, Si_b, S)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbd, iS_b, SI_d, S)*cc(bdi, bdd, bdb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - gd*SI_d
      + dd*IR_d
      + lm*Is_d
      + lm*Si_d
      - om*SI_d
      - al*pa(kid, sS_i, SI_d, S)*cc(idi, idd, idb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iS_i, SI_d, S)*cc(idi, idd, idb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rS_i, SI_d, S)*cc(idi, idd, idb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, SI_d, Is_i, I)*cc(dii, did, dib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, SI_d, Ir_i, I)*cc(dii, did, dib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IS_i, SI_d, S)*cc(idi, idd, idb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kdi, SI_d, II_i, I)*cc(dii, did, dib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iS_i, SI_d, S)*cc(idi, idd, idb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sS_b, SI_d, S)*cc(bdi, bdd, bdb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iS_b, SI_d, S)*cc(bdi, bdd, bdb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rS_b, SI_d, S)*cc(bdi, bdd, bdb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, SI_d, Is_b, I)*cc(dbi, dbd, dbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, SI_d, Ii_b, I)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, SI_d, Ir_b, I)*cc(dbi, dbd, dbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IS_b, SI_d, S)*cc(bdi, bdd, bdb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kdb, SI_d, II_b, I)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iS_b, SI_d, S)*cc(bdi, bdd, bdb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kdb, SI_d, Ii_b, I)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    SI_i_t = + bd*pa(kid, SS_i, SI_d, S)*cc(idi, idd, idb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bd*pa(kdi, IS_d, SI_i, S)*cc(dii, did, dib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + b1*pa(kid, SS_i, Si_d, S)*cc(idi, idd, idb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdi, iS_d, SI_i, S)*cc(dii, did, dib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kib, SS_i, SI_b, S)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbi, IS_b, SI_i, S)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kib, SS_i, Si_b, S)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbi, iS_b, SI_i, S)*cc(bii, bid, bib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - gd*SI_i
      + dd*IR_i
      + lm*Is_i
      + lm*Si_i
      - om*SI_i
      - al*pa(kii, sS_i, SI_i, S)*cc(iii, iid, iib, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rS_i, SI_i, S)*cc(iii, iid, iib, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, SI_i, Is_i, I)*cc(iii, iid, iib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, SI_i, Ir_i, I)*cc(iii, iid, iib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*SI_i
      - nu*pa(kii, IS_i, SI_i, S)*cc(iii, iid, iib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, SI_i, II_i, I)*cc(iii, iid, iib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sS_b, SI_i, S)*cc(bii, bid, bib, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iS_b, SI_i, S)*cc(bii, bid, bib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rS_b, SI_i, S)*cc(bii, bid, bib, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, SI_i, Is_b, I)*cc(ibi, ibd, ibb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, SI_i, Ii_b, I)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, SI_i, Ir_b, I)*cc(ibi, ibd, ibb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IS_b, SI_i, S)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kib, SI_i, II_b, I)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iS_b, SI_i, S)*cc(bii, bid, bib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kib, SI_i, Ii_b, I)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    SI_b_t = + bd*pa(kbd, SS_b, SI_d, S)*cc(bdi, bdd, bdb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bd*pa(kdb, IS_d, SI_b, S)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + b1*pa(kbd, SS_b, Si_d, S)*cc(bdi, bdd, bdb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdb, iS_d, SI_b, S)*cc(dbi, dbd, dbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - bdbr*SI_b
      + bdbr*pa(kbb, SS_b, SI_b, S)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbb, IS_b, SI_b, S)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbb, SS_b, Si_b, S)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbb, iS_b, SI_b, S)*cc(bbi, bbd, bbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - gd*SI_b
      + dd*IR_b
      + lm*Is_b
      + lm*Si_b
      - om*SI_b
      - al*pa(kib, sS_i, SI_b, S)*cc(ibi, ibd, ibb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iS_i, SI_b, S)*cc(ibi, ibd, ibb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rS_i, SI_b, S)*cc(ibi, ibd, ibb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, SI_b, Is_i, I)*cc(bii, bid, bib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, SI_b, Ii_i, I)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, SI_b, Ir_i, I)*cc(bii, bid, bib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IS_i, SI_b, S)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kbi, SI_b, II_i, I)*cc(bii, bid, bib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iS_i, SI_b, S)*cc(ibi, ibd, ibb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kbi, SI_b, Ii_i, I)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sS_b, SI_b, S)*cc(bbi, bbd, bbb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iS_b, SI_b, S)*cc(bbi, bbd, bbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rS_b, SI_b, S)*cc(bbi, bbd, bbb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, SI_b, Is_b, I)*cc(bbi, bbd, bbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, SI_b, Ii_b, I)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, SI_b, Ir_b, I)*cc(bbi, bbd, bbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*SI_b
      - nub*pa(kbb, IS_b, SI_b, S)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, SI_b, II_b, I)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iS_b, SI_b, S)*cc(bbi, bbd, bbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, SI_b, Ii_b, I)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    SR_d_t = - bd*pa(kdd, IS_d, SR_d, S)*cc(ddi, ddd, ddb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdd, iS_d, SR_d, S)*cc(ddi, ddd, ddb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbd, IS_b, SR_d, S)*cc(bdi, bdd, bdb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbd, iS_b, SR_d, S)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + gd*SI_d
      - dd*SR_d
      + dd*RR_d
      + lm*Rs_d
      + lm*Sr_d
      - al*pa(kid, sS_i, SR_d, S)*cc(idi, idd, idb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iS_i, SR_d, S)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rS_i, SR_d, S)*cc(idi, idd, idb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, SR_d, Rs_i, R)*cc(dii, did, dib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, SR_d, Rr_i, R)*cc(dii, did, dib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IS_i, SR_d, S)*cc(idi, idd, idb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kdi, SR_d, RI_i, R)*cc(dii, did, dib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iS_i, SR_d, S)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sS_b, SR_d, S)*cc(bdi, bdd, bdb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iS_b, SR_d, S)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rS_b, SR_d, S)*cc(bdi, bdd, bdb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, SR_d, Rs_b, R)*cc(dbi, dbd, dbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, SR_d, Ri_b, R)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, SR_d, Rr_b, R)*cc(dbi, dbd, dbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IS_b, SR_d, S)*cc(bdi, bdd, bdb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kdb, SR_d, RI_b, R)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iS_b, SR_d, S)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kdb, SR_d, Ri_b, R)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    SR_i_t = - bd*pa(kdi, IS_d, SR_i, S)*cc(dii, did, dib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdi, iS_d, SR_i, S)*cc(dii, did, dib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbi, IS_b, SR_i, S)*cc(bii, bid, bib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbi, iS_b, SR_i, S)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + gd*SI_i
      - dd*SR_i
      + dd*RR_i
      + lm*Rs_i
      + lm*Sr_i
      - al*pa(kii, sS_i, SR_i, S)*cc(iii, iid, iib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rS_i, SR_i, S)*cc(iii, iid, iib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, SR_i, Rs_i, R)*cc(iii, iid, iib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, SR_i, Rr_i, R)*cc(iii, iid, iib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IS_i, SR_i, S)*cc(iii, iid, iib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, SR_i, RI_i, R)*cc(iii, iid, iib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sS_b, SR_i, S)*cc(bii, bid, bib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iS_b, SR_i, S)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rS_b, SR_i, S)*cc(bii, bid, bib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, SR_i, Rs_b, R)*cc(ibi, ibd, ibb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, SR_i, Ri_b, R)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, SR_i, Rr_b, R)*cc(ibi, ibd, ibb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IS_b, SR_i, S)*cc(bii, bid, bib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kib, SR_i, RI_b, R)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iS_b, SR_i, S)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kib, SR_i, Ri_b, R)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    SR_b_t = - bd*pa(kdb, IS_d, SR_b, S)*cc(dbi, dbd, dbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdb, iS_d, SR_b, S)*cc(dbi, dbd, dbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbb, IS_b, SR_b, S)*cc(bbi, bbd, bbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbb, iS_b, SR_b, S)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + gd*SI_b
      - dd*SR_b
      + dd*RR_b
      + lm*Rs_b
      + lm*Sr_b
      - al*pa(kib, sS_i, SR_b, S)*cc(ibi, ibd, ibb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iS_i, SR_b, S)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rS_i, SR_b, S)*cc(ibi, ibd, ibb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, SR_b, Rs_i, R)*cc(bii, bid, bib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, SR_b, Ri_i, R)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, SR_b, Rr_i, R)*cc(bii, bid, bib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IS_i, SR_b, S)*cc(ibi, ibd, ibb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kbi, SR_b, RI_i, R)*cc(bii, bid, bib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iS_i, SR_b, S)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kbi, SR_b, Ri_i, R)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sS_b, SR_b, S)*cc(bbi, bbd, bbb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iS_b, SR_b, S)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rS_b, SR_b, S)*cc(bbi, bbd, bbb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, SR_b, Rs_b, R)*cc(bbi, bbd, bbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, SR_b, Ri_b, R)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, SR_b, Rr_b, R)*cc(bbi, bbd, bbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IS_b, SR_b, S)*cc(bbi, bbd, bbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, SR_b, RI_b, R)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iS_b, SR_b, S)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, SR_b, Ri_b, R)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    Ss_d_t = - bd*pa(kdd, IS_d, Ss_d, S)*cc(ddi, ddd, ddb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bi*pa(kdd, Ss_d, si_d, s)*cc(ddi, ddd, ddb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdd, iS_d, Ss_d, S)*cc(ddi, ddd, ddb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdd, Ss_d, sI_d, s)*cc(ddi, ddd, ddb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbd, IS_b, Ss_d, S)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kdb, Ss_d, si_b, s)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbd, iS_b, Ss_d, S)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kdb, Ss_d, sI_b, s)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      + dd*Rs_d
      + di*Sr_d
      - lm*Ss_d
      + lm*ss_d
      + al*pa(kdi, SS_d, Ss_i, S)*cc(dii, did, dib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, sS_i, Ss_d, S)*cc(idi, idd, idb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SS_d, Si_i, S)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SS_d, Sr_i, S)*cc(dii, did, dib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rS_i, Ss_d, S)*cc(idi, idd, idb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, SS_d, SI_i, S)*cc(dii, did, dib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IS_i, Ss_d, S)*cc(idi, idd, idb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, SS_d, Si_i, S)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SS_d, Ss_b, S)*cc(dbi, dbd, dbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sS_b, Ss_d, S)*cc(bdi, bdd, bdb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SS_d, Si_b, S)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iS_b, Ss_d, S)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SS_d, Sr_b, S)*cc(dbi, dbd, dbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rS_b, Ss_d, S)*cc(bdi, bdd, bdb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, SS_d, SI_b, S)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IS_b, Ss_d, S)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, SS_d, Si_b, S)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iS_b, Ss_d, S)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Ss_i_t = - bd*pa(kdi, IS_d, Ss_i, S)*cc(dii, did, dib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bi*pa(kid, Ss_i, si_d, s)*cc(idi, idd, idb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdi, iS_d, Ss_i, S)*cc(dii, did, dib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2*pa(kid, Ss_i, sI_d, s)*cc(idi, idd, idb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbi, IS_b, Ss_i, S)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kib, Ss_i, si_b, s)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbi, iS_b, Ss_i, S)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kib, Ss_i, sI_b, s)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      + dd*Rs_i
      + di*Sr_i
      - lm*Ss_i
      + lm*ss_i
      - al*Ss_i
      + al*pa(kii, SS_i, Ss_i, S)*cc(iii, iid, iib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, sS_i, Ss_i, S)*cc(iii, iid, iib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SS_i, Si_i, S)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SS_i, Sr_i, S)*cc(iii, iid, iib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rS_i, Ss_i, S)*cc(iii, iid, iib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, SS_i, SI_i, S)*cc(iii, iid, iib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IS_i, Ss_i, S)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, SS_i, Si_i, S)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SS_i, Ss_b, S)*cc(ibi, ibd, ibb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sS_b, Ss_i, S)*cc(bii, bid, bib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SS_i, Si_b, S)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iS_b, Ss_i, S)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SS_i, Sr_b, S)*cc(ibi, ibd, ibb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rS_b, Ss_i, S)*cc(bii, bid, bib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, SS_i, SI_b, S)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IS_b, Ss_i, S)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, SS_i, Si_b, S)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iS_b, Ss_i, S)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Ss_b_t = - bd*pa(kdb, IS_d, Ss_b, S)*cc(dbi, dbd, dbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bi*pa(kbd, Ss_b, si_d, s)*cc(bdi, bdd, bdb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdb, iS_d, Ss_b, S)*cc(dbi, dbd, dbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2*pa(kbd, Ss_b, sI_d, s)*cc(bdi, bdd, bdb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbb, IS_b, Ss_b, S)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbb, Ss_b, si_b, s)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbb, iS_b, Ss_b, S)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbb, Ss_b, sI_b, s)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      + dd*Rs_b
      + di*Sr_b
      - lm*Ss_b
      + lm*ss_b
      + al*pa(kbi, SS_b, Ss_i, S)*cc(bii, bid, bib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, sS_i, Ss_b, S)*cc(ibi, ibd, ibb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SS_b, Si_i, S)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iS_i, Ss_b, S)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SS_b, Sr_i, S)*cc(bii, bid, bib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rS_i, Ss_b, S)*cc(ibi, ibd, ibb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, SS_b, SI_i, S)*cc(bii, bid, bib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IS_i, Ss_b, S)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, SS_b, Si_i, S)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iS_i, Ss_b, S)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*Ss_b
      + alb*pa(kbb, SS_b, Ss_b, S)*cc(bbi, bbd, bbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sS_b, Ss_b, S)*cc(bbi, bbd, bbb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SS_b, Si_b, S)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iS_b, Ss_b, S)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SS_b, Sr_b, S)*cc(bbi, bbd, bbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rS_b, Ss_b, S)*cc(bbi, bbd, bbb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, SS_b, SI_b, S)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IS_b, Ss_b, S)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, SS_b, Si_b, S)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iS_b, Ss_b, S)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Si_d_t = - bd*pa(kdd, IS_d, Si_d, S)*cc(ddi, ddd, ddb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bi*pa(kdd, Ss_d, si_d, s)*cc(ddi, ddd, ddb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*Si_d
      - b1*pa(kdd, iS_d, Si_d, S)*cc(ddi, ddd, ddb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kdd, Ss_d, sI_d, s)*cc(ddi, ddd, ddb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbd, IS_b, Si_d, S)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kdb, Ss_d, si_b, s)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbd, iS_b, Si_d, S)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kdb, Ss_d, sI_b, s)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - gi*Si_d
      + dd*Ri_d
      + lm*si_d
      - lm*Si_d
      + om*SI_d
      - al*pa(kid, sS_i, Si_d, S)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iS_i, Si_d, S)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rS_i, Si_d, S)*cc(idi, idd, idb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SI_d, Is_i, I)*cc(dii, did, dib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SI_d, Ir_i, I)*cc(dii, did, dib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IS_i, Si_d, S)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, SI_d, II_i, I)*cc(dii, did, dib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iS_i, Si_d, S)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, SI_d, Ii_i, I)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sS_b, Si_d, S)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iS_b, Si_d, S)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rS_b, Si_d, S)*cc(bdi, bdd, bdb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SI_d, Is_b, I)*cc(dbi, dbd, dbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SI_d, Ii_b, I)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SI_d, Ir_b, I)*cc(dbi, dbd, dbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IS_b, Si_d, S)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, SI_d, II_b, I)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iS_b, Si_d, S)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, SI_d, Ii_b, I)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    Si_i_t = - bd*pa(kdi, IS_d, Si_i, S)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bi*pa(kid, Ss_i, si_d, s)*cc(idi, idd, idb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdi, iS_d, Si_i, S)*cc(dii, did, dib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kid, Ss_i, sI_d, s)*cc(idi, idd, idb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbi, IS_b, Si_i, S)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kib, Ss_i, si_b, s)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbi, iS_b, Si_i, S)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kib, Ss_i, sI_b, s)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - gi*Si_i
      + dd*Ri_i
      + lm*si_i
      - lm*Si_i
      + om*SI_i
      - al*pa(kii, sS_i, Si_i, S)*cc(iii, iid, iib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*Si_i
      - al*pa(kii, iS_i, Si_i, S)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rS_i, Si_i, S)*cc(iii, iid, iib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SI_i, Is_i, I)*cc(iii, iid, iib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SI_i, Ir_i, I)*cc(iii, iid, iib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IS_i, Si_i, S)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, SI_i, II_i, I)*cc(iii, iid, iib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*Si_i
      - nu*pa(kii, iS_i, Si_i, S)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, SI_i, Ii_i, I)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sS_b, Si_i, S)*cc(bii, bid, bib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iS_b, Si_i, S)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rS_b, Si_i, S)*cc(bii, bid, bib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SI_i, Is_b, I)*cc(ibi, ibd, ibb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SI_i, Ii_b, I)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SI_i, Ir_b, I)*cc(ibi, ibd, ibb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IS_b, Si_i, S)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, SI_i, II_b, I)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iS_b, Si_i, S)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, SI_i, Ii_b, I)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    Si_b_t = - bd*pa(kdb, IS_d, Si_b, S)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bi*pa(kbd, Ss_b, si_d, s)*cc(bdi, bdd, bdb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdb, iS_d, Si_b, S)*cc(dbi, dbd, dbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kbd, Ss_b, sI_d, s)*cc(bdi, bdd, bdb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbb, IS_b, Si_b, S)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kbb, Ss_b, si_b, s)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - b1b*Si_b
      - b1b*pa(kbb, iS_b, Si_b, S)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kbb, Ss_b, sI_b, s)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - gi*Si_b
      + dd*Ri_b
      + lm*si_b
      - lm*Si_b
      + om*SI_b
      - al*pa(kib, sS_i, Si_b, S)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iS_i, Si_b, S)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rS_i, Si_b, S)*cc(ibi, ibd, ibb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SI_b, Is_i, I)*cc(bii, bid, bib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SI_b, Ii_i, I)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SI_b, Ir_i, I)*cc(bii, bid, bib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IS_i, Si_b, S)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, SI_b, II_i, I)*cc(bii, bid, bib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iS_i, Si_b, S)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, SI_b, Ii_i, I)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sS_b, Si_b, S)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*Si_b
      - alb*pa(kbb, iS_b, Si_b, S)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rS_b, Si_b, S)*cc(bbi, bbd, bbb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SI_b, Is_b, I)*cc(bbi, bbd, bbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SI_b, Ii_b, I)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SI_b, Ir_b, I)*cc(bbi, bbd, bbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IS_b, Si_b, S)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, SI_b, II_b, I)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*Si_b
      - nub*pa(kbb, iS_b, Si_b, S)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, SI_b, Ii_b, I)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    Sr_d_t = - bd*pa(kdd, IS_d, Sr_d, S)*cc(ddi, ddd, ddb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdd, iS_d, Sr_d, S)*cc(ddi, ddd, ddb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbd, IS_b, Sr_d, S)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbd, iS_b, Sr_d, S)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + gi*Si_d
      + dd*Rr_d
      - di*Sr_d
      + lm*sr_d
      - lm*Sr_d
      - al*pa(kid, sS_i, Sr_d, S)*cc(idi, idd, idb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rS_i, Sr_d, S)*cc(idi, idd, idb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SR_d, Rs_i, R)*cc(dii, did, dib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, SR_d, Rr_i, R)*cc(dii, did, dib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IS_i, Sr_d, S)*cc(idi, idd, idb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, SR_d, RI_i, R)*cc(dii, did, dib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, SR_d, Ri_i, R)*cc(dii, did, dib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sS_b, Sr_d, S)*cc(bdi, bdd, bdb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iS_b, Sr_d, S)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rS_b, Sr_d, S)*cc(bdi, bdd, bdb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SR_d, Rs_b, R)*cc(dbi, dbd, dbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SR_d, Ri_b, R)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, SR_d, Rr_b, R)*cc(dbi, dbd, dbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IS_b, Sr_d, S)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, SR_d, RI_b, R)*cc(dbi, dbd, dbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iS_b, Sr_d, S)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, SR_d, Ri_b, R)*cc(dbi, dbd, dbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    Sr_i_t = - bd*pa(kdi, IS_d, Sr_i, S)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdi, iS_d, Sr_i, S)*cc(dii, did, dib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbi, IS_b, Sr_i, S)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbi, iS_b, Sr_i, S)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + gi*Si_i
      + dd*Rr_i
      - di*Sr_i
      + lm*sr_i
      - lm*Sr_i
      - al*pa(kii, sS_i, Sr_i, S)*cc(iii, iid, iib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*Sr_i
      - al*pa(kii, rS_i, Sr_i, S)*cc(iii, iid, iib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SR_i, Rs_i, R)*cc(iii, iid, iib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, SR_i, Rr_i, R)*cc(iii, iid, iib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IS_i, Sr_i, S)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, SR_i, RI_i, R)*cc(iii, iid, iib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, SR_i, Ri_i, R)*cc(iii, iid, iib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sS_b, Sr_i, S)*cc(bii, bid, bib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iS_b, Sr_i, S)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rS_b, Sr_i, S)*cc(bii, bid, bib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SR_i, Rs_b, R)*cc(ibi, ibd, ibb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SR_i, Ri_b, R)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, SR_i, Rr_b, R)*cc(ibi, ibd, ibb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IS_b, Sr_i, S)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, SR_i, RI_b, R)*cc(ibi, ibd, ibb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iS_b, Sr_i, S)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, SR_i, Ri_b, R)*cc(ibi, ibd, ibb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    Sr_b_t = - bd*pa(kdb, IS_d, Sr_b, S)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - b1*pa(kdb, iS_d, Sr_b, S)*cc(dbi, dbd, dbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - bdbr*pa(kbb, IS_b, Sr_b, S)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - b1b*pa(kbb, iS_b, Sr_b, S)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + gi*Si_b
      + dd*Rr_b
      - di*Sr_b
      + lm*sr_b
      - lm*Sr_b
      - al*pa(kib, sS_i, Sr_b, S)*cc(ibi, ibd, ibb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iS_i, Sr_b, S)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rS_i, Sr_b, S)*cc(ibi, ibd, ibb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SR_b, Rs_i, R)*cc(bii, bid, bib, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SR_b, Ri_i, R)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, SR_b, Rr_i, R)*cc(bii, bid, bib, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IS_i, Sr_b, S)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, SR_b, RI_i, R)*cc(bii, bid, bib, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iS_i, Sr_b, S)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, SR_b, Ri_i, R)*cc(bii, bid, bib, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sS_b, Sr_b, S)*cc(bbi, bbd, bbb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iS_b, Sr_b, S)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*Sr_b
      - alb*pa(kbb, rS_b, Sr_b, S)*cc(bbi, bbd, bbb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SR_b, Rs_b, R)*cc(bbi, bbd, bbb, Ss_i, Ss_d, Ss_b, S*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SR_b, Ri_b, R)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, SR_b, Rr_b, R)*cc(bbi, bbd, bbb, Sr_i, Sr_d, Sr_b, S*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IS_b, Sr_b, S)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, SR_b, RI_b, R)*cc(bbi, bbd, bbb, SI_i, SI_d, SI_b, S*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iS_b, Sr_b, S)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, SR_b, Ri_b, R)*cc(bbi, bbd, bbb, Si_i, Si_d, Si_b, S*i, N_Qd, N_Qi, N_Qb);

    II_d_t = 2*(
                + bd*SI_d
                + bd*pa(kdd, IS_d, SI_d, S)*cc(ddi, ddd, ddb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                + b1*pa(kdd, iS_d, SI_d, S)*cc(ddi, ddd, ddb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
                + bdbr*pa(kbd, IS_b, SI_d, S)*cc(bdi, bdd, bdb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                + b1b*pa(kbd, iS_b, SI_d, S)*cc(bdi, bdd, bdb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
                - gd*II_d
                + lm*Ii_d
                - om*II_d
                - al*pa(kdi, II_d, Is_i, I)*cc(dii, did, dib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kdi, II_d, Ii_i, I)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kdi, II_d, Ir_i, I)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kdi, II_d, II_i, I)*cc(dii, did, dib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kdi, II_d, Ii_i, I)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, II_d, Is_b, I)*cc(dbi, dbd, dbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, II_d, Ii_b, I)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, II_d, Ir_b, I)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kdb, II_d, II_b, I)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kdb, II_d, Ii_b, I)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                );

    II_i_t = 2*(
                + bd*pa(kdi, IS_d, SI_i, S)*cc(dii, did, dib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                + b1*pa(kdi, iS_d, SI_i, S)*cc(dii, did, dib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
                + bdbr*pa(kbi, IS_b, SI_i, S)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                + b1b*pa(kbi, iS_b, SI_i, S)*cc(bii, bid, bib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
                - gd*II_i
                + lm*Ii_i
                - om*II_i
                - al*pa(kii, II_i, Is_i, I)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kii, II_i, Ii_i, I)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kii, II_i, Ir_i, I)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                - nu*II_i
                - nu*pa(kii, II_i, II_i, I)*cc(iii, iid, iib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, II_i, Is_b, I)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, II_i, Ii_b, I)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, II_i, Ir_b, I)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kib, II_i, II_b, I)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kib, II_i, Ii_b, I)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                );

    II_b_t = 2*(
                + bd*pa(kdb, IS_d, SI_b, S)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                + b1*pa(kdb, iS_d, SI_b, S)*cc(dbi, dbd, dbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
                + bdbr*SI_b
                + bdbr*pa(kbb, IS_b, SI_b, S)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                + b1b*pa(kbb, iS_b, SI_b, S)*cc(bbi, bbd, bbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
                - gd*II_b
                + lm*Ii_b
                - om*II_b
                - al*pa(kbi, II_b, Is_i, I)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kbi, II_b, Ii_i, I)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kbi, II_b, Ir_i, I)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kbi, II_b, II_i, I)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kbi, II_b, Ii_i, I)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, II_b, Is_b, I)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, II_b, Ii_b, I)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, II_b, Ir_b, I)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                - nub*II_b
                - nub*pa(kbb, II_b, II_b, I)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kbb, II_b, Ii_b, I)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                );

    IR_d_t = + bd*pa(kdd, IS_d, SR_d, S)*cc(ddi, ddd, ddb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdd, iS_d, SR_d, S)*cc(ddi, ddd, ddb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbd, IS_b, SR_d, S)*cc(bdi, bdd, bdb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbd, iS_b, SR_d, S)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + gd*II_d
      - gd*IR_d
      - dd*IR_d
      + lm*Ri_d
      + lm*Ir_d
      - om*IR_d
      - al*pa(kid, sI_i, IR_d, I)*cc(idi, idd, idb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iI_i, IR_d, I)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rI_i, IR_d, I)*cc(idi, idd, idb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, IR_d, Rs_i, R)*cc(dii, did, dib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kdi, IR_d, Rr_i, R)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, II_i, IR_d, I)*cc(idi, idd, idb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kdi, IR_d, RI_i, R)*cc(dii, did, dib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iI_i, IR_d, I)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sI_b, IR_d, I)*cc(bdi, bdd, bdb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iI_b, IR_d, I)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rI_b, IR_d, I)*cc(bdi, bdd, bdb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, IR_d, Rs_b, R)*cc(dbi, dbd, dbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, IR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kdb, IR_d, Rr_b, R)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, II_b, IR_d, I)*cc(bdi, bdd, bdb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kdb, IR_d, RI_b, R)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iI_b, IR_d, I)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kdb, IR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb);

    IR_i_t = + bd*pa(kdi, IS_d, SR_i, S)*cc(dii, did, dib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdi, iS_d, SR_i, S)*cc(dii, did, dib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbi, IS_b, SR_i, S)*cc(bii, bid, bib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbi, iS_b, SR_i, S)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + gd*II_i
      - gd*IR_i
      - dd*IR_i
      + lm*Ri_i
      + lm*Ir_i
      - om*IR_i
      - al*pa(kii, sI_i, IR_i, I)*cc(iii, iid, iib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rI_i, IR_i, I)*cc(iii, iid, iib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, II_i, IR_i, I)*cc(iii, iid, iib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*IR_i
      - nu*pa(kii, IR_i, RI_i, R)*cc(iii, iid, iib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sI_b, IR_i, I)*cc(bii, bid, bib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iI_b, IR_i, I)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rI_b, IR_i, I)*cc(bii, bid, bib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, IR_i, Rs_b, R)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, IR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kib, IR_i, Rr_b, R)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, II_b, IR_i, I)*cc(bii, bid, bib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kib, IR_i, RI_b, R)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iI_b, IR_i, I)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kib, IR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb);

    IR_b_t = + bd*pa(kdb, IS_d, SR_b, S)*cc(dbi, dbd, dbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdb, iS_d, SR_b, S)*cc(dbi, dbd, dbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbb, IS_b, SR_b, S)*cc(bbi, bbd, bbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbb, iS_b, SR_b, S)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + gd*II_b
      - gd*IR_b
      - dd*IR_b
      + lm*Ri_b
      + lm*Ir_b
      - om*IR_b
      - al*pa(kib, sI_i, IR_b, I)*cc(ibi, ibd, ibb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iI_i, IR_b, I)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rI_i, IR_b, I)*cc(ibi, ibd, ibb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, IR_b, Rs_i, R)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, IR_b, Ri_i, R)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kbi, IR_b, Rr_i, R)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, II_i, IR_b, I)*cc(ibi, ibd, ibb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kbi, IR_b, RI_i, R)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iI_i, IR_b, I)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kbi, IR_b, Ri_i, R)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sI_b, IR_b, I)*cc(bbi, bbd, bbb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iI_b, IR_b, I)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rI_b, IR_b, I)*cc(bbi, bbd, bbb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, IR_b, Rs_b, R)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, IR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, IR_b, Rr_b, R)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, II_b, IR_b, I)*cc(bbi, bbd, bbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*IR_b
      - nub*pa(kbb, IR_b, RI_b, R)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iI_b, IR_b, I)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb);

    Is_d_t = + bd*pa(kdd, IS_d, Ss_d, S)*cc(ddi, ddd, ddb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bi*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, ddb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdd, iS_d, Ss_d, S)*cc(ddi, ddd, ddb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2*Is_d
      - b2*pa(kdd, Is_d, sI_d, s)*cc(ddi, ddd, ddb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbd, IS_b, Ss_d, S)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kdb, Is_d, si_b, s)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbd, iS_b, Ss_d, S)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kdb, Is_d, sI_b, s)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - gd*Is_d
      + di*Ir_d
      - lm*Is_d
      + lm*si_d
      - om*Is_d
      + al*pa(kid, sS_i, SI_d, S)*cc(idi, idd, idb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iS_i, SI_d, S)*cc(idi, idd, idb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rS_i, SI_d, S)*cc(idi, idd, idb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, sI_i, Is_d, I)*cc(idi, idd, idb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iI_i, Is_d, I)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rI_i, Is_d, I)*cc(idi, idd, idb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, IS_i, SI_d, S)*cc(idi, idd, idb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, II_i, Is_d, I)*cc(idi, idd, idb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iS_i, SI_d, S)*cc(idi, idd, idb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iI_i, Is_d, I)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sS_b, SI_d, S)*cc(bdi, bdd, bdb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iS_b, SI_d, S)*cc(bdi, bdd, bdb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rS_b, SI_d, S)*cc(bdi, bdd, bdb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sI_b, Is_d, I)*cc(bdi, bdd, bdb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iI_b, Is_d, I)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rI_b, Is_d, I)*cc(bdi, bdd, bdb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, IS_b, SI_d, S)*cc(bdi, bdd, bdb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, II_b, Is_d, I)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iS_b, SI_d, S)*cc(bdi, bdd, bdb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iI_b, Is_d, I)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Is_i_t = + bd*pa(kdi, IS_d, Ss_i, S)*cc(dii, did, dib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bi*pa(kid, Is_i, si_d, s)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdi, iS_d, Ss_i, S)*cc(dii, did, dib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2*pa(kid, Is_i, sI_d, s)*cc(idi, idd, idb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbi, IS_b, Ss_i, S)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kib, Is_i, si_b, s)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbi, iS_b, Ss_i, S)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kib, Is_i, sI_b, s)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - gd*Is_i
      + di*Ir_i
      - lm*Is_i
      + lm*si_i
      - om*Is_i
      + al*pa(kii, sS_i, SI_i, S)*cc(iii, iid, iib, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rS_i, SI_i, S)*cc(iii, iid, iib, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - al*Is_i
      - al*pa(kii, sI_i, Is_i, I)*cc(iii, iid, iib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iI_i, Is_i, I)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rI_i, Is_i, I)*cc(iii, iid, iib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*SI_i
      + nu*pa(kii, IS_i, SI_i, S)*cc(iii, iid, iib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, II_i, Is_i, I)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iS_i, SI_i, S)*cc(iii, iid, iib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iI_i, Is_i, I)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sS_b, SI_i, S)*cc(bii, bid, bib, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iS_b, SI_i, S)*cc(bii, bid, bib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rS_b, SI_i, S)*cc(bii, bid, bib, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sI_b, Is_i, I)*cc(bii, bid, bib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iI_b, Is_i, I)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rI_b, Is_i, I)*cc(bii, bid, bib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, IS_b, SI_i, S)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, II_b, Is_i, I)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iS_b, SI_i, S)*cc(bii, bid, bib, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iI_b, Is_i, I)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Is_b_t = + bd*pa(kdb, IS_d, Ss_b, S)*cc(dbi, dbd, dbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bi*pa(kbd, Is_b, si_d, s)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdb, iS_d, Ss_b, S)*cc(dbi, dbd, dbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2*pa(kbd, Is_b, sI_d, s)*cc(bdi, bdd, bdb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbb, IS_b, Ss_b, S)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbb, Is_b, si_b, s)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbb, iS_b, Ss_b, S)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - b2b*Is_b
      - b2b*pa(kbb, Is_b, sI_b, s)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - gd*Is_b
      + di*Ir_b
      - lm*Is_b
      + lm*si_b
      - om*Is_b
      + al*pa(kib, sS_i, SI_b, S)*cc(ibi, ibd, ibb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iS_i, SI_b, S)*cc(ibi, ibd, ibb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rS_i, SI_b, S)*cc(ibi, ibd, ibb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, sI_i, Is_b, I)*cc(ibi, ibd, ibb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iI_i, Is_b, I)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rI_i, Is_b, I)*cc(ibi, ibd, ibb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, IS_i, SI_b, S)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, II_i, Is_b, I)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iS_i, SI_b, S)*cc(ibi, ibd, ibb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iI_i, Is_b, I)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sS_b, SI_b, S)*cc(bbi, bbd, bbb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iS_b, SI_b, S)*cc(bbi, bbd, bbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rS_b, SI_b, S)*cc(bbi, bbd, bbb, rI_i, rI_d, rI_b, r*I, N_Qd, N_Qi, N_Qb)
      - alb*Is_b
      - alb*pa(kbb, sI_b, Is_b, I)*cc(bbi, bbd, bbb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iI_b, Is_b, I)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rI_b, Is_b, I)*cc(bbi, bbd, bbb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*SI_b
      + nub*pa(kbb, IS_b, SI_b, S)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, II_b, Is_b, I)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iS_b, SI_b, S)*cc(bbi, bbd, bbb, iI_i, iI_d, iI_b, i*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iI_b, Is_b, I)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Ii_d_t = + bd*pa(kdd, IS_d, Si_d, S)*cc(ddi, ddd, ddb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bi*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, ddb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1*Si_d
      + b1*pa(kdd, iS_d, Si_d, S)*cc(ddi, ddd, ddb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*Is_d
      + b2*pa(kdd, Is_d, sI_d, s)*cc(ddi, ddd, ddb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbd, IS_b, Si_d, S)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kdb, Is_d, si_b, s)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbd, iS_b, Si_d, S)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kdb, Is_d, sI_b, s)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - gd*Ii_d
      - gi*Ii_d
      - lm*Ii_d
      + lm*ii_d
      + om*II_d
      - om*Ii_d
      + al*pa(kdi, II_d, Is_i, I)*cc(dii, did, dib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, sI_i, Ii_d, I)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, II_d, Ii_i, I)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, II_d, Ir_i, I)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rI_i, Ii_d, I)*cc(idi, idd, idb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, II_d, II_i, I)*cc(dii, did, dib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, II_i, Ii_d, I)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, II_d, Ii_i, I)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, II_d, Is_b, I)*cc(dbi, dbd, dbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sI_b, Ii_d, I)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, II_d, Ii_b, I)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iI_b, Ii_d, I)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, II_d, Ir_b, I)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rI_b, Ii_d, I)*cc(bdi, bdd, bdb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, II_d, II_b, I)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, II_b, Ii_d, I)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, II_d, Ii_b, I)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iI_b, Ii_d, I)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    Ii_i_t = + bd*pa(kdi, IS_d, Si_i, S)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bi*pa(kid, Is_i, si_d, s)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdi, iS_d, Si_i, S)*cc(dii, did, dib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kid, Is_i, sI_d, s)*cc(idi, idd, idb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbi, IS_b, Si_i, S)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kib, Is_i, si_b, s)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbi, iS_b, Si_i, S)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kib, Is_i, sI_b, s)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - gd*Ii_i
      - gi*Ii_i
      - lm*Ii_i
      + lm*ii_i
      + om*II_i
      - om*Ii_i
      + al*pa(kii, II_i, Is_i, I)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, sI_i, Ii_i, I)*cc(iii, iid, iib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*Ii_i
      + al*pa(kii, II_i, Ii_i, I)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, II_i, Ir_i, I)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rI_i, Ii_i, I)*cc(iii, iid, iib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*II_i
      + nu*pa(kii, II_i, II_i, I)*cc(iii, iid, iib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nu*Ii_i
      + nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, II_i, Is_b, I)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sI_b, Ii_i, I)*cc(bii, bid, bib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, II_i, Ii_b, I)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iI_b, Ii_i, I)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, II_i, Ir_b, I)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rI_b, Ii_i, I)*cc(bii, bid, bib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, II_i, II_b, I)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, II_b, Ii_i, I)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, II_i, Ii_b, I)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iI_b, Ii_i, I)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    Ii_b_t = + bd*pa(kdb, IS_d, Si_b, S)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bi*pa(kbd, Is_b, si_d, s)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdb, iS_d, Si_b, S)*cc(dbi, dbd, dbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kbd, Is_b, sI_d, s)*cc(bdi, bdd, bdb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbb, IS_b, Si_b, S)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kbb, Is_b, si_b, s)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + b1b*Si_b
      + b1b*pa(kbb, iS_b, Si_b, S)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*Is_b
      + b2b*pa(kbb, Is_b, sI_b, s)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - gd*Ii_b
      - gi*Ii_b
      - lm*Ii_b
      + lm*ii_b
      + om*II_b
      - om*Ii_b
      + al*pa(kbi, II_b, Is_i, I)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, sI_i, Ii_b, I)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, II_b, Ii_i, I)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iI_i, Ii_b, I)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, II_b, Ir_i, I)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rI_i, Ii_b, I)*cc(ibi, ibd, ibb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, II_b, II_i, I)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, II_i, Ii_b, I)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, II_b, Ii_i, I)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iI_i, Ii_b, I)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, II_b, Is_b, I)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sI_b, Ii_b, I)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*Ii_b
      + alb*pa(kbb, II_b, Ii_b, I)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iI_b, Ii_b, I)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, II_b, Ir_b, I)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rI_b, Ii_b, I)*cc(bbi, bbd, bbb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*II_b
      + nub*pa(kbb, II_b, II_b, I)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, II_b, Ii_b, I)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nub*Ii_b
      + nub*pa(kbb, II_b, Ii_b, I)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iI_b, Ii_b, I)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    Ir_d_t = + bd*pa(kdd, IS_d, Sr_d, S)*cc(ddi, ddd, ddb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdd, iS_d, Sr_d, S)*cc(ddi, ddd, ddb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbd, IS_b, Sr_d, S)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbd, iS_b, Sr_d, S)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - gd*Ir_d
      + gi*Ii_d
      - di*Ir_d
      + lm*ir_d
      - lm*Ir_d
      - om*Ir_d
      - al*pa(kid, sI_i, Ir_d, I)*cc(idi, idd, idb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rI_i, Ir_d, I)*cc(idi, idd, idb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, IR_d, Rs_i, R)*cc(dii, did, dib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, IR_d, Rr_i, R)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, II_i, Ir_d, I)*cc(idi, idd, idb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, IR_d, RI_i, R)*cc(dii, did, dib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, IR_d, Ri_i, R)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sI_b, Ir_d, I)*cc(bdi, bdd, bdb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iI_b, Ir_d, I)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rI_b, Ir_d, I)*cc(bdi, bdd, bdb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, IR_d, Rs_b, R)*cc(dbi, dbd, dbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, IR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, IR_d, Rr_b, R)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, II_b, Ir_d, I)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, IR_d, RI_b, R)*cc(dbi, dbd, dbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iI_b, Ir_d, I)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, IR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb);

    Ir_i_t = + bd*pa(kdi, IS_d, Sr_i, S)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdi, iS_d, Sr_i, S)*cc(dii, did, dib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbi, IS_b, Sr_i, S)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbi, iS_b, Sr_i, S)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - gd*Ir_i
      + gi*Ii_i
      - di*Ir_i
      + lm*ir_i
      - lm*Ir_i
      - om*Ir_i
      - al*pa(kii, sI_i, Ir_i, I)*cc(iii, iid, iib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*Ir_i
      - al*pa(kii, rI_i, Ir_i, I)*cc(iii, iid, iib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, II_i, Ir_i, I)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*IR_i
      + nu*pa(kii, IR_i, RI_i, R)*cc(iii, iid, iib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sI_b, Ir_i, I)*cc(bii, bid, bib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iI_b, Ir_i, I)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rI_b, Ir_i, I)*cc(bii, bid, bib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, IR_i, Rs_b, R)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, IR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, IR_i, Rr_b, R)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, II_b, Ir_i, I)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, IR_i, RI_b, R)*cc(ibi, ibd, ibb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iI_b, Ir_i, I)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, IR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb);

    Ir_b_t = + bd*pa(kdb, IS_d, Sr_b, S)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + b1*pa(kdb, iS_d, Sr_b, S)*cc(dbi, dbd, dbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + bdbr*pa(kbb, IS_b, Sr_b, S)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + b1b*pa(kbb, iS_b, Sr_b, S)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - gd*Ir_b
      + gi*Ii_b
      - di*Ir_b
      + lm*ir_b
      - lm*Ir_b
      - om*Ir_b
      - al*pa(kib, sI_i, Ir_b, I)*cc(ibi, ibd, ibb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iI_i, Ir_b, I)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rI_i, Ir_b, I)*cc(ibi, ibd, ibb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, IR_b, Rs_i, R)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, IR_b, Ri_i, R)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, IR_b, Rr_i, R)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, II_i, Ir_b, I)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, IR_b, RI_i, R)*cc(bii, bid, bib, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iI_i, Ir_b, I)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, IR_b, Ri_i, R)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sI_b, Ir_b, I)*cc(bbi, bbd, bbb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iI_b, Ir_b, I)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*Ir_b
      - alb*pa(kbb, rI_b, Ir_b, I)*cc(bbi, bbd, bbb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, IR_b, Rs_b, R)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, IR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, IR_b, Rr_b, R)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, II_b, Ir_b, I)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*IR_b
      + nub*pa(kbb, IR_b, RI_b, R)*cc(bbi, bbd, bbb, II_i, II_d, II_b, I*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iI_b, Ir_b, I)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, IR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb);

    RR_d_t = 2*(
                + gd*IR_d
                - dd*RR_d
                + lm*Rr_d
                - al*pa(kdi, RR_d, Rs_i, R)*cc(dii, did, dib, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, dib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kdi, RR_d, Rr_i, R)*cc(dii, did, dib, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kdi, RR_d, RI_i, R)*cc(dii, did, dib, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, dib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, RR_d, Rs_b, R)*cc(dbi, dbd, dbb, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, RR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kdb, RR_d, Rr_b, R)*cc(dbi, dbd, dbb, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kdb, RR_d, RI_b, R)*cc(dbi, dbd, dbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kdb, RR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                );

    RR_i_t = 2*(
                + gd*IR_i
                - dd*RR_i
                + lm*Rr_i
                - al*pa(kii, RR_i, Rs_i, R)*cc(iii, iid, iib, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, iib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kii, RR_i, Rr_i, R)*cc(iii, iid, iib, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kii, RR_i, RI_i, R)*cc(iii, iid, iib, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, iib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, RR_i, Rs_b, R)*cc(ibi, ibd, ibb, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, RR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kib, RR_i, Rr_b, R)*cc(ibi, ibd, ibb, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kib, RR_i, RI_b, R)*cc(ibi, ibd, ibb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kib, RR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                );

    RR_b_t = 2*(
                + gd*IR_b
                - dd*RR_b
                + lm*Rr_b
                - al*pa(kbi, RR_b, Rs_i, R)*cc(bii, bid, bib, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
                - al*pa(kbi, RR_b, Ri_i, R)*cc(bii, bid, bib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - al*pa(kbi, RR_b, Rr_i, R)*cc(bii, bid, bib, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
                - nu*pa(kbi, RR_b, RI_i, R)*cc(bii, bid, bib, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
                - nu*pa(kbi, RR_b, Ri_i, R)*cc(bii, bid, bib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, RR_b, Rs_b, R)*cc(bbi, bbd, bbb, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, RR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                - alb*pa(kbb, RR_b, Rr_b, R)*cc(bbi, bbd, bbb, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
                - nub*pa(kbb, RR_b, RI_b, R)*cc(bbi, bbd, bbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
                - nub*pa(kbb, RR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
                );

    Rs_d_t = - bi*pa(kdd, Rs_d, si_d, s)*cc(ddi, ddd, ddb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdd, Rs_d, sI_d, s)*cc(ddi, ddd, ddb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kdb, Rs_d, si_b, s)*cc(dbi, dbd, dbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kdb, Rs_d, sI_b, s)*cc(dbi, dbd, dbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + gd*Is_d
      - dd*Rs_d
      + di*Rr_d
      - lm*Rs_d
      + lm*sr_d
      + al*pa(kid, sS_i, SR_d, S)*cc(idi, idd, idb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iS_i, SR_d, S)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rS_i, SR_d, S)*cc(idi, idd, idb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, sR_i, Rs_d, R)*cc(idi, idd, idb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rR_i, Rs_d, R)*cc(idi, idd, idb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, IS_i, SR_d, S)*cc(idi, idd, idb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IR_i, Rs_d, R)*cc(idi, idd, idb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iS_i, SR_d, S)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sS_b, SR_d, S)*cc(bdi, bdd, bdb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iS_b, SR_d, S)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rS_b, SR_d, S)*cc(bdi, bdd, bdb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sR_b, Rs_d, R)*cc(bdi, bdd, bdb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iR_b, Rs_d, R)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rR_b, Rs_d, R)*cc(bdi, bdd, bdb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, IS_b, SR_d, S)*cc(bdi, bdd, bdb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IR_b, Rs_d, R)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iS_b, SR_d, S)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iR_b, Rs_d, R)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Rs_i_t = - bi*pa(kid, Rs_i, si_d, s)*cc(idi, idd, idb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - b2*pa(kid, Rs_i, sI_d, s)*cc(idi, idd, idb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kib, Rs_i, si_b, s)*cc(ibi, ibd, ibb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kib, Rs_i, sI_b, s)*cc(ibi, ibd, ibb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + gd*Is_i
      - dd*Rs_i
      + di*Rr_i
      - lm*Rs_i
      + lm*sr_i
      + al*pa(kii, sS_i, SR_i, S)*cc(iii, iid, iib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rS_i, SR_i, S)*cc(iii, iid, iib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*Rs_i
      - al*pa(kii, sR_i, Rs_i, R)*cc(iii, iid, iib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rR_i, Rs_i, R)*cc(iii, iid, iib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, IS_i, SR_i, S)*cc(iii, iid, iib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iS_i, SR_i, S)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sS_b, SR_i, S)*cc(bii, bid, bib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iS_b, SR_i, S)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rS_b, SR_i, S)*cc(bii, bid, bib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sR_b, Rs_i, R)*cc(bii, bid, bib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iR_b, Rs_i, R)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rR_b, Rs_i, R)*cc(bii, bid, bib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, IS_b, SR_i, S)*cc(bii, bid, bib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IR_b, Rs_i, R)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iS_b, SR_i, S)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iR_b, Rs_i, R)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Rs_b_t = - bi*pa(kbd, Rs_b, si_d, s)*cc(bdi, bdd, bdb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - b2*pa(kbd, Rs_b, sI_d, s)*cc(bdi, bdd, bdb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbb, Rs_b, si_b, s)*cc(bbi, bbd, bbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbb, Rs_b, sI_b, s)*cc(bbi, bbd, bbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + gd*Is_b
      - dd*Rs_b
      + di*Rr_b
      - lm*Rs_b
      + lm*sr_b
      + al*pa(kib, sS_i, SR_b, S)*cc(ibi, ibd, ibb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iS_i, SR_b, S)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rS_i, SR_b, S)*cc(ibi, ibd, ibb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, sR_i, Rs_b, R)*cc(ibi, ibd, ibb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iR_i, Rs_b, R)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rR_i, Rs_b, R)*cc(ibi, ibd, ibb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, IS_i, SR_b, S)*cc(ibi, ibd, ibb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IR_i, Rs_b, R)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iS_i, SR_b, S)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iR_i, Rs_b, R)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sS_b, SR_b, S)*cc(bbi, bbd, bbb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iS_b, SR_b, S)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rS_b, SR_b, S)*cc(bbi, bbd, bbb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*Rs_b
      - alb*pa(kbb, sR_b, Rs_b, R)*cc(bbi, bbd, bbb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iR_b, Rs_b, R)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rR_b, Rs_b, R)*cc(bbi, bbd, bbb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, IS_b, SR_b, S)*cc(bbi, bbd, bbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IR_b, Rs_b, R)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iS_b, SR_b, S)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iR_b, Rs_b, R)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    Ri_d_t = + bi*pa(kdd, Rs_d, si_d, s)*cc(ddi, ddd, ddb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kdd, Rs_d, sI_d, s)*cc(ddi, ddd, ddb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kdb, Rs_d, si_b, s)*cc(dbi, dbd, dbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kdb, Rs_d, sI_b, s)*cc(dbi, dbd, dbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + gd*Ii_d
      - gi*Ri_d
      - dd*Ri_d
      - lm*Ri_d
      + lm*ir_d
      + om*IR_d
      + al*pa(kid, sI_i, IR_d, I)*cc(idi, idd, idb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iI_i, IR_d, I)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rI_i, IR_d, I)*cc(idi, idd, idb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, sR_i, Ri_d, R)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rR_i, Ri_d, R)*cc(idi, idd, idb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, II_i, IR_d, I)*cc(idi, idd, idb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IR_i, Ri_d, R)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iI_i, IR_d, I)*cc(idi, idd, idb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sI_b, IR_d, I)*cc(bdi, bdd, bdb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iI_b, IR_d, I)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rI_b, IR_d, I)*cc(bdi, bdd, bdb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sR_b, Ri_d, R)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iR_b, Ri_d, R)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rR_b, Ri_d, R)*cc(bdi, bdd, bdb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, II_b, IR_d, I)*cc(bdi, bdd, bdb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IR_b, Ri_d, R)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iI_b, IR_d, I)*cc(bdi, bdd, bdb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iR_b, Ri_d, R)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    Ri_i_t = + bi*pa(kid, Rs_i, si_d, s)*cc(idi, idd, idb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kid, Rs_i, sI_d, s)*cc(idi, idd, idb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kib, Rs_i, si_b, s)*cc(ibi, ibd, ibb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kib, Rs_i, sI_b, s)*cc(ibi, ibd, ibb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + gd*Ii_i
      - gi*Ri_i
      - dd*Ri_i
      - lm*Ri_i
      + lm*ir_i
      + om*IR_i
      + al*pa(kii, sI_i, IR_i, I)*cc(iii, iid, iib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rI_i, IR_i, I)*cc(iii, iid, iib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, sR_i, Ri_i, R)*cc(iii, iid, iib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*Ri_i
      - al*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rR_i, Ri_i, R)*cc(iii, iid, iib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, II_i, IR_i, I)*cc(iii, iid, iib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iI_i, IR_i, I)*cc(iii, iid, iib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*Ri_i
      - nu*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sI_b, IR_i, I)*cc(bii, bid, bib, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iI_b, IR_i, I)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rI_b, IR_i, I)*cc(bii, bid, bib, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sR_b, Ri_i, R)*cc(bii, bid, bib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iR_b, Ri_i, R)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rR_b, Ri_i, R)*cc(bii, bid, bib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, II_b, IR_i, I)*cc(bii, bid, bib, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IR_b, Ri_i, R)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iI_b, IR_i, I)*cc(bii, bid, bib, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iR_b, Ri_i, R)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    Ri_b_t = + bi*pa(kbd, Rs_b, si_d, s)*cc(bdi, bdd, bdb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kbd, Rs_b, sI_d, s)*cc(bdi, bdd, bdb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kbb, Rs_b, si_b, s)*cc(bbi, bbd, bbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kbb, Rs_b, sI_b, s)*cc(bbi, bbd, bbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      + gd*Ii_b
      - gi*Ri_b
      - dd*Ri_b
      - lm*Ri_b
      + lm*ir_b
      + om*IR_b
      + al*pa(kib, sI_i, IR_b, I)*cc(ibi, ibd, ibb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iI_i, IR_b, I)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rI_i, IR_b, I)*cc(ibi, ibd, ibb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, sR_i, Ri_b, R)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iR_i, Ri_b, R)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rR_i, Ri_b, R)*cc(ibi, ibd, ibb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, II_i, IR_b, I)*cc(ibi, ibd, ibb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iI_i, IR_b, I)*cc(ibi, ibd, ibb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iR_i, Ri_b, R)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sI_b, IR_b, I)*cc(bbi, bbd, bbb, sR_i, sR_d, sR_b, s*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iI_b, IR_b, I)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rI_b, IR_b, I)*cc(bbi, bbd, bbb, rR_i, rR_d, rR_b, r*R, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sR_b, Ri_b, R)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - alb*Ri_b
      - alb*pa(kbb, iR_b, Ri_b, R)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rR_b, Ri_b, R)*cc(bbi, bbd, bbb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, II_b, IR_b, I)*cc(bbi, bbd, bbb, IR_i, IR_d, IR_b, I*R, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iI_b, IR_b, I)*cc(bbi, bbd, bbb, iR_i, iR_d, iR_b, i*R, N_Qd, N_Qi, N_Qb)
      - nub*Ri_b
      - nub*pa(kbb, iR_b, Ri_b, R)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    Rr_d_t = + gd*Ir_d
      + gi*Ri_d
      - dd*Rr_d
      - di*Rr_d
      - lm*Rr_d
      + lm*rr_d
      + al*pa(kdi, RR_d, Rs_i, R)*cc(dii, did, dib, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, sR_i, Rr_d, R)*cc(idi, idd, idb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, dib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kdi, RR_d, Rr_i, R)*cc(dii, did, dib, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kid, rR_i, Rr_d, R)*cc(idi, idd, idb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, RR_d, RI_i, R)*cc(dii, did, dib, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, IR_i, Rr_d, R)*cc(idi, idd, idb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kdi, RR_d, Ri_i, R)*cc(dii, did, dib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, RR_d, Rs_b, R)*cc(dbi, dbd, dbb, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, sR_b, Rr_d, R)*cc(bdi, bdd, bdb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, RR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, iR_b, Rr_d, R)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kdb, RR_d, Rr_b, R)*cc(dbi, dbd, dbb, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbd, rR_b, Rr_d, R)*cc(bdi, bdd, bdb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, RR_d, RI_b, R)*cc(dbi, dbd, dbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, IR_b, Rr_d, R)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kdb, RR_d, Ri_b, R)*cc(dbi, dbd, dbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbd, iR_b, Rr_d, R)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb);

    Rr_i_t = + gd*Ir_i
      + gi*Ri_i
      - dd*Rr_i
      - di*Rr_i
      - lm*Rr_i
      + lm*rr_i
      + al*pa(kii, RR_i, Rs_i, R)*cc(iii, iid, iib, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, sR_i, Rr_i, R)*cc(iii, iid, iib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, iib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - al*Rr_i
      + al*pa(kii, RR_i, Rr_i, R)*cc(iii, iid, iib, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kii, rR_i, Rr_i, R)*cc(iii, iid, iib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, RR_i, RI_i, R)*cc(iii, iid, iib, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, RR_i, Ri_i, R)*cc(iii, iid, iib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, RR_i, Rs_b, R)*cc(ibi, ibd, ibb, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, sR_b, Rr_i, R)*cc(bii, bid, bib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, RR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, iR_b, Rr_i, R)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kib, RR_i, Rr_b, R)*cc(ibi, ibd, ibb, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbi, rR_b, Rr_i, R)*cc(bii, bid, bib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, RR_i, RI_b, R)*cc(ibi, ibd, ibb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, IR_b, Rr_i, R)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kib, RR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbi, iR_b, Rr_i, R)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb);

    Rr_b_t = + gd*Ir_b
      + gi*Ri_b
      - dd*Rr_b
      - di*Rr_b
      - lm*Rr_b
      + lm*rr_b
      + al*pa(kbi, RR_b, Rs_i, R)*cc(bii, bid, bib, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, sR_i, Rr_b, R)*cc(ibi, ibd, ibb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, RR_b, Ri_i, R)*cc(bii, bid, bib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, iR_i, Rr_b, R)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kbi, RR_b, Rr_i, R)*cc(bii, bid, bib, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
      - al*pa(kib, rR_i, Rr_b, R)*cc(ibi, ibd, ibb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, RR_b, RI_i, R)*cc(bii, bid, bib, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, IR_i, Rr_b, R)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kbi, RR_b, Ri_i, R)*cc(bii, bid, bib, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - nu*pa(kib, iR_i, Rr_b, R)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, RR_b, Rs_b, R)*cc(bbi, bbd, bbb, Rs_i, Rs_d, Rs_b, R*s, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, sR_b, Rr_b, R)*cc(bbi, bbd, bbb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, RR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, iR_b, Rr_b, R)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - alb*Rr_b
      + alb*pa(kbb, RR_b, Rr_b, R)*cc(bbi, bbd, bbb, Rr_i, Rr_d, Rr_b, R*r, N_Qd, N_Qi, N_Qb)
      - alb*pa(kbb, rR_b, Rr_b, R)*cc(bbi, bbd, bbb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, RR_b, RI_b, R)*cc(bbi, bbd, bbb, RI_i, RI_d, RI_b, R*I, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, IR_b, Rr_b, R)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, RR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ri_i, Ri_d, Ri_b, R*i, N_Qd, N_Qi, N_Qb)
      - nub*pa(kbb, iR_b, Rr_b, R)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb);

    ss_d_t = 2*(
                - bi*pa(kdd, ss_d, si_d, s)*cc(ddi, ddd, ddb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                - b2*pa(kdd, ss_d, sI_d, s)*cc(ddi, ddd, ddb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
                - bibr*pa(kdb, ss_d, si_b, s)*cc(dbi, dbd, dbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                - b2b*pa(kdb, ss_d, sI_b, s)*cc(dbi, dbd, dbb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
                + di*sr_d
                - lm*ss_d
                + al*pa(kid, sS_i, Ss_d, S)*cc(idi, idd, idb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
                + al*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + al*pa(kid, rS_i, Ss_d, S)*cc(idi, idd, idb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
                + nu*pa(kid, IS_i, Ss_d, S)*cc(idi, idd, idb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                + nu*pa(kid, iS_i, Ss_d, S)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, sS_b, Ss_d, S)*cc(bdi, bdd, bdb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, iS_b, Ss_d, S)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, rS_b, Ss_d, S)*cc(bdi, bdd, bdb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbd, IS_b, Ss_d, S)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbd, iS_b, Ss_d, S)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                );

    ss_i_t = 2*(
                - bi*pa(kid, ss_i, si_d, s)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                - b2*pa(kid, ss_i, sI_d, s)*cc(idi, idd, idb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
                - bibr*pa(kib, ss_i, si_b, s)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                - b2b*pa(kib, ss_i, sI_b, s)*cc(ibi, ibd, ibb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
                + di*sr_i
                - lm*ss_i
                + al*Ss_i
                + al*pa(kii, sS_i, Ss_i, S)*cc(iii, iid, iib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
                + al*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + al*pa(kii, rS_i, Ss_i, S)*cc(iii, iid, iib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
                + nu*pa(kii, IS_i, Ss_i, S)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                + nu*pa(kii, iS_i, Ss_i, S)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, sS_b, Ss_i, S)*cc(bii, bid, bib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, iS_b, Ss_i, S)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, rS_b, Ss_i, S)*cc(bii, bid, bib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbi, IS_b, Ss_i, S)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbi, iS_b, Ss_i, S)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                );

    ss_b_t = 2*(
                - bi*pa(kbd, ss_b, si_d, s)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                - b2*pa(kbd, ss_b, sI_d, s)*cc(bdi, bdd, bdb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
                - bibr*pa(kbb, ss_b, si_b, s)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                - b2b*pa(kbb, ss_b, sI_b, s)*cc(bbi, bbd, bbb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
                + di*sr_b
                - lm*ss_b
                + al*pa(kib, sS_i, Ss_b, S)*cc(ibi, ibd, ibb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
                + al*pa(kib, iS_i, Ss_b, S)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + al*pa(kib, rS_i, Ss_b, S)*cc(ibi, ibd, ibb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
                + nu*pa(kib, IS_i, Ss_b, S)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                + nu*pa(kib, iS_i, Ss_b, S)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + alb*Ss_b
                + alb*pa(kbb, sS_b, Ss_b, S)*cc(bbi, bbd, bbb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbb, iS_b, Ss_b, S)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbb, rS_b, Ss_b, S)*cc(bbi, bbd, bbb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbb, IS_b, Ss_b, S)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbb, iS_b, Ss_b, S)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
                );

    si_d_t = - bi*si_d
      + bi*pa(kdd, ss_d, si_d, s)*cc(ddi, ddd, ddb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - bi*pa(kdd, is_d, si_d, s)*cc(ddi, ddd, ddb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kdd, ss_d, sI_d, s)*cc(ddi, ddd, ddb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, ddb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kdb, ss_d, si_b, s)*cc(dbi, dbd, dbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbd, is_b, si_d, s)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kdb, ss_d, sI_b, s)*cc(dbi, dbd, dbb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbd, Is_b, si_d, s)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - gi*si_d
      + di*ir_d
      - lm*si_d
      - lm*si_d
      + om*Is_d
      + al*pa(kid, sS_i, Si_d, S)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iS_i, Si_d, S)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rS_i, Si_d, S)*cc(idi, idd, idb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, sI_i, Is_d, I)*cc(idi, idd, idb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iI_i, Is_d, I)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rI_i, Is_d, I)*cc(idi, idd, idb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, IS_i, Si_d, S)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, II_i, Is_d, I)*cc(idi, idd, idb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iS_i, Si_d, S)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iI_i, Is_d, I)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sS_b, Si_d, S)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iS_b, Si_d, S)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rS_b, Si_d, S)*cc(bdi, bdd, bdb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sI_b, Is_d, I)*cc(bdi, bdd, bdb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iI_b, Is_d, I)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rI_b, Is_d, I)*cc(bdi, bdd, bdb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, IS_b, Si_d, S)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, II_b, Is_d, I)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iS_b, Si_d, S)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iI_b, Is_d, I)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    si_i_t = + bi*pa(kid, ss_i, si_d, s)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - bi*pa(kdi, is_d, si_i, s)*cc(dii, did, dib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kid, ss_i, sI_d, s)*cc(idi, idd, idb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdi, Is_d, si_i, s)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kib, ss_i, si_b, s)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbi, is_b, si_i, s)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kib, ss_i, sI_b, s)*cc(ibi, ibd, ibb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbi, Is_b, si_i, s)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - gi*si_i
      + di*ir_i
      - lm*si_i
      - lm*si_i
      + om*Is_i
      + al*pa(kii, sS_i, Si_i, S)*cc(iii, iid, iib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*Si_i
      + al*pa(kii, iS_i, Si_i, S)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rS_i, Si_i, S)*cc(iii, iid, iib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + al*Is_i
      + al*pa(kii, sI_i, Is_i, I)*cc(iii, iid, iib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iI_i, Is_i, I)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rI_i, Is_i, I)*cc(iii, iid, iib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, IS_i, Si_i, S)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, II_i, Is_i, I)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*Si_i
      + nu*pa(kii, iS_i, Si_i, S)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iI_i, Is_i, I)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sS_b, Si_i, S)*cc(bii, bid, bib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iS_b, Si_i, S)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rS_b, Si_i, S)*cc(bii, bid, bib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sI_b, Is_i, I)*cc(bii, bid, bib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iI_b, Is_i, I)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rI_b, Is_i, I)*cc(bii, bid, bib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, IS_b, Si_i, S)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, II_b, Is_i, I)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iS_b, Si_i, S)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iI_b, Is_i, I)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    si_b_t = + bi*pa(kbd, ss_b, si_d, s)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - bi*pa(kdb, is_d, si_b, s)*cc(dbi, dbd, dbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2*pa(kbd, ss_b, sI_d, s)*cc(bdi, bdd, bdb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdb, Is_d, si_b, s)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - bibr*si_b
      + bibr*pa(kbb, ss_b, si_b, s)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbb, is_b, si_b, s)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kbb, ss_b, sI_b, s)*cc(bbi, bbd, bbb, sI_i, sI_d, sI_b, s*I, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbb, Is_b, si_b, s)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      - gi*si_b
      + di*ir_b
      - lm*si_b
      - lm*si_b
      + om*Is_b
      + al*pa(kib, sS_i, Si_b, S)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iS_i, Si_b, S)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rS_i, Si_b, S)*cc(ibi, ibd, ibb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, sI_i, Is_b, I)*cc(ibi, ibd, ibb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iI_i, Is_b, I)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rI_i, Is_b, I)*cc(ibi, ibd, ibb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, IS_i, Si_b, S)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, II_i, Is_b, I)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iS_i, Si_b, S)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iI_i, Is_b, I)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sS_b, Si_b, S)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*Si_b
      + alb*pa(kbb, iS_b, Si_b, S)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rS_b, Si_b, S)*cc(bbi, bbd, bbb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + alb*Is_b
      + alb*pa(kbb, sI_b, Is_b, I)*cc(bbi, bbd, bbb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iI_b, Is_b, I)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rI_b, Is_b, I)*cc(bbi, bbd, bbb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, IS_b, Si_b, S)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, II_b, Is_b, I)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*Si_b
      + nub*pa(kbb, iS_b, Si_b, S)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iI_b, Is_b, I)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    sr_d_t = - bi*pa(kdd, is_d, sr_d, s)*cc(ddi, ddd, ddb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdd, Is_d, sr_d, s)*cc(ddi, ddd, ddb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbd, is_b, sr_d, s)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbd, Is_b, sr_d, s)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + gi*si_d
      - di*sr_d
      + di*rr_d
      - lm*sr_d
      - lm*sr_d
      + al*pa(kid, sS_i, Sr_d, S)*cc(idi, idd, idb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rS_i, Sr_d, S)*cc(idi, idd, idb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, sR_i, Rs_d, R)*cc(idi, idd, idb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rR_i, Rs_d, R)*cc(idi, idd, idb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, IS_i, Sr_d, S)*cc(idi, idd, idb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, IR_i, Rs_d, R)*cc(idi, idd, idb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iS_i, Sr_d, S)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iR_i, Rs_d, R)*cc(idi, idd, idb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sS_b, Sr_d, S)*cc(bdi, bdd, bdb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iS_b, Sr_d, S)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rS_b, Sr_d, S)*cc(bdi, bdd, bdb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sR_b, Rs_d, R)*cc(bdi, bdd, bdb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iR_b, Rs_d, R)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rR_b, Rs_d, R)*cc(bdi, bdd, bdb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, IS_b, Sr_d, S)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, IR_b, Rs_d, R)*cc(bdi, bdd, bdb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iS_b, Sr_d, S)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iR_b, Rs_d, R)*cc(bdi, bdd, bdb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    sr_i_t = - bi*pa(kdi, is_d, sr_i, s)*cc(dii, did, dib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdi, Is_d, sr_i, s)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbi, is_b, sr_i, s)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbi, Is_b, sr_i, s)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + gi*si_i
      - di*sr_i
      + di*rr_i
      - lm*sr_i
      - lm*sr_i
      + al*pa(kii, sS_i, Sr_i, S)*cc(iii, iid, iib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*Sr_i
      + al*pa(kii, rS_i, Sr_i, S)*cc(iii, iid, iib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*Rs_i
      + al*pa(kii, sR_i, Rs_i, R)*cc(iii, iid, iib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rR_i, Rs_i, R)*cc(iii, iid, iib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, IS_i, Sr_i, S)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, IR_i, Rs_i, R)*cc(iii, iid, iib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iS_i, Sr_i, S)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iR_i, Rs_i, R)*cc(iii, iid, iib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sS_b, Sr_i, S)*cc(bii, bid, bib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iS_b, Sr_i, S)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rS_b, Sr_i, S)*cc(bii, bid, bib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sR_b, Rs_i, R)*cc(bii, bid, bib, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iR_b, Rs_i, R)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rR_b, Rs_i, R)*cc(bii, bid, bib, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, IS_b, Sr_i, S)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, IR_b, Rs_i, R)*cc(bii, bid, bib, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iS_b, Sr_i, S)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iR_b, Rs_i, R)*cc(bii, bid, bib, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    sr_b_t = - bi*pa(kdb, is_d, sr_b, s)*cc(dbi, dbd, dbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - b2*pa(kdb, Is_d, sr_b, s)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      - bibr*pa(kbb, is_b, sr_b, s)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      - b2b*pa(kbb, Is_b, sr_b, s)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + gi*si_b
      - di*sr_b
      + di*rr_b
      - lm*sr_b
      - lm*sr_b
      + al*pa(kib, sS_i, Sr_b, S)*cc(ibi, ibd, ibb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iS_i, Sr_b, S)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rS_i, Sr_b, S)*cc(ibi, ibd, ibb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, sR_i, Rs_b, R)*cc(ibi, ibd, ibb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iR_i, Rs_b, R)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rR_i, Rs_b, R)*cc(ibi, ibd, ibb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, IS_i, Sr_b, S)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, IR_i, Rs_b, R)*cc(ibi, ibd, ibb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iS_i, Sr_b, S)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iR_i, Rs_b, R)*cc(ibi, ibd, ibb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sS_b, Sr_b, S)*cc(bbi, bbd, bbb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iS_b, Sr_b, S)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*Sr_b
      + alb*pa(kbb, rS_b, Sr_b, S)*cc(bbi, bbd, bbb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*Rs_b
      + alb*pa(kbb, sR_b, Rs_b, R)*cc(bbi, bbd, bbb, ss_i, ss_d, ss_b, s*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iR_b, Rs_b, R)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rR_b, Rs_b, R)*cc(bbi, bbd, bbb, rs_i, rs_d, rs_b, r*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, IS_b, Sr_b, S)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, IR_b, Rs_b, R)*cc(bbi, bbd, bbb, Is_i, Is_d, Is_b, I*s, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iS_b, Sr_b, S)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iR_b, Rs_b, R)*cc(bbi, bbd, bbb, is_i, is_d, is_b, i*s, N_Qd, N_Qi, N_Qb);

    ii_d_t = 2*(
                + bi*si_d
                + bi*pa(kdd, is_d, si_d, s)*cc(ddi, ddd, ddb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + b2*pa(kdd, Is_d, si_d, s)*cc(ddi, ddd, ddb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + bibr*pa(kbd, is_b, si_d, s)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + b2b*pa(kbd, Is_b, si_d, s)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - gi*ii_d
                - lm*ii_d
                + om*Ii_d
                + al*pa(kid, sI_i, Ii_d, I)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                + al*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + al*pa(kid, rI_i, Ii_d, I)*cc(idi, idd, idb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
                + nu*pa(kid, II_i, Ii_d, I)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + nu*pa(kid, iI_i, Ii_d, I)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, sI_b, Ii_d, I)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, iI_b, Ii_d, I)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, rI_b, Ii_d, I)*cc(bdi, bdd, bdb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbd, II_b, Ii_d, I)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbd, iI_b, Ii_d, I)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                );

    ii_i_t = 2*(
                + bi*pa(kdi, is_d, si_i, s)*cc(dii, did, dib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + b2*pa(kdi, Is_d, si_i, s)*cc(dii, did, dib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + bibr*pa(kbi, is_b, si_i, s)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + b2b*pa(kbi, Is_b, si_i, s)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - gi*ii_i
                - lm*ii_i
                + om*Ii_i
                + al*pa(kii, sI_i, Ii_i, I)*cc(iii, iid, iib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                + al*Ii_i
                + al*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + al*pa(kii, rI_i, Ii_i, I)*cc(iii, iid, iib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
                + nu*pa(kii, II_i, Ii_i, I)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + nu*Ii_i
                + nu*pa(kii, iI_i, Ii_i, I)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, sI_b, Ii_i, I)*cc(bii, bid, bib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, iI_b, Ii_i, I)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, rI_b, Ii_i, I)*cc(bii, bid, bib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbi, II_b, Ii_i, I)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbi, iI_b, Ii_i, I)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                );

    ii_b_t = 2*(
                + bi*pa(kdb, is_d, si_b, s)*cc(dbi, dbd, dbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + b2*pa(kdb, Is_d, si_b, s)*cc(dbi, dbd, dbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + bibr*si_b
                + bibr*pa(kbb, is_b, si_b, s)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + b2b*pa(kbb, Is_b, si_b, s)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                - gi*ii_b
                - lm*ii_b
                + om*Ii_b
                + al*pa(kib, sI_i, Ii_b, I)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                + al*pa(kib, iI_i, Ii_b, I)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + al*pa(kib, rI_i, Ii_b, I)*cc(ibi, ibd, ibb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
                + nu*pa(kib, II_i, Ii_b, I)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + nu*pa(kib, iI_i, Ii_b, I)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbb, sI_b, Ii_b, I)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
                + alb*Ii_b
                + alb*pa(kbb, iI_b, Ii_b, I)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbb, rI_b, Ii_b, I)*cc(bbi, bbd, bbb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbb, II_b, Ii_b, I)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
                + nub*Ii_b
                + nub*pa(kbb, iI_b, Ii_b, I)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
                );

    ir_d_t = + bi*pa(kdd, is_d, sr_d, s)*cc(ddi, ddd, ddb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + b2*pa(kdd, Is_d, sr_d, s)*cc(ddi, ddd, ddb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kbd, is_b, sr_d, s)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kbd, Is_b, sr_d, s)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + gi*ii_d
      - gi*ir_d
      - di*ir_d
      - lm*ir_d
      - lm*ir_d
      + om*Ir_d
      + al*pa(kid, sI_i, Ir_d, I)*cc(idi, idd, idb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rI_i, Ir_d, I)*cc(idi, idd, idb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, sR_i, Ri_d, R)*cc(idi, idd, idb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kid, rR_i, Ri_d, R)*cc(idi, idd, idb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, II_i, Ir_d, I)*cc(idi, idd, idb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, IR_i, Ri_d, R)*cc(idi, idd, idb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iI_i, Ir_d, I)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kid, iR_i, Ri_d, R)*cc(idi, idd, idb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sI_b, Ir_d, I)*cc(bdi, bdd, bdb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iI_b, Ir_d, I)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rI_b, Ir_d, I)*cc(bdi, bdd, bdb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, sR_b, Ri_d, R)*cc(bdi, bdd, bdb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, iR_b, Ri_d, R)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbd, rR_b, Ri_d, R)*cc(bdi, bdd, bdb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, II_b, Ir_d, I)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, IR_b, Ri_d, R)*cc(bdi, bdd, bdb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iI_b, Ir_d, I)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbd, iR_b, Ri_d, R)*cc(bdi, bdd, bdb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    ir_i_t = + bi*pa(kdi, is_d, sr_i, s)*cc(dii, did, dib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + b2*pa(kdi, Is_d, sr_i, s)*cc(dii, did, dib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kbi, is_b, sr_i, s)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kbi, Is_b, sr_i, s)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + gi*ii_i
      - gi*ir_i
      - di*ir_i
      - lm*ir_i
      - lm*ir_i
      + om*Ir_i
      + al*pa(kii, sI_i, Ir_i, I)*cc(iii, iid, iib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*Ir_i
      + al*pa(kii, rI_i, Ir_i, I)*cc(iii, iid, iib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, sR_i, Ri_i, R)*cc(iii, iid, iib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*Ri_i
      + al*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kii, rR_i, Ri_i, R)*cc(iii, iid, iib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, II_i, Ir_i, I)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, IR_i, Ri_i, R)*cc(iii, iid, iib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kii, iI_i, Ir_i, I)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*Ri_i
      + nu*pa(kii, iR_i, Ri_i, R)*cc(iii, iid, iib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sI_b, Ir_i, I)*cc(bii, bid, bib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iI_b, Ir_i, I)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rI_b, Ir_i, I)*cc(bii, bid, bib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, sR_b, Ri_i, R)*cc(bii, bid, bib, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, iR_b, Ri_i, R)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbi, rR_b, Ri_i, R)*cc(bii, bid, bib, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, II_b, Ir_i, I)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, IR_b, Ri_i, R)*cc(bii, bid, bib, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iI_b, Ir_i, I)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbi, iR_b, Ri_i, R)*cc(bii, bid, bib, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    ir_b_t = + bi*pa(kdb, is_d, sr_b, s)*cc(dbi, dbd, dbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + b2*pa(kdb, Is_d, sr_b, s)*cc(dbi, dbd, dbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + bibr*pa(kbb, is_b, sr_b, s)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + b2b*pa(kbb, Is_b, sr_b, s)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + gi*ii_b
      - gi*ir_b
      - di*ir_b
      - lm*ir_b
      - lm*ir_b
      + om*Ir_b
      + al*pa(kib, sI_i, Ir_b, I)*cc(ibi, ibd, ibb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iI_i, Ir_b, I)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rI_i, Ir_b, I)*cc(ibi, ibd, ibb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, sR_i, Ri_b, R)*cc(ibi, ibd, ibb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, iR_i, Ri_b, R)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + al*pa(kib, rR_i, Ri_b, R)*cc(ibi, ibd, ibb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, II_i, Ir_b, I)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, IR_i, Ri_b, R)*cc(ibi, ibd, ibb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iI_i, Ir_b, I)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nu*pa(kib, iR_i, Ri_b, R)*cc(ibi, ibd, ibb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sI_b, Ir_b, I)*cc(bbi, bbd, bbb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, iI_b, Ir_b, I)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + alb*Ir_b
      + alb*pa(kbb, rI_b, Ir_b, I)*cc(bbi, bbd, bbb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, sR_b, Ri_b, R)*cc(bbi, bbd, bbb, si_i, si_d, si_b, s*i, N_Qd, N_Qi, N_Qb)
      + alb*Ri_b
      + alb*pa(kbb, iR_b, Ri_b, R)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb)
      + alb*pa(kbb, rR_b, Ri_b, R)*cc(bbi, bbd, bbb, ri_i, ri_d, ri_b, r*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, II_b, Ir_b, I)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, IR_b, Ri_b, R)*cc(bbi, bbd, bbb, Ii_i, Ii_d, Ii_b, I*i, N_Qd, N_Qi, N_Qb)
      + nub*pa(kbb, iI_b, Ir_b, I)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
      + nub*Ri_b
      + nub*pa(kbb, iR_b, Ri_b, R)*cc(bbi, bbd, bbb, ii_i, ii_d, ii_b, i*i, N_Qd, N_Qi, N_Qb);

    rr_d_t = 2*(
                + gi*ir_d
                - di*rr_d
                - lm*rr_d
                + al*pa(kid, sR_i, Rr_d, R)*cc(idi, idd, idb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
                + al*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + al*pa(kid, rR_i, Rr_d, R)*cc(idi, idd, idb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
                + nu*pa(kid, IR_i, Rr_d, R)*cc(idi, idd, idb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                + nu*pa(kid, iR_i, Rr_d, R)*cc(idi, idd, idb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, sR_b, Rr_d, R)*cc(bdi, bdd, bdb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, iR_b, Rr_d, R)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbd, rR_b, Rr_d, R)*cc(bdi, bdd, bdb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbd, IR_b, Rr_d, R)*cc(bdi, bdd, bdb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbd, iR_b, Rr_d, R)*cc(bdi, bdd, bdb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                );

    rr_i_t = 2*(
                + gi*ir_i
                - di*rr_i
                - lm*rr_i
                + al*pa(kii, sR_i, Rr_i, R)*cc(iii, iid, iib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
                + al*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + al*Rr_i
                + al*pa(kii, rR_i, Rr_i, R)*cc(iii, iid, iib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
                + nu*pa(kii, IR_i, Rr_i, R)*cc(iii, iid, iib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                + nu*pa(kii, iR_i, Rr_i, R)*cc(iii, iid, iib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, sR_b, Rr_i, R)*cc(bii, bid, bib, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, iR_b, Rr_i, R)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbi, rR_b, Rr_i, R)*cc(bii, bid, bib, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbi, IR_b, Rr_i, R)*cc(bii, bid, bib, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbi, iR_b, Rr_i, R)*cc(bii, bid, bib, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                );

    rr_b_t = 2*(
                + gi*ir_b
                - di*rr_b
                - lm*rr_b
                + al*pa(kib, sR_i, Rr_b, R)*cc(ibi, ibd, ibb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
                + al*pa(kib, iR_i, Rr_b, R)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + al*pa(kib, rR_i, Rr_b, R)*cc(ibi, ibd, ibb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
                + nu*pa(kib, IR_i, Rr_b, R)*cc(ibi, ibd, ibb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                + nu*pa(kib, iR_i, Rr_b, R)*cc(ibi, ibd, ibb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbb, sR_b, Rr_b, R)*cc(bbi, bbd, bbb, sr_i, sr_d, sr_b, s*r, N_Qd, N_Qi, N_Qb)
                + alb*pa(kbb, iR_b, Rr_b, R)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                + alb*Rr_b
                + alb*pa(kbb, rR_b, Rr_b, R)*cc(bbi, bbd, bbb, rr_i, rr_d, rr_b, r*r, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbb, IR_b, Rr_b, R)*cc(bbi, bbd, bbb, Ir_i, Ir_d, Ir_b, I*r, N_Qd, N_Qi, N_Qb)
                + nub*pa(kbb, iR_b, Rr_b, R)*cc(bbi, bbd, bbb, ir_i, ir_d, ir_b, i*r, N_Qd, N_Qi, N_Qb)
                );

    return GSL_SUCCESS;         
  }
      
}; // FullModelPairApprox2

//------------------------------------------------------------
         
#endif /* FULL_MODEL_H */
