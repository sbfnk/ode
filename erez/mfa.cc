
/******************************************************************/

#include "mfa.hh"

using namespace std;

/******************************************************************/

MeanField::MeanField() : Ode()
{
   //OdePluginFuncs(&derivs,&jac);
   OdePluginFuncs(&derivs,NULL);
}

/******************************************************************/

MeanField::~MeanField()
{
   delete static_cast<ModelParams *>(GetModelParams());   
}

/******************************************************************/
//
// S = y[0]
// I = y[1]
// R = y[2]
// s = y[3]
// i = y[4]
// r = y[5]
//
int MeanField::derivs (double t, const double y[], double rhs[], void *params)
{   
   ModelParams p = *(ModelParams *)params;
   double bd=p.beta_d, gd=p.gamma_d, dd=p.delta_d;
   double bi=p.beta_i, gi=p.gamma_i, di=p.delta_i;
   double bm1=p.beta_m1, bm2=p.beta_m2, al=p.alpha, nu=p.nu, lm=p.lambda;
   double Qd=p.Qd, Qi=p.Qi, N=p.N;
   double Qd_N=Qd/N, Qi_N=Qi/N;
   
   rhs[0] = -(Qd_N) * ( bd*y[1] + bm1*y[4] ) * y[0] + dd * y[2]
      - al * (Qi_N) * ( y[3] + y[4] + y[5] ) * y[0] + lm * y[3]
      - nu * (Qi_N) * ( y[1] + y[4] ) * y[0];
   
   rhs[1] =  (Qd_N) * ( bd*y[1] + bm1*y[4] ) * y[0] - gd * y[1]
      - al * (Qi_N) * ( y[3] + y[4] + y[5] ) * y[1] + lm * y[4]
      - nu * (Qi_N) * ( y[1] + y[4] ) * y[1];
   
   rhs[2] =  gd * y[1] - dd * y[2] + lm * y[5]
      - al * (Qi_N) * ( y[3] + y[4] + y[5] ) * y[2]
      - nu * (Qi_N) * ( y[1] + y[4] ) * y[2];
   
   rhs[3] = -(Qd_N) * ( bi*y[4] + bm2*y[1] ) * y[3] + di * y[5]
      + al * (Qi_N) * ( y[3] + y[4] + y[5] ) * y[0] - lm * y[3]
      + nu * (Qi_N) * ( y[1] + y[4] ) * y[0];
   
   rhs[4] =  (Qd_N) * ( bi*y[4] + bm2*y[1] ) * y[3] - gi * y[4]
      + al * (Qi_N) * ( y[3] + y[4] + y[5] ) * y[1] - lm * y[4]
      + nu * (Qi_N) * ( y[1] + y[4] ) * y[1];
   
   rhs[5] =  gi * y[4] - di * y[5] - lm * y[5]
      + al * (Qi_N) * ( y[3] + y[4] + y[5] ) * y[2]
      + nu * (Qi_N) * ( y[1] + y[4] ) * y[2];
   
   return GSL_SUCCESS;
}

/******************************************************************/
/*
  int MeanField::jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
  {
  ModelParams prm = *(static_cast<ModelParams *>(params));
  double mu = prm.mu1+prm.mu2;
  int n=prm.njac;
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, n, n);
  gsl_matrix *m = &dfdy_mat.matrix; 
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, 1.0);
  gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
  gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  
  return GSL_SUCCESS;
  }
*/
/******************************************************************/ 

void MeanField::PrtModelPrms()
{
   PrtOdePrms();

   ModelParams *p=static_cast<ModelParams *>(GetModelParams());
   
   cout << "Model parameters:" << endl
        << "-----------------" << endl
        << "beta_d  = " << p->beta_d << endl
        << "gamma_d = " << p->gamma_d << endl
        << "delta_d = " << p->delta_d << endl
        << "beta_i  = " << p->beta_i << endl
        << "gamma_i = " << p->gamma_i << endl
        << "delta_i = " << p->delta_i << endl
        << "beta_m1 = " << p->beta_m1 << endl
        << "beta_m2 = " << p->beta_m2 << endl
        << "alpha   = " << p->alpha << endl
        << "nu      = " << p->nu << endl
        << "lambda  = " << p->lambda << endl
        << "Qd      = " << p->Qd << endl
        << "Qi      = " << p->Qi << endl
        << "N       = " << p->N << endl
        << endl;
}

/******************************************************************/


