/********************************************************************/

#include "van_der_pol.h"

using namespace std;

/********************************************************************/

VanDerPol::VanDerPol() : Ode()
{
   OdePluginFuncs(&derivs,&jac);
}

/********************************************************************/

VanDerPol::~VanDerPol()
{
   delete static_cast<ModelParams *>(GetModelParams());   
}

/********************************************************************/

int VanDerPol::derivs (double t, const double y[], double rhs[], void *params)
{   
   ModelParams prm = *(ModelParams *)params;
   double mu = prm.mu1+prm.mu2;
   
   rhs[0] = y[1];
   rhs[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
   
   return GSL_SUCCESS;
}

/********************************************************************/

int VanDerPol::jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
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

/********************************************************************/   
   
void VanDerPol::PrtModelPrms()
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
        << "beta_m  = " << p->beta_m << endl
        << "alpha   = " << p->alpha << endl
        << "nu      = " << p->nu << endl
        << "lambda  = " << p->lambda << endl
        << endl;
}

/********************************************************************/  


