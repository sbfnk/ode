/********************************************************************/

#include "van_der_pol.h"

using namespace std;

/********************************************************************/

VanDerPol::VanDerPol() : Ode()
{
   model_params = NULL;
   ode_plugin_funcs(&derivs,&jac);
}

/********************************************************************/

VanDerPol::~VanDerPol()
{
   delete model_params;
}

/********************************************************************/

void VanDerPol::set_model_params(ModelParams *model_params)
{
   VanDerPol::model_params=model_params;
}   

/********************************************************************/

void VanDerPol::plugin_model_params()
{
   set_params(static_cast<void *>(model_params));
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
   ModelParams prm = *(ModelParams *)params;
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
   
void VanDerPol::prt_model_prms()
{
   prt_ode_prms();

   cout << "Model parameters:" << endl
        << "-----------------" << endl
        << "beta_d  = " << get_model_params()->beta_d << endl
        << "gamma_d = " << get_model_params()->gamma_d << endl
        << "delta_d = " << get_model_params()->delta_d << endl
        << "beta_i  = " << get_model_params()->beta_i << endl
        << "gamma_i = " << get_model_params()->gamma_i << endl
        << "delta_i = " << get_model_params()->delta_i << endl
        << "beta_m  = " << get_model_params()->beta_m << endl
        << "alpha   = " << get_model_params()->alpha << endl
        << "nu      = " << get_model_params()->nu << endl
        << "lambda  = " << get_model_params()->lambda << endl
        << endl;
}

/********************************************************************/  


