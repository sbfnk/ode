
/******************************************************************/

#include "model_ode.hh"
#include "pa_macros.hh"
#include <iomanip>

using namespace std;

/******************************************************************/

ModelOde::ModelOde() : Ode() {}

/******************************************************************/

ModelOde::~ModelOde()
{
   delete static_cast<ModelParams *>(GetModelParams());   
}

/******************************************************************/

void ModelOde::PrtModelPrms() const
{
   ModelParams *p=static_cast<ModelParams *>(GetModelParams());
   
   cout << "Model parameters:" << endl
        << "-----------------" << endl
        << "beta--  = " << p->beta[0][0] << endl
        << "beta+-  = " << p->beta[1][0] << endl
        << "beta-+  = " << p->beta[0][1] << endl
        << "beta++  = " << p->beta[1][1] << endl
        << "gamma-  = " << p->gamma[0] << endl
        << "gamma+  = " << p->gamma[1] << endl
        << "delta-  = " << p->delta[0] << endl
        << "delta+  = " << p->delta[1] << endl
        << "alpha   = " << p->alpha << endl
        << "nu      = " << p->nu << endl
        << "lambda+ = " << p->lambda << endl
        << "omega-  = " << p->omega << endl
        << "Qd      = " << p->Qd << endl
        << "Qi      = " << p->Qi << endl
        << "N       = " << p->N << endl
        << "R_0 d   = " << (p->beta[0][0])*(p->Qd)/(p->gamma[0]) << endl
        << "R_0 i   = " << (p->alpha)*(p->Qi)/(p->lambda) << endl
        << endl;
}

/******************************************************************/ 

void ModelOde::PrtGraphPrms() const
{
   ModelParams *p=static_cast<ModelParams *>(GetModelParams());
   
   cout << "Graph parameters:" << endl
        << "-----------------" << endl
        << "C_ddi  = " << p->Cddi << endl
        << "C_ddd  = " << p->Cddd << endl
        << "C_dii  = " << p->Cdii << endl
        << "C_did  = " << p->Cdid << endl
        << "C_idi  = " << p->Cidi << endl
        << "C_idd  = " << p->Cidd << endl
        << "C_iii  = " << p->Ciii << endl
        << "C_iid  = " << p->Ciid << endl
        << endl;
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
int ModelOde::MFderivs (double t, const double y[], double rhs[], void *params)
{   
   ModelParams p = *(ModelParams *)params;
   double bd=p.beta[0][0], gd=p.gamma[0], dd=p.delta[0];
   double bi=p.beta[1][1], gi=p.gamma[1], di=p.delta[1];
   double bm1=p.beta[0][1], bm2=p.beta[1][0], al=p.alpha, nu=p.nu, lm=p.lambda;
   double Qd=p.Qd, Qi=p.Qi, N=p.N;
   //double Qd_N=Qd/N, Qi_N=Qi/N;
   double Qd_N=1./(double)N, Qi_N=1./double(N);
   
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
  int ModelOde::MFjac (double t, const double y[], double *dfdy, double dfdt[], void *params)
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

// Pair Approximation
double  ModelOde::PA(double kk, double ij, double jk, double j)
{
   //const double eps = 1e-12;
   
   //if(j<eps)
   if(j == 0.0)
      return 0.0;
   else
      return kk*((ij*jk)/j);
}

/******************************************************************/

// Clustering (Coefficient) Correction
double  ModelOde::CC(double C_i, double C_d, double ki_i, double ki_d,
                      double ki, double N_Qd, double N_Qi)
{
   const double eps = 1e-12;

   if(ki<eps)
      return ((1.0 - C_i - C_d));
   else
      return ((1.0 - C_i - C_d) + C_i*N_Qi*ki_i/ki + C_d*N_Qd*ki_d/ki);   
}

/******************************************************************/
//
// S = 0 , SS_d = 6  , ss_d = 12 , Ss_d = 18 , 
// I = 1 , SI_d = 7  , si_d = 13 , Si_d = 19 , sI_d = 24 ,
// R = 2 , SR_d = 8  , sr_d = 14 , Sr_d = 20 , sR_d = 25 ,
// s = 3 , II_d = 9  , ii_d = 15 , Ii_d = 21 ,
// i = 4 , IR_d = 10 , ir_d = 16 , Ir_d = 22 , iR_d = 26 ,
// r = 5 , RR_d = 11 , rr_d = 17 , Rr_d = 23 ,
//
//         SS_i = 27 , ss_i = 33 , Ss_i = 39 , 
//         SI_i = 28 , si_i = 34 , Si_i = 40 , sI_i = 45 ,
//         SR_i = 29 , sr_i = 35 , Sr_i = 41 , sR_i = 46 ,
//         II_i = 30 , ii_i = 36 , Ii_i = 42 ,
//         IR_i = 31 , ir_i = 37 , Ir_i = 43 , iR_i = 47 .
//         RR_i = 32 , rr_i = 38 , Rr_i = 44 ,
//
int ModelOde::PAderivs (double t, const double y[], double rhs[], void *params)
{
   ModelParams p = *(static_cast<ModelParams *>(params));   
   
   // local readable short variables
   double Qd=p.Qd, Qi=p.Qi, N=p.N;
   double bd=p.beta[0][0], gd=p.gamma[0], dd=p.delta[0];
   double bi=p.beta[1][1], gi=p.gamma[1], di=p.delta[1];
   double b1=p.beta[0][1], b2=p.beta[1][0];
   double al=p.alpha, nu=p.nu, lm=p.lambda;
   double N_Qd=N/Qd, N_Qi=N/Qi;

   //bd/=Qd; bi/=Qd; b1/=Qd; b2/=Qd;
   //al/=Qd; nu/=Qd;
     
   // Clustering corrections
   double idi=p.Cidi, idd=p.Cidd;
   double dii=p.Cdii, did=p.Cdid;
   double ddi=p.Cddi, ddd=p.Cddd;
   double iii=p.Ciii, iid=p.Ciid;
   double kdd=(Qd-1.0)/Qd, kii=(Qi-1.0)/Qi;
   double kdi=1.0, kid=1.0;

   // CHECK si_d, si_i info. loss TERMS FOR FACTOR 2 !!!!!!!!!!!!!!!!!!!!!


   
   // read to InitLocalParams() from the outside  - Just Once !!!

   //--------------------------------------------------------------
   // singlet equations
   //--------------------------------------------------------------
   
   S_t =
      - bd*SI_d - b1*Si_d + dd*R + lm*s
      - al*(Ss_i + Si_i + Sr_i) - nu*(SI_i + Si_i);
   
   I_t =
      + bd*SI_d + b1*Si_d - gd*I + lm*i
      - al*(Is_i + Ii_i + Ir_i) - nu*(II_i + Ii_i);
   
   R_t =
      + gd*I - dd*R + lm*r
      - al*(Rs_i + Ri_i +Rr_i) - nu*(RI_i + Ri_i);
   
   s_t =
      - bi*si_d - b2*sI_d + di*r - lm*s
      + al*(Ss_i + Si_i + Sr_i) + nu*(SI_i + Si_i);

   i_t =
      + bi*si_d + b2*sI_d - gi*i - lm*i
      + al*(Is_i + Ii_i + Ir_i) + nu*(II_i + Ii_i);
      
   r_t =
      + gi*i - di*r - lm*r
      + al*(Rs_i + Ri_i +Rr_i) + nu*(RI_i + Ri_i);

   // pair equations

   //--------------------------------------------------------------
   // un-informed d-edges dynamics
   //--------------------------------------------------------------
   
//----------------
   SS_d_t =

      2.0*(
      
      // infection
      - bd*PA(kdd, SS_d, SI_d,S)*CC(ddi, ddd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - b1*PA(kdd, SS_d, Si_d,S)*CC(ddi, ddd, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // susceptibility
      + dd*SR_d
      
      // info. trans.
      - al*PA(kdi, SS_d, Ss_i,S)*CC(dii, did, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      - al*PA(kdi, SS_d, Si_i,S)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      - al*PA(kdi, SS_d, Sr_i,S)*CC(dii, did, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kdi, SS_d, SI_i,S)*CC(dii, did, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - nu*PA(kdi, SS_d, Si_i,S)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*sS_d   );

//----------------
   SI_d_t =
      
      // infection
      - bd*SI_d

      + bd*PA(kdd, SS_d, SI_d,S)*CC(ddi, ddd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + b1*PA(kdd, SS_d, Si_d,S)*CC(ddi, ddd, Si_i, Si_d,S*i, N_Qd, N_Qi)

      - bd*PA(kdd, IS_d, SI_d,S)*CC(ddi, ddd, II_i, II_d,I*I, N_Qd, N_Qi)
      - b1*PA(kdd, iS_d, SI_d,S)*CC(ddi, ddd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      // recovery
      - gd*SI_d

      // susceptibility
      + dd*RI_d
      
      // info. trans.
      - al*PA(kid, sS_i, SI_d,S)*CC(idi, idd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - al*PA(kid, iS_i, SI_d,S)*CC(idi, idd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      - al*PA(kid, rS_i, SI_d,S)*CC(idi, idd, rI_i, rI_d,r*I, N_Qd, N_Qi)
      
      - al*PA(kdi, SI_d, Is_i,I)*CC(dii, did, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      - al*PA(kdi, SI_d, Ii_i,I)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      - al*PA(kdi, SI_d, Ir_i,I)*CC(dii, did, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kid, IS_i, SI_d,S)*CC(idi, idd, II_i, II_d,I*I, N_Qd, N_Qi)
      - nu*PA(kid, iS_i, SI_d,S)*CC(idi, idd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      - nu*PA(kdi, SI_d, II_i,I)*CC(dii, did, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - nu*PA(kdi, SI_d, Ii_i,I)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*sI_d + lm*Si_d;
   
//----------------   
   SR_d_t =

      // infection
      - bd*PA(kdd, IS_d, SR_d,S)*CC(ddi, ddd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - b1*PA(kdd, iS_d, SR_d,S)*CC(ddi, ddd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // recovery
      + gd*SI_d

      // susceptibility
      - dd*SR_d + dd*RR_d
      
      // info. trans.
      - al*PA(kid, sS_i, SR_d,S)*CC(idi, idd, sR_i, sR_d,s*R, N_Qd, N_Qi)
      - al*PA(kid, iS_i, SR_d,S)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      - al*PA(kid, rS_i, SR_d,S)*CC(idi ,idd, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kdi, SR_d, Rs_i,R)*CC(dii, did, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      - al*PA(kdi, SR_d, Ri_i,R)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      - al*PA(kdi, SR_d, Rr_i,R)*CC(dii, did, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kid, IS_i, SR_d,S)*CC(idi, idd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - nu*PA(kid, iS_i, SR_d,S)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      - nu*PA(kdi, SR_d, RI_i,R)*CC(dii, did, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - nu*PA(kdi, SR_d, Ri_i,R)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*sR_d + lm*Sr_d;      

//----------------
   II_d_t =

      2.0*(
      
      // infection
      + bd*SI_d
      
      + bd*PA(kdd, IS_d, SI_d,S)*CC(ddi, ddd, II_i, II_d,I*I, N_Qd, N_Qi)
      + b1*PA(kdd, IS_d, Si_d,S)*CC(ddi, ddd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // recovery
      - gd*II_d
      
      // info. trans.
      - al*PA(kdi, II_d, Is_i,I)*CC(dii, did, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - al*PA(kdi, II_d, Ii_i,I)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - al*PA(kdi, II_d, Ir_i,I)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kdi, II_d, II_i,I)*CC(dii, did, II_i, II_d,I*I, N_Qd, N_Qi)
      - nu*PA(kdi, II_d, Ii_i,I)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*Ii_d );       

//----------------
   IR_d_t =
      
      // infection
      + bd*PA(kdd, IS_d, SR_d,S)*CC(ddi, ddd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + b1*PA(kdd, iS_d, SR_d,S)*CC(ddi, ddd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // recovery
      + gd*II_d - gd*IR_d
      
      // susceptibility
      - dd*IR_d
      
      // info. trans.
      - al*PA(kid, sI_i, IR_d,I)*CC(idi, idd, sR_i, sR_d,s*R, N_Qd, N_Qi)
      - al*PA(kid, iI_i, IR_d,I)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      - al*PA(kid, rI_i, IR_d,I)*CC(idi, idd, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kdi, IR_d, Rs_i,R)*CC(dii, did, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - al*PA(kdi, IR_d, Ri_i,R)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - al*PA(kdi, IR_d, Rr_i,R)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kid, II_i, IR_d,I)*CC(idi, idd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - nu*PA(kid, iI_i, IR_d,I)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      - nu*PA(kdi, IR_d, RI_i,R)*CC(dii, did, II_i, II_d,I*I, N_Qd, N_Qi)
      - nu*PA(kdi, IR_d, Ri_i,R)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*iR_d + lm*Ir_d;        

//----------------
   RR_d_t =

      2.0*(
      
      // recovery
      + gd*IR_d
      
      // susceptibility
      - dd*RR_d
      
      // info. trans.
      - al*PA(kdi, RR_d, Rs_i,R)*CC(dii, did, Rs_i, Rs_d,R*s, N_Qd, N_Qi)
      - al*PA(kdi, RR_d, Ri_i,R)*CC(dii, did, Ri_i, Ri_d,R*i, N_Qd, N_Qi)
      - al*PA(kdi, RR_d, Rr_i,R)*CC(dii, did, Rr_i, Rr_d,R*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kdi, RR_d, RI_i,R)*CC(dii, did, RI_i, RI_d,R*I, N_Qd, N_Qi)
      - nu*PA(kdi, RR_d, Ri_i,R)*CC(dii, did, Ri_i, Ri_d,R*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*Rr_d );
   
   //--------------------------------------------------------------
   // informed d-edges dynamics
   //--------------------------------------------------------------
   
//----------------
   ss_d_t =

      2.0*(
      
      // infection
      - bi*PA(kdd, ss_d, si_d,s)*CC(ddi, ddd, si_i, si_d,s*i, N_Qd, N_Qi)
      - b2*PA(kdd, ss_d, sI_d,s)*CC(ddi, ddd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      
      // susceptibility
      + di*sr_d
      
      // info. trans.
      + al*PA(kdi, sS_d, Ss_i,S)*CC(dii, did, ss_i, ss_d,s*s, N_Qd, N_Qi)
      + al*PA(kdi, sS_d, Si_i,S)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kdi, sS_d, Sr_i,S)*CC(dii, did, sr_i, sr_d,s*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kdi, sS_d, SI_i,S)*CC(dii, did, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + nu*PA(kdi, sS_d, Si_i,S)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      
      // info. loss
      - lm*ss_d );   

//----------------
   si_d_t =
      
      // infection
      - bi*si_d
      
      + bi*PA(kdd, ss_d, si_d,s)*CC(ddi, ddd, si_i, si_d,s*i, N_Qd, N_Qi)
      + b2*PA(kdd, ss_d, sI_d,s)*CC(ddi, ddd, sI_i, sI_d,s*I, N_Qd, N_Qi)

      - bi*PA(kdd, is_d, si_d,s)*CC(ddi, ddd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - b2*PA(kdd, Is_d, si_d,s)*CC(ddi, ddd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // recovery
      - gi*si_d

      // susceptibility
      + di*ri_d
      
      // info. trans.
      + al*PA(kdi, sI_d, Is_i,I)*CC(dii, did, ss_i, ss_d,s*s, N_Qd, N_Qi)
      + al*PA(kdi, sI_d, Ii_i,I)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kdi, sI_d, Ir_i,I)*CC(dii, did, sr_i, sr_d,s*r, N_Qd, N_Qi)
      
      + al*PA(kid, sS_i, Si_d,S)*CC(idi, idd, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kid, iS_i, Si_d,S)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + al*PA(kid, rS_i, Si_d,S)*CC(idi, idd, ri_i, ri_d,r*i, N_Qd, N_Qi)
            
      // info. source
      + nu*PA(kid, iS_i, Si_d,S)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + nu*PA(kid, IS_i, Si_d,S)*CC(idi, idd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      + nu*PA(kdi, sI_d, II_i,I)*CC(dii, did, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + nu*PA(kdi, sI_d, Ii_i,I)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      
      // info. loss
      - 2.0 * lm*si_d;
   //- lm*si_d;

//----------------
   sr_d_t = 
      
      // infection
      - bi*PA(kdd, is_d, sr_d,s)*CC(ddi, ddd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      - b2*PA(kdd, Is_d, sr_d,s)*CC(ddi, ddd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // recovery
      + gi*si_d

      // susceptibility
      - di*sr_d + di*rr_d
      
      // info. trans.
      + al*PA(kdi, sR_d, Rs_i,R)*CC(dii, did, ss_i, ss_d,s*s, N_Qd, N_Qi)
      + al*PA(kdi, sR_d, Ri_i,R)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kdi, sR_d, Rr_i,R)*CC(dii, did, sr_i, sr_d,s*r, N_Qd, N_Qi)

      + al*PA(kid, sS_i, Sr_d,S)*CC(idi, idd, sr_i, sr_d,s*r, N_Qd, N_Qi)
      + al*PA(kid, iS_i, Sr_d,S)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + al*PA(kid, rS_i, Sr_d,S)*CC(idi, idd, rr_i, rr_d,r*r, N_Qd, N_Qi)      
      
      // info. source
      + nu*PA(kid, iS_i, Sr_d,S)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + nu*PA(kid, IS_i, Sr_d,S)*CC(idi, idd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      + nu*PA(kdi, sR_d, RI_i,R)*CC(dii, did, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + nu*PA(kdi, sR_d, Ri_i,R)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      
      // info. loss
      - 2.0 * lm*sr_d;
   //- lm*sr_d; 
   
//----------------
   ii_d_t =

      2.0*(
      
      // infection
      + bi*si_d
      
      + bi*PA(kdd, is_d, si_d,s)*CC(ddi, ddd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + b2*PA(kdd, is_d, sI_d,s)*CC(ddi, ddd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      // recovery
      - gi*ii_d
      
      // info. trans.
      + al*PA(kdi, iI_d, Is_i,I)*CC(dii, did, is_i, is_d,i*s, N_Qd, N_Qi)
      + al*PA(kdi, iI_d, Ii_i,I)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + al*PA(kdi, iI_d, Ir_i,I)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kdi, iI_d, II_i,I)*CC(dii, did, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + nu*PA(kdi, iI_d, Ii_i,I)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // info. loss
      - lm*ii_d );       

//----------------
   ir_d_t =  

      // infection
      + bi*PA(kdd, is_d, sr_d,s)*CC(ddi, ddd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + b2*PA(kdd, Is_d, sr_d,s)*CC(ddi, ddd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // recovery
      + gi*ii_d - gi*ir_d
      
      // susceptibility
      - di*ir_d
      
      // info. trans.
      + al*PA(kdi, iR_d, Rs_i,R)*CC(dii, did, is_i, is_d,i*s, N_Qd, N_Qi)
      + al*PA(kdi, iR_d, Ri_i,R)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + al*PA(kdi, iR_d, Rr_i,R)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      + al*PA(kid, sI_i, Ir_d,I)*CC(idi, idd, sr_i, sr_d,s*r, N_Qd, N_Qi)
      + al*PA(kid, iI_i, Ir_d,I)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + al*PA(kid, rI_i, Ir_d,I)*CC(idi, idd, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kid, II_i, Ir_d,I)*CC(idi, idd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      + nu*PA(kid, iI_i, Ir_d,I)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      + nu*PA(kdi, iR_d, RI_i,R)*CC(dii, did, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + nu*PA(kdi, iR_d, Ri_i,R)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // info. loss
      - 2.0 * lm*ir_d;
      //- lm*ir_d;

//----------------
   rr_d_t =

      2.0*(
      
      // recovery
      + gi*ir_d
      
      // susceptibility
      - di*rr_d
      
      // info. trans.
      + al*PA(kdi, rR_d, Rs_i,R)*CC(dii, did, rs_i, rs_d,r*s, N_Qd, N_Qi)
      + al*PA(kdi, rR_d, Ri_i,R)*CC(dii, did, ri_i, ri_d,r*i, N_Qd, N_Qi)
      + al*PA(kdi, rR_d, Rr_i,R)*CC(dii, did, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kdi, rR_d, RI_i,R)*CC(dii, did, rI_i, rI_d,r*I, N_Qd, N_Qi)
      + nu*PA(kdi, rR_d, Ri_i,R)*CC(dii, did, ri_i, ri_d,r*i, N_Qd, N_Qi)
      
      // info. loss
      - lm*rr_d );
   
   //--------------------------------------------------------------
   // mix d-edges dynamics      
   //--------------------------------------------------------------
   
//----------------
   Ss_d_t =
      
      // infection
      - b2*PA(kdd, Ss_d, sI_d,s)*CC(ddi, ddd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - bi*PA(kdd, Ss_d, si_d,s)*CC(ddi, ddd, Si_i, Si_d,S*i, N_Qd, N_Qi)

      - bd*PA(kdd, IS_d, Ss_d,S)*CC(ddi, ddd, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - b1*PA(kdd, iS_d, Ss_d,S)*CC(ddi, ddd, is_i, is_d,i*s, N_Qd, N_Qi)
      
      // susceptibility
      + di*Sr_d + dd*sR_d
      
      // info. trans.
      + al*PA(kdi, SS_d, Ss_i,S)*CC(dii, did, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      + al*PA(kdi, SS_d, Si_i,S)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      + al*PA(kdi, SS_d, Sr_i,S)*CC(dii, did, Sr_i, Sr_d,S*r, N_Qd, N_Qi)

      - al*PA(kid, sS_i, Ss_d,S)*CC(idi, idd, ss_i, ss_d,s*s, N_Qd, N_Qi)
      - al*PA(kid, iS_i, Ss_d,S)*CC(idi, idd, is_i, is_d,i*s, N_Qd, N_Qi)
      - al*PA(kid, rS_i, Ss_d,S)*CC(idi, idd, rs_i, rs_d,r*s, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kdi, SS_d, SI_i,S)*CC(dii, did, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + nu*PA(kdi, SS_d, Si_i,S)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)

      - nu*PA(kid, IS_i, Ss_d,S)*CC(idi, idd, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - nu*PA(kid, iS_i, Ss_d,S)*CC(idi, idd, is_i, is_d,i*s, N_Qd, N_Qi)
      
      // info. loss
      - lm*Ss_d + lm*ss_d;

//----------------
   Si_d_t =

      // infection
      - b1*Si_d
      
      + b2*PA(kdd, Ss_d, sI_d,s)*CC(ddi, ddd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + bi*PA(kdd, Ss_d, si_d,s)*CC(ddi, ddd, Si_i, Si_d,S*i, N_Qd, N_Qi)      

      - bd*PA(kdd, IS_d, Si_d,S)*CC(ddi, ddd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - b1*PA(kdd, iS_d, Si_d,S)*CC(ddi, ddd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // recovery
      - gi*Si_d
      
      // susceptibility
      + dd*Ri_d
      
      // info. trans.
      + al*PA(kdi, SI_d, Is_i,I)*CC(dii, did, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      + al*PA(kdi, SI_d, Ii_i,I)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      + al*PA(kdi, SI_d, Ir_i,I)*CC(dii, did, Sr_i, Sr_d,S*r, N_Qd, N_Qi)

      - al*PA(kid, sS_i, Si_d,S)*CC(idi, idd, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kid, iS_i, Si_d,S)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - al*PA(kid, rS_i, Si_d,S)*CC(idi, idd, ri_i, ri_d,r*i, N_Qd, N_Qi)      

      // info. source
      - nu*PA(kid, IS_i, Si_d,S)*CC(idi, idd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - nu*PA(kid, iS_i, Si_d,S)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      + nu*PA(kdi, SI_d, II_i,I)*CC(dii, did, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + nu*PA(kdi, SI_d, Ii_i,I)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*si_d - lm*Si_d;   

//----------------   
   Sr_d_t =
      
      // infection
      - bd*PA(kdd, IS_d, Sr_d,S)*CC(ddi, ddd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      - b1*PA(kdd, iS_d, Sr_d,S)*CC(ddi, ddd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // recovery
      + gi*Si_d

      // susceptibility
      - di*Sr_d + dd*Rr_d
      
      // info. trans.
      + al*PA(kdi, SR_d, Rs_i,R)*CC(dii, did, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      + al*PA(kdi, SR_d, Ri_i,R)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      + al*PA(kdi, SR_d, Rr_i,R)*CC(dii, did, Sr_i, Sr_d,S*r, N_Qd, N_Qi)

      - al*PA(kid, sS_i, Sr_d,S)*CC(idi, idd, sr_i, sr_d,s*r, N_Qd, N_Qi)
      - al*PA(kid, iS_i, Sr_d,S)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      - al*PA(kid, rS_i, Sr_d,S)*CC(idi ,idd, rr_i, rr_d,r*r, N_Qd, N_Qi)      
      
      // info. source
      + nu*PA(kdi, SR_d, RI_i,R)*CC(dii, did, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + nu*PA(kdi, SR_d, Ri_i,R)*CC(dii, did, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      - nu*PA(kid, IS_i, Sr_d,S)*CC(idi, idd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      - nu*PA(kid, iS_i, Sr_d,S)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. loss
      + lm*sr_d - lm*Sr_d;      

//----------------      
   sI_d_t =

      // infection
      - b2*sI_d

      + bd*PA(kdd, sS_d, SI_d,S)*CC(ddi, ddd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + b1*PA(kdd, sS_d, Si_d,S)*CC(ddi, ddd, si_i, si_d,s*i, N_Qd, N_Qi)

      - b2*PA(kdd, Is_d, sI_d,s)*CC(ddi, ddd, II_i, II_d,I*I, N_Qd, N_Qi)
      - bi*PA(kdd, is_d, sI_d,s)*CC(ddi, ddd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      // recovery
      - gd*sI_d
      
      // susceptibility
      + di*rI_d
      
      // info. trans.
      - al*PA(kdi, sI_d, Is_i,I)*CC(dii, did, ss_i, ss_d,s*s, N_Qd, N_Qi)
      - al*PA(kdi, sI_d, Ii_i,I)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kdi, sI_d, Ir_i,I)*CC(dii, did, sr_i, sr_d,s*r, N_Qd, N_Qi)

      + al*PA(kid, sS_i, SI_d,S)*CC(idi, idd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + al*PA(kid, iS_i, SI_d,S)*CC(idi, idd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + al*PA(kid, rS_i, SI_d,S)*CC(idi, idd, rI_i, rI_d,r*I, N_Qd, N_Qi)      
      
      // info. source
      - nu*PA(kdi, sI_d, II_i,I)*CC(dii, did, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - nu*PA(kdi, sI_d, Ii_i,I)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      
      + nu*PA(kid, iS_i, SI_d,S)*CC(idi, idd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + nu*PA(kid, IS_i, SI_d,S)*CC(idi, idd, II_i, II_d,I*I, N_Qd, N_Qi)
      
      // info. loss
      - lm*sI_d + lm*si_d;   

//----------------      
   Ii_d_t =
      
      // infection
      + b2*sI_d
      + b1*Si_d

      + b2*PA(kdd, Is_d, sI_d,s)*CC(ddi, ddd, II_i, II_d,I*I, N_Qd, N_Qi)
      + bi*PA(kdd, Is_d, si_d,s)*CC(ddi, ddd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      + bd*PA(kdd, IS_d, Si_d,S)*CC(ddi, ddd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      + b1*PA(kdd, iS_d, Si_d,S)*CC(ddi, ddd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // recovery
      - gd*Ii_d - gi*Ii_d
      
      // info. trans.
      + al*PA(kdi, II_d, Is_i,I)*CC(dii, did, Is_i, Is_d,I*s, N_Qd, N_Qi)
      + al*PA(kdi, II_d, Ii_i,I)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      + al*PA(kdi, II_d, Ir_i,I)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)

      - al*PA(kid, sI_i, Ii_d,I)*CC(idi, idd, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kid, iI_i, Ii_d,I)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - al*PA(kid, rI_i, Ii_d,I)*CC(idi, idd, ri_i, ri_d,r*i, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kdi, II_d, II_i,I)*CC(dii, did, II_i, II_d,I*I, N_Qd, N_Qi)
      + nu*PA(kdi, II_d, Ii_i,I)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      - nu*PA(kid, II_i, Ii_d,I)*CC(idi, idd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - nu*PA(kid, iI_i, Ii_d,I)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*ii_d - lm*Ii_d;       

//----------------      
   Ir_d_t =

      // infection
      + bd*PA(kdd, IS_d, Sr_d,S)*CC(ddi, ddd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      + b1*PA(kdd, iS_d, Sr_d,S)*CC(ddi, ddd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // recovery
      + gi*Ii_d - gd*Ir_d
      
      // susceptibility
      - di*Ir_d
      
      // info. trans.
      + al*PA(kdi, IR_d, Rs_i,R)*CC(dii, did, Is_i, Is_d,I*s, N_Qd, N_Qi)
      + al*PA(kdi, IR_d, Ri_i,R)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      + al*PA(kdi, IR_d, Rr_i,R)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      - al*PA(kid, sI_i, Ir_d,I)*CC(idi, idd, sr_i, sr_d,s*r, N_Qd, N_Qi)
      - al*PA(kid, iI_i, Ir_d,I)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      - al*PA(kid, rI_i, Ir_d,I)*CC(idi, idd, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kdi, IR_d, RI_i,R)*CC(dii, did, II_i, II_d,I*I, N_Qd, N_Qi)
      + nu*PA(kdi, IR_d, Ri_i,R)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)

      - nu*PA(kid, II_i, Ir_d,I)*CC(idi, idd, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      - nu*PA(kid, iI_i, Ir_d,I)*CC(idi, idd, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. loss
      + lm*ir_d - lm*Ir_d;        

//----------------
   sR_d_t =

      // infection
      - b2*PA(kdd, Is_d, sR_d,s)*CC(ddi, ddd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - bi*PA(kdd, is_d, sR_d,s)*CC(ddi, ddd, iR_i, iR_d,i*R, N_Qd, N_Qi)      
      
      // recovery
      + gd*sI_d

      // susceptibility
      - dd*sR_d + di*rR_d
      
      // info. trans.
      + al*PA(kid, sS_i, SR_d,S)*CC(idi, idd, sR_i, sR_d,s*R, N_Qd, N_Qi)
      + al*PA(kid, iS_i, SR_d,S)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      + al*PA(kid, rS_i, SR_d,S)*CC(idi ,idd, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kdi, sR_d, Rs_i,R)*CC(dii, did, ss_i, ss_d,s*s, N_Qd, N_Qi)
      - al*PA(kdi, sR_d, Ri_i,R)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kdi, sR_d, Rr_i,R)*CC(dii, did, sr_i, sr_d,s*r, N_Qd, N_Qi)
      
      // info. source      
      - nu*PA(kdi, sR_d, RI_i,R)*CC(dii, did, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - nu*PA(kdi, sR_d, Ri_i,R)*CC(dii, did, si_i, si_d,s*i, N_Qd, N_Qi)

      + nu*PA(kid, IS_i, SR_d,S)*CC(idi, idd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + nu*PA(kid, iS_i, SR_d,S)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // info. loss
      + lm*sr_d - lm*sR_d;      

//----------------      
   iR_d_t =
      
      // infection
      - b2*PA(kdd, Is_d, sR_d,s)*CC(ddi, ddd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - bi*PA(kdd, is_d, sR_d,s)*CC(ddi, ddd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // recovery
      + gd*iI_d - gi*iR_d
      
      // susceptibility
      - dd*iR_d
      
      // info. trans.
      + al*PA(kid, sI_i, IR_d,I)*CC(idi, idd, sR_i, sR_d,s*R, N_Qd, N_Qi)
      + al*PA(kid, iI_i, IR_d,I)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      + al*PA(kid, rI_i, IR_d,I)*CC(idi, idd, rR_i, rR_d,r*R, N_Qd, N_Qi)

      - al*PA(kdi, iR_d, Rs_i,R)*CC(dii, did, is_i, is_d,i*s, N_Qd, N_Qi)
      - al*PA(kdi, iR_d, Ri_i,R)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - al*PA(kdi, iR_d, Rr_i,R)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kid, II_i, IR_d,I)*CC(idi, idd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + nu*PA(kid, iI_i, IR_d,I)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)

      - nu*PA(kdi, iR_d, RI_i,R)*CC(dii, did, iI_i, iI_d,i*I, N_Qd, N_Qi)
      - nu*PA(kdi, iR_d, Ri_i,R)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)

      // info. loss
      + lm*ir_d - lm*iR_d;        

//----------------
   rR_d_t =
      
      // recovery
      + gi*iR_d + gd*rI_d
      
      // susceptibility
      - di*rR_d - dd*rR_d
      
      // info. trans.
      + al*PA(kid, sR_i, RR_d,R)*CC(idi, idd, sR_i, sR_d,s*R, N_Qd, N_Qi)
      + al*PA(kid, iR_i, RR_d,R)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)
      + al*PA(kid, rR_i, RR_d,R)*CC(idi, idd, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kdi, rR_d, Rs_i,R)*CC(dii, did, rs_i, rs_d,r*s, N_Qd, N_Qi)
      - al*PA(kdi, rR_d, Ri_i,R)*CC(dii, did, ri_i, ri_d,r*i, N_Qd, N_Qi)
      - al*PA(kdi, rR_d, Rr_i,R)*CC(dii, did, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kid, IR_i, RR_d,R)*CC(idi, idd, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + nu*PA(kid, iR_i, RR_d,R)*CC(idi, idd, iR_i, iR_d,i*R, N_Qd, N_Qi)

      - nu*PA(kdi, rR_d, RI_i,R)*CC(dii, did, rI_i, rI_d,r*I, N_Qd, N_Qi)
      - nu*PA(kdi, rR_d, Ri_i,R)*CC(dii, did, ri_i, ri_d,r*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*rr_d - lm*rR_d ;


   //--------------------------------------------------------------
   // un-informed i-edges dynamics
   //--------------------------------------------------------------
   
//----------------
   SS_i_t =

      2.0*(
      
      // infection
      - bd*PA(kid, SS_i, SI_d,S)*CC(idi, idd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - b1*PA(kid, SS_i, Si_d,S)*CC(idi, idd, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // susceptibility
      + dd*SR_i
      
      // info. trans.
      - al*PA(kii, SS_i, Ss_i,S)*CC(iii, iid, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      - al*PA(kii, SS_i, Si_i,S)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      - al*PA(kii, SS_i, Sr_i,S)*CC(iii, iid, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kii, SS_i, SI_i,S)*CC(iii, iid, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - nu*PA(kii, SS_i, Si_i,S)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*sS_i );
   
//----------------   
   SI_i_t =
      
      // infection
      + bd*PA(kid, SS_i, SI_d,S)*CC(idi, idd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + b1*PA(kid, SS_i, Si_d,S)*CC(idi, idd, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      - bd*PA(kdi, IS_d, SI_i,S)*CC(dii, did, II_i, II_d,I*I, N_Qd, N_Qi)
      - b1*PA(kdi, iS_d, SI_i,S)*CC(dii, did, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      // recovery
      - gd*SI_i
      
      // susceptibility
      + dd*RI_i
      
      // info. trans.
      - al*PA(kii, sS_i, SI_i,S)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - al*PA(kii, iS_i, SI_i,S)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      - al*PA(kii, rS_i, SI_i,S)*CC(iii, iid, rI_i, rI_d,r*I, N_Qd, N_Qi)
      
      - al*PA(kii, SI_i, Is_i,I)*CC(iii, iid, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      - al*PA(kii, SI_i, Ii_i,I)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      - al*PA(kii, SI_i, Ir_i,I)*CC(iii, iid, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      // info. source
      - nu*SI_i
      
      - nu*PA(kii, IS_i, SI_i,S)*CC(iii, iid, II_i, II_d,I*I, N_Qd, N_Qi)
      - nu*PA(kii, iS_i, SI_i,S)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      - nu*PA(kii, SI_i, II_i,I)*CC(iii, iid, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - nu*PA(kii, SI_i, Ii_i,I)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*sI_i + lm*Si_i;
   
//----------------
   SR_i_t =

      // infection
      - bd*PA(kdi, IS_d, SR_i,S)*CC(dii, did, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - b1*PA(kdi, iS_d, SR_i,S)*CC(dii, did, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // recovery
      + gd*SI_i

      // susceptibility
      - dd*SR_i + dd*RR_i
      
      // info. trans.
      - al*PA(kii, sS_i, SR_i,S)*CC(iii, iid, sR_i, sR_d,s*R, N_Qd, N_Qi)
      - al*PA(kii, iS_i, SR_i,S)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      - al*PA(kii, rS_i, SR_i,S)*CC(iii, iid, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kii, SR_i, Rs_i,R)*CC(iii, iid, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      - al*PA(kii, SR_i, Ri_i,R)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      - al*PA(kii, SR_i, Rr_i,R)*CC(iii, iid, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kii, IS_i, SR_i,S)*CC(iii, iid, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - nu*PA(kii, iS_i, SR_i,S)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      - nu*PA(kii, SR_i, RI_i,R)*CC(iii, iid, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - nu*PA(kii, SR_i, Ri_i,R)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*sR_i + lm*Sr_i;

//----------------
   II_i_t =

      2.0*(
      
      // infection
      + bd*PA(kid, IS_i, SI_d,S)*CC(idi, idd, II_i, II_d,I*I, N_Qd, N_Qi)
      + b1*PA(kid, IS_i, Si_d,S)*CC(idi, idd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // recovery
      - gd*II_i
      
      // info. trans.
      - al*PA(kii, II_i, Is_i,I)*CC(iii, iid, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - al*PA(kii, II_i, Ii_i,I)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - al*PA(kii, II_i, Ir_i,I)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // info. source
      - nu*II_i
      
      - nu*PA(kii, II_i, II_i,I)*CC(iii, iid, II_i, II_d,I*I, N_Qd, N_Qi)
      - nu*PA(kii, II_i, Ii_i,I)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*Ii_i );        

//----------------
   IR_i_t =
            
      // infection
      + bd*PA(kdi, IS_d, SR_i,S)*CC(dii, did, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + b1*PA(kdi, iS_d, SR_i,S)*CC(dii, did, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // recovery
      + gd*II_i - gd*IR_i
      
      // susceptibility
      - dd*IR_i
      
      // info. trans.
      - al*PA(kii, sI_i, IR_i,I)*CC(iii, iid, sR_i, sR_d,s*R, N_Qd, N_Qi)
      - al*PA(kii, iI_i, IR_i,I)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      - al*PA(kii, rI_i, IR_i,I)*CC(iii, iid, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kii, IR_i, Rs_i,R)*CC(iii, iid, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - al*PA(kii, IR_i, Ri_i,R)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - al*PA(kii, IR_i, Rr_i,R)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // info. source
      - nu*IR_i
      
      - nu*PA(kii, II_i, IR_i,I)*CC(iii, iid, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - nu*PA(kii, iI_i, IR_i,I)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      - nu*PA(kii, IR_i, RI_i,R)*CC(iii, iid, II_i, II_d,I*I, N_Qd, N_Qi)
      - nu*PA(kii, IR_i, Ri_i,R)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*iR_i + lm*Ir_i;      

//----------------
   RR_i_t =

      2.0*(
      
      // recovery
      + gd*IR_i
      
      // susceptibility
      - dd*RR_i
      
      // info. trans.
      - al*PA(kii, RR_i, Rs_i,R)*CC(iii, iid, Rs_i, Rs_d,R*s, N_Qd, N_Qi)
      - al*PA(kii, RR_i, Ri_i,R)*CC(iii, iid, Ri_i, Ri_d,R*i, N_Qd, N_Qi)
      - al*PA(kii, RR_i, Rr_i,R)*CC(iii, iid, Rr_i, Rr_d,R*r, N_Qd, N_Qi)
      
      // info. source
      - nu*PA(kii, RR_i, RI_i,R)*CC(iii, iid, RI_i, RI_d,R*I, N_Qd, N_Qi)
      - nu*PA(kii, RR_i, Ri_i,R)*CC(iii, iid, Ri_i, Ri_d,R*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*Rr_i );
      
   //--------------------------------------------------------------
   // informed i-edges dynamics
   //--------------------------------------------------------------
   
//----------------
   ss_i_t =

      2.0*(
      
      // infection
      - bi*PA(kid, ss_i, si_d,s)*CC(idi, idd, si_i, si_d,s*i, N_Qd, N_Qi)
      - b2*PA(kid, ss_i, sI_d,s)*CC(idi, idd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      
      // susceptibility
      + di*sr_i
      
      // info. trans.
      + al*Ss_i
      
      + al*PA(kii, sS_i, Ss_i,S)*CC(iii, iid, ss_i, ss_d,s*s, N_Qd, N_Qi)
      + al*PA(kii, sS_i, Si_i,S)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kii, sS_i, Sr_i,S)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kii, sS_i, SI_i,S)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + nu*PA(kii, sS_i, Si_i,S)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      
      // info. loss
      - lm*ss_i );

//----------------
   si_i_t =
      
      // infection
      + bi*PA(kid, ss_i, si_d,s)*CC(idi, idd, si_i, si_d,s*i, N_Qd, N_Qi)
      + b2*PA(kid, ss_i, sI_d,s)*CC(idi, idd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - bi*PA(kdi, is_d, si_i,s)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - b2*PA(kdi, Is_d, si_i,s)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      // recovery
      - gi*si_i

      // susceptibility
      + di*ri_i
      
      // info. trans.
      + al*sI_i + al*Si_i
      
      + al*PA(kii, sI_i, Is_i,I)*CC(iii, iid, ss_i, ss_d,s*s, N_Qd, N_Qi)
      + al*PA(kii, sI_i, Ii_i,I)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kii, sI_i, Ir_i,I)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)

      + al*PA(kii, sS_i, Si_i,S)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kii, iS_i, Si_i,S)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + al*PA(kii, rS_i, Si_i,S)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)      
      
      // info. source
      + nu*Si_i

      + nu*PA(kii, iS_i, Si_i,S)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + nu*PA(kii, IS_i, Si_i,S)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      + nu*PA(kii, sI_i, II_i,I)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + nu*PA(kii, sI_i, Ii_i,I)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      
      // info. loss
      - 2.0 * lm*si_i;
      //- lm*si_i; 

//----------------
   sr_i_t =

      // infection
      - bi*PA(kdi, is_d, sr_i,s)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      - b2*PA(kdi, Is_d, sr_i,s)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // recovery
      + gi*si_i

      // susceptibility
      - di*sr_i + di*rr_i
      
      // info. trans.
      + al*sR_i + al*Sr_i
      
      + al*PA(kii, sR_i, Rs_i,R)*CC(iii, iid, ss_i, ss_d,s*s, N_Qd, N_Qi)
      + al*PA(kii, sR_i, Ri_i,R)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      + al*PA(kii, sR_i, Rr_i,R)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)

      + al*PA(kii, sS_i, Sr_i,S)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)
      + al*PA(kii, iS_i, Sr_i,S)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + al*PA(kii, rS_i, Sr_i,S)*CC(iii, iid, rr_i, rr_d,r*r, N_Qd, N_Qi)      
      
      // info. source
      + nu*PA(kii, iS_i, Sr_i,S)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + nu*PA(kii, IS_i, Sr_i,S)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      + nu*PA(kii, sR_i, RI_i,R)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + nu*PA(kii, sR_i, Ri_i,R)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      
      // info. loss
      - 2.0 * lm*sr_i;
      //- lm*sr_i; 
   
//----------------
   ii_i_t =

      2.0*(
      
      // infection
      + bi*PA(kid, is_i, si_d,s)*CC(idi, idd, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + b2*PA(kid, is_i, sI_d,s)*CC(idi, idd, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      // recovery
      - gi*ii_i
      
      // info. trans.
      + al*Ii_i
      
      + al*PA(kii, iI_i, Is_i,I)*CC(iii, iid, is_i, is_d,i*s, N_Qd, N_Qi)
      + al*PA(kii, iI_i, Ii_i,I)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + al*PA(kii, iI_i, Ir_i,I)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. source
      + nu*Ii_i
      
      + nu*PA(kii, iI_i, II_i,I)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + nu*PA(kii, iI_i, Ii_i,I)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // info. loss
      - lm*ii_i );       

//----------------
   ir_i_t =

      // infection
      + bi*PA(kdi, is_d, sr_i,s)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + b2*PA(kdi, Is_d, sr_i,s)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      // recovery
      + gi*ii_i - gi*ir_i
      
      // susceptibility
      - di*ir_i
      
      // info. trans.
      + al*iR_i + al*Ir_i
      
      + al*PA(kii, iR_i, Rs_i,R)*CC(iii, iid, is_i, is_d,i*s, N_Qd, N_Qi)
      + al*PA(kii, iR_i, Ri_i,R)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      + al*PA(kii, iR_i, Rr_i,R)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)

      + al*PA(kii, sI_i, Ir_i,I)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)
      + al*PA(kii, iI_i, Ir_i,I)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      + al*PA(kii, rI_i, Ir_i,I)*CC(iii, iid, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*iR_i
      
      + nu*PA(kii, II_i, Ir_i,I)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      + nu*PA(kii, iI_i, Ir_i,I)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      + nu*PA(kii, iR_i, RI_i,R)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + nu*PA(kii, iR_i, Ri_i,R)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // info. loss
      - 2 * lm*ir_i;
      //- lm*ir_i; 
   
//----------------
   rr_i_t =

      2.0*(
      
      // recovery
      + gi*ir_i
      
      // susceptibility
      - di*rr_i
      
      // info. trans.
      + al*rR_i
      
      + al*PA(kii, rR_i, Rs_i,R)*CC(iii, iid, rs_i, rs_d,r*s, N_Qd, N_Qi)
      + al*PA(kii, rR_i, Ri_i,R)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)
      + al*PA(kii, rR_i, Rr_i,R)*CC(iii, iid, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kii, rR_i, RI_i,R)*CC(iii, iid, rI_i, rI_d,r*I, N_Qd, N_Qi)
      + nu*PA(kii, rR_i, Ri_i,R)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)
      
      // info. loss
      - lm*rr_i );
   
   //--------------------------------------------------------------
   // mix i-edges dynamics
   //--------------------------------------------------------------
   
//----------------
   Ss_i_t =
      // infection
      - b2*PA(kid, Ss_i, sI_d,s)*CC(idi, idd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      - bi*PA(kid, Ss_i, si_d,s)*CC(idi, idd, Si_i, Si_d,S*i, N_Qd, N_Qi)

      - bd*PA(kdi, IS_d, Ss_i,S)*CC(dii, did, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - b1*PA(kdi, iS_d, Ss_i,S)*CC(dii, did, is_i, is_d,i*s, N_Qd, N_Qi)
      
      // susceptibility
      + di*Sr_i + dd*sR_i
      
      // info. trans.
      - al*Ss_i
      
      + al*PA(kii, SS_i, Ss_i,S)*CC(iii, iid, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      + al*PA(kii, SS_i, Si_i,S)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      + al*PA(kii, SS_i, Sr_i,S)*CC(iii, iid, Sr_i, Sr_d,S*r, N_Qd, N_Qi)
      
      - al*PA(kii, sS_i, Ss_i,S)*CC(iii, iid, ss_i, ss_d,s*s, N_Qd, N_Qi)
      - al*PA(kii, iS_i, Ss_i,S)*CC(iii, iid, is_i, is_d,i*s, N_Qd, N_Qi)
      - al*PA(kii, rS_i, Ss_i,S)*CC(iii, iid, rs_i, rs_d,r*s, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kii, SS_i, SI_i,S)*CC(iii, iid, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + nu*PA(kii, SS_i, Si_i,S)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      - nu*PA(kii, IS_i, Ss_i,S)*CC(iii, iid, Is_i, Is_d,I*s, N_Qd, N_Qi)
      - nu*PA(kii, iS_i, Ss_i,S)*CC(iii, iid, is_i, is_d,i*s, N_Qd, N_Qi)
      
      // info. loss
      - lm*Ss_i + lm*ss_i;
   
//----------------
   Si_i_t =
      
      // infection
      + b2*PA(kid, Ss_i, sI_d,s)*CC(idi, idd, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + bi*PA(kid, Ss_i, si_d,s)*CC(idi, idd, Si_i, Si_d,S*i, N_Qd, N_Qi)      

      - bd*PA(kdi, IS_d, Si_i,S)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - b1*PA(kdi, iS_d, Si_i,S)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // recovery
      - gi*Si_i
      
      // susceptibility
      + dd*Ri_i
      
      // info. trans.
      - al*Si_i
      
      + al*PA(kii, SI_i, Is_i,I)*CC(iii, iid, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      + al*PA(kii, SI_i, Ii_i,I)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      + al*PA(kii, SI_i, Ir_i,I)*CC(iii, iid, Sr_i, Sr_d,S*r, N_Qd, N_Qi)

      - al*PA(kii, sS_i, Si_i,S)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kii, iS_i, Si_i,S)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - al*PA(kii, rS_i, Si_i,S)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)      

      // info. source
      - nu*Si_i
      
      - nu*PA(kii, IS_i, Si_i,S)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - nu*PA(kii, iS_i, Si_i,S)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      + nu*PA(kii, SI_i, II_i,I)*CC(iii, iid, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + nu*PA(kii, SI_i, Ii_i,I)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*si_i - lm*Si_i;   

//----------------
   Sr_i_t =
      
      // infection
      - bd*PA(kdi, IS_d, Sr_i,S)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      - b1*PA(kdi, iS_d, Sr_i,S)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // recovery
      + gi*Si_i
      
      // susceptibility
      - di*Sr_i + dd*Rr_i
      
      // info. trans.
      - al*Sr_i
      
      + al*PA(kii, SR_i, Rs_i,R)*CC(iii, iid, Ss_i, Ss_d,S*s, N_Qd, N_Qi)
      + al*PA(kii, SR_i, Ri_i,R)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      + al*PA(kii, SR_i, Rr_i,R)*CC(iii, iid, Sr_i, Sr_d,S*r, N_Qd, N_Qi)

      - al*PA(kii, sS_i, Sr_i,S)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)
      - al*PA(kii, iS_i, Sr_i,S)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      - al*PA(kii, rS_i, Sr_i,S)*CC(iii ,iid, rr_i, rr_d,r*r, N_Qd, N_Qi)      
      
      // info. source
      + nu*PA(kii, SR_i, RI_i,R)*CC(iii, iid, SI_i, SI_d,S*I, N_Qd, N_Qi)
      + nu*PA(kii, SR_i, Ri_i,R)*CC(iii, iid, Si_i, Si_d,S*i, N_Qd, N_Qi)
      
      - nu*PA(kii, IS_i, Sr_i,S)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      - nu*PA(kii, iS_i, Sr_i,S)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. loss
      + lm*sr_i - lm*Sr_i;      

//----------------
   sI_i_t =
      
      // infection
      + bd*PA(kid, sS_i, SI_d,S)*CC(idi, idd, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + b1*PA(kid, sS_i, Si_d,S)*CC(idi, idd, si_i, si_d,s*i, N_Qd, N_Qi)

      - b2*PA(kdi, Is_d, sI_i,s)*CC(dii, did, II_i, II_d,I*I, N_Qd, N_Qi)
      - bi*PA(kdi, is_d, sI_i,s)*CC(dii, did, iI_i, iI_d,i*I, N_Qd, N_Qi)
      
      // recovery
      - gd*sI_i
      
      // susceptibility
      + di*rI_i
      
      // info. trans.
      - al*sI_i
      
      - al*PA(kii, sI_i, Is_i,I)*CC(iii, iid, ss_i, ss_d,s*s, N_Qd, N_Qi)
      - al*PA(kii, sI_i, Ii_i,I)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kii, sI_i, Ir_i,I)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)

      + al*PA(kii, sS_i, SI_i,S)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      + al*PA(kii, iS_i, SI_i,S)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + al*PA(kii, rS_i, SI_i,S)*CC(iii, iid, rI_i, rI_d,r*I, N_Qd, N_Qi)      
      
      // info. source
      + nu*SI_i
      
      - nu*PA(kii, sI_i, II_i,I)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - nu*PA(kii, sI_i, Ii_i,I)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      
      + nu*PA(kii, iS_i, SI_i,S)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      + nu*PA(kii, IS_i, SI_i,S)*CC(iii, iid, II_i, II_d,I*I, N_Qd, N_Qi)
      
      // info. loss
      - lm*sI_i + lm*si_i;

//----------------   
   Ii_i_t =
      
      // infection
      + b2*PA(kid, Is_i, sI_d,s)*CC(idi, idd, II_i, II_d,I*I, N_Qd, N_Qi)
      + bi*PA(kid, Is_i, si_d,s)*CC(idi, idd, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      + bd*PA(kdi, IS_d, Si_i,S)*CC(dii, did, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      + b1*PA(kdi, iS_d, Si_i,S)*CC(dii, did, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // recovery
      - gd*Ii_i - gi*Ii_i
      
      // info. trans.
      - al*Ii_i

      + al*PA(kii, II_i, Is_i,I)*CC(iii, iid, Is_i, Is_d,I*s, N_Qd, N_Qi)
      + al*PA(kii, II_i, Ii_i,I)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      + al*PA(kii, II_i, Ir_i,I)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)

      - al*PA(kii, sI_i, Ii_i,I)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kii, iI_i, Ii_i,I)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - al*PA(kii, rI_i, Ii_i,I)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)
      
      // info. source
      + nu*II_i - nu*Ii_i
      
      + nu*PA(kii, II_i, II_i,I)*CC(iii, iid, II_i, II_d,I*I, N_Qd, N_Qi)
      + nu*PA(kii, II_i, Ii_i,I)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      
      - nu*PA(kii, II_i, Ii_i,I)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      - nu*PA(kii, iI_i, Ii_i,I)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*ii_i - lm*Ii_i;       

//----------------
   Ir_i_t =
      
      // infection
      + bd*PA(kdi, IS_d, Sr_i,S)*CC(dii, did, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      + b1*PA(kdi, iS_d, Sr_i,S)*CC(dii, did, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // recovery
      + gi*Ii_i - gd*Ir_i
      
      // susceptibility
      - di*Ir_i
      
      // info. trans.
      - al*Ir_i
      
      + al*PA(kii, IR_i, Rs_i,R)*CC(iii, iid, Is_i, Is_d,I*s, N_Qd, N_Qi)
      + al*PA(kii, IR_i, Ri_i,R)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)
      + al*PA(kii, IR_i, Rr_i,R)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      
      - al*PA(kii, sI_i, Ir_i,I)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)
      - al*PA(kii, iI_i, Ir_i,I)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      - al*PA(kii, rI_i, Ir_i,I)*CC(iii, iid, rr_i, rr_d,r*r, N_Qd, N_Qi)
            
      // info. source
      + nu*IR_i
      
      + nu*PA(kii, IR_i, RI_i,R)*CC(iii, iid, II_i, II_d,I*I, N_Qd, N_Qi)
      + nu*PA(kii, IR_i, Ri_i,R)*CC(iii, iid, Ii_i, Ii_d,I*i, N_Qd, N_Qi)

      - nu*PA(kii, II_i, Ir_i,I)*CC(iii, iid, Ir_i, Ir_d,I*r, N_Qd, N_Qi)
      - nu*PA(kii, iI_i, Ir_i,I)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. loss
      + lm*ir_i - lm*Ir_i;        

//----------------
   sR_i_t =      
      
      // infection
      - b2*PA(kdi, Is_d, sR_i,s)*CC(dii, did, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - bi*PA(kdi, is_d, sR_i,s)*CC(dii, did, iR_i, iR_d,i*R, N_Qd, N_Qi)      
      
      // recovery
      + gd*sI_i
      
      // susceptibility
      - dd*sR_i + di*rR_i
      
      // info. trans.
      - al*sR_i

      + al*PA(kii, sS_i, SR_i,S)*CC(iii, iid, sR_i, sR_d,s*R, N_Qd, N_Qi)
      + al*PA(kii, iS_i, SR_i,S)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      + al*PA(kii, rS_i, SR_i,S)*CC(iii ,iid, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kii, sR_i, Rs_i,R)*CC(iii, iid, ss_i, ss_d,s*s, N_Qd, N_Qi)
      - al*PA(kii, sR_i, Ri_i,R)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)
      - al*PA(kii, sR_i, Rr_i,R)*CC(iii, iid, sr_i, sr_d,s*r, N_Qd, N_Qi)
      
      // info. source      
      - nu*PA(kii, sR_i, RI_i,R)*CC(iii, iid, sI_i, sI_d,s*I, N_Qd, N_Qi)
      - nu*PA(kii, sR_i, Ri_i,R)*CC(iii, iid, si_i, si_d,s*i, N_Qd, N_Qi)

      + nu*PA(kii, IS_i, SR_i,S)*CC(iii, iid, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + nu*PA(kii, iS_i, SR_i,S)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // info. loss
      + lm*sr_i - lm*sR_i;      
      
//----------------
   iR_i_t =
      
      // infection
      - b2*PA(kdi, Is_d, sR_i,s)*CC(dii, did, IR_i, IR_d,I*R, N_Qd, N_Qi)
      - bi*PA(kdi, is_d, sR_i,s)*CC(dii, did, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      // recovery
      + gd*iI_i - gi*iR_i
      
      // susceptibility
      - dd*iR_i
      
      // info. trans.
      -al*iR_i
      
      + al*PA(kii, sI_i, IR_i,I)*CC(iii, iid, sR_i, sR_d,s*R, N_Qd, N_Qi)
      + al*PA(kii, iI_i, IR_i,I)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      + al*PA(kii, rI_i, IR_i,I)*CC(iii, iid, rR_i, rR_d,r*R, N_Qd, N_Qi)

      - al*PA(kii, iR_i, Rs_i,R)*CC(iii, iid, is_i, is_d,i*s, N_Qd, N_Qi)
      - al*PA(kii, iR_i, Ri_i,R)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)
      - al*PA(kii, iR_i, Rr_i,R)*CC(iii, iid, ir_i, ir_d,i*r, N_Qd, N_Qi)
      
      // info. source
      - nu*iR_i
      
      + nu*PA(kii, II_i, IR_i,I)*CC(iii, iid, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + nu*PA(kii, iI_i, IR_i,I)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)

      - nu*PA(kii, iR_i, RI_i,R)*CC(iii, iid, iI_i, iI_d,i*I, N_Qd, N_Qi)
      - nu*PA(kii, iR_i, Ri_i,R)*CC(iii, iid, ii_i, ii_d,i*i, N_Qd, N_Qi)

      // info. loss
      + lm*ir_i - lm*iR_i;        
   
//----------------
   rR_i_t =
      
      // recovery
      + gi*iR_i + gd*rI_i
      
      // susceptibility
      - di*rR_i - dd*rR_i
      
      // info. trans.
      - al*rR_i
      
      + al*PA(kii, sR_i, RR_i,R)*CC(iii, iid, sR_i, sR_d,s*R, N_Qd, N_Qi)
      + al*PA(kii, iR_i, RR_i,R)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      + al*PA(kii, rR_i, RR_i,R)*CC(iii, iid, rR_i, rR_d,r*R, N_Qd, N_Qi)
      
      - al*PA(kii, rR_i, Rs_i,R)*CC(iii, iid, rs_i, rs_d,r*s, N_Qd, N_Qi)
      - al*PA(kii, rR_i, Ri_i,R)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)
      - al*PA(kii, rR_i, Rr_i,R)*CC(iii, iid, rr_i, rr_d,r*r, N_Qd, N_Qi)
      
      // info. source
      + nu*PA(kii, IR_i, RR_i,R)*CC(iii, iid, IR_i, IR_d,I*R, N_Qd, N_Qi)
      + nu*PA(kii, iR_i, RR_i,R)*CC(iii, iid, iR_i, iR_d,i*R, N_Qd, N_Qi)
      
      - nu*PA(kii, rR_i, RI_i,R)*CC(iii, iid, rI_i, rI_d,r*I, N_Qd, N_Qi)
      - nu*PA(kii, rR_i, Ri_i,R)*CC(iii, iid, ri_i, ri_d,r*i, N_Qd, N_Qi)
      
      // info. loss
      + lm*rr_i - lm*rR_i ;
   
   
#ifdef DEBUG   

   double tmp1=0;
   double tmp2=0;
   for (int jj=0; jj<6; jj++)
   {
      tmp1+=y[jj];
      tmp2+=rhs[jj];
   }

   double tmp3=0;
   double tmp4=0;
   for (jj=6; jj<48; jj++)
   {
      tmp3+=y[jj];
      tmp4+=rhs[jj];
   }
   
   //cout.setf(ios::showpos | ios:: showpoint | ios::left);
   //cout.setf(ios::scientific | ios::uppercase);
   cout.precision(5);
   cout.setf (ios_base::left, ios_base::basefield );
   cout.setf (ios_base::fixed, ios_base::floatfield );
   cout.width(12);
   
   cout << "single: y= " << setw(12) << tmp1
        << "   rhs= " << setw(12) << tmp2
        << "   pair: y= " << setw(12) << tmp3
        << "   rhs= " << setw(12) << tmp4  << endl;

#endif
   return GSL_SUCCESS;
}

/******************************************************************/

void ModelOde::InitParameters()
{
  // calculate Qd and Qi in pair approximation
  if (GetNvars() > 6) {
    
    // calculate total number of disease pairs
    double dPairs =
      SS_d_t + SI_d_t + SR_d_t + Ss_d_t + Si_d_t + Sr_d_t +
      II_d_t + IR_d_t + Is_d_t + Ii_d_t + Ir_d_t +
      RR_d_t + Rs_d_t + Ri_d_t + Rr_d_t +
      ss_d_t + si_d_t + sr_d_t +
      ii_d_t + ir_d_t +
      rr_d_t;
    
    SetQd(2*dPairs / GetN());
    
    // calculate total number of disease pairs
    double iPairs =
      SS_i_t + SI_i_t + SR_i_t + Ss_i_t + Si_i_t + Sr_i_t +
      II_i_t + IR_i_t + Is_i_t + Ii_i_t + Ir_i_t +
      RR_i_t + Rs_i_t + Ri_i_t + Rr_i_t +
      ss_i_t + si_i_t + sr_i_t +
      ii_i_t + ir_i_t +
      rr_i_t;
    
    SetQi(2*iPairs / GetN());

  } else {
    SetQd(0);
    SetQi(0);
  }

  SetQdi(Q_di);
   
  SetCddi(0.0);
  SetCddd(0.0);
  SetCdii(0.0);
  SetCdid(0.0);
  SetCidi(0.0);
  SetCidd(0.0);
  SetCiii(0.0);
  SetCiid(0.0);
   
}
    

/******************************************************************/
