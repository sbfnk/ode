/********************************************************************/

#include <iostream>
#include "ode.h"
#include "van_der_pol.h"
#include "io_utils.h"

using namespace std;

/********************************************************************/

int main(int argc, char *argv[])
{
   VanDerPol test;

   /**********/
   
   read_ode_params("init.dat", &test);
   read_model_params("init.dat", &test);
   
   test.plugin_model_params();
   
   test.prt_model_prms();

   /**********/

   test.ode_solve();
   
   return 0;
}

   
   
