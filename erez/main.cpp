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
   
   ReadOdeParams("init.dat", &test);
   ReadModelParams("init.dat", &test);
   
   test.PluginModelParams();
   
   test.PrtModelPrms();

   /**********/

   test.OdeSolve();
   
   return 0;
}

   
   
