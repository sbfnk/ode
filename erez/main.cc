/******************************************************************/

#include <iostream>
#include <cstdlib>
#include "ode.hh"
#include "van_der_pol.hh"
#include "ode_io_utils.hh"

using namespace std;

/******************************************************************/

int main(int argc, char *argv[])
{
   VanDerPol test;

   /********/

   if(argc<2)
   {
      cout << "Usage: a.out [parameters file]\n"
           << "             e.g. file.prm\n";
      exit(1);
   }

   cout << endl
        << "Starting ODE solver\n"
        << "-------------------\n";

   // reading ode parameters 
   ReadOdeParams(argv[1], test);

   // reading model parameters 
   ReadModelParams(argv[1], test);

   // printing ode + model parameters 
   test.PrtModelPrms();
   
   // solve the system    
   test.OdeSolve();
   
   return 0;
}

   
   
