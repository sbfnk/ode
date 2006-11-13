/******************************************************************/

#include <iostream>
#include <cstdlib>
#include "ode.hh"
#include "mfa.hh"
#include "ode_io_utils.hh"

using namespace std;

/******************************************************************/

int main(int argc, char *argv[])
{
   MeanField Mfa;
   
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
   ReadOdeParams(argv[1], Mfa);

   // reading model parameters 
   ReadModelParams(argv[1], Mfa);

   // printing ode + model parameters 
   Mfa.PrtModelPrms();
   
   // solve the system    
   Mfa.OdeSolve();
   
   return 0;
}

   
   
