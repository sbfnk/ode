/*******************************************************************/
//
// File ode_io_utils
// -----------------
//
// This file containes utility I/O functions for reading parameters
// from files and initialize Ode and derived classes objects. The
// main reading functions uses the Ode mutators extensively. In
// addition, a few more helping function are provided, mainly to
// help converting input strings into numbers.
//
// Functions:
// ----------
// ReadOdeParams - This function reads the Ode related parameters from
//                 a given file. The inputs for the function are the
//                 input file name and an Ode object. The function
//                 opens the input file, looks for the string
//                 "Ode parameters", and then strarts reading parameters
//                 by order ! Each parameter is then given to the object
//                 using the class mutators. When all parameters were
//                 read, the function closes the input file.
//
// ReadModelParams -  This function reads the Model related parameters
//                    from a given file. The inputs for the function are
//                    the input file name and a Model object. The function
//                    opens the input file, looks for the string
//                    "Model parameters", and then strarts reading
//                    parameters by order !
//                    The parameters structure for the Model object is
//                    actually allocated in this function as a local
//                    variable. The parameters that are read are stored
//                    in this local structure, and once all parameters
//                    were read, a pointer to this allocated structure
//                    is passed to the object, using the member function
//                    Ode::SetModelParams. Since the GSL functions
//                    require a void *, the relevant cast is done.
//                    This allocated memory is freed by the object's
//                    destructor. The function then closes the input file.
//
// read_..._val - These functions receive as an input a reference to an
//                input stream (ifstream) and read the next line. The
//                required value is extracted by finding the "=" sign,
//                which separates the definition from the value it self
//                in the input file. All parameters therefor should be
//                included in the input files as follows
//
//                step_type    = rkf45
//                abs_tol      = 1e-5
//
//                with nothing after the parameter value.
//
// Special notices:
// ----------------
// - The user must include the header file of the derived Ode class,
//   for which the parameters should be read.
// - The Read...Params functions receive as input a REFERENCE to an
//   object.
//
/******************************************************************/

#ifndef ODE_IO_UTILS_H
#define ODE_IO_UTILS_H

/******************************************************************/

#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include "ode.hh"
#include "mfa.hh"

using namespace std;

/******************************************************************/

void ReadOdeParams(const char *ifile_name, Ode& obj);
void ReadModelParams(const char *ifile_name, MeanField& obj);

int read_int_val(ifstream& ifile);
double read_dbl_val(ifstream& ifile);
const char *read_str_val(ifstream& ifile);

/******************************************************************/

#endif // ODE_IO_UTILS_H
