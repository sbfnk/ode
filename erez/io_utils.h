/***********************************************************************

iostream - I/O to std output
fstream  - I/O to files
cstdlib  - exit()
string   - string class 
cstring  - strcpy,...
vector   - vector class

************************************************************************/

#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <iostream> 
#include <fstream>
#include <cstdlib>
#include <string>
#include <cstring>
#include "ode.h"
#include "van_der_pol.h"

using namespace std;

/********************************************************************/

void ReadOdeParams(const char *ifile_name, Ode *obj);
void ReadModelParams(const char *ifile_name, VanDerPol *obj);

int read_int_val(ifstream *ifile);
double read_dbl_val(ifstream *ifile);
const char *read_str_val(ifstream *ifile);

/********************************************************************/

#endif // IO_UTILS_H
