/********************************************************************/  

#include "io_utils.h"
#include "ode.h"
#include "van_der_pol.h"

using namespace std;

/********************************************************************/  

void read_ode_params(const char *ifile_name, Ode *obj)
{
   string line;
   ifstream ifile;
   
   /*****/
   
   // open init file
   cout << "... opening ode input file: " << ifile_name << endl;
   
   ifile.open (ifile_name, ios::in);  
   if(ifile.fail())
   {
      cout << "... unable to open " << ifile_name << endl;
      exit(1); // uses stdlib.h
   }
   else
   {
      cout << "... done\n\n";
   }
   
   /*****/
   
   // reading Ode content of init.dat
   cout << "... reading ode parameters\n\n"; 
   
   int i=-1;
   while (i<0)
   {
      getline(ifile,line); // reading a line
      i=line.find("Ode parameters"); // looking for "Ode ... %%%"
   }
   
   /* ode output file name */
   obj->set_ofile_name(read_str_val(&ifile));
  
   /* tmax */
   obj->set_tmax(read_dbl_val(&ifile));
  
   /* dt */
   obj->set_dt(read_dbl_val(&ifile));
  
   /* nsave */
   obj->set_nsave(read_int_val(&ifile));
  
   /* nvars */
   obj->set_nvars(read_int_val(&ifile));

   /* step_algo */
   obj->set_step_algo(read_str_val(&ifile));

   /* abs_tol, rel_tol */
   obj->set_tol(read_dbl_val(&ifile),read_dbl_val(&ifile));
   
   /* ode i.c. file name */
   obj->set_ic_file_name(read_str_val(&ifile));

   /*****/
  
   // closing init.dat
   cout << "... closing ode input file: " << ifile_name << endl;
  
   ifile.close(); 
   if(!ifile.is_open())
      cout << "... done\n\n";
}

/********************************************************************/  

void read_model_params(const char *ifile_name, VanDerPol *obj)
{
   string line;
   ifstream ifile;
   ModelParams *model_params = new ModelParams;
   
   /*****/
   
   // open init file
   cout << "... opening model input file: " << ifile_name << endl;
   
   ifile.open (ifile_name, ios::in); 
   if(ifile.fail())
   {
      cout << "... unable to open " << ifile_name << endl;
      exit(1); // uses stdlib.h
   }
   else
   {
      cout << "... done\n\n";
   }
   
   /*****/
   
   // reading Model content of init.dat   
   cout << "... reading model parameters\n\n"; 

   int i=-1;
   while (i<0)
   {
      getline(ifile,line); // reading a line
      i=line.find("Model parameters"); // looking for "Model ... %%%"
   }
  
   /* beta_d */   
   model_params->beta_d=read_dbl_val(&ifile);
  
   /* gamma_d */
   model_params->gamma_d=read_dbl_val(&ifile);
  
   /* delta_d */
   model_params->delta_d=read_dbl_val(&ifile);
      
   /* beta_i */
   model_params->beta_i=read_dbl_val(&ifile);
   
   /* gamma_i */
   model_params->gamma_i=read_dbl_val(&ifile);
   
   /* delta_i */
   model_params->delta_i=read_dbl_val(&ifile);
   
   /* beta_m */
   model_params->beta_m=read_dbl_val(&ifile);
   
   /* alpha */
   model_params->alpha=read_dbl_val(&ifile);
   
   /* nu */
   model_params->nu=read_dbl_val(&ifile);
   
   /* lambda */
   model_params->lambda=read_dbl_val(&ifile);
   
   /* mu1 */
   model_params->mu1=read_dbl_val(&ifile);
   
   /* mu2 */
   model_params->mu2=read_dbl_val(&ifile);

   /* njac */
   model_params->njac=obj->get_nvars();
   
   /*****/
   
   // plugin parameters to class object
   obj->set_model_params(model_params); // model_params will be delete
                                        // by object destructor
   
   
   // closing init.dat
   cout << "... closing model input file: " << ifile_name << endl;
   
   ifile.close(); 
   if(!ifile.is_open())
      cout << "... done\n\n";
}

/********************************************************************/  

int read_int_val(ifstream *ifile)
{
   string line, str_tmp;
   
   getline(*ifile,line); 
   int i=line.find("="); // looking for the =  
   if(i>=0) // check if = was found 
   {
      str_tmp=line.substr(i+2,line.length()-i); // extracting substring
   }
   else
   {
      cout << "... problem reading params, no = was found\n";
      exit(1);
   }

   return atoi(str_tmp.c_str());
}

/********************************************************************/  

double read_dbl_val(ifstream *ifile)
{
   string line, str_tmp;
   
   getline(*ifile,line); 
   int i=line.find("="); // looking for the =  
   if(i>=0) // check if = was found 
   {
      str_tmp=line.substr(i+2,line.length()-i); // extracting substring
   }
   else
   {
      cout << "... problem reading params, no = was found\n";
      exit(1);
   }

   return atof(str_tmp.c_str());

}

/********************************************************************/

const char *read_str_val(ifstream *ifile)
{
   string line, str_tmp;
   char *cstr=new char[64];
   
   getline(*ifile,line); 
   int i=line.find("="); // looking for the =    
   if(i>=0) // check if = was found
   {
      str_tmp=line.substr(i+2,line.length()-i);
   }
   else
   {
      cout << "... problem reading params, no = was found\n";
      exit(1);
   }

   strcpy(cstr, str_tmp.c_str());
   return cstr;
}

/********************************************************************/








