/******************************************************************/  

#include "ode_io_utils.hh"
#include "ode.hh"
#include "mfa.hh"

using namespace std;

/******************************************************************/  

void ReadOdeParams(const char *ifile_name, Ode& obj)
{
   string line;
   ifstream ifile;
   
   // open init file
   cout << "... opening ode input file: " << ifile_name;
   
   ifile.open (ifile_name, ios::in);  
   if(ifile.fail())
   {
      cout << "... unable to open " << ifile_name << endl;
      exit(1); // uses stdlib.h
   }
   else
   {
      cout << " ... done\n";
   }
   
   // reading Ode content of init.dat
   cout << "... reading ode parameters"; 
   
   int i=-1;
   while (i<0)
   {
      getline(ifile,line); // reading a line
      i=line.find("Ode parameters"); // looking for "Ode ... %%%"
   }
   
   // ode output file name 
   obj.SetoFileName(read_str_val(ifile));
  
   // tmax 
   obj.SetTmax(read_dbl_val(ifile));
   
   // dt 
   obj.SetDt(read_dbl_val(ifile));
   
   // nsave 
   obj.SetNsave(read_int_val(ifile));
   
   // nvars 
   obj.SetNvars(read_int_val(ifile));

   // step_algo 
   obj.SetStepAlgo(read_str_val(ifile));

   // abs_tol 
   obj.SetAbsTol(read_dbl_val(ifile));

   // rel_tol 
   obj.SetRelTol(read_dbl_val(ifile));
   
   // ode i.c. file name 
   obj.SeticFileName(read_str_val(ifile));

   cout << " ... done\n";
   
   // closing init.dat
   cout << "... closing ode input file: " << ifile_name;
  
   ifile.close(); 
   if(!ifile.is_open())
      cout << " ... done\n";
}

/******************************************************************/  

void ReadModelParams(const char *ifile_name, MeanField& obj)
{
   string line;
   ifstream ifile;
   ModelParams *model_params = new ModelParams;
   
   // open init file
   cout << "... opening model input file: " << ifile_name;
   
   ifile.open (ifile_name, ios::in); 
   if(ifile.fail())
   {
      cout << "... unable to open " << ifile_name << endl;
      exit(1); // uses stdlib.h
   }
   else
   {
      cout << " ... done\n";
   }
   
   // reading Model content of init.dat   
   cout << "... reading model parameters"; 
   
   int i=-1;
   while (i<0)
   {
      getline(ifile,line); // reading a line
      i=line.find("Model parameters"); // looking for "Model ... %%%"
   }
   
   // beta_d    
   model_params->beta_d=read_dbl_val(ifile);
  
   // gamma_d 
   model_params->gamma_d=read_dbl_val(ifile);
  
   // delta_d 
   model_params->delta_d=read_dbl_val(ifile);
      
   // beta_i 
   model_params->beta_i=read_dbl_val(ifile);
   
   // gamma_i 
   model_params->gamma_i=read_dbl_val(ifile);
   
   // delta_i 
   model_params->delta_i=read_dbl_val(ifile);
   
   // beta_m1 
   model_params->beta_m1=read_dbl_val(ifile);
   
   // beta_m2 
   model_params->beta_m2=read_dbl_val(ifile);
   
   // alpha 
   model_params->alpha=read_dbl_val(ifile);
   
   // nu 
   model_params->nu=read_dbl_val(ifile);
   
   // lambda 
   model_params->lambda=read_dbl_val(ifile);
   
   // Qd 
   model_params->Qd=read_dbl_val(ifile);
   
   // Qi 
   model_params->Qi=read_dbl_val(ifile);
   
   // N
   model_params->N=read_dbl_val(ifile);
   
   // njac 
   model_params->njac=obj.GetNvars();

   cout << " ... done\n";
   
   // passing a pointer to model_params to Ode class member _params 
   obj.SetModelParams(static_cast<void *>(model_params));
   
   // closing init.dat
   cout << "... closing model input file: " << ifile_name;
   
   ifile.close(); 
   if(!ifile.is_open())
      cout << " ... done\n";
}

/******************************************************************/  

int read_int_val(ifstream& ifile)
{
   string line, str_tmp;
   
   getline(ifile,line); 
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

/******************************************************************/  

double read_dbl_val(ifstream& ifile)
{
   string line, str_tmp;
   
   getline(ifile,line); 
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

/******************************************************************/

const char *read_str_val(ifstream& ifile)
{
   string line, str_tmp;
   char *cstr=new char[64];
   
   getline(ifile,line); 
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

/******************************************************************/








