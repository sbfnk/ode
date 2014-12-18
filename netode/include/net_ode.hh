/*! \file net_ode.hh
*/
#ifndef NET_ODE_HH
#define NET_ODE_HH

bool read_single_degree_dist(std::vector<double>& degreeDist,
                             std::string fileName)
{
  degreeDist.clear();
  std::ifstream degreeFile;
  try {
    degreeFile.open(fileName.c_str(),
                    std::ios::in);
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    return false;
  }
  std::string line;
  if (degreeFile.is_open()) {
    getline(degreeFile, line);
    while(!degreeFile.eof()) {
      //read line
      if (line.size() > 0) {
        std::istringstream linestream(line);
        unsigned int degree = 0;
        double prob = 0.;
        linestream >> degree >> prob;
        if (degree + 1 > degreeDist.size() && prob > 0) {
          degreeDist.resize(degree+1, 0.);
          degreeDist[degree] = prob;
        }
      }
      getline(degreeFile, line);
    }
  } else {
    return false;
  }
  
  return true;
}

//------------------------------------------------------------

#endif // NET_ODE_HH
