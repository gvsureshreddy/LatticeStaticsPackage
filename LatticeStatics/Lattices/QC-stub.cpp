#include <iostream>

using std::cerr;

extern "C" void get_qc_(int* mode,int& nfree,double* u,double& t,double& E,double* Eu,
                        double* Euu,double* Eut)
{
   cerr << "Error QC object is not designed for stand-alone use.\n";
   exit(-1);
}
