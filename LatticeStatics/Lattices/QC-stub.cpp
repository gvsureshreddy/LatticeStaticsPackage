#include <iostream>

using std::cerr;

extern "C" void qcbfb_energy_(int* mode,int& nfree,double* u,double& t,double& E,double* Eu,
                              double* Euu,double* Eut)
{
   cerr << "Error QC object is not designed for stand-alone use.\n";
   exit(-1);
}

extern "C" void qcbfb_restart_(char* filename)
{
   cerr << "Error QC object is not designed for stand-alone use.\n";
   exit(-1);
}
