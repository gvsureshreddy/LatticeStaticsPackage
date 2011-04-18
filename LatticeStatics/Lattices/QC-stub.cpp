#include <iostream>
#include <cstdlib>

using std::cerr;

extern "C" void qcbfb_energy_(int* mode, int& nfree, double* u, double& t, double& E, double* Eu,
                              double* Euu, double* Eut)
{
   cerr << "Error QC object is not designed for stand-alone use.\n";
   exit(-1);
}

extern "C" void qcbfb_restart_(char* filename, int const& n)
{
   cerr << "Error QC object is not designed for stand-alone use.\n";
   exit(-1);
}

extern "C" void qcbfb_output_(int& nfree, double* u, double& prop, int& nint, int* intdata, int& ndouble, double* doubledata)
{
   cerr << "Error QC object is not designed for stand-alone use.\n";
   exit(-1);
}

extern "C" void qcbfb_scan_(int& nfree, double* u, double& prop, double& dprop, double* tangent, double& alphamin, double& alphamax, int& steps)
{
   cerr << "erro not standalone\n";
   exit(-1);
}

