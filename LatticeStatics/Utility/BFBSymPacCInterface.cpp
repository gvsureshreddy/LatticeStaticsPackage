#include <iostream>
#include "Vector.h"
using std::cout;

extern "C" void get_qc_(int* mode,int& nfree,double* u,double& t,double& E,double* Eu,
                        double* Euu,double* Eut);

extern "C" void bfb_wrapper_(int& nfree)
{
   cout << "inside ryan_wrapper_()\n";
   cout << "-->" << nfree << "<--\n";
   Vector u(nfree,0.0);
   double t=0.0;

   double E=0.0;
   Vector Eu(nfree);

   get_qc_(0,nfree,&(u[0]),t,E,&(Eu[0]),0,0);

   cout << "E = " << E << "\n";
   cout << "Eu.Norm() = " << Eu.Norm() << "\n";

}

// pgif90 -o xxxxx xxxxx.o -lBFBSymPac -lstdc++ `perl =MExtUtils::Embed -e ldopts`
