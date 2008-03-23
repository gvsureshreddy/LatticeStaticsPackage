#include "KnownLattices.h"
#include "PerlInput.h"
#include <fstream>

char *builddate();

using namespace std;

void GetMainSettings(int &Width,int &Presision,PerlInput &Input);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
           << " ParamFile cutoff <Echo>" << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n";
      exit(-1);
   }
   
   char *datafile = argv[1];
   int cutoff = atoi(argv[2]);

   PerlInput Input;
   Input.Readfile(datafile,"Input File:");
   
   Lattice *Lat;
   
   int Width,Precision,Echo=0;
   
   if (argc == 4) Echo = 1;
   
   GetMainSettings(Width,Precision,Input);
   
   Lat = InitializeLattice(Input,Echo);
   
   cout  << setiosflags(ios::fixed) << setprecision(Precision);
   
   Vector DOF((Lat->DOF()).Dim());
   double Temp;
   
   while (cin >> Temp)
   {
      cin >> DOF;
      
      Lat->SetTemp(Temp);
      Lat->SetDOF(DOF);
      cout << setw(Width);
      Lat->NeighborDistances(cutoff,cout);
   }
   
   return 1;
}

void GetMainSettings(int &Width, int &Precision,PerlInput &Input)
{
   Width = Input.getInt("Main","FieldWidth");
   Precision = Input.getInt("Main","Precision");
}
