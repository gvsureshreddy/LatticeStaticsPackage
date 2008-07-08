#include "KnownLattices.h"
#include "PerlInput.h"
#include <fstream>

char *builddate();

using namespace std;

void GetMainSettings(int& Width,int& Presision,PerlInput const& Input);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
           << " ParamFile GridSize Echo" << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n";
      exit(-1);
   }
   
   char *datafile = argv[1];
   int const GridSize = atoi(argv[2]);

   PerlInput Input;
   Input.Readfile(datafile,"Input File:");
   
   Lattice *Lat;
   
   int Width,Precision,Echo=0;
   
   if (argc == 4) Echo = 1;
   
   GetMainSettings(Width,Precision,Input);
   
   Lat = InitializeLattice(Input,Echo);
   Lat->SetGridSize(GridSize);
   
   cout  << setiosflags(ios::fixed) << setprecision(Precision);
   
   Vector DOF((Lat->DOF()).Dim()),K(3);
   int BlochWaveStable;
   double Temp,Lambda;
   
   while (cin >> Temp)
   {
      cin >> Lambda;
      cin >> DOF;
      
      Lat->SetTemp(Temp);
      Lat->SetLambda(Lambda);
      Lat->SetDOF(DOF);
      
      BlochWaveStable = Lat->BlochWave(K);
      
      cout << "BlochWave Stability (GridSize=" << GridSize << "):"
           << setw(Width) << BlochWaveStable << ", "
           << setw(Width) << K << "\n" << flush;
   }
   
   return 1;
}

void GetMainSettings(int& Width,int& Precision,PerlInput const& Input)
{
   Width = Input.getInt("Main","FieldWidth");
   Precision = Input.getInt("Main","Precision");
}
