#include "KnownLattices.h"
#include "UtilityFunctions.h"
#include <fstream>

using namespace std;

void GetMainSettings(int &Width,int &Presision,char *datafile,const char *prefix);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile GridSize Echo" << endl;
      cerr << "Built on:               " << builddate() << endl
	   << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
	   << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }

   char *datafile = argv[1],
      prefix[LINELENGTH];
   int GridSize = atoi(argv[2]);

   strcpy(prefix,"^Input File:");
   
   Lattice *Lat;

   int Width,Precision,Echo=0;

   if (argc == 4) Echo = 1;

   GetMainSettings(Width,Precision,datafile,prefix);

   Lat = InitializeLattice(datafile,prefix,Echo);
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
	   << setw(Width) << K << endl << flush;
   }
   
   return 1;
}

void GetMainSettings(int &Width, int &Precision,char *datafile,const char *prefix)
{
   if(!GetParameter(prefix,"MainFieldWidth",datafile,"%d",&Width)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,"%d",&Precision)) exit(-1);   
}
