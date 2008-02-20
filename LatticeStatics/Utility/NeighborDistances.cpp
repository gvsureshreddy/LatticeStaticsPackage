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
           << " ParamFile cutoff <Echo>" << endl;
      cerr << "Built on:               " << builddate() << endl
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
           << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }
   
   char *datafile = argv[1],
      prefix[LINELENGTH];
   int cutoff = atoi(argv[2]);
   
   strcpy(prefix,"^Input File:");
   
   Lattice *Lat;
   
   int Width,Precision,Echo=0;
   
   if (argc == 4) Echo = 1;
   
   GetMainSettings(Width,Precision,datafile,prefix);
   
   Lat = InitializeLattice(datafile,prefix,Echo);
   
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

void GetMainSettings(int &Width, int &Precision,char *datafile,const char *prefix)
{
   if(!GetParameter(prefix,"MainFieldWidth",datafile,'i',&Width)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,'i',&Precision)) exit(-1);
}
