#include "KnownLattices.h"
#include "UtilityFunctions.h"
#include <fstream.h>

void GetMainSettings(int &Width,int &Presision,char *datafile);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,
			  Lattice *Lat,int Precision,int Width);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile OutputFile" << endl;
      cerr << "Built on:               " << builddate() << endl
	   << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
	   << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }

   char *datafile = argv[1],
      *outputfile = argv[2];

   Lattice *Lat;

   int Width,Precision;

   GetMainSettings(Width,Precision,datafile);

   Lat = InitializeLattice(datafile);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,Lat,Precision,Width);

   Lat->DispersionCurves(datafile,out);

   out.close();
   return 1;
}
   



void GetMainSettings(int &Width, int &Precision,char *datafile)
{
   if(!GetParameter("^MainFieldWidth",datafile,"%d",&Width)) exit(-1);
   if(!GetParameter("^MainPrecision",datafile,"%d",&Precision)) exit(-1);   
}

void InitializeOutputFile(fstream &out,char *outfile,char *datafile,
			  Lattice *Lat,int Precision,int Width)
{
   fstream input,devnull;
   char dataline[LINELENGTH];

   input.open(datafile,ios::in);
   if (input.fail())
   {
      cerr << "Error: Unable to open file : " << datafile << " for read"
	   << endl;
      exit(-1);
   }
   out.open(outfile,ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << outfile << " for write"
	   << endl;
      exit(-1);
   }
   devnull.open("/dev/null",ios::out);

   while (!input.eof())
   {
      input.getline(dataline,LINELENGTH-1);
      out << "#Input File: " << dataline << endl;
   }

   input.close();

   cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);

   cout << "Built on:               " << builddate() << endl
	<< "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
	<< "MyMath Built on:        " << MyMathBuildDate() << endl
	<< setw(Width);
   out << "#Built on:               " << builddate() << endl
       << "#LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
       << "#MyMath Built on:        " << MyMathBuildDate() << endl
       << setw(Width);

   devnull << setw(Width);
   Lat->Print(devnull,Lattice::PrintLong);

   devnull.close();
}
