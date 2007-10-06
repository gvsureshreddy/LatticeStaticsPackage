#include "KnownLattices.h"
#include "KnownModes.h"
#include "KnownSolutionMethods.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"

#include <fstream>

using namespace std;

enum YN {No,Yes};
void GetMainSettings(int &Width, int &Precision,YN &BisectCP,int &Echo,char *datafile);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,char *startfile,
			  int Precision,int Width,int Echo);

int main(int argc, char *argv[])
{
   // Check commandline arguments
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " [--debug]"
	   << " ParamFile OutputFile <StartData>" << endl;
      cerr << "Built on:               " << builddate() << endl
	   << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
	   << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }
   int Debug;
   if (!strcmp(argv[1],"--debug"))
      Debug=1;
   else
      Debug=0;

   char *datafile = argv[1+Debug],
      *outputfile = argv[2+Debug],
      *startfile;
   
   if (argc == 4+Debug)
      startfile = argv[3+Debug];
   else
      startfile = NULL;   
   

   Lattice *Lat;
   LatticeMode *Mode;
   SolutionMethod *SolveMe;

   int Width,Precision,Echo;
   YN BisectCP;

   if (GetParameter("^","CommandFile",datafile,'s',UTILITYechocommand,0))
   {
      UTILITYechocmd = UTILITYechocommand;
   }
   GetMainSettings(Width,Precision,BisectCP,Echo,datafile);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,startfile,Precision,Width,Echo);
   
   Lat = InitializeLattice(datafile,"^",Echo,Width,Debug);
   Lat->Print(out,Lattice::PrintLong);

   Mode = InitializeMode(Lat,datafile,"^");
   SolveMe = InitializeSolution(Mode,datafile,startfile,Lat,out,Width,Echo);

   out << "Mode: " << Mode->ModeName() << endl;
   if (Echo) cout << "Mode: " << Mode->ModeName() << endl;
   
   int success = 1;
   int Nulity=-1,
      OldNulity=Nulity;

   
   while (!SolveMe->AllSolutionsFound())
   {
      SolveMe->FindNextSolution(success);

      if (success)
      {	 
	 // Check for Critical Point Crossing
	 OldNulity = Nulity;
	 Nulity = Lat->StiffnessNulity();
	 if ((OldNulity != Nulity) && (BisectCP == Yes) && (OldNulity != -1))
	    SolveMe->BisectAlert(Lat,datafile,"^",Width,out);
	 
	 // Send Output
	 out << setw(Width) << Lat << "Success = 1" << endl;
      }
   }

   out.close();
   delete SolveMe;
   delete Mode;
   delete Lat;

   return 0;
}




void GetMainSettings(int &Width, int &Precision,YN &BisectCP,int &Echo,char *datafile)
{
   char bisect[LINELENGTH];

   if(!GetParameter("^","MainFieldWidth",datafile,'i',&Width)) exit(-1);
   if(!GetParameter("^","MainPrecision",datafile,'i',&Precision)) exit(-1);
   if(!GetParameter("^","MainEcho",datafile,'i',&Echo,0)) Echo = 1;
   if(!GetParameter("^","MainBisectCP",datafile,'s',bisect)) exit(-1);
   if ((!strcmp("Yes",bisect)) || (!strcmp("yes",bisect)))
      BisectCP = Yes;
   else if ((!strcmp("No",bisect)) || (!strcmp("no",bisect)))
      BisectCP = No;
   else
   {
      cerr << "Unknown BisectCP option : " << bisect << endl;
      exit(-1);
   }
}


void InitializeOutputFile(fstream &out,char *outfile,char *datafile,char *startfile,
			  int Precision,int Width,int Echo)
{
   fstream input,start;
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

   while (!input.eof())
   {
      input.getline(dataline,LINELENGTH-1);
      out << "Input File:" << dataline << endl;
   }

   input.close();

   if (startfile != NULL)
   {
      start.open(startfile,ios::in);
      if (start.fail())
      {
	 cerr << "Error: Unable to open file : " << startfile << " for read"
	      << endl;
	 exit(-1);
      }

      while (!start.eof())
      {
	 start.getline(dataline,LINELENGTH-1);
	 out << "Start File:" << dataline << endl;
      }

      start.close();
   }
   
   if (Echo) cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);

   if (Echo) cout << "Built on:               " << builddate() << endl
		  << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
		  << "MyMath Built on:        " << MyMathBuildDate() << endl
		  << setw(Width);
   out << "Built on:               " << builddate() << endl
       << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
       << "MyMath Built on:        " << MyMathBuildDate() << endl
       << setw(Width);
}
