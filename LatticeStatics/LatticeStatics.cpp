#include "KnownLattices.h"
#include "KnownModes.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"

#include <fstream.h>

enum YN {No,Yes};
void GetMainSettings(int &Width, int &Precision,YN &BisectCP,int &Echo,char *datafile);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,char *startfile,
			  Lattice *Lat,int Precision,int Width,int Echo);
SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out,int Width,
				   int Echo);

int main(int argc, char *argv[])
{
   // Check commandline arguments
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile OutputFile <StartData>" << endl;
      cerr << "Built on:               " << builddate() << endl
	   << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
	   << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }
   char *datafile = argv[1],
      *outputfile = argv[2],
      *startfile;
   if (argc == 4)
      startfile = argv[3];
   else
      startfile = NULL;   
   

   Lattice *Lat;
   LatticeMode *Mode;
   SolutionMethod *SolveMe;

   int Width,Precision,Echo;
   YN BisectCP;

   
   GetMainSettings(Width,Precision,BisectCP,Echo,datafile);

   Lat = InitializeLattice(datafile,"^",Echo);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,startfile,Lat,Precision,Width,Echo);
   
   Mode = InitializeMode(Lat,datafile,"^");
   out << "Mode: " << Mode->ModeName() << endl;
   if (Echo) cout << "Mode: " << Mode->ModeName() << endl;
   
   SolveMe = InitializeSolution(Mode,datafile,startfile,Lat,out,Width,Echo);


   double uncertainty;
   int success = 1;
   int Nulity=-1,
      OldNulity=Nulity;

   
   while (!SolveMe->AllSolutionsFound())
   {
      uncertainty = SolveMe->FindNextSolution(success);

      if (success)
      {	 
	 // Check for Critical Point Crossing
	 OldNulity = Nulity;
	 Nulity = Lat->StiffnessNulity();
	 if ((OldNulity != Nulity) && (BisectCP == Yes) && (OldNulity != -1))
	    SolveMe->BisectAlert(Lat,datafile,"^",Width,out);
	 
	 // Send Output
	 out << setw(Width) << Lat
	     << "Uncertainty = " << setw(Width) << uncertainty << endl
	     << "Success = 1" << endl;
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

   if(!GetParameter("^","MainFieldWidth",datafile,"%d",&Width)) exit(-1);
   if(!GetParameter("^","MainPrecision",datafile,"%d",&Precision)) exit(-1);
   if(!GetParameter("^","MainEcho",datafile,"%d",&Echo,0)) Echo = 1;
   if(!GetParameter("^","MainBisectCP",datafile,"%s",bisect)) exit(-1);
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
			  Lattice *Lat,int Precision,int Width,int Echo)
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
   
   Lat->Print(out,Lattice::PrintLong);
}

SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out,int Width,
				   int Echo)
{
   char slvmthd[LINELENGTH];

   enum solution {Scanning,ArcLen};
   solution solu;

   if(!GetParameter("^","MainSolutionMethod",datafile,"%s",slvmthd)) exit(-1);
   if ((!strcmp("Scanning",slvmthd))
       || (!strcmp("scanning",slvmthd)))
      solu = Scanning;
   else if ((!strcmp("ArcLength",slvmthd)) || (!strcmp("arclength",slvmthd)))
      solu = ArcLen;
   else
   {
      cerr << "Unknown SolutionMethod : " << slvmthd << endl;
      exit(-1);
   }

   switch (solu)
   {
      case Scanning:
      {
	 return new ScanningSolution(Mode,datafile,"^",Echo);
      }
      break;
      case ArcLen:
      {
	 int good;
	 Vector One = Mode->ArcLenDef(),
	    Two = Mode->ArcLenDef();

	 if ( startfile == NULL)
	 {
	    ScanningSolution ScanMe(Mode,datafile,"^",Echo);
	    
	    while (!ScanMe.AllSolutionsFound())
	    {
	       One = Two;
	       ScanMe.FindNextSolution(good);
	       if (good)
	       {
		  out << setw(Width) << Lat;
		  Two = Mode->ArcLenDef();
	       }
	    }

	    return new ArcLengthSolution(Mode,datafile,"^",One,Two,Echo);
	 }
	 else
	 {
	    return new ArcLengthSolution(Mode,datafile,"^",startfile,out,Echo);
	 }
      }
      break;
   }

   return NULL;
}
