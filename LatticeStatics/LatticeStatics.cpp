#include "KnownLattices.h"
#include "KnownModes.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"

#include <fstream.h>

enum YN {No,Yes};
void GetMainSettings(int &Width, int &Precision,YN &BisectCP,char *datafile);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,char *startfile,
			  Lattice *Lat,int Precision,int Width);
SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out, int Width);

int main(int argc, char *argv[])
{
   // Check commandline arguments
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile OutputFile <StartData>" << endl;
#ifdef BUILD_DATE
      cerr << "Built on: " << BUILD_DATE << endl;
#endif
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

   int Width,Precision;
   YN BisectCP;

   
   GetMainSettings(Width,Precision,BisectCP,datafile);

   Lat = InitializeLattice(datafile);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,startfile,Lat,Precision,Width);
   
   Mode = InitializeMode(Lat,datafile);
   out << "Mode: " << Mode->ModeName() << endl;
   cout << "Mode: " << Mode->ModeName() << endl;
   
   SolveMe = InitializeSolution(Mode,datafile,startfile,Lat,out,Width);


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
	    SolveMe->BisectAlert(Lat,Width,out);
	 
	 // Send Output
	 cout << setw(Width) << Lat
	      << "Uncertainty = " << setw(Width) << uncertainty << endl
	      << "Success = 1" << endl;
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




void GetMainSettings(int &Width, int &Precision,YN &BisectCP,char *datafile)
{
   char bisect[LINELENGTH];

   GetParameter("^MainFieldWidth",datafile,"%d",&Width);
   GetParameter("^MainPrecision",datafile,"%d",&Precision);
   
   GetParameter("^MainBisectCP",datafile,"%s",bisect);
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
			  Lattice *Lat,int Precision,int Width)
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
      out << "Input File: " << dataline << endl;
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
	 out << "Start File: " << dataline << endl;
      }

      start.close();
   }
   
   cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);

#ifdef BUILD_DATE
   cout << "Built on: " << BUILD_DATE << endl << setw(Width);
   out << "Built on: " << BUILD_DATE << endl << setw(Width);
#endif
   
   Lat->Print(cout,Lattice::PrintLong);
   Lat->Print(out,Lattice::PrintLong);
}



SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out,int Width)
{
   char slvmthd[LINELENGTH];

   enum solution {Scanning,ArcLen};
   solution solu;

   GetParameter("^MainSolutionMethod",datafile,"%s",slvmthd);
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
	 return new ScanningSolution(Mode,datafile);
      }
      break;
      case ArcLen:
      {
	 int good;
	 Vector One = Mode->ArcLenDef(),
	    Two = Mode->ArcLenDef();

	 if ( startfile == NULL)
	 {
	    ScanningSolution ScanMe(Mode,datafile);
	    
	    while (!ScanMe.AllSolutionsFound())
	    {
	       One = Two;
	       ScanMe.FindNextSolution(good);
	       if (good)
	       {
		  cout << setw(Width) << Lat;
		  Two = Mode->ArcLenDef();
	       }
	    }

	    return new ArcLengthSolution(Mode,datafile,One,Two);
	 }
	 else
	 {
	    return new ArcLengthSolution(Mode,datafile,startfile,out);
	 }
      }
      break;
   }

   return NULL;
}
