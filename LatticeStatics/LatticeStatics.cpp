#include "KnownLattices.h"
#include "KnownModes.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"
#include <string.h>

#include <fstream.h>

enum YN {No,Yes};
void GetMainSettings(int &Width, int &Precision,YN &BisectCP,char *datafile);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,
			  Lattice *Lat,int Precision,int Width);
SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,int Width);

int main(int argc, char *argv[])
{
   // Check commandline arguments
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile OutputFile <StartData>" << endl;
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
   InitializeOutputFile(out,outputfile,datafile,Lat,Precision,Width);
   
   Mode = InitializeMode(Lat,datafile);
   out << "Mode: " << Mode->ModeName() << endl;
   cout << "Mode: " << Mode->ModeName() << endl;
   
   SolveMe = InitializeSolution(Mode,datafile,startfile,Lat,Width);


   
   int success;
   int Nulity=Lat->StiffnessNulity(),
      OldNulity=Nulity;

   
   while (!SolveMe->AllSolutionsFound())
   {
      
      success=SolveMe->FindNextSolution();

      if (success)
      {	 
	 // Check for Critical Point Crossing
	 OldNulity = Nulity;
	 Nulity = Lat->StiffnessNulity();
	 if ((OldNulity != Nulity) && (BisectCP == Yes))
	    SolveMe->BisectAlert(Lat,Width,out);
	 
	 // Send Output
	 cout << setw(Width) << Lat << "Success = 1" << endl;
	 out << setw(Width) << Lat << "Success = 1" << endl;
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
   FILE *pipe;
   char command[LINELENGTH];

   char width[]="^MainFieldWidth";
   SetPerlCommand(command,datafile,width);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%d",&Width);
   if (pclose(pipe)) Errfun(width);

   char prec[]="^MainPrecision";
   SetPerlCommand(command,datafile,prec);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%d",&Precision);
   if (pclose(pipe)) Errfun(prec);

   char bisect[]="^MainBisectCP";
   SetPerlCommand(command,datafile,bisect);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%s",command);
   if (pclose(pipe)) Errfun(bisect);
   if ((!strcmp("Yes",command)) || (!strcmp("yes",command)))
      BisectCP = Yes;
   else if ((!strcmp("No",command)) || (!strcmp("no",command)))
      BisectCP = No;
   else
   {
      cerr << "Unknown BisectCP option : " << command << endl;
      exit(-1);
   }
}


void InitializeOutputFile(fstream &out,char *outfile,char *datafile,
			  Lattice *Lat,int Precision,int Width)
{
   fstream input;
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

   
   cout << setiosflags(ios::fixed) << setprecision(Precision) << setw(Width);
   out  << setiosflags(ios::fixed) << setprecision(Precision) << setw(Width);

   Lat->Print(cout,Lattice::PrintLong);
   Lat->Print(out,Lattice::PrintLong);
}



SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,int Width)
{
   FILE *pipe;
   char command[LINELENGTH];

   enum solution {Scanning,ArcLen};
   solution solu;

   
   char solv[]="^MainSolutionMethod";
   SetPerlCommand(command,datafile,solv);
   pipe=OpenPipe(command,"r");
   fscanf(pipe,"%s",command);
   if (pclose(pipe)) Errfun(solv);
   if ((!strcmp("Scanning",command))
       || (!strcmp("scanning",command)))
      solu = Scanning;
   else if ((!strcmp("ArcLength",command)) || (!strcmp("arclength",command)))
      solu = ArcLen;
   else
   {
      cerr << "Unknown SolutionMethod : " << command << endl;
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
	 ScanningSolution ScanMe(Mode,datafile);
	 Vector One = Mode->ArcLenDef(),
	    Two = Mode->ArcLenDef();

	 while (!ScanMe.AllSolutionsFound())
	 {
	    One = Two;
	    good = ScanMe.FindNextSolution();
	    if (good)
	    {
	       cout << setw(Width) << Lat;
	       Two = Mode->ArcLenDef();
	    }
	 }

	 return new ArcLengthSolution(Mode,datafile,One,Two);
      }
      break;
   }

   return NULL;
}
