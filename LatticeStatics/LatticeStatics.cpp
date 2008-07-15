#include "PerlInput.h"
#include "KnownLattices.h"
#include "KnownModes.h"
#include "KnownSolutionMethods.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"
#define LINELENGTH 600

#include <fstream>

char *builddate();

using namespace std;

enum YN {No,Yes};
void GetMainSettings(int& Width,int& Precision,YN& BisectCP,int& Echo,PerlInput const& Input);
void InitializeOutputFile(fstream& out,char const* const outfile,char const* const datafile,
                          char const* const startfile,int const& Precision,int const& Width,
                          int const& Echo);

int main(int argc, char *argv[])
{
   // Check commandline arguments
   if ((argc < 3) || ((argc <4) && (!strcmp(argv[1],"--debug"))))
   {
      cerr << "Usage: " << argv[0]
           << " [--debug]"
           << " ParamFile OutputFile <StartData>" << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n";
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

   PerlInput Input(datafile);

   if (argc == 4+Debug)
   {
      startfile = argv[3+Debug];
      Input.Readfile(startfile);
   }
   else
   {
      startfile = 0;
   }
   
   Lattice *Lat;
   LatticeMode *Mode;
   SolutionMethod *SolveMe;
   
   int Width,Precision,Echo;
   YN BisectCP;

   GetMainSettings(Width,Precision,BisectCP,Echo,Input);
   
   fstream out;
   InitializeOutputFile(out,outputfile,datafile,startfile,Precision,Width,Echo);
   
   Lat = InitializeLattice(Input,Echo,Width,Debug);
   Lat->Print(out,Lattice::PrintLong);
   
   Mode = InitializeMode(Lat,Input);

   out << "Mode: " << Mode->ModeName() << "\n";
   if (Echo) cout << "Mode: " << Mode->ModeName() << "\n";
   
   SolveMe = InitializeSolution(Mode,Input,Lat,out,Width,Echo);

   int success = 1;
   int TestValue=-1,
      OldTestValue=TestValue;
   Vector EigenValues(Lat->DOF().Dim());
   
   while (!SolveMe->AllSolutionsFound())
   {
      success=SolveMe->FindNextSolution();
      
      if (success)
      {
         // Check for Critical Point Crossing
         OldTestValue = TestValue;
         TestValue = Lat->TestFunctions(EigenValues);
         if ((OldTestValue != TestValue) && (BisectCP == Yes) && (OldTestValue != -1))
            SolveMe->FindCriticalPoint(Lat,Input,Width,out);
         
         // Send Output
         out << setw(Width) << *Lat << "Success = 1" << "\n";
      }
   }
   
   out.close();
   delete SolveMe;
   delete Mode;
   delete Lat;
   
   return 0;
}

void GetMainSettings(int& Width,int& Precision,YN& BisectCP,int& Echo,PerlInput const& Input)
{
   char bisect[LINELENGTH];

   Width = Input.getInt("Main","FieldWidth");
   Precision = Input.getInt("Main","Precision");
   if (Input.ParameterOK("Main","Echo"))
   {
      Echo = Input.getInt("Main","Echo");
   }
   else
   {
      Echo = Input.useInt(1,"Main","Echo"); // Default value
   }
   char const* const bisectcp = Input.getString("Main","BisectCP");
   if (!strcmp("Yes",bisectcp))
      BisectCP = Yes;
   else if (!strcmp("No",bisectcp))
      BisectCP = No;
   else
   {
      cerr << "Unknown BisectCP option : " << bisect << "\n";
      exit(-1);
   }
   Input.EndofInputSection();
}

void InitializeOutputFile(fstream& out,char const* const outfile,char const* const datafile,
                          char const* const startfile,int const& Precision,int const& Width,
                          int const& Echo)
{
   fstream input,start;
   char dataline[LINELENGTH];
   
   input.open(datafile,ios::in);
   if (input.fail())
   {
      cerr << "Error: Unable to open file : " << datafile << " for read"
           << "\n";
      exit(-1);
   }
   out.open(outfile,ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << outfile << " for write"
           << "\n";
      exit(-1);
   }
   
   while (!input.eof())
   {
      input.getline(dataline,LINELENGTH-1);
      out << "Input File:" << dataline << "\n";
   }
   
   input.close();
   
   if (startfile != 0)
   {
      start.open(startfile,ios::in);
      if (start.fail())
      {
         cerr << "Error: Unable to open file : " << startfile << " for read"
              << "\n";
         exit(-1);
      }
      
      while (!start.eof())
      {
         start.getline(dataline,LINELENGTH-1);
         out << "Start File:" << dataline << "\n";
      }
      
      start.close();
   }
   
   if (Echo) cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);
   
   if (Echo) cout << "Built on:               " << builddate() << "\n"
                  << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
                  << "MyMath Built on:        " << MyMathBuildDate() << "\n"
                  << setw(Width);
   out << "Built on:               " << builddate() << "\n"
       << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
       << "MyMath Built on:        " << MyMathBuildDate() << "\n"
       << setw(Width);
}
