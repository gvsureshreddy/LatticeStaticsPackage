#include "PerlInput.h"
#include "KnownLattices.h"
#include "KnownRestrictions.h"
#include "KnownSolutionMethods.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"

#include <string>
#include <fstream>

char* builddate();

using namespace std;

void GetMainSettings(int& Width, int& Precision, int& Echo, PerlInput const& Input);
void InitializeOutputFile(fstream& out, char const* const outfile, char const* const datafile,
                          char const* const startfile, int const& Precision, int const& Width,
                          int const& Echo);

// define as global to allow outatexit to flush
fstream out;
void outatexit(void)
{
   out.flush();
   cout.flush();
   cerr.flush();
   cout << "outatexit called" << endl;
}


int main(int argc, char* argv[])
{
   // Check commandline arguments
   if ((argc < 3) || ((argc < 4) && (!strcmp(argv[1], "--debug"))))
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
   if (!strcmp(argv[1], "--debug"))
   {
      Debug = 1;
   }
   else
   {
      Debug = 0;
   }

   char* datafile = argv[1 + Debug],
   * outputfile = argv[2 + Debug],
   * startfile;

   PerlInput Input(datafile);

   if (argc == 4 + Debug)
   {
      startfile = argv[3 + Debug];
      Input.Readfile(startfile);
   }
   else
   {
      startfile = 0;
   }

   Lattice* Lat;
   Restriction* Restrict;
   SolutionMethod* SolveMe;

   int Width, Precision, Echo;

   GetMainSettings(Width, Precision, Echo, Input);

   // out declared at global level
   InitializeOutputFile(out, outputfile, datafile, startfile, Precision, Width, Echo);
   // flush buffers when exit() is called
   atexit(outatexit);

   Lat = InitializeLattice(Input, Echo, Width, Debug);
   Lat->Print(out, Lattice::PrintLong, Lattice::NotSolutionPt);

   Restrict = InitializeRestriction(Lat, Input);

   out << "Restriction: " << Restrict->Name() << "\n";
   if (Echo)
   {
      cout << "Restriction: " << Restrict->Name() << "\n";
   }

   SolveMe = InitializeSolution(Restrict, Input, Lat, out, Width, Echo);

   int success = 1;

   while ((success) && (!SolveMe->AllSolutionsFound()))
   {
      success = SolveMe->FindNextSolution(Input, Width, out);
   }
   if (!success)
   {
      cout << "Error encountered. Exiting." << endl;
   }

   out.close();
   delete SolveMe;
   delete Restrict;
   delete Lat;

   return 0;
}

void GetMainSettings(int& Width, int& Precision, int& Echo, PerlInput const& Input)
{
   Width = Input.getInt("Main", "FieldWidth");
   Precision = Input.getInt("Main", "Precision");
   if (Input.ParameterOK("Main", "Echo"))
   {
      Echo = Input.getInt("Main", "Echo");
   }
   else
   {
      Echo = Input.useInt(1, "Main", "Echo"); // Default value
   }
   Input.EndofInputSection();
}

void InitializeOutputFile(fstream& out, char const* const outfile, char const* const datafile,
                          char const* const startfile, int const& Precision, int const& Width,
                          int const& Echo)
{
   fstream input, start;
   string dataline;

   input.open(datafile, ios::in);
   if (input.fail())
   {
      cerr << "Error: Unable to open file : " << datafile << " for read"
           << "\n";
      exit(-1);
   }
   out.open(outfile, ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << outfile << " for write"
           << "\n";
      exit(-1);
   }

   while (!input.eof())
   {
      getline(input, dataline);
      out << "Input File:" << dataline << "\n";
   }

   input.close();

   if (startfile != 0)
   {
      start.open(startfile, ios::in);
      if (start.fail())
      {
         cerr << "Error: Unable to open file : " << startfile << " for read"
              << "\n";
         exit(-1);
      }

      while (!start.eof())
      {
         getline(start, dataline);
         out << "Start File:" << dataline << "\n";
      }

      start.close();
   }

   if (Echo)
   {
      cout << setiosflags(ios::fixed) << setprecision(Precision);
   }
   out << setiosflags(ios::fixed) << setprecision(Precision);

   if (Echo)
   {
      cout << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n"
           << setw(Width);
   }
   out << "Built on:               " << builddate() << "\n"
       << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
       << "MyMath Built on:        " << MyMathBuildDate() << "\n"
       << setw(Width);
}
