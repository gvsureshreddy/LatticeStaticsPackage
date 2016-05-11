#include "TwoBarTruss.h"
#include "PerlInput.h"
#include <fstream>
#include <cstdlib>

using namespace std;

char* builddate();

void GetMainSettings(int& Width, int& Presision, int& Echo, PerlInput const& Input);
void InitializeOutputFile(fstream& out, char const* const outfile, char const* const datafile,
                          int const& argc, int const& Precision);

int main(int argc, char* argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
           << " ParamFile OutputFile" << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n";
      exit(-1);
   }

   char const* const datafile = argv[1];
   char const* const outputfile = argv[2];

   PerlInput Input(datafile);

   int Width, Precision, Echo;

   GetMainSettings(Width, Precision, Echo, Input);

   TwoBarTruss TBT(Input, Echo, Width);

   fstream out;
   InitializeOutputFile(out, outputfile, datafile, argc, Precision);

   Vector DOF(2);
   double lambda;

   Input.getVector(DOF, "TwoBarTrussModel", "DOF");
   lambda = Input.getDouble("TwoBarTrussModel", "Lambda");

   TBT.SetDOF(DOF);
   TBT.SetLambda(lambda);

   out << setw(Width) << TBT.E0() << "\n";
   out << setw(Width) << TBT.E1() << "\n";
   out << setw(Width) << TBT.E1DLoad() << "\n";
   out << setw(Width) << TBT.E2() << "\n";

   out.close();

   return 1;
}

void GetMainSettings(int& Width, int& Precision, int& Echo, PerlInput const& Input)
{
   Width = Input.getPosInt("Main", "FieldWidth");
   Precision = Input.getPosInt("Main", "Precision");
   Echo = Input.usePosInt(0, "Main", "Echo"); // Default Value
   Input.EndofInputSection();
}

void InitializeOutputFile(fstream& out, char const* const outfile, char const* const datafile,
                          int const& argc, int const& Precision)
{
   out.open(outfile, ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << outfile << " for write"
           << "\n";
      exit(-1);
   }

   out << setiosflags(ios::scientific) << setprecision(Precision);
}
