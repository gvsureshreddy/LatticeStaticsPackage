#include <iostream>
#include <fstream>
#include "PerlInput.h"
#include "KnownLattices.h"
#include "KnownRestrictions.h"
#include "KnownSolutionMethods.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"

#include <string>
#include <sstream>

char* builddate();

using namespace std;

void GetMainSettings(int& Width, int& Precision, int& Echo, PerlInput const& Input);
void InitializeOutputFile(char const* const datafile, char const* const startfile,
                          int const& Precision, int const& Width, int const& Echo);

PerlInput Input;
Lattice* Lat;
Restriction* Restrict;
SolutionMethod* SolveMe;
int success = 1;
int TestValue;
int StableValue = -1;
Vector EigenValues;
int Width, Precision, Echo;

extern "C" void bfb_gettang_(double* tang)
{
   Vector tan1tmp(Input.getArrayLength("StartType","Tangent"));
   Input.getVector(tan1tmp,"StartType","Tangent");

   for (int i=0;i<tan1tmp.Dim()-1;++i)
   {
      tang[i] = tan1tmp[i];
   }

   return;
}

extern "C" void bfb_init_wrapper_(int& nfree, double* ufree_init, double& t, char* bfbfile)
{
   Vector utmp(nfree);
   // Vector sol(nfree+1);
   for (int i = 0; i < nfree; ++i)
   {
      utmp[i] = ufree_init[i];
      // sol[i] = ufree_init[i];
   }
   // sol[nfree] = t;

   Input.Readfile(bfbfile);
   ostringstream tmp;
   tmp << "$Main{Echo} = 0;";
   Input.EvaluateString(tmp.str().c_str());
   tmp.str("");
   tmp << "$Lattice{QC}{DOFS} = " << nfree << ";";
   Input.EvaluateString(tmp.str().c_str());
   tmp.str("");
   tmp << "$Lattice{NewtonPCSolution}{NumSolutions} = 100000;";
   Input.EvaluateString(tmp.str().c_str());
   GetMainSettings(Width, Precision, Echo, Input);

   InitializeOutputFile(bfbfile, 0, Precision, Width, Echo);

   Lat = InitializeLattice(Input, Echo, Width, 0);
   Lat->SetDOF(utmp);
   Lat->SetLambda(t);

   Restrict = InitializeRestriction(Lat, Input);

   cout << "Restriction: " << Restrict->Name() << "\n";

   // Let RefineEqbmSolution take care of giving the first solution, utmp, to NewtonPCSolution
   // tmp.str("");
   // tmp << scientific << setprecision(Precision);
   // tmp << "$StartType{Solution} = [";
   // tmp << setw(Width) << sol[0];
   // for (int i=1;i<sol.Dim();++i)
   //   tmp << "," << setw(Width) << sol[i];
   // tmp << "];";
   // Input.EvaluateString(tmp.str().c_str());

   SolveMe = InitializeSolution(Restrict, Input, Lat, cout, Width, Echo);
   EigenValues.Resize(Lat->NumTestFunctions());

   cout << setw(Width);
   Lat->Print(cout, Lattice::PrintLong);
}

// bfbreturncode: 0-regular point, 1-terminate, 2-critical point
extern "C" void bfb_wrapper_(int& bfbstable, int& bfbreturncode)
{
   success = SolveMe->FindNextSolution(Input, Width, cout);
   // always returns 1

   TestValue = Lat->TestFunctions(EigenValues);
   StableValue = 0;
   for (int i = Lat->DOF().Dim(); i < EigenValues.Dim(); ++i)
   {
      if (EigenValues[i] < 0.0)
      {
         ++StableValue;
      }
   }
   if (0 == StableValue)
   {
      bfbstable = 0;
   }
   else
   {
      bfbstable = 1;
   }

   if (success == 2)
   {
      bfbreturncode = 2;
   }
   else
   {
      bfbreturncode = 0;
   }

   if ((!success) || (SolveMe->AllSolutionsFound()))
   {
      if (!success)
      {
         cout << "Error encountered. Exiting." << "\n";
      }
      bfbreturncode = 1;
   }
}

extern "C" void bfb_term_wrapper_()
{
   delete SolveMe;
   delete Restrict;
   delete Lat;
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

void InitializeOutputFile(char const* const datafile, char const* const startfile,
                          int const& Precision, int const& Width, int const& Echo)
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

   while (!input.eof())
   {
      getline(input, dataline);
      cout << "Input File:" << dataline << "\n";
   }

   input.close();

   cout << setiosflags(ios::fixed) << setprecision(Precision);

   cout << "Built on:               " << builddate() << "\n"
        << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
        << "MyMath Built on:        " << MyMathBuildDate() << "\n"
        << setw(Width);
}
