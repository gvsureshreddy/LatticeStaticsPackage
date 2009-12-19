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

char *builddate();

using namespace std;

enum YN {No,Yes};
void GetMainSettings(int& Width,int& Precision,YN& BisectCP,int& Echo,PerlInput const& Input);
void InitializeOutputFile(char const* const datafile,char const* const startfile,
                          int const& Precision,int const& Width,int const& Echo);

PerlInput Input;
Lattice *Lat;
Restriction *Restrict;
SolutionMethod *SolveMe;
int success = 1;
int* TotalNumCPs = 0;
int TestValue=-1;
int StableValue=-1;
int OldTestValue=TestValue;
Vector EigenValues;
int Width,Precision,Echo;
YN BisectCP;

extern "C" void bfb_init_wrapper_(int& nfree,double* ufree_init,double& t,char* bfbfile)
{
   Vector utmp(nfree);
   //Vector sol(nfree+1);
   for (int i=0;i<nfree;++i)
   {
      utmp[i] = ufree_init[i];
      //sol[i] = ufree_init[i];
   }
   //sol[nfree] = t;
   
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
   GetMainSettings(Width,Precision,BisectCP,Echo,Input);
   
   InitializeOutputFile(bfbfile,0,Precision,Width,Echo);
   
   Lat = InitializeLattice(Input,Echo,Width,0);
   Lat->SetDOF(utmp);
   Lat->SetLambda(t);

   TotalNumCPs = new int[Lat->NumTestFunctions()];
   for (int i=0;i<Lat->NumTestFunctions();++i) TotalNumCPs[i] = 0;
   
   Restrict = InitializeRestriction(Lat,Input);
   
   cout << "Restriction: " << Restrict->Name() << "\n";

   // Let RefineEqbmSolution take care of giving the first solution, utmp, to NewtonPCSolution
   //tmp.str("");
   //tmp << scientific << setprecision(Precision);
   //tmp << "$StartType{Solution} = [";
   //tmp << setw(Width) << sol[0];
   //for (int i=1;i<sol.Dim();++i)
   //   tmp << "," << setw(Width) << sol[i];
   //tmp << "];";
   //Input.EvaluateString(tmp.str().c_str());
   
   SolveMe = InitializeSolution(Restrict,Input,Lat,cout,Width,Echo);
   EigenValues.Resize(Lat->NumTestFunctions());

   cout << setw(Width);
   Lat->Print(cout,Lattice::PrintLong);
}

// bfbreturncode: 0-regular point, 1-terminate, 2-critical point
extern "C" void bfb_wrapper_(int& bfbstable,int& bfbreturncode)
{
   success=SolveMe->FindNextSolution();
   // always returns 1
   
   // Check for Critical Point Crossing
   OldTestValue = TestValue;
   TestValue = Lat->TestFunctions(EigenValues);

   StableValue = 0;
   for (int i=Lat->DOF().Dim();i<EigenValues.Dim();++i)
   {
      if (EigenValues[i] < 0.0) ++StableValue;
   }
   if (0 == StableValue)
   {
      bfbstable = 0;
   }
   else
   {
      bfbstable = 1;
   }
   
   if ((OldTestValue != TestValue) && (BisectCP == Yes) && (OldTestValue != -1))
   {
      bfbreturncode = 2;
      SolveMe->FindCriticalPoint(Lat,TotalNumCPs,Input,Width,cout);
   }
   else
   {
      bfbreturncode = 0;
   }
   // Send Output
   cout << setw(Width) << *Lat << "Success = 1" << "\n";
   
   if (SolveMe->AllSolutionsFound())
   {
      bfbreturncode=1;
   }
}

extern "C" void bfb_term_wrapper_()
{
   delete SolveMe;
   delete Restrict;
   delete Lat;
   delete [] TotalNumCPs;
}

void GetMainSettings(int& Width,int& Precision,YN& BisectCP,int& Echo,PerlInput const& Input)
{
   string bisect;
   
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

void InitializeOutputFile(char const* const datafile,char const* const startfile,
                          int const& Precision,int const& Width,int const& Echo)
{
   fstream input,start;
   string dataline;
   
   input.open(datafile,ios::in);
   if (input.fail())
   {
      cerr << "Error: Unable to open file : " << datafile << " for read"
           << "\n";
      exit(-1);
   }
   
   while (!input.eof())
   {
      getline(input,dataline);
      cout << "Input File:" << dataline << "\n";
   }
   
   input.close();
   
   cout  << setiosflags(ios::fixed) << setprecision(Precision);
   
   cout << "Built on:               " << builddate() << "\n"
        << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
        << "MyMath Built on:        " << MyMathBuildDate() << "\n"
        << setw(Width);
}
