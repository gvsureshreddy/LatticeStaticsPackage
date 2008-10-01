#include <iostream>
#include "PerlInput.h"
#include "KnownLattices.h"
#include "KnownRestrictions.h"
#include "KnownSolutionMethods.h"
#include "ScanningSolution.h"
#include "ArcLengthSolution.h"

#include "UtilityFunctions.h"

#include <string>
#include <sstream>
#include <fstream>

char *builddate();

using namespace std;

enum YN {No,Yes};
void GetMainSettings(int& Width,int& Precision,YN& BisectCP,int& Echo,PerlInput const& Input);
void InitializeOutputFile(fstream& out,char const* const outfile,char const* const datafile,
                          char const* const startfile,int const& Precision,int const& Width,
                          int const& Echo);

extern "C" void get_qc_(int& mode,int& nfree,double* u,double& t,double& E,double* Eu,
                        double* Euu,double* Eut);

PerlInput Input;
Lattice *Lat;
Restriction *Restrict;
SolutionMethod *SolveMe;
int success = 1;
int TotalNumCPs = 0;
int TestValue=-1;
int OldTestValue=TestValue;
Vector EigenValues;
int Width,Precision,Echo;
YN BisectCP;
fstream out;

extern "C" void bfb_init_wrapper_(int& nfree,double* ufree_init,double& t,char* bfbfile)
{
  ostringstream outfile;
  outfile << bfbfile << ".out";

  Vector utmp(nfree);
  Vector sol(nfree+1);
  for (int i=0;i<nfree;++i)
  {
    utmp[i] = ufree_init[i];
    sol[i] = ufree_init[i];
  }
  sol[nfree] = t;


     Input.Readfile(bfbfile);
     ostringstream tmp;
     tmp << "$Lattice{QC}{DOFS} = " << nfree << ";";
     Input.EvaluateString(tmp.str().c_str());
     GetMainSettings(Width,Precision,BisectCP,Echo,Input);
     
     InitializeOutputFile(out,outfile.str().c_str(),bfbfile,0,Precision,Width,Echo);

     Lat = InitializeLattice(Input,Echo,Width,0);
     Lat->SetDOF(utmp);
     Lat->SetLambda(t);
     //Lat->Print(out,Lattice::PrintLong);
     
     Restrict = InitializeRestriction(Lat,Input);
     
     out << "Restriction: " << Restrict->Name() << "\n";
     if (Echo) cout << "Restriction: " << Restrict->Name() << "\n";
     
     tmp.str("");
     tmp << scientific << setprecision(Precision);
     tmp << "$StartType{Solution} = [";
     tmp << setw(Width) << sol[0];
     for (int i=1;i<sol.Dim();++i)
        tmp << "," << setw(Width) << sol[i];
     tmp << "];";
     Input.EvaluateString(tmp.str().c_str());

     SolveMe = InitializeSolution(Restrict,Input,Lat,out,Width,Echo);
     EigenValues.Resize(Lat->DOF().Dim());
}

// bfbreturncode: 0-regular point, 1-terminate, 2-critical point
extern "C" void bfb_wrapper_(int& bfbreturncode)
{
   success=SolveMe->FindNextSolution();
   // always returns 1
   
   // Check for Critical Point Crossing
   OldTestValue = TestValue;
   TestValue = Lat->TestFunctions(EigenValues);
   if ((OldTestValue != TestValue) && (BisectCP == Yes) && (OldTestValue != -1))
   {
      bfbreturncode = 2;
      SolveMe->FindCriticalPoint(Lat,TotalNumCPs,Input,Width,out);
   }
   else
   {
      bfbreturncode = 0;
   }
   // Send Output
   out << setw(Width) << *Lat << "Success = 1" << "\n";

   if (SolveMe->AllSolutionsFound())
   {
     out.close();
     delete SolveMe;
     delete Restrict;
     delete Lat;

     bfbreturncode=1;
   }
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

void InitializeOutputFile(fstream& out,char const* const outfile,char const* const datafile,
                          char const* const startfile,int const& Precision,int const& Width,
                          int const& Echo)
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
   out.open(outfile,ios::out);
   if (out.fail())
   {
      cerr << "Error: Unable to open file : " << outfile << " for write"
           << "\n";
      exit(-1);
   }
   
   while (!input.eof())
   {
      getline(input,dataline);
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
         getline(start,dataline);
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
