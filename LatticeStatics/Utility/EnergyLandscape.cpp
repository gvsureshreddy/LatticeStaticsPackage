#include "KnownLattices.h"
#include "KnownModes.h"
#include "PerlInput.h"
#include <fstream>

char *builddate();

using namespace std;

void GetMainSettings(int& Width,int& Presision,int& Echo,PerlInput const& Input);
void InitializeOutputFile(fstream& out,char const* const outfile,char const* const datafile,
                          int const& argc,Lattice const* const Lat,int const& Precision,
                          int const& Width,int const& Echo);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
           << " ParamFile OutputFile" << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n";
      exit(-1);
   }
   
   char const* const datafile = argv[1];
   char const* const outputfile = argv[2];

   PerlInput Input(datafile);
   
   LatticeMode *Mode;
   
   int Width,Precision,Echo;
   
   GetMainSettings(Width,Precision,Echo,Input);
   
   Lattice* const Lat = InitializeLattice(Input,Echo);
   
   fstream out;
   InitializeOutputFile(out,outputfile,datafile,argc,Lat,Precision,Width,Echo);
   
   Mode = InitializeMode(Lat,Input);
   
   int NoDims;
   NoDims = Input.getPosInt("EnergyLandscape","Directions");
   char tmp[LINELENGTH];
   int* const divs = new int[NoDims+1];
   Input.getIntVector(divs,NoDims,"EnergyLandscape","DirectionDivisions");
   divs[NoDims]=1;
   
   Vector* const Corners = new Vector[NoDims+1];
   
   Corners[0].Resize(Mode->ModeDOF().Dim());
   Input.getVector(Corners[0],"EnergyLandscape","Origin");
   for (int i=1;i<=NoDims;++i)
   {
      Corners[i].Resize(Mode->ModeDOF().Dim());
      sprintf(tmp,"Corner_%u",i);
      Input.getVector(Corners[i],"EnergyLandscape",tmp);
   }
   Input.EndofInputSection();
   
   Vector* const Directions = new Vector[NoDims];
   for (int i=0;i<NoDims;++i)
   {
      Directions[i].Resize(Corners[0].Dim());
      Directions[i] = (Corners[i+1] - Corners[0])/(divs[i]);
   }
   
   int* const counter = new int[NoDims+1];
   for (int i=0;i<=NoDims;++i) counter[i]=0;
   
   Vector state(Mode->ModeDOF().Dim());
   int i;
   do
   {
      state=Corners[0];
      for (i=0;i<NoDims;++i) state += counter[i]*Directions[i];
      //cout << setw(Width) << state << "\n";
      Mode->SetModeDOF(state);
      
      for (i=0;i<NoDims;++i) out << setw(Width) << counter[i];
      out << setw(Width) << Mode->ModeEnergy() << "\n";
      
      ++(counter[0]);
      i=0;
      while (counter[i]>divs[i])
      {
         if (i<2) out << "\n";
         counter[i]=0;
         if (i < NoDims) ++(counter[++i]);
      }
   }
   while(counter[NoDims]==0);
   
   delete [] Corners;
   delete [] divs;
   delete [] Directions;
   delete [] counter;
   
   delete Mode;
   delete Lat;
   
   out.close();
   return 1;
}

void GetMainSettings(int& Width,int& Precision,int& Echo,PerlInput const& Input)
{
   Width = Input.getPosInt("Main","FieldWidth");
   Precision = Input.getPosInt("Main","Precision");
   if (Input.ParameterOK("Main","Echo"))
   {
      Echo = Input.getPosInt("Main","Echo");
   }
   else
   {
      Echo = Input.usePosInt(1,"Main","Echo"); // Default Value
   }
   Input.EndofInputSection();
}

void InitializeOutputFile(fstream& out,char const* const outfile,char const* const datafile,
                          int const& argc,Lattice const* const Lat,int const& Precision,
                          int const& Width,int const& Echo)
{
   fstream input;
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
   
   if (argc >= 4)
   {
      while (!input.eof())
      {
         input.getline(dataline,LINELENGTH-1);
         out << "# Input File:" << dataline << "\n";
      }
   }
   else
   {
      input.getline(dataline,LINELENGTH-1);
      while ((strstr(dataline,"Input File:") != 0) ||
             (strstr(dataline,"Start File:") != 0))
      {
         out << "# " << dataline << "\n";
         input.getline(dataline,LINELENGTH-1);
      }
   }
   
   input.close();
   
   cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);
   
   cout << "Built on:               " << builddate() << "\n"
        << "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
        << "MyMath Built on:        " << MyMathBuildDate() << "\n";
   out << "# Built on:               " << builddate() << "\n"
       << "# LinearAlgebra Build on: " << LinearAlgebraBuildDate() << "\n"
       << "# MyMath Built on:        " << MyMathBuildDate() << "\n";
}
