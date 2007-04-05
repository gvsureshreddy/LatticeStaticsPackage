#include "KnownLattices.h"
#include "KnownModes.h"
#include "UtilityFunctions.h"
#include <fstream>

using namespace std;

void GetMainSettings(int &Width,int &Presision,int &Echo,char *datafile,
		     const char *prefix);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,const char *prefix,
			  Lattice *Lat,int Precision,int Width,int Echo);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile OutputFile" << endl;
      cerr << "Built on:               " << builddate() << endl
	   << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
	   << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }

   char *datafile = argv[1],
      *outputfile = argv[2],
      prefix[LINELENGTH];

   if (GetParameter("","^Input File:",datafile,"%s",prefix,0))
   {
      strcpy(prefix,"^Input File:");
   }
   else
   {
      strcpy(prefix,"^");
   }
   
   Lattice *Lat;
   LatticeMode *Mode;

   int Width,Precision,Echo;

   GetMainSettings(Width,Precision,Echo,datafile,prefix);

   Lat = InitializeLattice(datafile,prefix,Echo);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,prefix,Lat,Precision,Width,Echo);

   Mode = InitializeMode(Lat,datafile,prefix);

   int NoDims;
   if(!GetParameter(prefix,"Directions",datafile,"%u",&NoDims)) exit(-1);
   char tmp[LINELENGTH];
   int *divs;
   divs = new int[NoDims+1];
   for (int i=0;i<NoDims;++i)
   {
      sprintf(tmp,"DirectionDivisions_%u",i);
      if(!GetParameter(prefix,tmp,datafile,"%u",&(divs[i]))) exit(-1);
   }
   divs[NoDims]=1;
	 
   Vector *Corners;
   Corners = new Vector[NoDims+1];

   Corners[0].Resize(Mode->ModeDOF().Dim());
   if(!GetVectorParameter(prefix,"Origin",datafile,&(Corners[0]))) exit(-1);
   for (int i=1;i<=NoDims;++i)
   {
      Corners[i].Resize(Mode->ModeDOF().Dim());
      sprintf(tmp,"Corner_%u",i);
      if(!GetVectorParameter(prefix,tmp,datafile,&(Corners[i]))) exit(-1);
   }

   Vector *Directions;
   Directions = new Vector[NoDims];
   for (int i=0;i<NoDims;++i)
   {
      Directions[i].Resize(Corners[0].Dim());
      Directions[i] = (Corners[i+1] - Corners[0])/(divs[i]);
   }

   int *counter;
   counter = new int[NoDims+1];
   for (int i=0;i<=NoDims;++i) counter[i]=0;

   Vector state(Mode->ModeDOF().Dim());
   int i;
   do
   {
      state=Corners[0];
      for (i=0;i<NoDims;++i) state += counter[i]*Directions[i];
      //cout << setw(Width) << state << endl;
      Mode->SetModeDOF(state);

      for (i=0;i<NoDims;++i) out << setw(Width) << counter[i];
      out << setw(Width) << Mode->ModeEnergy() << endl;

      ++(counter[0]);
      i=0;
      while (counter[i]>divs[i])
      {
	 if (i<2) out << endl;
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




void GetMainSettings(int &Width,int &Precision,int &Echo,char *datafile,
		     const char *prefix)
{
   if(!GetParameter(prefix,"MainFieldWidth",datafile,"%d",&Width)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,"%d",&Precision)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,"%d",&Echo,0)) Echo=1;   
}

void InitializeOutputFile(fstream &out,char *outfile,char *datafile,const char *prefix,
			  Lattice *Lat,int Precision,int Width,int Echo)
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

   if (!strcmp("^",prefix))
   {
      while (!input.eof())
      {
	 input.getline(dataline,LINELENGTH-1);
	 out << "# Input File:" << dataline << endl;
      }
   }
   else
   {
      input.getline(dataline,LINELENGTH-1);
      while ((strstr(dataline,"Input File:") != NULL) ||
	     (strstr(dataline,"Start File:") != NULL))
      {
	 out << "# " << dataline << endl;
	 input.getline(dataline,LINELENGTH-1);
      }
   }

   input.close();

   cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);

   cout << "Built on:               " << builddate() << endl
	<< "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
	<< "MyMath Built on:        " << MyMathBuildDate() << endl;
   out << "# Built on:               " << builddate() << endl
       << "# LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
       << "# MyMath Built on:        " << MyMathBuildDate() << endl;
}

