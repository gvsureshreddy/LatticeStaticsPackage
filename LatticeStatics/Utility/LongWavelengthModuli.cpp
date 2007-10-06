#include "KnownLattices.h"
#include "UtilityFunctions.h"
#include <fstream>

using namespace std;

void GetMainSettings(int &Width,int &Presision,int &Echo,char *datafile,
		     const char *prefix);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,const char *prefix,
			  Lattice *Lat,int Precision,int Width,int Echo);

void RefineEquilibrium(Lattice *Lat,double Tol,int Width,int Echo);

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

   if (GetParameter("","^Input File:",datafile,'s',prefix,0))
   {
      strcpy(prefix,"^Input File:");
   }
   else
   {
      strcpy(prefix,"^");
   }
   
   Lattice *Lat;

   int Width,Precision,Echo;

   GetMainSettings(Width,Precision,Echo,datafile,prefix);

   Lat = InitializeLattice(datafile,prefix,Echo);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,prefix,Lat,Precision,Width,Echo);

   int gridsize;
   double dk;
   if(!GetParameter(prefix,"LongWavelengthGridSize",datafile,'u',&gridsize)) exit(-1);
   if(!GetParameter(prefix,"LongWavelengthDK",datafile,'l',&dk)) exit(-1);

   if (strcmp(prefix,"^Input File:"))
   {
      Lat->LongWavelengthModuli(dk,gridsize,"",out);
      out << endl;
      if (Echo) cout << endl;
   }
   else
   {
      FILE *pipe;
      char format[]=
      {"perl -e '$_=<>;while(! m/^Mode:/){$_=<>;}while(<>){if(/^Temperature/){"\
       "@fld=split(/:/,$_);print $fld[1];}if(/^DOF/){$_=<>;print $_;}}print"\
       "\"DONE\\n\";' %s"};
      
      char strng[LINELENGTH];
      char tmp[LINELENGTH];
      sprintf(strng,format,datafile);

      pipe = popen(strng,"r");
      
      double temp;
      Vector DOF((Lat->DOF()).Dim());
      fscanf(pipe,"%s",tmp);
      while (strcmp("DONE",tmp))
      {
	 temp = atof(tmp);
	 for (int j=0;j<DOF.Dim();++j)
	 {
	    fscanf(pipe,"%lf",&(DOF[j]));
	 }
	 
	 Lat->SetTemp(temp);
	 Lat->SetDOF(DOF);

	 out << "#" << setw(Width) << temp << endl
	     << "#" << setw(Width) << DOF << endl << setw(Width);
	 if (Echo) cout << "#" << setw(Width) << temp << endl
			<< "#" << setw(Width) << DOF << endl;

	 Lat->LongWavelengthModuli(dk,gridsize,"",out);
	 out << endl;
	 if (Echo) cout << endl;
	 
	 fscanf(pipe,"%s",tmp);
      }
      pclose(pipe);
      
      out << endl;
      if (Echo) cout << endl;
   }
      
   out.close();
   return 1;
}




void GetMainSettings(int &Width,int &Precision,int &Echo,char *datafile,
		     const char *prefix)
{
   if(!GetParameter(prefix,"MainFieldWidth",datafile,'i',&Width)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,'i',&Precision)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,'i',&Echo,0)) Echo=1;   
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

void RefineEquilibrium(Lattice *Lat,double Tol,int Width,int Echo)
{
   Vector S = -Lat->E1();
   Matrix K = Lat->E2();
   Vector Du(S.Dim(),0.0);
   Vector u = Lat->DOF();

   while (S.Norm() > Tol)
   {
      Du = SolveSVD(K,S,MAXCONDITION,Echo);
      u += Du;
      
      Lat->SetDOF(u);
      S = -Lat->E1();
      K = Lat->E2();
      if (Echo)
      {
	 cout << setw(Width) << S.Norm() << endl;
      }
   }
}
