#include "KnownLattices.h"
#include "UtilityFunctions.h"
#include <fstream.h>

void GetMainSettings(int &Width,int &Presision,char *datafile,const char *prefix);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,const char *prefix,
			  Lattice *Lat,int Precision,int Width);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
	   << " ParamFile OutputFile <PathFile>" << endl;
      cerr << "Built on:               " << builddate() << endl
	   << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << endl
	   << "MyMath Built on:        " << MyMathBuildDate() << endl;
      exit(-1);
   }

   char *datafile = argv[1],
      *outputfile = argv[2],
      prefix[]="^";

   GenericLat *Lat;

   int Width,Precision;

   GetMainSettings(Width,Precision,datafile,prefix);

   Lat = (GenericLat*) InitializeLattice(datafile,prefix);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,prefix,Lat,Precision,Width);

   if (argc<4)
   {
      Lat->DispersionCurves(datafile,prefix,out);
   }
   else
   {
      FILE *pipe;
      char format[]=
      {"perl -e '$_=<>;while(! m/^Mode:/){$_=<>;}while(<>){if(/^Temperature/){"\
       "@fld=split(/:/,$_);print $fld[1];}if(/^DOF/){$_=<>;print $_;}}print"\
       "\"DONE\\n\";' %s"};
      
      char strng[LINELENGTH];
      sprintf(strng,format,argv[3]);
      pipe = popen(strng,"r");

      double temp;
      Vector DOF((Lat->DOF()).Dim());
      fscanf(pipe,"%s",strng);
      while (strcmp("DONE",strng))
      {
	 temp = atof(strng);
	 for (int i=0;i<DOF.Dim();++i)
	 {
	    fscanf(pipe,"%lf",&(DOF[i]));
	 }

	 Lat->SetTemp(temp);
	 Lat->SetDOF(DOF);

	 out << setw(0);
	 out << "#" << setw(Width) << temp << endl
	     << "#" << setw(Width) << DOF << endl << setw(Width);
	 cout << setw(0);
	 cout << "#" << setw(Width) << temp << endl
	      << "#" << setw(Width) << DOF << endl;;
	 
	 Lat->DispersionCurves(datafile,prefix,out);

	 fscanf(pipe,"%s",strng);
      }
      pclose(pipe);
   }

   out.close();
   return 1;
}




void GetMainSettings(int &Width, int &Precision,char *datafile,const char *prefix)
{
   if(!GetParameter(prefix,"MainFieldWidth",datafile,"%d",&Width)) exit(-1);
   if(!GetParameter(prefix,"MainPrecision",datafile,"%d",&Precision)) exit(-1);   
}

void InitializeOutputFile(fstream &out,char *outfile,char *datafile,const char *prefix,
			  Lattice *Lat,int Precision,int Width)
{
   fstream input,devnull;
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
   devnull.open("/dev/null",ios::out);

   if (!strcmp("^",prefix))
   {
      while (!input.eof())
      {
	 input.getline(dataline,LINELENGTH-1);
	 out << "#Input File:" << dataline << endl;
      }
   }
   else
   {
      input.getline(dataline,LINELENGTH-1);
      while ((strstr(dataline,"Input File:") != NULL) ||
	     (strstr(dataline,"Start File:") != NULL))
      {
	 out << "#" << dataline << endl;
      }
   }

   input.close();

   cout << setiosflags(ios::fixed) << setprecision(Precision);
   out  << setiosflags(ios::fixed) << setprecision(Precision);

   cout << "Built on:               " << builddate() << endl
	<< "LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
	<< "MyMath Built on:        " << MyMathBuildDate() << endl
	<< setw(Width);
   out << "#Built on:               " << builddate() << endl
       << "#LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
       << "#MyMath Built on:        " << MyMathBuildDate() << endl
       << setw(Width);

   devnull << setw(Width);
   Lat->Print(devnull,Lattice::PrintLong);

   devnull.close();
}
