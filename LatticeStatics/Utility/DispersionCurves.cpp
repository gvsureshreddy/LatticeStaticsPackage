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

   int NoDirs,NoPTS;
   if(!GetParameter(prefix,"DispersionDirections",datafile,"%u",&NoDirs)) exit(-1);
   if(!GetParameter(prefix,"DispersionPoints",datafile,"%u",&NoPTS)) exit(-1);

   Vector *Direction;
   Direction = new Vector[NoDirs];
   char tmp[LINELENGTH];

   for (int i=0;i<NoDirs;++i)
   {
      Direction[i].Resize(3);
      sprintf(tmp,"DispersionDir_%u",i);
      if(!GetVectorParameter(prefix,tmp,datafile,&(Direction[i]))) exit(-1);
   }

   if (argc<4)
   {
      for (int i=0;i<NoDirs;++i)
      {
	 out << "#" << setw(Width) << Direction[i] << endl << setw(Width);
	 cout << "#" << setw(Width) << Direction[i] << endl << setw(Width);
	 Lat->DispersionCurves(Direction[i],NoPTS,"",out);
	 out << endl << endl;
	 cout << endl << endl;
      }
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
      sprintf(strng,format,argv[3]);

      for (int i=0;i<NoDirs;++i)
      {
	 out << "#" << setw(Width) << Direction[i] << endl;
	 cout << "#" << setw(Width) << Direction[i] << endl;
	 
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
	    cout << "#" << setw(Width) << temp << endl
		 << "#" << setw(Width) << DOF << endl;;
	    
	    Lat->DispersionCurves(Direction[i],NoPTS,"",out);
	    out << endl;
	    cout << endl;
	    
	    fscanf(pipe,"%s",tmp);
	 }
	 pclose(pipe);
	 
	 out << endl;
	 cout << endl;
      }
   }  
      
   delete [] Direction;
   
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
	<< "MyMath Built on:        " << MyMathBuildDate() << endl;
   out << "#Built on:               " << builddate() << endl
       << "#LinearAlgebra Build on: " << LinearAlgebraBuildDate() << endl
       << "#MyMath Built on:        " << MyMathBuildDate() << endl;

   devnull << setw(Width);
   Lat->Print(devnull,Lattice::PrintLong);

   devnull.close();
}

