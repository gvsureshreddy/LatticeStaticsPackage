#include "KnownLattices.h"
#include "PerlInput.h"
#include <fstream>

char *builddate();

using namespace std;

void GetMainSettings(int &Width,int &Presision,int &Echo,PerlInput &Input);
void InitializeOutputFile(fstream &out,char *outfile,char *datafile,int argc,
                          Lattice *Lat,int Precision,int Width,int Echo);

int main(int argc,char *argv[])
{
   // Check commandline args
   if (argc < 3)
   {
      cerr << "Usage: " << argv[0]
           << " ParamFile OutputFile [1 - regular inputfile (not an output of LatticeStatics)]"
           << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n";
      exit(-1);
   }
   
   char *datafile = argv[1],
      *outputfile = argv[2];

   PerlInput Input;
   if (argc >= 4)
   {
      Input.Readfile(datafile);
   }
   else
   {
      Input.Readfile(datafile,"Input File:");
   }

   Lattice *Lat;
   
   int Width,Precision,Echo;
   
   GetMainSettings(Width,Precision,Echo,Input);

   Lat = InitializeLattice(Input,Echo);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,argc,Lat,Precision,Width,Echo);
   
   unsigned NoLines,NoPTS;
   NoLines = Input.getUnsigned("DispersionCurves","NoLines");
   NoPTS = Input.getUnsigned("DispersionCurves","Points");
   
   Vector *Line;
   Line = new Vector[NoLines];

   for (unsigned i=0;i<NoLines;++i)
   {
      Line[i].Resize(6);
      Input.getVector(Line[i],"DispersionCurves","Lines",i);
   }

   if (argc >= 4)
   {
      for (unsigned i=0;i<NoLines;++i)
      {
         out << "#" << setw(Width) << Line[i] << "\n" << setw(Width);
         if (Echo) cout << "#" << setw(Width) << Line[i] << "\n" << setw(Width);
         Lat->DispersionCurves(Line[i],NoPTS,"",out);
         out << "\n" << "\n";
         if (Echo) cout << "\n" << "\n";
      }
   }
   else
   {
      FILE *pipe;
      char format[]=
         {"perl -e '$_=<>;while(! m/^Mode:/){$_=<>;}while(<>){if(/^Temperature/){"\
          "@fld=split(/:/,$_);print $fld[1];}if(/^Lambda/){@fld=split(/:/,$_);print $fld[1];}"\
          "if(/^DOF/){$_=<>;print $_;}}print \"DONE\\n\";' %s"};
      
      char strng[LINELENGTH];
      char tmp[LINELENGTH];
      sprintf(strng,format,datafile);
      
      for (unsigned i=0;i<NoLines;++i)
      {
         out << "#" << setw(Width) << Line[i] << "\n";
         if (Echo) cout << "#" << setw(Width) << Line[i] << "\n";
         
         pipe = popen(strng,"r");
         
         double temp,lambda;
         Vector DOF((Lat->DOF()).Dim());
         fscanf(pipe,"%s",tmp);
         while (strcmp("DONE",tmp))
         {
            temp = atof(tmp);
            fscanf(pipe,"%lf",&lambda);
            for (unsigned j=0;j<DOF.Dim();++j)
            {
               fscanf(pipe,"%lf",&(DOF[j]));
            }
            
            Lat->SetTemp(temp);
            Lat->SetLambda(lambda);
            Lat->SetDOF(DOF);
            
            out << "# Temp= " << setw(Width) << temp << " Lambda= "
                << setw(Width) << lambda << "\n"
                << "#" << setw(Width) << DOF << "\n" << setw(Width);
            if (Echo) cout << "# Temp= " << setw(Width) << temp << " Lambda= "
                           << setw(Width) << lambda << "\n"
                           << "#" << setw(Width) << DOF << "\n";;
            
            Lat->DispersionCurves(Line[i],NoPTS,"",out);
            out << "\n" << "\n";
            if (Echo) cout << "\n" << "\n";
            
            fscanf(pipe,"%s",tmp);
         }
         pclose(pipe);
         
         out << "\n";
         if (Echo) cout << "\n";
      }
   }
   
   delete [] Line;
   
   delete Lat;
   
   out.close();
   return 1;
}




void GetMainSettings(int &Width,int &Precision,int &Echo,PerlInput &Input)
{
   Width = Input.getInt("Main","FieldWidth");
   Precision = Input.getInt("Main","Precision");
   if (Input.ParameterOK("Main","Echo"))
   {
      Echo = Input.getInt("Main","Echo");
   }
   else
   {
      Echo = 1;
   }
}

void InitializeOutputFile(fstream &out,char *outfile,char *datafile,int argc,
                          Lattice *Lat,int Precision,int Width,int Echo)
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
      while ((strstr(dataline,"Input File:") != NULL) ||
             (strstr(dataline,"Start File:") != NULL))
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

