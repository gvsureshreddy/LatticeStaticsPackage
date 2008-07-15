#include "KnownLattices.h"
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
           << " ParamFile OutputFile [1 - regular inputfile (not an output of LatticeStatics)]"
           << "\n";
      cerr << "Built on:               " << builddate() << "\n"
           << "LinearAlgebra Built on: " << LinearAlgebraBuildDate() << "\n"
           << "MyMath Built on:        " << MyMathBuildDate() << "\n";
      exit(-1);
   }
   
   char const* const datafile = argv[1];
   char const* const outputfile = argv[2];

   PerlInput Input;
   if (argc >= 4)
   {
      Input.Readfile(datafile);
   }
   else
   {
      Input.Readfile(datafile,"Input File:");
   }

   int Width,Precision,Echo;
   
   GetMainSettings(Width,Precision,Echo,Input);

   Lattice* const Lat = InitializeLattice(Input,Echo);

   fstream out;
   InitializeOutputFile(out,outputfile,datafile,argc,Lat,Precision,Width,Echo);
   
   int NoLines,NoPTS;
   NoLines = Input.getPosInt("DispersionCurves","NoLines");
   NoPTS = Input.getPosInt("DispersionCurves","Points");
   
   Vector* const Line = new Vector[NoLines];

   for (int i=0;i<NoLines;++i)
   {
      Line[i].Resize(6);
      Input.getVector(Line[i],"DispersionCurves","Lines",i);
   }
   Input.EndofInputSection();

   if (argc >= 4)
   {
      for (int i=0;i<NoLines;++i)
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
      FILE* pipe;
      char format[]=
         {"perl -e '$_=<>;while(! m/^Mode:/){$_=<>;}while(<>){if(/^Temperature/){"\
          "@fld=split(/:/,$_);print $fld[1];}if(/^Lambda/){@fld=split(/:/,$_);print $fld[1];}"\
          "if(/^DOF/){$_=<>;print $_;}}print \"DONE\\n\";' %s"};
      
      char strng[LINELENGTH];
      char tmp[LINELENGTH];
      sprintf(strng,format,datafile);
      
      for (int i=0;i<NoLines;++i)
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
            for (int j=0;j<DOF.Dim();++j)
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

