#include "KnownSolutionMethods.h"

#include "UtilityFunctions.h"

SolutionMethod *InitializeSolution(LatticeMode *Mode,char *datafile,
				   char *startfile,Lattice *Lat,fstream &out,int Width,
				   int Echo)
{
   char slvmthd[LINELENGTH];
   
   enum solution {Scanning,ArcLen,NewtonPC,NewtonUpdatePC};
   solution solu;

   if(!GetParameter("^","MainSolutionMethod",datafile,'s',slvmthd)) exit(-1);
   if ((!strcmp("Scanning",slvmthd))
       || (!strcmp("scanning",slvmthd)))
      solu = Scanning;
   else if ((!strcmp("ArcLength",slvmthd)) || (!strcmp("arclength",slvmthd)))
      solu = ArcLen;
   else if ((!strcmp("NewtonPC",slvmthd)) || (!strcmp("newtonpc",slvmthd)))
      solu = NewtonPC;
   else if ((!strcmp("NewtonUpdatePC",slvmthd)) || (!strcmp("newtonupdatepc",slvmthd)))
      solu = NewtonUpdatePC;  
   else
   {
      cerr << "Unknown SolutionMethod : " << slvmthd << endl;
      exit(-1);
   }

   switch (solu)
   {
      case Scanning:
      {
	 return new ScanningSolution(Mode,datafile,"^",Echo);
      }
      case ArcLen:
      {
	 int good = 1;
	 int count = 0;
	 Vector One = Mode->ArcLenDef(),
	    Two = Mode->ArcLenDef();

	 if ( startfile == NULL)
	 {
	    ScanningSolution ScanMe(Mode,datafile,"^",Echo);
	    
	    while (!ScanMe.AllSolutionsFound())
	    {
	       One = Two;
	       ScanMe.FindNextSolution(good);
	       if (good)
	       {
		  count++;
		  out << setw(Width) << Lat << "Success = 1" << endl;
		  Two = Mode->ArcLenDef();
	       }
	    }
	    
	    if (count < 2)
	    {
	       cout << "Did not find two solutions with Scanning Solutions with "
		    << "which to initialize ArcLengthSolution." << endl;
	       exit(-55);
	    }
	    else
	    {
	       return new ArcLengthSolution(Mode,datafile,"^",One,Two,Echo);
	    }
	 }
	 else
	 {
	    return new ArcLengthSolution(Mode,datafile,"^",startfile,out,Echo);
	 }
      }
      case NewtonPC:
      {
	 int good=1;
	 int count=0;
	 Vector One(Mode->ModeDOF().Dim());
	 if ( startfile == NULL)
	 {
	    ScanningSolution ScanMe(Mode,datafile,"^",Echo);
	    
	    ScanMe.FindNextSolution(good);
	    if (good)
	    {
	       count++;
	       out << setw(Width) << Lat << "Success = 1" << endl;
	       One = Mode->ArcLenDef();
	       return new NewtonPCSolution(Mode,datafile,"^", One,Echo);
	    }
	 }
	 else
	 {
	    return new NewtonPCSolution(Mode,datafile,"^",startfile,out,Echo);
	 }
      }
      case NewtonUpdatePC:
      {
	 int good=1;
	 int count=0;
	 Vector One(Mode->ModeDOF().Dim());
	 if ( startfile == NULL)
	 {
	    ScanningSolution ScanMe(Mode,datafile,"^",Echo);
	    
	    ScanMe.FindNextSolution(good);
	    if (good)
	    {
	       count++;
	       out << setw(Width) << Lat << "Success = 1" << endl;
	       One = Mode->ArcLenDef();
	       return new NewtonUpdatePCSolution(Mode,datafile,"^", One,Echo);
	    }
	 }
	 else
	 {
	    return new NewtonUpdatePCSolution(Mode,datafile,"^",startfile,out,Echo);
	 }
      }
   }
   
   return NULL;
}
