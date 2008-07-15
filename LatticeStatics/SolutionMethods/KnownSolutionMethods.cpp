#include "KnownSolutionMethods.h"

SolutionMethod* InitializeSolution(LatticeMode* const Mode,PerlInput const& Input,
                                   Lattice* const Lat,fstream& out,int const& Width,
                                   int const& Echo)
{
   enum solution {Scanning,ArcLen,NewtonPC};
   solution solu;
   
   const char *slvmthd = Input.getString("SolutionMethod","Type");
   if (!strcmp("ScanningSolution",slvmthd))
      solu = Scanning;
   else if (!strcmp("ArcLengthSolution",slvmthd))
      solu = ArcLen;
   else if (!strcmp("NewtonPCSolution",slvmthd))
      solu = NewtonPC;
   else
   {
      cerr << "Unknown SolutionMethod : " << slvmthd << "\n";
      exit(-1);
   }
   
   switch (solu)
   {
      case Scanning:
      {
         return new ScanningSolution(Mode,Input,Echo);
      }
      case ArcLen:
      {
         int good = 1;
         int count = 0;
         Vector One = Mode->ModeDOF(),
            Two = Mode->ModeDOF();

         if (Input.HashOK("StartType"))
         {
            return new ArcLengthSolution(Mode,Input,Echo);
         }
         else
         {
            ScanningSolution ScanMe(Mode,Input,Echo);

            while (!ScanMe.AllSolutionsFound())
            {
               One = Two;
               good = ScanMe.FindNextSolution();
               if (good)
               {
                  count++;
                  out << setw(Width) << *Lat << "Success = 1" << "\n";
                  Two = Mode->ModeDOF();
               }
            }
            
            if (count < 2)
            {
               cout << "Did not find two solutions with Scanning Solutions with "
                    << "which to initialize ArcLengthSolution." << "\n";
               exit(-55);
            }
            else
            {
               return new ArcLengthSolution(Mode,Input,One,Two,Echo);
            }
         }
      }
      case NewtonPC:
      {
         int good=1;
         int count=0;
         Vector One(Mode->ModeDOF().Dim());
         if (Input.HashOK("StartType"))
         {
            return new NewtonPCSolution(Mode,Input,Echo);
         }
         else
         {
            ScanningSolution ScanMe(Mode,Input,Echo);

            good = ScanMe.FindNextSolution();
            if (good)
            {
               count++;
               out << setw(Width) << *Lat << "Success = 1" << "\n";
               One = Mode->ModeDOF();
               return new NewtonPCSolution(Mode,Input,One,Echo);
            }
         }
      }
   }
   
   return 0;
}
