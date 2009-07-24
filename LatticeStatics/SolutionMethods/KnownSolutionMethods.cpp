#include "KnownSolutionMethods.h"

SolutionMethod* InitializeSolution(Restriction* const Restrict,PerlInput const& Input,
                                   Lattice* const Lat,ostream& out,int const& Width,
                                   int const& Echo)
{
   enum solution {RefineEqbm,Scanning,ArcLen,NewtonPC};
   solution solu;
   
   const char *slvmthd = Input.getString("SolutionMethod","Type");
   if (!strcmp("RefineEqbmSolution",slvmthd))
      solu = RefineEqbm;
   else if (!strcmp("ScanningSolution",slvmthd))
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
      case RefineEqbm:
      {
         return new RefineEqbmSolution(Restrict,Input,0);
      }
      case Scanning:
      {
         return new ScanningSolution(Restrict,Input,Echo);
      }
      case ArcLen:
      {
         int good = 1;
         int count = 0;
         Vector One = Restrict->DOF(),
            Two = Restrict->DOF();

         if (Input.HashOK("StartType"))
         {
            return new ArcLengthSolution(Restrict,Input,Echo);
         }
         else
         {
            ScanningSolution ScanMe(Restrict,Input,Echo);

            while (!ScanMe.AllSolutionsFound())
            {
               One = Two;
               good = ScanMe.FindNextSolution();
               if (good)
               {
                  count++;
                  out << setw(Width) << *Lat << "Success = 1" << "\n";
                  Two = Restrict->DOF();
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
               return new ArcLengthSolution(Restrict,Input,One,Two,Echo);
            }
         }
      }
      case NewtonPC:
      {
         int good=1;
         int count=0;
         Vector One(Restrict->DOF().Dim());
         if (Input.HashOK("StartType"))
         {
            return new NewtonPCSolution(Restrict,Input,Echo);
         }
         else
         {
            // find and setup first solution point method
            SolutionMethod* slnmthd;
            PerlInput::HashStruct Hash = Input.getHash("SolutionMethod","NewtonPCSolution");
            if (Input.ParameterOK(Hash,"FirstPointMethod"))
            {
               const char *frstpt = Input.getString(Hash,"FirstPointMethod");
               if (!strcmp("RefineEqbmSolution",frstpt))
               {
                  slnmthd = new RefineEqbmSolution(Restrict,Input,0);
               }
               else if (!strcmp("ScanningSolution",frstpt))
               {
                  slnmthd = new ScanningSolution(Restrict,Input,Echo);
               }
               else
               {
                  cerr << "Unknown FirstPointMehtod: " << frstpt << "\nExiting!\n";
                  exit(-21);
               }
            }
            else
            {
               // default to RefineEqbmSolution
               slnmthd = new RefineEqbmSolution(Restrict,Input,0);
            }

            // Find first solution point
            while (!slnmthd->AllSolutionsFound())
            {
               good = slnmthd->FindNextSolution();
               if (good)
               {
                  count++;
                  out << setw(Width) << *Lat << "Success = 1" << "\n";
               }
            }
            // teardown slnmthd
            delete slnmthd;
            
            One = Restrict->DOF();
            return new NewtonPCSolution(Restrict,Input,One,Echo);
         }
      }
   }
   
   return 0;
}
