#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "UtilityFunctions.h"
#include <cmath>

using namespace std;


NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   const Vector &one,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   // get needed parameters
   
   // Set Lattice to solution "one"
   Mode_->SetModeDOF(one);
}

int NewtonPCSolution::AllSolutionsFound()
{
   if (CurrentSolution_ < NumSolutions_)
   { 
      return 0;
   }
   else
   {
      return 1;
   }
   
}

NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   char *startfile,fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{

   // get needed parameters
   

   const char *StartType[] = {"Bifurcation","Continuation","ConsistencyCheck"};
   switch (GetStringParameter(prefix,"StartType",startfile,StartType,3))
   {
      case -1:
      {
         cerr << "Unknown StartType!" << endl;
         exit(-1);
         break;
      }
      case 0:
      {
	 // Bifurcation

	 // read in bifurcation point and tangent then proceed as in case 1
         break;
      }
      case 1:
      {
	 // Continuation
	 
         // Get solution1
         Vector one(Mode_->ArcLenDef().Dim());
         if(!GetVectorParameter(prefix,"Solution1",startfile,&one)) exit(-1);
         Mode_->SetModeDOF(one);

         break;
      }
      case 2:
      {
	 // ConsistencyCheck

	 // do nothing for now
         break;
      }
   }
}


double NewtonPCSolution::FindNextSolution(int &good)
{
   //Finds the next solution

   //perform predictor and corrector steps with adaptive step length
}


int NewtonPCSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,int Width,
				  fstream &out)
{
   // for the moment, do nothing
   return 1;
}
