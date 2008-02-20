#ifndef __NewtonPCSolution
#define __NewtonPCSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

class Lattice;

class NewtonPCSolution : public SolutionMethod
{
private:
   LatticeMode *Mode_;
   
   int CurrentSolution_;
   int Echo_;
   int NumSolutions_;
   
   double MaxDS_;
   double CurrentDS_;		//initial Steplength h > 0
   double cont_rate_nom_;	//Nominal contraction rate
   double delta_nom_;		//Nominal distance to (from predicted to corrected point) curve
   double alpha_nom_;		//Nominal angle to curve
   double Converge_;		//Convergence criteria
   double MinDSRatio_;		// Minimum Stepsize ratio
   int ClosedLoopStart_;	//Closed loop test variable
   int Direction_;			//Direction of tangent
   
   Vector FirstSolution_;	//Initial point on curve
   Vector Tangent1_;		//Tangent vector of ith point
   Vector Tangent2_;		//Tangent Vector of ith + 1 point
   
   void MoorePenrose(const Matrix &Q, const Matrix &R, const Vector &Force, Vector &Corrector);
   
public:
   Vector Previous_Solution_;
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,const Vector &one,
		    int Echo=1,  int Direction=0);
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,char *startfile,
		    fstream &out,int Echo);
   ~NewtonPCSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int FindCriticalPoint(Lattice *Lat, char *datafile,const char *prefix,int Width,
                                 fstream &out);
};

#endif
