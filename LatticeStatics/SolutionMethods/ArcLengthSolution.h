#ifndef __ArcLengthSolution
#define __ArcLengthSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

class Lattice;

class ArcLengthSolution : public SolutionMethod
{
private:
   LatticeMode *Mode_;
   unsigned MaxIter_;
   double Tolerance_;

   double DSMax_;
   double DSMin_;
   double CurrentDS_;
   double AngleCutoff_;
   double AngleIncrease_;
   double Aspect_;

   unsigned NumSolutions_;
   unsigned CurrentSolution_;

   Vector Difference_;

   // Consistency Check data
   double ConsistencyEpsilon_;

   double ArcLengthNewton(int &good);
   int ConsistencyCheck();

public:
   ArcLengthSolution(LatticeMode *Mode,char *datafile,
		     const Vector &one,const Vector &two);
   ArcLengthSolution(LatticeMode *Mode,char *datafile,char *startfile,fstream &out);
   ~ArcLengthSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int BisectAlert(Lattice *Lat,char *datafile,int Width,fstream &out);
   
};

#endif
