#ifndef __ArcLengthSolution
#define __ArcLengthSolution

#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

class Lattice;

class ArcLengthSolution : public SolutionMethod
{
private:
   int Echo_;
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
   int ClosedLoopStart_;
   Vector FirstSolution_;
   
   Vector Difference_;

   // Consistency Check data
   double ConsistencyEpsilon_;

   double ArcLengthNewton(int &good);
   int ConsistencyCheck();

public:
   ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
		     const Vector &one,const Vector &two,int Echo=1);
   ArcLengthSolution(LatticeMode *Mode,char *datafile,const char *prefix,
		     char *startfile,fstream &out,int Echo=1);
   ~ArcLengthSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
			   int Width,fstream &out);
   
};

#endif
