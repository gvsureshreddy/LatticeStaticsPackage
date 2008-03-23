#ifndef __NewtonQRUpdatePCSolution
#define __NewtonQRUpdatePCSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

#define CLOSEDDEFAULT 30

class Lattice;

class NewtonQRUpdatePCSolution : public SolutionMethod
{
private:
   LatticeMode *Mode_;
   
   int Echo_;
   unsigned CurrentSolution_;
   unsigned NumSolutions_;
   
   double MaxDS_;
   double CurrentDS_;           //initial Steplength h > 0
   double cont_rate_nom_;       //Nominal contraction rate
   double delta_nom_;           //Nominal distance to (from predicted to corrected point) curve
   double alpha_nom_;           //Nominal angle to curve
   double Converge_;            //Convergence criteria
   double MinDSRatio_;          // Minimum Stepsize ratio
   unsigned ClosedLoopStart_;   //Closed loop test variable
   int Direction_;              //Direction of tangent
   
   Vector FirstSolution_;       //Initial point on curve
   Vector Tangent1_;            //Tangent vector of ith point
   Vector Tangent2_;            //Tangent Vector of ith + 1 point
   
   void MoorePenrose(const Matrix& Q,const Matrix& R,const Vector& Force,Vector& Corrector);
   void QRUpdate(const Vector& Force,const Vector& difference,Matrix& Q,Matrix& R);
   
public:
   Vector Previous_Solution_;
   NewtonQRUpdatePCSolution(LatticeMode *Mode,const Vector &one,
                            unsigned CurrentSolution,unsigned NumSolutions,double MaxDS,
                            double CurrentDS,double cont_rate_nom,double delta_nom,
                            double alpha_nom,double Converge,double MinDSRatio,
                            const Vector &FirstSolution,int Direction=1,
                            unsigned ClosedLoopStart=CLOSEDDEFAULT,int Echo=1);
   NewtonQRUpdatePCSolution(LatticeMode *Mode,PerlInput &Input,const Vector &one,int Echo=1);
   NewtonQRUpdatePCSolution(LatticeMode *Mode,PerlInput &Input,int Echo);
   ~NewtonQRUpdatePCSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual int FindNextSolution();
   virtual int FindCriticalPoint(Lattice *Lat,PerlInput &Input,int Width,fstream &out);
};

#endif
