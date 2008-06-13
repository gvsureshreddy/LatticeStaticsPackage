#ifndef __NewtonPCSolution
#define __NewtonPCSolution

#include "PerlInput.h"
#include "SolutionMethod.h"
#include "LatticeMode.h"

using namespace std;

#define CLOSEDDEFAULT 30

class Lattice;

class NewtonPCSolution : public SolutionMethod
{
private:
   LatticeMode *Mode_;
   
   int Echo_;
   int CurrentSolution_;
   int UpdateType_;             // 0-QR update (default), 1-Stiffness update, 2-none
   int NumSolutions_;
   
   double MaxDS_;
   double CurrentDS_;           //Initial Steplength h > 0
   double cont_rate_nom_;       //Nominal contraction rate
   double delta_nom_;           //Nominal distance to (from predicted to corrected point) curve
   double alpha_nom_;           //Nominal angle to curve
   double Converge_;            //Convergence criteria
   double MinDSRatio_;          //Minimum Stepsize ratio
   int ClosedLoopStart_;        //Closed loop test variable
   int StopAtCPNum_;            //Stop at critical point test flag
   int Direction_;              //Direction of tangent
   
   Vector FirstSolution_;       //Initial point on curve
   Vector Tangent1_;            //Tangent vector of ith point
   Vector Tangent2_;            //Tangent Vector of ith + 1 point
   
   void GetQR(Vector& Force,Vector& diff,Matrix& Q,Matrix& R);
   void MoorePenrose(const Matrix& Q,const Matrix& R,const Vector& Force,Vector& Corrector);
   void QRUpdate(const Vector& Force,const Vector& difference,Matrix& Q,Matrix& R);
   
public:
   Vector Previous_Solution_;
   NewtonPCSolution(LatticeMode *Mode,const Vector &one,
                    int CurrentSolution,int UpdateType,int NumSolutions,double MaxDS,
                    double CurrentDS,double cont_rate_nom,double delta_nom,
                    double alpha_nom,double Converge,double MinDSRatio,
                    const Vector &FirstSolution,int Direction=1,
                    int ClosedLoopStart=CLOSEDDEFAULT,int StopAtCPNum=-1,int Echo=1);
   NewtonPCSolution(LatticeMode *Mode,PerlInput &Input,const Vector &one,int Echo=1);
   NewtonPCSolution(LatticeMode *Mode,PerlInput &Input,int Echo);
   ~NewtonPCSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual int FindNextSolution();
   virtual int FindCriticalPoint(Lattice *Lat,PerlInput &Input,int Width,fstream &out);
   
private:
   // "static" member variables
   // FindNextSolution
   Vector v_static;
   Vector w_static;
   Vector Force_static;
   Vector Corrector_static;
   Vector difference_static;
   Matrix Q_static;
   Matrix R_static;
   // GetQR
   Matrix Stiff_static;
   // MoorePenrose
   Vector y_static;
   // QRUpdate
   Vector u_static;
   Vector a_static;
   Vector e_static;
};

#endif
