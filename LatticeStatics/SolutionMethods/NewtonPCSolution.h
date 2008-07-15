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
public:
   enum UpdateType {NoUpdate,QRUpdate,StiffnessUpdate,Exact};
private:
   LatticeMode *Mode_;
   
   int Echo_;
   int CurrentSolution_;
   UpdateType UpdateType_;             // 0-QR update (default), 1-Stiffness update, 2-none
   int NumSolutions_;
   
   double MaxDS_;
   double CurrentDS_;           //Initial Steplength h > 0
   double MinDS_;               //Minimum Stepsize
   double cont_rate_max_;       //Max contraction rate
   double delta_max_;           //Max distance to (from predicted to corrected point) curve
   double alpha_max_;           //Max angle to curve (must be less than pi/2)
   double Converge_;            //Convergence criteria
   int ClosedLoopStart_;        //Closed loop test variable
   int StopAtCPNum_;            //Stop at critical point test flag
   int Direction_;              //Direction of tangent
   double Omega_;               //Multiplier to help traverse bifurcation points
   double accel_max_;           //Max acceleration rate
   
   Vector FirstSolution_;       //Initial point on curve
   Vector PreviousSolution_;    //Previous point on curve
   Vector Tangent1_;            //Tangent vector of ith point
   Vector Tangent2_;            //Tangent Vector of ith + 1 point
   
   void GetQR(Vector const& Force,Vector const& diff,Matrix& Q,Matrix& R) const;
   void MoorePenrose(Matrix const& Q,Matrix const& R,Vector const& Force,Vector& Corrector)
      const;
   void UpdateQR(Vector const& Force,Vector const& difference,Matrix& Q,Matrix& R)
      const;
   
public:
   NewtonPCSolution(LatticeMode* const Mode,Vector const& one,
                    int const& CurrentSolution,UpdateType const& Type,int const& NumSolutions,
                    double const& MaxDS,double const& CurrentDS,double const& MinDS,
                    double const& cont_rate_max,double const& delta_max,double const& alpha_max,
                    double const& Converge,Vector const& FirstSolution,int const& Direction=1,
                    double const& accel_max=2.0,int const& ClosedLoopStart=CLOSEDDEFAULT,
                    int const& StopAtCPNum=-1,int const& Echo=1);
   NewtonPCSolution(LatticeMode* const Mode,PerlInput const& Input,Vector const& one,
                    int const& Echo=1);
   NewtonPCSolution(LatticeMode* const Mode,PerlInput const& Input,int const& Echo);
   ~NewtonPCSolution() {}
   
   // Functions required by SolutionMethod
   virtual int AllSolutionsFound() const;
   virtual int FindNextSolution();
   virtual int FindCriticalPoint(Lattice* const Lat,PerlInput const& Input,int const& Width,
                                 fstream& out);
   
private:
   // "static" member variables
   // FindNextSolution
   mutable Vector v_static;
   mutable Vector w_static;
   mutable Vector Force_static;
   mutable Vector Corrector_static;
   mutable Vector difference_static;
   mutable Matrix Q_static;
   mutable Matrix R_static;
   // GetQR
   mutable Matrix Stiff_static;
   // MoorePenrose
   mutable Vector y_static;
   // QRUpdate
   mutable Vector u_static;
   mutable Vector a_static;
   mutable Vector e_static;
};

#endif
