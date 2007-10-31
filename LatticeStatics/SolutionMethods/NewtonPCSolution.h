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

   // other parameters?
	double CurrentDS_;					//initial Steplength h > 0
	double cont_rate_nom_;				//Nominal contraction rate
	double delta_nom_;					//Nominal distance to (from predicted to corrected point) curve 
	double alpha_nom_;					//Nominal angle to curve
	double Converge_;					//Convergence criteria
   
	int ClosedLoopStart_;				//Closed loop test variable
   
	Vector FirstSolution_;				//Initial point on curve
	Vector Tangent1_;					//Tangent vector of ith point
	Vector Tangent2_;					//Tangent Vector of ith + 1 point
	
	Matrix Stiffness_;					//Stiffness Matrix of point
	//Vector Force1_;
	//Vector Force2_;
	
   
public:
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,const Vector &one,
		    int Echo=1);
   NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,char *startfile,
		    fstream &out,int Echo);
   ~NewtonPCSolution() {}

   // Functions required by SolutionMethod
   virtual int AllSolutionsFound();
   virtual double FindNextSolution(int &good);
   virtual int BisectAlert(Lattice *Lat,char *datafile,const char *prefix,
			   int Width,fstream &out);
   virtual Vector MoorePenrose(const Matrix& Q, const Matrix& R, const Vector& Force);
};

#endif
