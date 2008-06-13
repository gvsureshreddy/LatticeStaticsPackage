#ifndef __TwoBarTruss
#define __TwoBarTruss

#include "PerlInput.h"
#include "Lattice.h"

using namespace std;

class TwoBarTruss : public Lattice
{
private:
   int DOFS_;
   
   // DOF[i] = [u v]
   Vector DOF_;
   enum LDeriv {L0,DL};
   double Lambda_;
   double Theta_;
   double COSTheta_;
   double SINTheta_;

   int Echo_;
   int Width_;

   int Caching_;
   int Cached_[3];
   double E0CachedValue_;
   Matrix E1CachedValue_;
   Matrix E2CachedValue_;
   int CallCount_[3];
   
public:
   // Functions provided by TwoBarTruss
   TwoBarTruss(PerlInput& Input,int Echo=1,int Width=20);
   ~TwoBarTruss();
   
   // Virtual Functions required by Lattice
   Vector DOF() {return DOF_;}
   void SetDOF(const Vector &dof)
   {DOF_ = dof; if (Caching_) for (int i=0;i<3;++i) Cached_[i]=0;}

   double Lambda() {return Lambda_;}
   void SetLambda(const double &lambda)
   {Lambda_ = lambda; if (Caching_) for (int i=0;i<3;++i) Cached_[i]=0;}
   
   virtual double E0();
   virtual Matrix E1();
   virtual Matrix E1DLoad();
   virtual Matrix E2();
   virtual void Print(ostream &out,PrintDetail flag);

   friend ostream &operator<<(ostream &out,TwoBarTruss &A);
   
   // ignore these
   virtual void CriticalPointInfo(const Vector &DrDt,int NumZeroEigenVals,
                                  double Tolerance,int Width,ostream &out) {}
   double Entropy() {return 0.0;}
   double HeatCapacity() {return 0.0;}
   Matrix StressDT() {return Matrix();}
   Matrix StiffnessDT() {return Matrix();}
   double Temp() {return 0.0;}
   void SetTemp(const double &Ntemp) {}
   Matrix StressDL() {return Matrix();}
   Matrix StiffnessDL() {return Matrix();}
   virtual Matrix E3() {return Matrix();}
   virtual Matrix E4() {return Matrix();}
   virtual void SetParameters(double *Vals,int ResetRef = 1) {}
   virtual void SetGridSize(int Grid) {}
   
private:
   
};

#endif
