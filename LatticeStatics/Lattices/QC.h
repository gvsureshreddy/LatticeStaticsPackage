#ifndef RSE__QC
#define RSE__QC

#include "PerlInput.h"
#include "Lattice.h"

using namespace std;

class QC : public Lattice
{
private:
   mutable int DOFS_;
   
   mutable Vector DOF_;
   mutable double Lambda_;

   int Echo_;
   int Width_;
   double Tolerance_;

   enum UpdateFlag {NoStiffness=0,NeedStiffness=1};
   void UpdateValues(UpdateFlag flag) const;

   static const int cachesize = 4;
   mutable int Cached_[cachesize];
   mutable double E0CachedValue_;
   mutable Vector E1CachedValue_;
   mutable Vector E1DLoadCachedValue_;
   mutable Matrix E2CachedValue_;
   mutable int CallCount_[cachesize];
   
public:
   // Functions provided by QC
   QC(PerlInput const& Input,int const& Echo=1,int const& Width=20);
   ~QC();
   
   // Virtual Functions required by Lattice
   Vector const& DOF() const {return DOF_;}
   void SetDOF(Vector const& dof)
   {DOF_ = dof; for (int i=0;i<cachesize;++i) Cached_[i]=0;}

   double Lambda() const {return Lambda_;}
   void SetLambda(double const& lambda)
   {Lambda_ = lambda; for (int i=0;i<cachesize;++i) Cached_[i]=0;}
   
   virtual double E0() const;
   virtual Vector const& E1() const;
   virtual Vector const& E1DLoad() const;
   virtual Vector const& StressDL() const {return E1DLoad();}
   virtual Matrix const& E2() const;
   virtual Matrix const& StiffnessDL() const;
   virtual Matrix const& E3() const;
   virtual void Print(ostream& out,PrintDetail const& flag);

   friend ostream& operator<<(ostream& out,QC& A);
   
   virtual int CriticalPointInfo(int const& CPCrossingNum,char const& CPSubNum,
                                 Vector const& DrDt,int const& NumZeroEigenVals,
                                 double const& Tolerance,int const& Width,
                                 PerlInput const& Input,ostream& out);
   // ignore these
   double Entropy() const {return 0.0;}
   double HeatCapacity() const {return 0.0;}
   Vector const& StressDT() const {return EmptyV_;}
   Matrix const& StiffnessDT() const {return EmptyM_;}
   double Temp() const {return 0.0;}
   void SetTemp(double const& Ntemp) {}
   virtual Matrix const& E4() const
   {cerr << "QC::E4() Not Programmed\n"; exit(-1); return EmptyM_;}
   virtual void SetParameters(double const* const Vals,int const& ResetRef = 1) {}
   virtual void SetGridSize(int const& Grid) {}
   
private:
   // statice for StiffnessDL and E3
   mutable Matrix stiffdl_static;
   mutable Matrix E3_static;
   // place holder
   Vector EmptyV_;
   Matrix EmptyM_;
};

#endif
