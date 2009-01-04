#ifndef RSE__DFTExternal
#define RSE__DFTExternal

#include "PerlInput.h"
#include "Lattice.h"
#include <fstream>

using namespace std;

class DFTExternal : public Lattice
{
private:
   int DOFS_;
   
   // DOF[i] = [B11, B22, B33, B23, B31, B12, S11, S12, S13, S21,...,SN3]
   Vector DOF_;
   enum LDeriv {L0,DL};
   double Lambda_;

   int Width_;

   enum UpdateFlag {NoStiffness,NeedStiffness};
   void UpdateValues(UpdateFlag flag) const;

   static const int cachesize = 4;
   mutable int Cached_[cachesize];
   mutable double E0CachedValue_;
   mutable Vector E1CachedValue_;
   mutable Vector E1DLoadCachedValue_;
   mutable Matrix E2CachedValue_;
   mutable int CallCount_[cachesize];
   
public:
   // Functions provided by DFTExternal
   DFTExternal(PerlInput const& Input,int const& Echo=1,int const& Width=20);
   ~DFTExternal();
   
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
   virtual Matrix const& E2() const;
   virtual char const* const Type() const {return "DFTExternal";}
   virtual void Print(ostream& out,PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);

   friend ostream& operator<<(ostream& out,DFTExternal& A);
   
   // ignore these
   double Entropy() const {return 0.0;}
   double HeatCapacity() const {return 0.0;}
   Vector const& StressDT() const {return EmptyV_;}
   Matrix const& StiffnessDT() const {return EmptyM_;}
   double Temp() const {return 0.0;}
   void SetTemp(double const& Ntemp) {}
   Vector const& StressDL() const {return E1DLoad();}
   Matrix const& StiffnessDL() const
   {cerr << "DFTExternal::StiffnessDL() Not Programmed\n"; return EmptyM_;}
   virtual Matrix const& E3() const
   {cerr << "DFTExternal::E3() Not Programmed\n"; exit(-1); return EmptyM_;}
   virtual Matrix const& E4() const
   {cerr << "DFTExternal::E4() Not Programmed\n"; exit(-1); return EmptyM_;}
   virtual void SetParameters(double const* const Vals,int const& ResetRef = 1) {}
   virtual void SetGridSize(int const& Grid) {}
   
private:
   // place holder
   Vector EmptyV_;
   Matrix EmptyM_;

   static const double Alt[3][3][3];
   static const double Del[3][3];
   static fstream dbug;
};

#endif
