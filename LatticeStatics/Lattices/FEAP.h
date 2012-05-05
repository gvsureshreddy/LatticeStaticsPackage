#ifndef RSE__FEAP
#define RSE__FEAP

#include <string>
#include <sstream>
#include "PerlInput.h"
#include "Lattice.h"

using namespace std;

class FEAP : public Lattice
{
private:
   mutable int DOFS_;

   mutable Vector DOF_;
   mutable double Lambda_;

   char const* ffin_; // FEAP input file name
   int ndf_; // number of DOFs per node
   int ndm_; // number of spatial dimensions
   int numnp_; // number of nodes in mesh
   int neq_; // number of reduced equations
   int* eqnID_; // equation number ID array
   int* bcID_; // displacement boundary condition ID array
   Vector X_; // Reference coordinates of nodes

   int Width_;
   double Tolerance_;

   enum UpdateFlag {NoStiffness = 0, NeedStiffness = 1};
   void UpdateValues(UpdateFlag flag) const;

   static const int cachesize = 4;
   mutable int Cached_[cachesize];
   mutable double E0CachedValue_;
   mutable Vector E1CachedValue_;
   mutable Vector E1DLoadCachedValue_;
   mutable Matrix E2CachedValue_;
   mutable int EvaluationCount_[2];
   mutable int CallCount_[cachesize];

public:
   // Functions provided by FEAP
   FEAP(PerlInput const& Input, int const& Echo = 1, int const& Width = 20);
   ~FEAP();

   // Virtual Functions required by Lattice
   Vector const& DOF() const
   {
      return DOF_;
   }

   void SetDOF(Vector const& dof)
   {
      DOF_ = dof;
      for (int i = 0; i < cachesize; ++i)
      {
         Cached_[i] = 0;
      }
   }

   double Lambda() const
   {
      return Lambda_;
   }

   void SetLambda(double const& lambda)
   {
      Lambda_ = lambda;
      for (int i = 0; i < cachesize; ++i)
      {
         Cached_[i] = 0;
      }
   }

   virtual double E0() const;
   virtual Vector const& E1() const;
   virtual Vector const& E1DLoad() const;
   virtual Vector const& StressDL() const
   {
      return E1DLoad();
   }

   virtual Matrix const& E2() const;
   virtual Matrix const& StiffnessDL() const;
   virtual char const* const Type() const
   {
      return "FEAP";
   }

   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);

   friend ostream& operator<<(ostream& out, FEAP& A);

   // ignore these
   double Entropy() const
   {
      return 0.0;
   }

   double HeatCapacity() const
   {
      return 0.0;
   }

   Vector const& StressDT() const
   {
      return EmptyV_;
   }

   Matrix const& StiffnessDT() const
   {
      return EmptyM_;
   }

   double Temp() const
   {
      return 0.0;
   }

   void SetTemp(double const& Ntemp)
   {
   }

   virtual Matrix const& E3() const
   {
      cerr << "FEAP::E3() Not Programmed\n"; exit(-1); return EmptyM_;
   }

   virtual Matrix const& E4() const
   {
      cerr << "FEAP::E4() Not Programmed\n"; exit(-1); return EmptyM_;
   }

   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
   {
   }

   virtual void SetGridSize(int const& Grid)
   {
   }

private:
   // statice for StiffnessDL
   mutable Matrix stiffdl_static;

   // place holder
   Vector EmptyV_;
   Matrix EmptyM_;
};

#endif
