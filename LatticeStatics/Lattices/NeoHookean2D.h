#ifndef RSE__NeoHookean2D
#define RSE__NeoHookean2D

#include <string>
#include <sstream>
#include <cstdlib>
#include "Lattice.h"

class NeoHookean2D : public Lattice
{
private:
   mutable int DOFS_;

   mutable Vector DOF_;
   mutable Vector RHS_;
   mutable Matrix Stiff_;
   mutable double Lambda_;

   int Width_;
   std::size_t system_size_;

public:
   // Functions provided by NeoHookean2D
   NeoHookean2D(PerlInput const& Input, int const& Echo = 1,
                int const& Width = 20);
   ~NeoHookean2D();

   // Virtual Functions required by Lattice
   Vector const& DOF() const
   {
      return DOF_;
   }

   void SetDOF(Vector const& dof)
   {
      DOF_ = dof;
   }

   double Lambda() const
   {
      return Lambda_;
   }

   void SetLambda(double const& lambda)
   {
      Lambda_ = lambda;
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
   virtual void ExtraTestFunctions(Vector& TF) const;
   virtual char const* const Type() const
   {
      return "NeoHookean2D";
   }

   virtual void Print(ostream& out, PrintDetail const& flag,
                      PrintPathSolutionType const& SolType = RegularPt);

   friend ostream& operator<<(ostream& out, NeoHookean2D& A);

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
     cerr << "NeoHookean2D::E3() Not Programmed\n"; exit(-1); return EmptyM_;
   }

   virtual Matrix const& E4() const
   {
      cerr << "NeoHookean2D::E4() Not Programmed\n"; exit(-1); return EmptyM_;
   }

   virtual void SetParameters(double const* const Vals, int const& ResetRef = 1)
   {
   }

   virtual void SetGridSize(int const& Grid)
   {
   }

private:
   // statice for StiffnessDL and E3
   mutable Matrix stiffdl_static;
   mutable Matrix E3_static;

   // place holder
   Vector EmptyV_;
   Matrix EmptyM_;
};

#endif
