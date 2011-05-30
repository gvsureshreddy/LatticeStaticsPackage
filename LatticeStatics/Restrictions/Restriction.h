#ifndef RSE__Restriction
#define RSE__Restriction

#include "Lattice.h"
#include "SparseMatrix.h"

class Restriction
{
protected:
   SparseMatrix* SymmetryCheck_;
   int SymmetryCheckCount_;
   double SymmetryCheckTol_;

public:
   virtual ~Restriction();

   Restriction(PerlInput const& Input);

   virtual double Energy() const = 0;
   virtual Vector const& DrDt(Vector const& Diff) const = 0;

   virtual Vector const& Force() const = 0;
   virtual Matrix const& Stiffness() const = 0;
   virtual Vector const& DOF() const = 0;
   virtual int SymmetryOK() const;
   virtual Vector RestrictDOF(Vector const& dof) = 0;
   virtual Vector UnRestrictDOF(Vector const& dof) = 0;
   virtual Vector TransformVector(Vector const& T) = 0;
   virtual Vector UnTransformVector(Vector const& T) = 0;
   virtual void SetDOF(Vector const& dof) = 0;
   virtual void UpdateDOF(Vector const& dr) = 0;

   void ConsistencyCheck(Vector const& dof, double const& ConsistencyEpsilon, int const& Width,
                         ostream& out);

   virtual char const* const Name() const
   {
      return "Restriction";
   }

   friend ostream& operator<<(ostream& out, Restriction const& R)
   {
      out << R.Name(); return out;
   }
};

#endif
