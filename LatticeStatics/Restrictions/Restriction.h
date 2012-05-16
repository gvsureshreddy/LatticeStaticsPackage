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
   Lattice* Lattice_;

public:
   virtual ~Restriction();

   Restriction(PerlInput const& Input);

   Lattice* const Lat() const
   {
      return Lattice_;
   }
   
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

   int const NumTestFunctions() const
   {
      return Lattice_->NumTestFunctions();
   }
   
   int TestFunctions(Vector& TF1, Lattice::StateType const& State = Lattice::LHS,
                     Vector* const EV2 = 0) const
   {
      return Lattice_->TestFunctions(TF1, State, EV2);
   }

   Matrix const& RelativeEigVects() const
   {
      return Lattice_->RelativeEigVects();
   }

   void ConsistencyCheck(Vector const& dof, double const& ConsistencyEpsilon, int const& Width,
                         ostream& out);

   void ConsistencyCheckRidders(Vector const& dof, double const& ConsistencyEpsilon, int const& Width,
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
