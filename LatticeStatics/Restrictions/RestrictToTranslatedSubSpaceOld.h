#ifndef RSE__RestrictToTranslatedSubSpaceOld
#define RSE__RestrictToTranslatedSubSpaceOld

#include "PerlInput.h"
#include "Restriction.h"
#include "Lattice.h"

#include <cmath>
#define DOFMAX 256

using namespace std;

class RestrictToTranslatedSubSpaceOld : public Restriction
{
private:
   Lattice *Lattice_;
   
   int DOFS_;
   int DOFindlen_[DOFMAX];
   int DOFindex_[DOFMAX][DOFMAX];
   Vector DOFMult_[DOFMAX];
   Vector ReferenceState_;
   Vector DOF_;
   
   void UpdateLatticeState();
public:
   RestrictToTranslatedSubSpaceOld(Lattice* const M,PerlInput const& Input);
   
   ~RestrictToTranslatedSubSpaceOld() {}
   
   // Functions required by Restriction
   virtual double Energy() const {return Lattice_->E0();}
   virtual Vector const& DrDt(Vector const& Diff) const;
   
   //----------------------------------------------------------------
   virtual Vector const& Force() const;
   virtual Matrix const& Stiffness() const;
   virtual Vector const& DOF() const {return DOF_;}
   virtual Vector RestrictDOF(Vector const& dof);
   virtual Vector UnRestrictDOF(Vector const& dof);
   virtual Vector TransformVector(Vector const& T);
   virtual Vector UnTransformVector(Vector const& T);
   virtual void SetDOF(Vector const& dof);
   virtual void UpdateDOF(Vector const& dr);
   //----------------------------------------------------------------
   virtual char const* const Name() const {return "RestrictToTranslatedSubSpaceOld";}

private:
   // "static" member variables
   // UpdateLatticeState
   int size_static;
   Vector DOF_static;
   // DrDt
   mutable Vector ddt_static;
   // Force
   mutable Vector force_static;
   mutable Vector stress_static;
   // Stiffness
   mutable Matrix K_static;
   mutable Matrix Stiff_static;
   mutable Vector stressdt_static;
};

#endif
