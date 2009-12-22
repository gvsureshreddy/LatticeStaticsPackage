#ifndef RSE__RestrictToTranslatedSubSpace
#define RSE__RestrictToTranslatedSubSpace

#include "PerlInput.h"
#include "Restriction.h"
#include "Lattice.h"
#include "SparseMatrix.h"

using namespace std;

class RestrictToTranslatedSubSpace : public Restriction
{
private:
   Lattice *Lattice_;
   
   int DOFS_;
   SparseMatrix ForceProjectionMatrix_;
   SparseMatrix DOFProjectionMatrix_;
   SparseMatrix* SymmetryCheck_;
   int SymmetryCheckCount_;
   double SymmetryCheckTol_;
   SparseMatrix const& ForceProject_;
   SparseMatrix const& DOFProject_;
   Vector ReferenceState_;
   Vector DOF_;
   
   void UpdateLatticeState();
public:
   RestrictToTranslatedSubSpace(Lattice* const M,PerlInput const& Input);
   
   ~RestrictToTranslatedSubSpace();
   
   // Functions required by Restriction
   virtual double Energy() const {++counter_[1]; return Lattice_->E0();}
   virtual Vector const& DrDt(Vector const& Diff) const;
   
   //----------------------------------------------------------------
   virtual Vector const& Force() const;
   virtual Matrix const& Stiffness() const;
   virtual Vector const& DOF() const {++counter_[5]; return DOF_;}
   virtual int SymmetryOK() const;
   virtual Vector RestrictDOF(Vector const& dof);
   virtual Vector UnRestrictDOF(Vector const& dof);
   virtual Vector TransformVector(Vector const& T);
   virtual Vector UnTransformVector(Vector const& T);
   virtual void SetDOF(Vector const& dof);
   virtual void UpdateDOF(Vector const& dr);
   //----------------------------------------------------------------
   virtual char const* const Name() const {return "RestrictToTranslatedSubSpace";}

private:
   enum CounterID {UPDATE,ENERGY,DRDT,FORCE,STIFFNESS,DOFCount,SYMMETRY,RESTRICT,UNRESTRICT,
                   TRANSFORM,UNTRANSFORM,SETDOF,UPDATEDOF,ENDOFCOUNTERLIST};
   // "static" member variables
   // UpdateLatticeState
   int size_static;
   Vector Lat_DOF_static;
   Vector Rest_DOF_static;
   // DrDt
   mutable Vector Lat_ddt_static;
   mutable Vector Rest_ddt_static;
   // Force
   mutable Vector force_static;
   mutable Vector stress_static;
   // Stiffness
   mutable Matrix K_static;
   mutable Matrix E2_tmp_static;
   mutable Matrix Stiff_static;
   mutable Vector stressdt_static;

   // counters
   static int const nocounters_ = ENDOFCOUNTERLIST;
   mutable int counter_[nocounters_];
};

#endif
