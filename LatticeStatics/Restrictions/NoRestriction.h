#ifndef RSE__NoRestriction
#define RSE__NoRestriction

#include "Restriction.h"
#include "Lattice.h"

using namespace std;

class NoRestriction : public Restriction
{
private:
   void UpdateLatticeDOF()
   {
      Vector tmpdof(Lattice_->DOF().Dim());
      for (int i = 0; i < tmpdof.Dim(); ++i)
      {
         tmpdof[i] = dof_static[i];
      }
      Lattice_->SetDOF(tmpdof);

      if (Lattice_->LoadParameter() == Lattice::Temperature)
      {
         Lattice_->SetTemp(dof_static[tmpdof.Dim()]);
      }
      else if (Lattice_->LoadParameter() == Lattice::Load)
      {
         Lattice_->SetLambda(dof_static[tmpdof.Dim()]);
      }
      else
      {
         cerr << "Error, Unknown LoadParameter type! exiting...\n";
         exit(-15);
      }
   }

public:
   NoRestriction(Lattice* const M, PerlInput const& Input) :
      Restriction(Input), dof_static(M->DOF().Dim() + 1),
      ddt_static(M->DOF().Dim()),
      Stiff_static(M->DOF().Dim(), M->DOF().Dim()),
      stressdt_static(M->DOF().Dim())
   {
      Lattice_ = M;
      
      int DOFS = dof_static.Dim();
      for (int i = 0; i < SymmetryCheckCount_; ++i)
      {
         if (SymmetryCheck_[i].Cols() != DOFS)
         {
            cerr << "Error. " << Name()
                 << " Incorrect number of columns in SymmetryCheckProjectionMatrix"
                 << " number " << i << "\n";
            exit(-38);
         }
      }

      // Make sure everything is synchronized
      Vector t = Lattice_->DOF();
      Vector tt(t.Dim() + 1);
      for (int i = 0; i < t.Dim(); ++i)
      {
         tt[i] = t[i];
      }
      tt[t.Dim()] = ((Lattice::Temperature == Lattice_->LoadParameter()) ?
                     Lattice_->Temp() : Lattice_->Lambda());
      SetDOF(tt);
   }

   ~NoRestriction()
   {
   }

   // Functions required by Restriction
   virtual double Energy() const
   {
      return Lattice_->E0();
   }

   virtual Vector const& DrDt(Vector const& Diff) const
   {
      for (int i = 0; i < ddt_static.Dim(); ++i)
      {
         ddt_static[i] = Diff[i] / Diff[ddt_static.Dim()];
      }
      return ddt_static;
   }

   // ----------------------------------------------------------------
   virtual Vector const& Force() const
   {
      return Lattice_->E1();
   }

   virtual Matrix const& Stiffness() const
   {
      K_static.Resize(Lattice_->DOF().Dim(), Lattice_->DOF().Dim() + 1, 0.0);
      Stiff_static = Lattice_->E2();
      stressdt_static = Lattice_->E1DLoad();
      for (int i = 0; i < Stiff_static.Rows(); ++i)
      {
         for (int j = 0; j < Stiff_static.Cols(); ++j)
         {
            K_static[i][j] = Stiff_static[i][j];
         }
         K_static[i][K_static.Cols() - 1] = stressdt_static[i];
      }
      return K_static;
   }

   virtual Vector const& DOF() const
   {
      return dof_static;
   }

   virtual Vector RestrictDOF(Vector const& dof)
   {
      return dof;
   }

   virtual Vector UnRestrictDOF(Vector const& dof)
   {
      return dof;
   }

   virtual Vector TransformVector(Vector const& T)
   {
      return T;
   }

   virtual Vector UnTransformVector(Vector const& T)
   {
      return T;
   }

   virtual void SetDOF(Vector const& dof)
   {
      dof_static = dof; UpdateLatticeDOF();
   }

   virtual void UpdateDOF(Vector const& dr)
   {
      dof_static += dr; UpdateLatticeDOF();
   }

   // ----------------------------------------------------------------
   virtual char const* const Name() const
   {
      return "NoRestriction";
   }

private:
   // "static" member variables
   mutable Vector dof_static;

   // DrDt
   mutable Vector ddt_static;

   // Stiffness
   mutable Matrix K_static;
   mutable Matrix Stiff_static;
   mutable Vector stressdt_static;
};

#endif
