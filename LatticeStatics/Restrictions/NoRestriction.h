#ifndef RSE__NoRestriction
#define RSE__NoRestriction

#include "Restriction.h"
#include "Lattice.h"

using namespace std;

class NoRestriction : public Restriction
{
private:
   Lattice *Lattice_;
   
public:
   NoRestriction(Lattice* const M): Lattice_(M),
                                    ddt_static(Lattice_->DOF().Dim()),
                                    Stiff_static(Lattice_->DOF().Dim(),Lattice_->DOF().Dim()),
                                    stressdt_static(1,Lattice_->DOF().Dim())
   {}
   
   ~NoRestriction() {}
   
   // Functions required by Restriction
   virtual double Energy() const {return Lattice_->E0();}
   virtual Vector const& DrDt(Vector const& Diff) const
   {
      for (int i=0;i<ddt_static.Dim();++i) ddt_static[i] = Diff[i]/Diff[ddt_static.Dim()];
      return ddt_static;
   }
   
   //----------------------------------------------------------------
   virtual Vector const& Force() const {return Lattice_->E1();}
   virtual Matrix const& Stiffness() const
   {
      K_static.Resize(Lattice_->DOF().Dim(),Lattice_->DOF().Dim()+1,0.0);
      Stiff_static = Lattice_->E2();
      stressdt_static = Lattice_->E1DLoad();
      for (int i=0;i<Stiff_static.Rows();++i)
      {
         for (int j=0;j<Stiff_static.Cols();++j)
         {
            K_static[i][j] = Stiff_static[i][j];
         }
         K_static[i][K_static.Cols()-1] = stressdt_static[i];
      }
      return K_static;
   }
   virtual Vector const& DOF() const {return Lattice_->DOF();}
   virtual void SetDOF(Vector const& dof) {Lattice_->SetDOF(dof);}
   virtual void UpdateDOF(Vector const& dr) {Lattice_->SetDOF(Lattice_->DOF()+dr);}
   //----------------------------------------------------------------
   virtual char const* const Name() const {return "NoRestriction";}

private:
   // "static" member variables
   // DrDt
   mutable Vector ddt_static;
   // Stiffness
   mutable Matrix K_static;
   mutable Matrix Stiff_static;
   mutable Vector stressdt_static;
};

#endif
