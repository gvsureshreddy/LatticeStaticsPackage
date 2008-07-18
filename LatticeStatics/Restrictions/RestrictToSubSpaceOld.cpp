#include "RestrictToSubSpaceOld.h"

RestrictToSubSpaceOld::RestrictToSubSpaceOld(Lattice* const M,PerlInput const& Input)
{
   char tmp[LINELENGTH];

   Lattice_ = (Lattice *) M;

   PerlInput::HashStruct Hash = Input.getHash("Restriction","RestrictToSubSpaceOld");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   DOF_.Resize(DOFS_+1,0.0);

   int LatDOFS = Lattice_->DOF().Dim();
   
   for (int i=0;i<DOFS_;++i)
   {
      sprintf(tmp,"DOF_%u",i);
      DOFindlen_[i] = Input.getArrayLength(Hash,tmp,0);
      Input.getIntVector(DOFindex_[i],DOFindlen_[i],Hash,tmp,0);
      DOFMult_[i].Resize(DOFindlen_[i]);
      Input.getVector(DOFMult_[i],Hash,tmp,1);
   }
   
   //Baseline DOF Initialization
   BaselineDOF_.Resize(LatDOFS,0.0);

   char const* UseBaseLineState;
   if (Input.ParameterOK(Hash,"UseBaseLineState"))
   {
      UseBaseLineState = Input.getString(Hash,"UseBaseLineState");
   }
   else
   {
      UseBaseLineState = Input.useString("No",Hash,"UseBaseLineState");
   }
   
   if (!strcmp("Yes",UseBaseLineState))
   {
      int sz = Input.getArrayLength(Hash,"BaseLineState",0);
      int* pos = new int[sz];
      Input.getPosIntVector(pos,sz,Hash,"BaseLineState",0);
      Vector Vals(sz);
      Input.getVector(Vals,Hash,"BaseLineState",1);
      for (int i=0;i<sz;++i)
      {
         BaselineDOF_[pos[i]] = Vals[i];
      }
      delete [] pos;
   }
   Input.EndofInputSection();

   // intitalize "static" member variables
   size_static = LatDOFS;
   DOF_static.Resize(size_static);
   force_static.Resize(DOFS_);
   stress_static.Resize(size_static);
   K_static.Resize(DOFS_,DOFS_+1);
   Stiff_static.Resize(size_static,size_static);
   stressdt_static.Resize(size_static);

}

// Functions required by Restriction
Vector const& RestrictToSubSpaceOld::DrDt(Vector const& Diff) const
{
   ddt_static.Resize(Lattice_->DOF().Dim(),0.0);
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         ddt_static[DOFindex_[i][j]] += DOFMult_[i][j]*(Diff[i]/Diff[DOFS_]);
      }
   }
   
   return ddt_static;
}

//----------------------------------------------------------------
void RestrictToSubSpaceOld::UpdateLatticeState()
{
   for (int i=0;i<size_static;++i)
   {
      DOF_static[i] = BaselineDOF_[i];
   }
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         DOF_static[DOFindex_[i][j]] += DOFMult_[i][j]*DOF_[i];
      }
   }
   
   Lattice_->SetDOF(DOF_static);
   Lattice_->SetLoadParameter(DOF_[DOFS_]);
}

Vector const& RestrictToSubSpaceOld::Force() const
{
   stress_static = Lattice_->E1();
   
   for (int i=0;i<DOFS_;++i)
   {
      force_static[i] = 0.0;
      for (int j=0;j<DOFindlen_[i];j++)
      {
         force_static[i] += DOFMult_[i][j]*stress_static[DOFindex_[i][j]];
      }
   }
   
   return force_static;
}

Matrix const& RestrictToSubSpaceOld::Stiffness() const
{
   K_static.Resize(DOFS_,DOFS_+1,0.0);
   Stiff_static = Lattice_->E2();
   stressdt_static = Lattice_->E1DLoad();
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         for (int k=0;k<DOFS_;++k)
         {
            for (int l=0;l<DOFindlen_[k];++l)
            {
               K_static[i][k] +=
                  DOFMult_[i][j]*(Stiff_static[DOFindex_[i][j]][DOFindex_[k][l]])*DOFMult_[k][l];
            }
         }
         
         K_static[i][DOFS_] += DOFMult_[i][j]*stressdt_static[DOFindex_[i][j]];
      }
   }
   
   return K_static;
}

void RestrictToSubSpaceOld::SetDOF(Vector const& dof)
{
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] = dof[i];
   }
   
   UpdateLatticeState();
}

void RestrictToSubSpaceOld::UpdateDOF(Vector const& dr)
{
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] += dr[i];
   }
   
   UpdateLatticeState();
}
