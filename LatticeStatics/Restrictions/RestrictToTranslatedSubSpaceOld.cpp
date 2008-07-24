#include "RestrictToTranslatedSubSpaceOld.h"
#include <sstream>

RestrictToTranslatedSubSpaceOld::RestrictToTranslatedSubSpaceOld(Lattice* const M,
                                                                 PerlInput const& Input)
{
   stringstream tmp;

   Lattice_ = (Lattice *) M;

   PerlInput::HashStruct Hash = Input.getHash("Restriction","RestrictToTranslatedSubSpaceOld");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   DOF_.Resize(DOFS_+1,0.0);

   int LatDOFS = Lattice_->DOF().Dim();
   
   for (int i=0;i<DOFS_;++i)
   {
      tmp.str("");
      tmp << "DOF_" << i;
      DOFindlen_[i] = Input.getArrayLength(Hash,tmp.str().c_str(),0);
      Input.getIntVector(DOFindex_[i],DOFindlen_[i],Hash,tmp.str().c_str(),0);
      DOFMult_[i].Resize(DOFindlen_[i]);
      Input.getVector(DOFMult_[i],Hash,tmp.str().c_str(),1);
   }
   
   //ReferenceState DOF Initialization
   ReferenceState_.Resize(LatDOFS+1,0.0);

   char const* UseReferenceState;
   if (Input.ParameterOK(Hash,"UseReferenceState"))
   {
      UseReferenceState = Input.getString(Hash,"UseReferenceState");
   }
   else
   {
      UseReferenceState = Input.useString("No",Hash,"UseReferenceState");
   }
   
   if (!strcmp("Yes",UseReferenceState))
   {
      int sz = Input.getArrayLength(Hash,"ReferenceState",0);
      int* pos = new int[sz];
      Input.getPosIntVector(pos,sz,Hash,"ReferenceState",0);
      Vector Vals(sz);
      Input.getVector(Vals,Hash,"ReferenceState",1);
      for (int i=0;i<sz;++i)
      {
         ReferenceState_[pos[i]] = Vals[i];
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
Vector const& RestrictToTranslatedSubSpaceOld::DrDt(Vector const& Diff) const
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
void RestrictToTranslatedSubSpaceOld::UpdateLatticeState()
{
   for (int i=0;i<size_static;++i)
   {
      DOF_static[i] = ReferenceState_[i];
   }
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         DOF_static[DOFindex_[i][j]] += DOFMult_[i][j]*DOF_[i];
      }
   }
   
   Lattice_->SetDOF(DOF_static);
   Lattice_->SetLoadParameter(ReferenceState_[size_static]+DOF_[DOFS_]);
}

Vector const& RestrictToTranslatedSubSpaceOld::Force() const
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

Matrix const& RestrictToTranslatedSubSpaceOld::Stiffness() const
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

Vector RestrictToTranslatedSubSpaceOld::RestrictDOF(Vector const& dof)
{
   if (dof.Dim() == DOFS_+1)
   {
      return dof;
   }
   else
   {
      cerr << "Error. " << Name() << " does not know how to RestrictDOF()\n";
      exit(-22);
   }
}

Vector RestrictToTranslatedSubSpaceOld::UnRestrictDOF(Vector const& dof)
{
   Vector UnRestricted(size_static+1);
   for (int i=0;i<size_static+1;++i)
   {
      UnRestricted[i] = ReferenceState_[i];
   }
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         UnRestricted[DOFindex_[i][j]] += DOFMult_[i][j]*dof[i];
      }
   }

   return UnRestricted;
}

Vector RestrictToTranslatedSubSpaceOld::TransformVector(Vector const& T)
{
   if (T.Dim() == DOFS_+1)
   {
      return T;
   }
   else
   {
      cerr << "Error. " << Name() << " does not know how to TransformVector()\n";
      exit(-22);
   }
}

Vector RestrictToTranslatedSubSpaceOld::UnTransformVector(Vector const& T)
{
   Vector UnTransformed(size_static+1);
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         UnTransformed[DOFindex_[i][j]] += DOFMult_[i][j]*T[i];
      }
   }

   return UnTransformed;
}

void RestrictToTranslatedSubSpaceOld::SetDOF(Vector const& dof)
{
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] = dof[i];
   }
   
   UpdateLatticeState();
}

void RestrictToTranslatedSubSpaceOld::UpdateDOF(Vector const& dr)
{
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] += dr[i];
   }
   
   UpdateLatticeState();
}
