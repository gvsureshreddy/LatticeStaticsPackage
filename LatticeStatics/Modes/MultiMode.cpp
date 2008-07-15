#include "MultiMode.h"

MultiMode::MultiMode(Lattice* const M,PerlInput const& Input)
{
   char tmp[LINELENGTH];

   PerlInput::HashStruct Hash = Input.getHash("Mode","MultiMode");
   DOFS_ = Input.getPosInt(Hash,"DOFS");
   ModeDOF_.Resize(DOFS_+1,0.0);
   
   for (int i=0;i<DOFS_;++i)
   {
      sprintf(tmp,"DOF_%u",i);
      DOFindlen_[i] = Input.getArrayLength(Hash,tmp,0);
      Input.getIntVector(DOFindex_[i],DOFindlen_[i],Hash,tmp,0);
      DOFMult_[i].Resize(DOFindlen_[i]);
      Input.getVector(DOFMult_[i],Hash,tmp,1);
   }
   
   Lattice_ = (Lattice *) M;
   
   //Baseline DOF Initialization
   int temp1 = (Lattice_->DOF()).Dim();
   int temp2;
   double temp3;
   int Baseline_DOFS;
   if (Input.ParameterOK(Hash,"Baseline_DOFS"))
   {
      Baseline_DOFS = Input.getPosInt(Hash,"Baseline_DOFS");
   }
   else
   {
      Baseline_DOFS = Input.usePosInt(0,Hash,"Baseline_DOFS");
   }
   
   BaselineDOF_.Resize(temp1,0.0);
   
   for(int i=0; i < Baseline_DOFS ; ++i)
   {
      sprintf(tmp,"Baseline_DOF_Index_%u",i);
      temp2 = Input.getPosInt(Hash,tmp);
      sprintf(tmp,"Baseline_DOF_Value_%u",i);
      temp3 = Input.getDouble(Hash,tmp);
      
      BaselineDOF_[temp2] =  temp3;
   }
   Input.EndofInputSection();
   
   // intitalize "static" member variables
   size_static = (Lattice_->DOF()).Dim();
   DOF_static.Resize(size_static);
   force_static.Resize(DOFS_);
   stress_static.Resize(1,size_static);
   K_static.Resize(DOFS_,DOFS_+1);
   Stiff_static.Resize(size_static,size_static);
   stressdt_static.Resize(1,size_static);

}

// Functions required by LatticeMode
Vector const& MultiMode::DrDt(Vector const& Diff) const
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
void MultiMode::UpdateLatticeState()
{
   for (int i=0;i<size_static;++i)
   {
      DOF_static[i] = BaselineDOF_[i];
   }
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         DOF_static[DOFindex_[i][j]] += DOFMult_[i][j]*ModeDOF_[i];
      }
   }
   
   Lattice_->SetDOF(DOF_static);
   Lattice_->SetLoadParameter(ModeDOF_[DOFS_]);
}

Vector const& MultiMode::ModeForce() const
{
   stress_static = Lattice_->E1();
   
   for (int i=0;i<DOFS_;++i)
   {
      force_static[i] = 0.0;
      for (int j=0;j<DOFindlen_[i];j++)
      {
         force_static[i] += DOFMult_[i][j]*stress_static[0][DOFindex_[i][j]];
      }
   }
   
   return force_static;
}

Matrix const& MultiMode::ModeStiffness() const
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
         
         K_static[i][DOFS_] += DOFMult_[i][j]*stressdt_static[0][DOFindex_[i][j]];
      }
   }
   
   return K_static;
}

void MultiMode::SetModeDOF(Vector const& dof)
{
   for (int i=0;i<=DOFS_;++i)
   {
      ModeDOF_[i] = dof[i];
   }
   
   UpdateLatticeState();
}

void MultiMode::UpdateModeDOF(Vector const& dr)
{
   for (int i=0;i<=DOFS_;++i)
   {
      ModeDOF_[i] += dr[i];
   }
   
   UpdateLatticeState();
}
