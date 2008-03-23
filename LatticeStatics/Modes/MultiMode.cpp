#include "MultiMode.h"

MultiMode::MultiMode(Lattice *M,PerlInput &Input)
{
   char tmp[LINELENGTH];

   PerlInput::HashStruct Hash = Input.getHash("Mode","MultiMode");
   DOFS_ = Input.getUnsigned(Hash,"DOFS");
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
   unsigned temp2;
   double temp3;
   unsigned Baseline_DOFS;
   if (Input.ParameterOK(Hash,"Baseline_DOFS"))
   {
      Baseline_DOFS = Input.getUnsigned(Hash,"Baseline_DOFS");
   }
   else
   {
      Baseline_DOFS=0;
   }
   
   BaselineDOF_.Resize(temp1,0.0);
   
   for(unsigned i=0; i < Baseline_DOFS ; ++i)
   {
      sprintf(tmp,"Baseline_DOF_Index_%u",i);
      temp2 = Input.getUnsigned(Hash,tmp);
      sprintf(tmp,"Baseline_DOF_Value_%u",i);
      temp3 = Input.getDouble(Hash,tmp);
      
      BaselineDOF_[temp2] =  temp3;
   }
}

// Functions required by LatticeMode
Vector MultiMode::DrDt(const Vector &Diff)
{
   Vector ddt((Lattice_->DOF()).Dim(),0.0);
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         ddt[DOFindex_[i][j]] += DOFMult_[i][j]*(Diff[i]/Diff[DOFS_]);
      }
   }
   
   return ddt;
}

//----------------------------------------------------------------
void MultiMode::UpdateLatticeState()
{
   static int size = (Lattice_->DOF()).Dim();
   static Vector DOF(size);
   for (int i=0;i<size;++i)
   {
      DOF[i] = BaselineDOF_[i];
   }
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         DOF[DOFindex_[i][j]] += DOFMult_[i][j]*ModeDOF_[i];
      }
   }
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetLoadParameter(ModeDOF_[DOFS_]);
}

Vector MultiMode::ModeForce()
{
   static Vector force(DOFS_);
   static Matrix stress(1,(Lattice_->DOF()).Dim());
   
   stress = Lattice_->E1();
   
   for (int i=0;i<DOFS_;++i)
   {
      force[i] = 0.0;
      for (int j=0;j<DOFindlen_[i];j++)
      {
         force[i] += DOFMult_[i][j]*stress[0][DOFindex_[i][j]];
      }
   }
   
   return force;
}

Matrix MultiMode::ModeStiffness()
{
   static Matrix K(DOFS_,DOFS_+1);
   static Matrix Stiff((Lattice_->DOF()).Dim(),(Lattice_->DOF()).Dim());
   static Matrix stressdt(1,(Lattice_->DOF()).Dim());
   
   K.Resize(DOFS_,DOFS_+1,0.0);
   Stiff = Lattice_->E2();
   stressdt = Lattice_->E1DLoad();
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<DOFindlen_[i];++j)
      {
         for (int k=0;k<DOFS_;++k)
         {
            for (int l=0;l<DOFindlen_[k];++l)
            {
               K[i][k] +=
                  DOFMult_[i][j]*(Stiff[DOFindex_[i][j]][DOFindex_[k][l]])*DOFMult_[k][l];
            }
         }
         
         K[i][DOFS_] += DOFMult_[i][j]*stressdt[0][DOFindex_[i][j]];
      }
   }
   
   return K;
}

void MultiMode::SetModeDOF(const Vector &dof)
{
   for (int i=0;i<=DOFS_;++i)
   {
      ModeDOF_[i] = dof[i];
   }
   
   UpdateLatticeState();
}

void MultiMode::UpdateModeDOF(const Vector &dr)
{
   for (int i=0;i<=DOFS_;++i)
   {
      ModeDOF_[i] += dr[i];
   }
   
   UpdateLatticeState();
}
