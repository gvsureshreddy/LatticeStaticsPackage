#include "RestrictToTranslatedSubSpace.h"
#include <sstream>

RestrictToTranslatedSubSpace::~RestrictToTranslatedSubSpace()
{
   cout.width(0);
   cout << "RestrictToTranslatedSubSpace Function Calls:\n"
        << "\tUpdateLatticeState - " << counter_[0] << "\n"
        << "\tEnergy - " << counter_[1] << "\n"
        << "\tDrDt - " << counter_[2] << "\n"
        << "\tForce - " << counter_[3] << "\n"
        << "\tStiffness - " << counter_[4] << "\n"
        << "\tDOF - " << counter_[5] << "\n"
        << "\tRestrictDOF - " << counter_[6] << "\n"
        << "\tUnRestrictDOF - " << counter_[7] << "\n"
        << "\tTransformVector - " << counter_[8] << "\n"
        << "\tUnTransformVector - " << counter_[9] << "\n"
        << "\tSetDOF - " << counter_[10] << "\n"
        << "\tUpdateDOF - " << counter_[11] << "\n";
}


RestrictToTranslatedSubSpace::RestrictToTranslatedSubSpace(Lattice* const M,PerlInput const& Input):
   ForceProjectionMatrix_(),
   DOFProjectionMatrix_(),
   ForceProject_(ForceProjectionMatrix_),
   DOFProject_(DOFProjectionMatrix_)
{
   for (int i=0;i<nocounters_;++i) counter_[i] = 0;
   
   stringstream tmp;

   Lattice_ = (Lattice *) M;

   PerlInput::HashStruct Hash = Input.getHash("Restriction","RestrictToTranslatedSubSpace");

   int LatDOFS=Lattice_->DOF().Dim();
   if (Input.ParameterOK(Hash,"ProjectionMatrix"))
   {
      DOFS_ = Input.getArrayLength(Hash,"ProjectionMatrix");
      DOF_.Resize(DOFS_+1,0.0);

      if (Input.getArrayLength(Hash,"ProjectionMatrix",0) != LatDOFS)
      {
         cerr << "Error. " << Name() << " Incorrect number of columns in ProjectionMatrix\n";
         exit(-37);
      }
      
      Matrix PM(DOFS_,LatDOFS);
      Input.getMatrix(PM,Hash,"ProjectionMatrix");
      int nononzero = 0;
      for (int i=0;i<DOFS_;++i)
      {
         for (int j=0;j<LatDOFS;++j)
         {
            if (fabs(PM[i][j]) > 1.0e-15) ++nononzero;
         }
      }

      ForceProjectionMatrix_.Resize(DOFS_,LatDOFS,nononzero);
      DOFProjectionMatrix_.Resize(LatDOFS+1,DOFS_+1,nononzero+1);
      
      int count=0;
      for (int i=0;i<DOFS_;++i)
      {
         for (int j=0;j<LatDOFS;++j)
         {
            if (fabs(PM[i][j]) > 1.0e-15)
            {
               ForceProjectionMatrix_.SetNonZeroEntry(count,i,j,PM[i][j]);
               DOFProjectionMatrix_.SetNonZeroEntry(count,j,i,PM[i][j]);
               ++count;
            }
         }
      }
      DOFProjectionMatrix_.SetNonZeroEntry(nononzero,LatDOFS,DOFS_,1.0);
   }
   else
   {
      DOFS_ = Input.getPosInt(Hash,"DOFS");
      DOF_.Resize(DOFS_+1,0.0);
      
      Vector* const Values = new Vector[DOFS_];
      int** const Positions = new int*[DOFS_];
      int len;
      int nononzero=0;
      
      for (int i=0;i<DOFS_;++i)
      {
         tmp.str("");
         tmp << "DOF_" << i;
         len = Input.getArrayLength(Hash,tmp.str().c_str(),0);
         nononzero += len;
         Positions[i] = new int[len];
         Input.getIntVector(Positions[i],len,Hash,tmp.str().c_str(),0);
         Values[i].Resize(len);
         Input.getVector(Values[i],Hash,tmp.str().c_str(),1);
         Values[i] /= Values[i].Norm(); // normalize the rows
      }
      ForceProjectionMatrix_.Resize(DOFS_,LatDOFS,nononzero);
      DOFProjectionMatrix_.Resize(LatDOFS+1,DOFS_+1,nononzero+1);
      
      int count=0;
      for (int i=0;i<DOFS_;++i)
      {
         for (int j=0;j<Values[i].Dim();++j)
         {
            ForceProjectionMatrix_.SetNonZeroEntry(count,i,Positions[i][j],Values[i][j]);
            DOFProjectionMatrix_.SetNonZeroEntry(count,Positions[i][j],i,Values[i][j]);
            ++count;
         }
         delete [] Positions[i];
      }
      delete [] Positions;
      delete [] Values;
      DOFProjectionMatrix_.SetNonZeroEntry(nononzero,LatDOFS,DOFS_,1.0);
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
   Lat_DOF_static.Resize(size_static);
   Rest_DOF_static.Resize(DOFS_);
   force_static.Resize(DOFS_);
   stress_static.Resize(size_static);
   K_static.Resize(DOFS_,DOFS_+1);
   E2_tmp_static.Resize(size_static,size_static+1);
   Stiff_static.Resize(size_static,size_static);
   stressdt_static.Resize(size_static);
   Lat_ddt_static.Resize(size_static);
   Rest_ddt_static.Resize(DOFS_);

   // Make sure everything is synchronized
   Vector t = Lattice_->DOF();
   Vector tt(t.Dim()+1);
   for (int i=0;i<t.Dim();++i) tt[i] = t[i];
   tt[t.Dim()] = ((Lattice::Temperature == Lattice_->LoadParameter()) ?
                  Lattice_->Temp() : Lattice_->Lambda());
   SetDOF(RestrictDOF(tt));
}

// Functions required by Restriction
Vector const& RestrictToTranslatedSubSpace::DrDt(Vector const& Diff) const
{
   ++counter_[2];
   
   for (int i=0;i<DOFS_;++i)
   {
      Rest_ddt_static[i] = Diff[i]/Diff[DOFS_];
   }
   Multiply(Lat_ddt_static,Rest_ddt_static,ForceProject_);
   
   return Lat_ddt_static;
}

//----------------------------------------------------------------
void RestrictToTranslatedSubSpace::UpdateLatticeState()
{
   ++counter_[0];
   
   for (int i=0;i<DOFS_;++i)
   {
      Rest_DOF_static[i] = DOF_[i];
   }
   Multiply(Lat_DOF_static,Rest_DOF_static,ForceProject_);
   for (int i=0;i<size_static;++i)
   {
      Lat_DOF_static[i] += ReferenceState_[i];
   }
   
   Lattice_->SetDOF(Lat_DOF_static);
   Lattice_->SetLoadParameter(ReferenceState_[size_static]+DOF_[DOFS_]);
}

Vector const& RestrictToTranslatedSubSpace::Force() const
{
   ++counter_[3];
   
   stress_static = Lattice_->E1();
   
   Multiply(force_static,ForceProject_,stress_static);

   return force_static;
}

Matrix const& RestrictToTranslatedSubSpace::Stiffness() const
{
   ++counter_[4];
   
   Stiff_static = Lattice_->E2();
   stressdt_static = Lattice_->E1DLoad();
   for (int i=0;i<size_static;++i)
   {
      for (int j=0;j<size_static;++j)
      {
         E2_tmp_static[i][j] = Stiff_static[i][j];
      }
      E2_tmp_static[i][size_static] = stressdt_static[i];
   }

   Multiply(K_static,ForceProject_,E2_tmp_static,DOFProject_);
   
   return K_static;
}

Vector RestrictToTranslatedSubSpace::RestrictDOF(Vector const& dof)
{
   ++counter_[6];
   
   if (dof.Dim() == DOFS_+1)
   {
      return dof;
   }
   else if (dof.Dim() == size_static+1)
   {
      Vector Restricted(DOFS_+1);
      Vector Translated(size_static+1);
      for (int i=0;i<size_static+1;++i)
      {
         Translated[i] = dof[i] - ReferenceState_[i];
      }
      Multiply(Restricted,Translated,DOFProject_);
      return Restricted;
   }
   else
   {
      cerr
         << "Error. " << Name() << " Unhandled dof vector size in RestrictDOF()";
      exit(-31);
   }
}

Vector RestrictToTranslatedSubSpace::UnRestrictDOF(Vector const& dof)
{
   ++counter_[7];
   
   if (dof.Dim() == DOFS_+1)
   {
      Vector UnRestricted(size_static+1);
      Multiply(UnRestricted,DOFProject_,dof);

      for (int i=0;i<size_static+1;++i)
      {
         UnRestricted[i] += ReferenceState_[i];
      }

      return UnRestricted;
   }
   else if (dof.Dim() == size_static+1)
   {
      return dof;
   }
   else
   {
      cerr << "Error. " << Name() << " Unhandled size dof vector in UnRestrictDOF().\n";
      exit(-31);
   }
}

Vector RestrictToTranslatedSubSpace::TransformVector(Vector const& T)
{
   ++counter_[8];
   
   if (T.Dim() == DOFS_+1)
   {
      return T;
   }
   else if (T.Dim() == size_static+1)
   {
      Vector Transformed(DOFS_+1);
      Multiply(Transformed,T,DOFProject_);
      return Transformed;
   }
   else
   {
      cerr << "Error. " << Name() << " Unhandled size T vector in TransformVector().\n";
      exit(-31);
   }
}

Vector RestrictToTranslatedSubSpace::UnTransformVector(Vector const& T)
{
   ++counter_[9];
   
   if (T.Dim() == DOFS_+1)
   {
      Vector UnTransformed(size_static+1);
      Multiply(UnTransformed,DOFProject_,T);
      return UnTransformed;
   }
   else if (T.Dim() == size_static+1)
   {
      return T;
   }
   else
   {
      cerr << "Error. " << Name() << " Unhandled size T vector in UnTransforVector().\n";
      exit(-31);
   }
}

void RestrictToTranslatedSubSpace::SetDOF(Vector const& dof)
{
   ++counter_[10];
   
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] = dof[i];
   }
   
   UpdateLatticeState();
}

void RestrictToTranslatedSubSpace::UpdateDOF(Vector const& dr)
{
   ++counter_[11];
   
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] += dr[i];
   }
   
   UpdateLatticeState();
}
