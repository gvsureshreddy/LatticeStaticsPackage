#include "RestrictToTranslatedSubSpace.h"
#include <sstream>

RestrictToTranslatedSubSpace::~RestrictToTranslatedSubSpace()
{
   if (SymmetryCheckCount_ > 0) delete [] SymmetryCheck_;
   
   cout.width(0);
   cout << "RestrictToTranslatedSubSpace Function Calls:\n"
        << "\tUpdateLatticeState - " << counter_[UPDATE] << "\n"
        << "\tEnergy - " << counter_[ENERGY] << "\n"
        << "\tDrDt - " << counter_[DRDT] << "\n"
        << "\tForce - " << counter_[FORCE] << "\n"
        << "\tStiffness - " << counter_[STIFFNESS] << "\n"
        << "\tDOF - " << counter_[DOFCount] << "\n"
        << "\tSymmetryOK - " << counter_[SYMMETRY] << "\n"
        << "\tRestrictDOF - " << counter_[RESTRICT] << "\n"
        << "\tUnRestrictDOF - " << counter_[UNRESTRICT] << "\n"
        << "\tTransformVector - " << counter_[TRANSFORM] << "\n"
        << "\tUnTransformVector - " << counter_[UNTRANSFORM] << "\n"
        << "\tSetDOF - " << counter_[SETDOF] << "\n"
        << "\tUpdateDOF - " << counter_[UPDATEDOF] << "\n";
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

   if (Input.ParameterOK(Hash,"SymmetryCheckProjectionMatrices"))
   {
      SymmetryCheckCount_ = Input.getArrayLength(Hash,"SymmetryCheckProjectionMatrices");
      if (SymmetryCheckCount_ == 0)
      {
         cerr << "Error. " << Name()
              << " SymmetryCheckProjectionMatrices is empty\n";
         exit(-37);
      }
      
      SymmetryCheck_ = new SparseMatrix[SymmetryCheckCount_];

      for (int i=0;i<SymmetryCheckCount_;++i)
      {
         if (Input.getArrayLength(Hash,"SymmetryCheckProjectionMatrices",i,0) != DOFS_)
         {
            cerr << "Error. " << Name()
                 << " Incorrect number of columns in SymmetryCheckProjectionMatrix"
                 << " number " << i << "\n";
            exit(-38);
         }

         int Rows = Input.getArrayLength(Hash,"SymmetryCheckProjectionMatrices",i);
         Matrix SCPM(Rows,DOFS_);
         Input.getMatrix(SCPM,Hash,"SymmetryCheckProjectionMatrices",i);
         int nononzero = 0;
         for (int j=0;j<Rows;++j)
         {
            for (int k=0;k<DOFS_;++k)
            {
               if (fabs(SCPM[j][k]) > 1.0e-15) ++nononzero;
            }
         }
         
         SymmetryCheck_[i].Resize(Rows,DOFS_+1,nononzero); // DOFS_+1 to ignore load value
         
         int count=0;
         for (int j=0;j<Rows;++j)
         {
            for (int k=0;k<DOFS_;++k)
            {
               if (fabs(SCPM[j][k]) > 1.0e-15)
               {
                  SymmetryCheck_[i].SetNonZeroEntry(count,j,k,SCPM[j][k]);
                  ++count;
               }
            }
         }
      }
   }
   else
   {
      SymmetryCheckCount_ = 0;
   }

   if (Input.ParameterOK(Hash,"SymmetryCheckTolerance"))
   {
      SymmetryCheckTol_ = Input.getDouble(Hash,"SymmetryCheckTolerance");
   }
   else
   {
      // Default Value
      SymmetryCheckTol_ = Input.useDouble(1.0e-14,Hash,"SymmetryCheckTolerance");
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
   ++counter_[DRDT];
   
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
   ++counter_[UPDATE];
   
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
   ++counter_[FORCE];
   
   stress_static = Lattice_->E1();
   
   Multiply(force_static,ForceProject_,stress_static);

   return force_static;
}

Matrix const& RestrictToTranslatedSubSpace::Stiffness() const
{
   ++counter_[STIFFNESS];
   
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

int RestrictToTranslatedSubSpace::SymmetryOK() const
{
   ++counter_[SYMMETRY];

   int retval = SymmetryCheckCount_;
   for (int i=0;i<SymmetryCheckCount_;++i)
   {
      if ( !((SymmetryCheck_[i]*DOF()).Norm() < SymmetryCheckTol_) )
      {
         --retval;
      }
   }

   return !retval;
}

Vector RestrictToTranslatedSubSpace::RestrictDOF(Vector const& dof)
{
   ++counter_[RESTRICT];
   
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
   ++counter_[UNRESTRICT];
   
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
   ++counter_[TRANSFORM];
   
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
   ++counter_[UNTRANSFORM];
   
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
   ++counter_[SETDOF];
   
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] = dof[i];
   }
   
   UpdateLatticeState();
}

void RestrictToTranslatedSubSpace::UpdateDOF(Vector const& dr)
{
   ++counter_[UPDATE];
   
   for (int i=0;i<=DOFS_;++i)
   {
      DOF_[i] += dr[i];
   }
   
   UpdateLatticeState();
}
