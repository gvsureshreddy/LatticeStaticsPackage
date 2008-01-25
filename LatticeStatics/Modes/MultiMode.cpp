#include "MultiMode.h"
#include "UtilityFunctions.h"

MultiMode::MultiMode(Lattice *M,const char *datafile,const char *prefix)
{
   char tmp[LINELENGTH];
   
   if (!GetParameter(prefix,"MultiMode_DOFS",datafile,'u',&DOFS_)) exit(-1);
   ModeDOF_.Resize(DOFS_+1,0.0);
   
   for (int i=0;i<DOFS_;++i)
   {
      sprintf(tmp,"MultiMode_DOF_%u_Len",i);
      if (!GetParameter(prefix,tmp,datafile,'u',&(DOFindlen_[i]))) exit(-1);
      sprintf(tmp,"MultiMode_DOF_%u_Val",i);
      if (!GetIntVectorParameter(prefix,tmp,datafile,DOFindlen_[i],DOFindex_[i])) exit(-1);
      DOFMult_[i].Resize(DOFindlen_[i]);
      sprintf(tmp,"MultiMode_DOF_%u_Mul",i);
      if (!GetVectorParameter(prefix,tmp,datafile,&(DOFMult_[i]))) exit(-1);
   }

   if (!GetParameter(prefix,"MultiMode_ScnDefParam",datafile,'u',&ScnDefParam_))
      exit(-1);
   if ((ScnDefParam_ < 0) || (ScnDefParam_ >= DOFS_))
   {
      cerr << "MultiMode: ScnDefParam too small or too big!" << endl;
      exit(-1);
   }
   
   Lattice_ = (Lattice *) M;
   
   //Baseline DOF Initialization 
   int temp1 = (Lattice_->DOF()).Dim();
   int temp2;
   double temp3;
   int Baseline_DOFS;
   if (!GetParameter(prefix,"MultiMode_Baseline_DOFS",datafile,'u',&Baseline_DOFS)) exit(-1);
   BaselineDOF_.Resize(temp1,0.0);
   
   for(int i=0; i< Baseline_DOFS ; ++i)
   {
      sprintf(tmp,"MultiMode_Baseline_DOF_Index_%u",i);
      if (!GetParameter(prefix,tmp,datafile,'u',&temp2)) exit(-1);
      sprintf(tmp,"MultiMode_Baseline_DOF_Value_%u",i);
      if (!GetParameter(prefix,tmp,datafile,'l',&temp3)) exit(-1);
      
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

//----------------------------------------------------------------
Vector MultiMode::ArcLenForce(double DS,const Vector &Diff,
			      double Aspect)
{
   static Vector force(DOFS_+1);
   static Vector mdfc(DOFS_);
   mdfc = ModeForce();
   
   force[DOFS_] = DS*DS - Diff[DOFS_]*Diff[DOFS_]/(Aspect*Aspect);
   for (int i=0;i<DOFS_;++i)
   {
      force[i] = mdfc[i];
      force[DOFS_] -= Diff[i]*Diff[i];
   }

   return force;
}

Matrix MultiMode::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   static Matrix K(DOFS_+1,DOFS_+1);
   static Matrix ModeK(DOFS_,DOFS_+1);

   ModeK = ModeStiffness();
   
   for (int i=0;i<DOFS_;++i)
   {
      for (int j=0;j<=DOFS_;++j)
      {
	 K[i][j] = ModeK[i][j];
      }
      K[DOFS_][i] = -2.0*Diff[i];
   }
   K[DOFS_][DOFS_] = -2.0*Diff[DOFS_]/(Aspect*Aspect);
   
   return K;
}

double MultiMode::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[DOFS_] /= Aspect;
   New[DOFS_] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

//----------------------------------------------------------------
double MultiMode::ScanningDefParameter()
{
   return ModeDOF_[ScnDefParam_];
}

void MultiMode::ScanningDefParamSet(const double val)
{
   ModeDOF_[ScnDefParam_] = val;

   UpdateLatticeState();
}

void MultiMode::ScanningDefParamUpdate(const double newval)
{
   ModeDOF_[ScnDefParam_] += newval;

   UpdateLatticeState();
}

double MultiMode::ScanningLoadParameter()
{
   return ModeDOF_[DOFS_];
}

void MultiMode::ScanningLoadParamSet(const double val)
{
   ModeDOF_[DOFS_] = val;
   Lattice_->SetLoadParameter(ModeDOF_[DOFS_]);
}

void MultiMode::ScanningLoadParamUpdate(const double newval)
{
   ModeDOF_[DOFS_] += newval;
   Lattice_->SetLoadParameter(ModeDOF_[DOFS_]);
}

double MultiMode::ScanningStressParameter()
{
   double str = 0.0;
   Matrix stress = Lattice_->E1();

   for (int i=0;i<DOFindlen_[ScnDefParam_];++i)
   {
      str += DOFMult_[ScnDefParam_][i]*stress[0][DOFindex_[ScnDefParam_][i]];
   }
   
   return str;
}
   
Vector MultiMode::ScanningForce()
{
   Matrix stress = Lattice_->E1();
   Vector force(DOFS_-1,0.0);
   int a=0;

   if (DOFS_ != 1)
   {
      for (int i=0;i<DOFS_;++i)
      {
	 if (i != ScnDefParam_)
	 {
	    for (int j=0;j<DOFindlen_[i];++j)
	    {
	       force[a] += DOFMult_[i][j]*stress[0][DOFindex_[i][j]];
	    }
	    ++a;
	 }
      }

      return force;
   }
   else
      return Vector(1,0.0);

}

Vector MultiMode::ScanningDef()
{
   Vector DEF(DOFS_-1,0.0);
   int a=0;

   if (DOFS_ != 1)
   {
      for (int i=0;i<DOFS_;++i)
      {
	 if (i != ScnDefParam_)
	 {
	    DEF[a] = ModeDOF_[i];
	    ++a;
	 }
      }
      
      return DEF;
   }
   else
      return Vector(1,0.0);
}

void MultiMode::ScanningSet(const Vector &val)
{
   for (int i=0;i<DOFS_;++i)
   {
      if (i != ScnDefParam_)
      {
	 ModeDOF_[i] = val[i];
      }
   }

   UpdateLatticeState();
}

void MultiMode::ScanningUpdate(const Vector &newval)
{
   for (int i=0;i<DOFS_;++i)
   {
      if (i != ScnDefParam_)
      {
	 ModeDOF_[i] += newval[(i>ScnDefParam_)?i-1:i];
      }
   }
   
   UpdateLatticeState();
}

Matrix MultiMode::ScanningStiffness()
{
   Matrix stiff = Lattice_->E2();
   Matrix K(DOFS_-1,DOFS_-1,0.0);
   int a=0,b=0;

   if (DOFS_ != 1)
   {
      for (int i=0;i<DOFS_;++i)
      {
	 if (i != ScnDefParam_)
	 {
	    for (int j=0;j<DOFindlen_[i];++j)
	    {
	       b=0;
	       for (int k=0;k<DOFS_;++k)
	       {
		  if (k != ScnDefParam_)
		  {
		     for (int l=0;l<DOFindlen_[k];++l)
		     {
			K[a][b] += DOFMult_[i][j]
			   *(stiff[DOFindex_[i][j]][DOFindex_[k][l]])*DOFMult_[k][l];
		     }
		     ++b;
		  }
	       }
	    }
	    ++a;
	 }
      }
      return K;
   }
   else
      return Matrix(1,1,1.0);
}
