#include "MultiMode.h"
#include "UtilityFunctions.h"

MultiMode::MultiMode(Lattice *M,const char *datafile,const char *prefix)
{
   char tmp[LINELENGTH];
   
   if (!GetParameter(prefix,"MultiMode_DOFS",datafile,"%u",&DOFS)) exit(-1);
   
   for (int i=0;i<DOFS;++i)
   {
      sprintf(tmp,"MultiMode_DOF_%u_Len",i);
      if (!GetParameter(prefix,tmp,datafile,"%u",&(DOFindlen[i]))) exit(-1);
      sprintf(tmp,"MultiMode_DOF_%u_Val",i);
      if (!GetIntVectorParameter(prefix,tmp,datafile,DOFindlen[i],DOFindex[i])) exit(-1);
      DOFMult[i].Resize(DOFMAX);
      sprintf(tmp,"MultiMode_DOF_%u_Mul",i);
      if (!GetVectorParameter(prefix,tmp,datafile,&(DOFMult[i]))) exit(-1);
   }

   if (!GetParameter(prefix,"MultiMode_ScnDefParam",datafile,"%u",&ScnDefParam))
      exit(-1);
   if ((ScnDefParam < 0) || (ScnDefParam >= DOFS))
   {
      cerr << "MultiMode: ScnDefParam too small or too big!" << endl;
      exit(-1);
   }

   Lattice_ = (Lattice *) M;
}

// Functions required by LatticeMode
Vector MultiMode::ArcLenRHS(double DS,const Vector &Diff,
			  double Aspect)
{
   static Vector rhs(DOFS+1);
   static Matrix stress(1,(Lattice_->DOF()).Dim());

   stress = Lattice_->Stress();

   for (int i=0;i<DOFS;++i)
   {
      rhs[i] = 0.0;
      for (int j=0;j<DOFindlen[i];j++)
      {
	 rhs[i] += DOFMult[i][j]*stress[0][DOFindex[i][j]];
      }
   }
   
   rhs[DOFS] = DS*DS - Diff[DOFS]*Diff[DOFS]/(Aspect*Aspect);
   for (int i=0;i<DOFS;++i)
   {
      rhs[DOFS] -= Diff[i]*Diff[i];
   }

   return rhs;
}

Vector MultiMode::ArcLenDef()
{
   static Vector def(DOFS+1);
   static Vector dof((Lattice_->DOF()).Dim());
   dof = Lattice_->DOF();

   for (int i=0;i<DOFS;++i)
   {
      def[i] = dof[DOFindex[i][0]]/DOFMult[i][0];
   }

   if (Lattice_->LoadParameter()==Lattice::Temperature)
      def[DOFS] = Lattice_->Temp();
   else if (Lattice_->LoadParameter()==Lattice::Load)
      def[DOFS] = Lattice_->Lambda();

   
   return def;
}

void MultiMode::ArcLenSet(const Vector &val)
{
   Vector DOF((Lattice_->DOF()).Dim(),0.0);
   
   for (int i=0;i<DOFS;++i)
   {
      for (int j=0;j<DOFindlen[i];++j)
      {
	 DOF[DOFindex[i][j]] += DOFMult[i][j]*val[i];
      }
   }

   Lattice_->SetDOF(DOF);
   if (Lattice_->LoadParameter()==Lattice::Temperature)
      Lattice_->SetTemp(val[DOFS]);
   else if (Lattice_->LoadParameter()==Lattice::Load)
      Lattice_->SetLambda(val[DOFS]);
}
   
void MultiMode::ArcLenUpdate(const Vector &newval)
{
   static Vector DOF((Lattice_->DOF()).Dim());

   DOF = Lattice_->DOF();

   for (int i=0;i<DOFS;++i)
   {
      for (int j=0;j<DOFindlen[i];++j)
      {
	 DOF[DOFindex[i][j]] -= DOFMult[i][j]*newval[i];
      }
   }

   Lattice_->SetDOF(DOF);
   if (Lattice_->LoadParameter()==Lattice::Temperature)
      Lattice_->SetTemp(Lattice_->Temp() - newval[DOFS]);
   else if (Lattice_->LoadParameter()==Lattice::Load)
      Lattice_->SetLambda(Lattice_->Lambda() - newval[DOFS]);
}

double MultiMode::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[DOFS] /= Aspect;
   New[DOFS] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector MultiMode::DrDt(const Vector &Diff)
{
   Vector ddt((Lattice_->DOF()).Dim(),0.0);

   for (int i=0;i<DOFS;++i)
   {
      for (int j=0;j<DOFindlen[i];++j)
      {
	 ddt[DOFindex[i][j]] += DOFMult[i][j]*(Diff[i]/Diff[DOFS]);
      }
   }
   
   return ddt;
}

Matrix MultiMode::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   static Matrix K(DOFS+1,DOFS+1);
   static Matrix Stiff((Lattice_->DOF()).Dim(),(Lattice_->DOF()).Dim());
   static Matrix stressdt(1,(Lattice_->DOF()).Dim());

   K.Resize(DOFS+1,DOFS+1,0.0);
   Stiff = Lattice_->Stiffness();
   if (Lattice_->LoadParameter()==Lattice::Temperature)
      stressdt = Lattice_->StressDT();
   else if (Lattice_->LoadParameter()==Lattice::Load)
      stressdt = Lattice_->StressDL();

   for (int i=0;i<DOFS;++i)
   {
      for (int j=0;j<DOFindlen[i];++j)
      {
	 for (int k=0;k<DOFS;++k)
	 {
	    for (int l=0;l<DOFindlen[k];++l)
	    {
	       K[i][k] +=
		  DOFMult[i][j]*(Stiff[DOFindex[i][j]][DOFindex[k][l]])*DOFMult[k][l];
	    }
	 }

	 K[i][DOFS] += DOFMult[i][j]*stressdt[0][DOFindex[i][j]];
      }
   }

   for (int i=0;i<DOFS;++i)
   {
      K[DOFS][i] = -2.0*Diff[i];
   }
   K[DOFS][DOFS] = -2.0*Diff[DOFS]/(Aspect*Aspect);

   return K;
}

double MultiMode::ScanningDefParameter()
{
   Vector DOF = Lattice_->DOF();

   return DOF[DOFindex[ScnDefParam][0]]/DOFMult[ScnDefParam][0];
}

void MultiMode::ScanningDefParamSet(const double val)
{
   Vector DOF=Lattice_->DOF();

   for (int i=0;i<DOFindlen[ScnDefParam];++i)
   {
      DOF[DOFindex[ScnDefParam][i]] += DOFMult[ScnDefParam][i]*val;
   }

   Lattice_->SetDOF(DOF);
}

void MultiMode::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   for (int i=0;i<DOFindlen[ScnDefParam];++i)
   {
      DOF[DOFindex[ScnDefParam][i]] -= DOFMult[ScnDefParam][i]*newval;
   }

   Lattice_->SetDOF(DOF);
}

double MultiMode::ScanningLoadParameter()
{
   if (Lattice_->LoadParameter()==Lattice::Temperature)
      return Lattice_->Temp();
   else if (Lattice_->LoadParameter()==Lattice::Load)
      return Lattice_->Lambda();
}

void MultiMode::ScanningLoadParamSet(const double val)
{
   if (Lattice_->LoadParameter()==Lattice::Temperature)
      Lattice_->SetTemp(val);
   else if (Lattice_->LoadParameter()==Lattice::Load)
      Lattice_->SetLambda(val);
}

void MultiMode::ScanningLoadParamUpdate(const double newval)
{
   if (Lattice_->LoadParameter()==Lattice::Temperature)
      Lattice_->SetTemp(Lattice_->Temp() - newval);
   else if (Lattice_->LoadParameter()==Lattice::Load)
      Lattice_->SetLambda(Lattice_->Lambda() - newval);
}

double MultiMode::ScanningStressParameter()
{
   double str = 0.0;
   Matrix stress = Lattice_->Stress();

   for (int i=0;i<DOFindlen[ScnDefParam];++i)
   {
      str += DOFMult[ScnDefParam][i]*stress[0][DOFindex[ScnDefParam][i]];
   }
   
   return str;
}
   
Vector MultiMode::ScanningRHS()
{
   Matrix stress = Lattice_->Stress();
   Vector RHS(DOFS-1,0.0);
   int a=0;

   if (DOFS != 1)
   {
      for (int i=0;i<DOFS;++i)
      {
	 if (i != ScnDefParam)
	 {
	    for (int j=0;j<DOFindlen[i];++j)
	    {
	       RHS[a] += DOFMult[i][j]*stress[0][DOFindex[i][j]];
	    }
	    ++a;
	 }
      }

      return RHS;
   }
   else
      return Vector(1,0.0);

}

Vector MultiMode::ScanningDef()
{
   Vector DOF = Lattice_->DOF();
   Vector DEF(DOFS-1,0.0);
   int a=0;

   if (DOFS != 1)
   {
      for (int i=0;i<DOFS;++i)
      {
	 if (i != ScnDefParam)
	 {
	    DEF[a] = DOF[DOFindex[i][0]]/DOFMult[i][0];
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
   Vector DOF((Lattice_->DOF()).Dim(),0.0);
   
   for (int i=0;i<DOFS;++i)
   {
      if (i != ScnDefParam)
      {
	 for (int j=0;j<DOFindlen[i];++j)
	 {
	    DOF[DOFindex[i][j]] += DOFMult[i][j]*val[i];
	 }
      }
   }

   Lattice_->SetDOF(DOF);
}

void MultiMode::ScanningUpdate(const Vector &newval)
{
   Vector dof=Lattice_->DOF();

   for (int i=0;i<DOFS;++i)
   {
      if (i != ScnDefParam)
      {
	 for (int j=0;j<DOFindlen[i];++j)
	 {
	    dof[DOFindex[i][j]] -= DOFMult[i][j]*newval[i>ScnDefParam?i-1:i];
	 }
      }
   }
   
   Lattice_->SetDOF(dof);
}

Matrix MultiMode::ScanningStiffness()
{
   Matrix stiff = Lattice_->Stiffness();
   Matrix K(DOFS-1,DOFS-1,0.0);
   int a=0,b=0;

   if (DOFS != 1)
   {
      for (int i=0;i<DOFS;++i)
      {
	 if (i != ScnDefParam)
	 {
	    for (int j=0;j<DOFindlen[i];++j)
	    {
	       b=0;
	       for (int k=0;k<DOFS;++k)
	       {
		  if (k != ScnDefParam)
		  {
		     for (int l=0;l<DOFindlen[k];++l)
		     {
			K[a][b] += DOFMult[i][j]
			   *(stiff[DOFindex[i][j]][DOFindex[k][l]])*DOFMult[k][l];
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
