#include "ML_Expand.h"

ML_Expand::ML_Expand(Lattice *M)
{
   Lattice_ = (Lattice *) M;
}

// Functions required by LatticeMode
Vector ML_Expand::ArcLenRHS(double DS,const Vector &Diff,
			  double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - Diff[1]*Diff[1]/(Aspect*Aspect) - Diff[0]*Diff[0];

   return rhs;
}

Vector ML_Expand::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DOF()[0];
   def[1] = Lattice_->Temp();

   return def;
}

void ML_Expand::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval[0];
   for (int i=3;i<DOF.Dim();i++)
      DOF[i] = 0.0;
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double ML_Expand::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector ML_Expand::DrDt(const Vector &Diff)
{
   Vector ddt((Lattice_->DOF()).Dim(),0.0);

   ddt[0] = Diff[0]/Diff[1];
   ddt[1] = ddt[2] = ddt[0];

   return ddt;
}

Matrix ML_Expand::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];
   K[0][1] = Lattice_->StressDT()[0][0];
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double ML_Expand::ScanningDefParameter()
{
   return Lattice_->DOF()[0];
}

void ML_Expand::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval;

   Lattice_->SetDOF(DOF);
}

double ML_Expand::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void ML_Expand::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double ML_Expand::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector ML_Expand::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector ML_Expand::ScanningDef()
{
   return Vector(1,0.0);
}

void ML_Expand::ScanningUpdate(const Vector &newval) {}

Matrix ML_Expand::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
