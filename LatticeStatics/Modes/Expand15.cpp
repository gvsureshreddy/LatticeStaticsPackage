#include "Expand15.h"

#define DOFS 15

Expand15::Expand15(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector Expand15::ArcLenRHS(double DS,const Vector &Diff,
			   double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - Diff[1]*Diff[1]/(Aspect*Aspect) - Diff[0]*Diff[0];

   return rhs;
}

Vector Expand15::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DOF()[0];
   def[1] = Lattice_->Temp();

   return def;
}

void Expand15::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval[0];
   for (int i=3;i<15;i++)
      DOF[i] = 0.0;
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double Expand15::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector Expand15::DrDt(const Vector &Diff)
{
   Vector ddt(DOFS,0.0);

   ddt[0] = Diff[0]/Diff[1];
   ddt[1] = ddt[2] = ddt[0];

   return ddt;
}

Matrix Expand15::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];
   K[0][1] = Lattice_->StressDT()[0][0];
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double Expand15::ScanningDefParameter()
{
   return Lattice_->DOF()[0];
}

void Expand15::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval;

   Lattice_->SetDOF(DOF);
}

double Expand15::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void Expand15::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double Expand15::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector Expand15::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector Expand15::ScanningDef()
{
   return Vector(1,0.0);
}

void Expand15::ScanningUpdate(const Vector &newval) {}

Matrix Expand15::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
