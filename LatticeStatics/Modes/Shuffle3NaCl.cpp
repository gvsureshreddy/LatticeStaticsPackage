#include "Shuffle3NaCl.h"

Shuffle3NaCl::Shuffle3NaCl(Lattice *M)
{
   Lattice_ = (GenericLat*) M;
}

// Functions required by LatticeMode
Vector Shuffle3NaCl::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - (Diff[1]*Diff[1])/(Aspect*Aspect)
      - Diff[0]*Diff[0];

   return rhs;
}

Vector Shuffle3NaCl::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DOF()[0];
   def[1] = Lattice_->Temp();

   return def;
}

void Shuffle3NaCl::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] = DOF[2] = DOF[0];
   DOF[3] = DOF[0]/4.0;
   DOF[4] = DOF[5] = DOF[3];
   DOF[6]=DOF[7]=DOF[8] = 0.0;

   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double Shuffle3NaCl::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix Shuffle3NaCl::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2] +
      (Stiff[0][3] + Stiff[0][4] + Stiff[0][5])/4.0;
   K[0][1] = Lattice_->StressDT()[0][0];
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double Shuffle3NaCl::ScanningDefParameter()
{
   return Lattice_->DOF()[0];
}

void Shuffle3NaCl::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -=newval;
   DOF[1] = DOF[2] = DOF[0];
   DOF[3] = DOF[4] = DOF[5] = DOF[0]/4.0;

   Lattice_->SetDOF(DOF);
}

double Shuffle3NaCl::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void Shuffle3NaCl::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double Shuffle3NaCl::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector Shuffle3NaCl::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector Shuffle3NaCl::ScanningDef()
{
   return Vector(1,0.0);
}

void Shuffle3NaCl::ScanningUpdate(const Vector &newval) {}

Matrix Shuffle3NaCl::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
