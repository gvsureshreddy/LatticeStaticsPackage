#include "Shuffle3Rhombo.h"

Shuffle3Rhombo::Shuffle3Rhombo(Lattice *M)
{
   Lattice_ = (GenericLat*) M;
}

// Functions required by LatticeMode
Vector Shuffle3Rhombo::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(3);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[0][3];
   rhs[2] = DS*DS - (Diff[2]*Diff[2])/(Aspect*Aspect)
      - Diff[0]*Diff[0] - Diff[1]*Diff[1];

   return rhs;
}

Vector Shuffle3Rhombo::ArcLenDef()
{
   Vector def(3);

   def[0] = Lattice_->DOF()[0];
   def[1] = Lattice_->DOF()[3];
   def[2] = Lattice_->Temp();

   return def;
}

void Shuffle3Rhombo::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] = DOF[2] = DOF[0];
   DOF[3] -= newval[1];
   DOF[4] = DOF[5] = DOF[3];
   DOF[6] = DOF[7] = DOF[8] = 0.0;

   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[2]);
}

double Shuffle3Rhombo::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[2] /= Aspect;
   New[2] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix Shuffle3Rhombo::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(3,3);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];
   K[0][1] = 2.0*(Stiff[0][3] + Stiff[0][4] + Stiff[0][5]);
   K[0][2] = SDT[0][0];
   K[1][0] = Stiff[3][0] + Stiff[3][1] + Stiff[3][2];
   K[1][1] = 2.0*(Stiff[3][3] + Stiff[3][4] + Stiff[3][5]);
   K[1][2] = SDT[0][3];
   K[2][0] = -2.0*Diff[0];
   K[2][1] = -2.0*Diff[1];
   K[2][2] = -2.0*Diff[2]/(Aspect*Aspect);
   
   return K;
}

double Shuffle3Rhombo::ScanningDefParameter()
{
   return Lattice_->DOF()[3];
}

void Shuffle3Rhombo::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[3] = DOF[4] = DOF[5] -= newval;

   Lattice_->SetDOF(DOF);
}

double Shuffle3Rhombo::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void Shuffle3Rhombo::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double Shuffle3Rhombo::ScanningStressParameter()
{
   return Lattice_->Stress()[0][3];
}
   
Vector Shuffle3Rhombo::ScanningRHS()
{
   Vector rhs(1);

   rhs[0] = Lattice_->Stress()[0][0];

   return rhs;
}

Vector Shuffle3Rhombo::ScanningDef()
{
   Vector def(1);

   def[0] = Lattice_->DOF()[0];

   return def;
}

void Shuffle3Rhombo::ScanningUpdate(const Vector &newval)
{
   Vector DOF = Lattice_->DOF();

   DOF[0] = DOF[1] = DOF[2] -= newval[0];

   Lattice_->SetDOF(DOF);
}

Matrix Shuffle3Rhombo::ScanningStiffness()
{
   Matrix K(1,1),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];

   return K;
}
