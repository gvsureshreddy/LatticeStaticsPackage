#include "Shuffle3Mono3M.h"

Shuffle3Mono3M::Shuffle3Mono3M(Lattice *M)
{
   Lattice_ = (GenericLat*) M;
}

// Functions required by LatticeMode
Vector Shuffle3Mono3M::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(5);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[0][2];
   rhs[2] = S[0][3];
   rhs[3] = S[0][4];
   rhs[4] = DS*DS - (Diff[4]*Diff[4])/(Aspect*Aspect)
      - Diff[0]*Diff[0] -Diff[1]*Diff[1] - Diff[2]*Diff[2] - Diff[3]*Diff[3];

   return rhs;
}

Vector Shuffle3Mono3M::ArcLenDef()
{
   Vector def(5);
   Vector DOF=Lattice_->DOF();

   def[0] = DOF[0];
   def[1] = DOF[2];
   def[2] = DOF[3];
   def[3] = DOF[4];
   def[4] = Lattice_->Temp();

   return def;
}

void Shuffle3Mono3M::ArcLenUpdate(const Vector &newval)
{
   Vector dof=Lattice_->DOF();

   dof[0]=dof[1] -= newval[0];
   dof[2] -= newval[1];
   dof[3] -= newval[2];
   dof[4] -= newval[3];
   dof[5] = -dof[4];
   dof[6]=dof[7]=dof[8] = 0.0;
   
   Lattice_->SetDOF(dof);
   Lattice_->SetTemp(Lattice_->Temp() - newval[4]);
}

double Shuffle3Mono3M::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[4] /= Aspect;
   New[4] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix Shuffle3Mono3M::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(5,5);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[0][2] = 2.0*(Stiff[0][3]);
   K[0][3] = 2.0*(Stiff[0][4] - Stiff[0][5]);
   K[0][4] = SDT[0][0];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];
   K[1][2] = 2.0*(Stiff[2][3]);
   K[1][3] = 2.0*(Stiff[2][4] - Stiff[2][5]);
   K[1][4] = SDT[0][2];
   K[2][0] = Stiff[3][0] + Stiff[3][1];
   K[2][1] = Stiff[3][2];
   K[2][2] = 2.0*(Stiff[3][3]);
   K[2][3] = 2.0*(Stiff[3][4] - Stiff[3][5]);
   K[2][4] = SDT[0][3];
   K[3][0] = Stiff[4][0] + Stiff[4][1];
   K[3][1] = Stiff[4][2];
   K[3][2] = 2.0*(Stiff[4][3]);
   K[3][3] = 2.0*(Stiff[4][4] - Stiff[4][5]);
   K[3][4] = SDT[0][4];
   K[4][0] = -2.0*Diff[0];
   K[4][1] = -2.0*Diff[1];
   K[4][2] = -2.0*Diff[2];
   K[4][3] = -2.0*Diff[3];
   K[4][4] = -2.0*Diff[4]/(Aspect*Aspect);
   
   return K;
}

double Shuffle3Mono3M::ScanningDefParameter()
{
   return Lattice_->DOF()[4];
}

void Shuffle3Mono3M::ScanningDefParamUpdate(const double newval)
{
   Vector dof=Lattice_->DOF();

   dof[4] -= newval;
   dof[5] = -dof[4];

   Lattice_->SetDOF(dof);
}

double Shuffle3Mono3M::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void Shuffle3Mono3M::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double Shuffle3Mono3M::ScanningStressParameter()
{
   return Lattice_->Stress()[0][4];
}
   
Vector Shuffle3Mono3M::ScanningRHS()
{
   Vector rhs(3);
   Matrix S=Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[0][2];
   rhs[2] = S[0][3];

   return rhs;
}

Vector Shuffle3Mono3M::ScanningDef()
{
   Vector def(3),
      DOF=Lattice_->DOF();

   def[0] = DOF[0];
   def[1] = DOF[2];
   def[2] = DOF[3];

   return def;
}

void Shuffle3Mono3M::ScanningUpdate(const Vector &newval)
{
   Vector dof = Lattice_->DOF();

   dof[0] = dof[1] -= newval[0];
   dof[2] -= newval[1];
   dof[3] -= newval[2];
   
   Lattice_->SetDOF(dof);
}

Matrix Shuffle3Mono3M::ScanningStiffness()
{
   Matrix K(3,3),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[0][2] = 2.0*(Stiff[0][3]);
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];
   K[1][2] = 2.0*(Stiff[2][3]);
   K[2][0] = Stiff[3][0] + Stiff[3][1];
   K[2][1] = Stiff[3][2];
   K[2][2] = 2.0*(Stiff[3][3]);

   return K;
}
