#include "OrthoShuffle.h"

OrthoShuffle::OrthoShuffle(Lattice *M)
{
   Lattice_ = (NiTiShuffleTPPLat *) M;
}

// Functions required by LatticeMode
Vector OrthoShuffle::ArcLenRHS(double DS,const Vector &Diff,
				double Aspect)
{
   Vector rhs(4);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = Lattice_->Stress()[0][2];
   rhs[2] = Lattice_->Stress()[0][3];
   rhs[3] = DS*DS - Diff[3]*Diff[3]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1]
      - Diff[2]*Diff[2];

   return rhs;
}

Vector OrthoShuffle::ArcLenDef()
{
   Vector def(4),
      DOF = Lattice_->DOF();

   def[0] = DOF[0];
   def[1] = DOF[2];
   def[2] = DOF[3];
   def[3] = Lattice_->Temp();

   return def;
}

void OrthoShuffle::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] = DOF[0];
   DOF[2] -= newval[1];
   DOF[3] -= newval[2];
   DOF[4]=DOF[5]=DOF[6]=DOF[7] = 0.0;

   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[3]);
}

double OrthoShuffle::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[3] /= Aspect;
   New[3] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix OrthoShuffle::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(4,4);
   Matrix Stiff=Lattice_->Stiffness(),
      StressDT = Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[0][2] = 2.0*Stiff[0][3];
   K[0][3] = StressDT[0][0];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];
   K[1][2] = 2.0*Stiff[2][3];
   K[1][3] = StressDT[0][2];
   K[2][0] = Stiff[3][0] + Stiff[3][1];
   K[2][1] = Stiff[3][2];
   K[2][2] = 2.0*Stiff[3][3];
   K[2][3] = StressDT[0][3];
   K[3][0] = -2.0*Diff[0];
   K[3][1] = -2.0*Diff[1];
   K[3][2] = -2.0*Diff[2];
   K[3][3] = -2.0*Diff[3]/(Aspect*Aspect);
   
   return K;
}

double OrthoShuffle::ScanningDefParameter()
{
   return Lattice_->DOF()[3];
}

void OrthoShuffle::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[3] -= newval;

   Lattice_->SetDOF(DOF);
}

double OrthoShuffle::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void OrthoShuffle::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double OrthoShuffle::ScanningStressParameter()
{
   return Lattice_->Stress()[0][3];
}
   
Vector OrthoShuffle::ScanningRHS()
{
   Matrix Stress=Lattice_->Stress();
   Vector RHS(2);

   RHS[0] = Stress[0][0];
   RHS[1] = Stress[0][2];

   return RHS;
}

Vector OrthoShuffle::ScanningDef()
{
   Vector DOF=Lattice_->DOF();
   Vector Def(2);

   Def[0] = DOF[0];
   Def[1] = DOF[2];

   return Def;
}

void OrthoShuffle::ScanningUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] = DOF[0];
   DOF[2] -= newval[1];

   Lattice_->SetDOF(DOF);
}

Matrix OrthoShuffle::ScanningStiffness()
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];

   return K;
}
