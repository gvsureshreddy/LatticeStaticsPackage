#include "NaClShuffle.h"

NaClShuffle::NaClShuffle(Lattice *M)
{
   Lattice_ = (NiTiShuffleTPPLat*) M;
}

// Functions required by LatticeMode
Vector NaClShuffle::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - (Diff[1]*Diff[1])/(Aspect*Aspect)
      - Diff[0]*Diff[0];

   return rhs;
}

Vector NaClShuffle::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DOF()[0];
   def[1] = Lattice_->Temp();

   return def;
}

void NaClShuffle::ArcLenUpdate(const Vector &newval)
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

double NaClShuffle::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix NaClShuffle::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2] +
      2.0*(Stiff[0][3] + Stiff[0][4] + Stiff[0][5])/4.0;
   K[0][1] = Lattice_->StressDT()[0][0];
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double NaClShuffle::ScanningDefParameter()
{
   return Lattice_->DOF()[0];
}

void NaClShuffle::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -=newval;
   DOF[1] = DOF[2] = DOF[0];
   DOF[3] = DOF[4] = DOF[5] = DOF[0]/4.0;

   Lattice_->SetDOF(DOF);
}

double NaClShuffle::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void NaClShuffle::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double NaClShuffle::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector NaClShuffle::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector NaClShuffle::ScanningDef()
{
   return Vector(1,0.0);
}

void NaClShuffle::ScanningUpdate(const Vector &newval) {}

Matrix NaClShuffle::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
