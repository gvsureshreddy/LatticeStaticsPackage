#include "UniDefTemp3D1MMono.h"

UniDefTemp3D1MMono::UniDefTemp3D1MMono(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3D1MMono::UniDefTemp3D1MMono(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3D1MMono::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(5);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[1][1];
   rhs[2] = S[1][2];
   rhs[3] = S[0][1];
   rhs[4] = DS*DS - (Diff[4]*Diff[4])/(Aspect*Aspect)
      - Diff[0]*Diff[0] -Diff[1]*Diff[1] - Diff[2]*Diff[2] - Diff[3]*Diff[3];

   return rhs;
}

Vector UniDefTemp3D1MMono::ArcLenDef()
{
   Vector def(5);
   Matrix Def=Lattice_->DefGrad();

   def[0] = Def[0][0];
   def[1] = Def[1][1];
   def[2] = Def[1][2];
   def[3] = Def[0][1];
   def[4] = Lattice_->Temp();

   return def;
}

void UniDefTemp3D1MMono::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0] -= newval[0];
   U[1][1]=U[2][2] -= newval[1];
   U[1][2]=U[2][1] -= newval[2];
   U[0][1]=U[1][0] -= newval[3];
   U[0][2]=U[2][0] = -U[0][1];
   
   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[4]);
}

double UniDefTemp3D1MMono::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[4] /= Aspect;
   New[4] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3D1MMono::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(5,5);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1] + Stiff[0][2];
   K[0][2] = 2.0*(Stiff[0][5]);
   K[0][3] = 2.0*(Stiff[0][3] - Stiff[0][4]);
   K[0][4] = SDT[0][0];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1] + Stiff[1][2];
   K[1][2] = 2.0*(Stiff[1][5]);
   K[1][3] = 2.0*(Stiff[1][3] - Stiff[1][4]);
   K[1][4] = SDT[1][1];
   K[2][0] = Stiff[5][0];
   K[2][1] = Stiff[5][1] + Stiff[5][2];
   K[2][2] = 2.0*(Stiff[5][5]);
   K[2][3] = 2.0*(Stiff[5][3] - Stiff[5][4]);
   K[2][4] = SDT[1][2];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1] + Stiff[3][2];
   K[3][2] = 2.0*(Stiff[3][5]);
   K[3][3] = 2.0*(Stiff[3][3] - Stiff[3][4]);
   K[3][4] = SDT[0][1];
   K[4][0] = -2.0*Diff[0];
   K[4][1] = -2.0*Diff[1];
   K[4][2] = -2.0*Diff[2];
   K[4][3] = -2.0*Diff[3];
   K[4][4] = -2.0*Diff[4]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp3D1MMono::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][1];
}

void UniDefTemp3D1MMono::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][1]=U[1][0] += newval;
   U[0][2]=U[2][0] = -U[0][1];
   
   Lattice_->SetDefGrad(U);
}

double UniDefTemp3D1MMono::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3D1MMono::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() + newval);
}

double UniDefTemp3D1MMono::ScanningStressParameter()
{
   return Lattice_->Stress()[0][1];
}
   
Vector UniDefTemp3D1MMono::ScanningRHS()
{
   Vector rhs(3);
   Matrix S=Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[1][1];
   rhs[2] = S[1][2];

   return rhs;
}

Vector UniDefTemp3D1MMono::ScanningDef()
{
   Vector def(3);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[1][1];
   def[2] = Lattice_->DefGrad()[1][2];

   return def;
}

void UniDefTemp3D1MMono::ScanningUpdate(const Vector &newval)
{
   Matrix U = Lattice_->DefGrad();

   U[0][0] -= newval[0];
   U[1][1]=U[2][2] -= newval[1];
   U[1][2]=U[2][1] -= newval[2];
   
   Lattice_->SetDefGrad(U);
}

Matrix UniDefTemp3D1MMono::ScanningStiffness()
{
   Matrix K(3,3),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1] + Stiff[0][2];
   K[0][2] = 2.0*(Stiff[0][5]);
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1] + Stiff[1][2];
   K[1][2] = 2.0*(Stiff[1][5]);
   K[2][0] = Stiff[5][0];
   K[2][1] = Stiff[5][1] + Stiff[5][2];
   K[2][2] = 2.0*(Stiff[5][5]);

   return K;
}
