#include "UniDefTemp3D3MMono.h"

UniDefTemp3D3MMono::UniDefTemp3D3MMono(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3D3MMono::UniDefTemp3D3MMono(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3D3MMono::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(5);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[2][2];
   rhs[2] = S[0][1];
   rhs[3] = S[0][2];
   rhs[4] = DS*DS - (Diff[4]*Diff[4])/(Aspect*Aspect)
      - Diff[0]*Diff[0] -Diff[1]*Diff[1] - Diff[2]*Diff[2] - Diff[3]*Diff[3];

   return rhs;
}

Vector UniDefTemp3D3MMono::ArcLenDef()
{
   Vector def(5);
   Matrix Def=Lattice_->DefGrad();

   def[0] = Def[0][0];
   def[1] = Def[2][2];
   def[2] = Def[0][1];
   def[3] = Def[0][2];
   def[4] = Lattice_->Temp();

   return def;
}

void UniDefTemp3D3MMono::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1] -= newval[0];
   U[2][2] -= newval[1];
   U[0][1]=U[1][0] -= newval[2];
   U[0][2]=U[2][0] -= newval[3];
   U[1][2]=U[2][1] = -U[0][2];
   
   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[4]);
}

double UniDefTemp3D3MMono::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[4] /= Aspect;
   New[4] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3D3MMono::ArcLenStiffness(const Vector &Diff,double Aspect)
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
   K[1][4] = SDT[2][2];
   K[2][0] = Stiff[3][0] + Stiff[3][1];
   K[2][1] = Stiff[3][2];
   K[2][2] = 2.0*(Stiff[3][3]);
   K[2][3] = 2.0*(Stiff[3][4] - Stiff[3][5]);
   K[2][4] = SDT[0][1];
   K[3][0] = Stiff[4][0] + Stiff[4][1];
   K[3][1] = Stiff[4][2];
   K[3][2] = 2.0*(Stiff[4][3]);
   K[3][3] = 2.0*(Stiff[4][4] - Stiff[4][5]);
   K[3][4] = SDT[0][2];
   K[4][0] = -2.0*Diff[0];
   K[4][1] = -2.0*Diff[1];
   K[4][2] = -2.0*Diff[2];
   K[4][3] = -2.0*Diff[3];
   K[4][4] = -2.0*Diff[4]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp3D3MMono::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][2];
}

void UniDefTemp3D3MMono::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][2]=U[2][0] += newval;
   U[1][2]=U[2][1] = - U[0][2];
   
   Lattice_->SetDefGrad(U);
}

double UniDefTemp3D3MMono::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3D3MMono::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() + newval);
}

double UniDefTemp3D3MMono::ScanningStressParameter()
{
   return Lattice_->Stress()[0][2];
}
   
Vector UniDefTemp3D3MMono::ScanningRHS()
{
   Vector rhs(3);
   Matrix S=Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[2][2];
   rhs[2] = S[0][1];

   return rhs;
}

Vector UniDefTemp3D3MMono::ScanningDef()
{
   Vector def(3);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[2][2];
   def[2] = Lattice_->DefGrad()[0][1];

   return def;
}

void UniDefTemp3D3MMono::ScanningUpdate(const Vector &newval)
{
   Matrix U = Lattice_->DefGrad();

   U[0][0] = U[1][1] -= newval[0];
   U[2][2] -= newval[1];
   U[0][1]=U[1][0] -= newval[2];
   
   Lattice_->SetDefGrad(U);
}

Matrix UniDefTemp3D3MMono::ScanningStiffness()
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
