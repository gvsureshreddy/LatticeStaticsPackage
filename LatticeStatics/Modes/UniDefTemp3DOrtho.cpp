#include "UniDefTemp3DOrtho.h"

UniDefTemp3DOrtho::UniDefTemp3DOrtho(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3DOrtho::UniDefTemp3DOrtho(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3DOrtho::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(4);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[2][2];
   rhs[2] = S[0][1];
   rhs[3] = DS*DS - (Diff[3]*Diff[3])/(Aspect*Aspect)
      - Diff[0]*Diff[0] -Diff[1]*Diff[1] - Diff[2]*Diff[2];

   return rhs;
}

Vector UniDefTemp3DOrtho::ArcLenDef()
{
   Vector def(4);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[2][2];
   def[2] = Lattice_->DefGrad()[0][1];
   def[3] = Lattice_->Temp();

   return def;
}

void UniDefTemp3DOrtho::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1] -= newval[0];
   U[2][2] -= newval[1];
   U[0][1]=U[1][0] -= newval[2];
   U[0][2]=U[2][0]=U[1][2]=U[2][1] = 0.0;
   
   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[3]);
}

double UniDefTemp3DOrtho::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[3] /= Aspect;
   New[3] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3DOrtho::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(4,4);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[0][2] = 2.0*(Stiff[0][3]);
   K[0][3] = SDT[0][0];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];
   K[1][2] = 2.0*(Stiff[2][3]);
   K[1][3] = SDT[2][2];
   K[2][0] = Stiff[3][0] + Stiff[3][1];
   K[2][1] = Stiff[3][2];
   K[2][2] = 2.0*(Stiff[3][3]);
   K[2][3] = SDT[0][1];
   K[3][0] = -2.0*Diff[0];
   K[3][1] = -2.0*Diff[1];
   K[3][2] = -2.0*Diff[2];
   K[3][3] = -2.0*Diff[3]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp3DOrtho::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][1];
}

void UniDefTemp3DOrtho::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][1]=U[1][0] -= newval;

   Lattice_->SetDefGrad(U);
}

double UniDefTemp3DOrtho::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3DOrtho::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double UniDefTemp3DOrtho::ScanningStressParameter()
{
   return Lattice_->Stress()[0][1];
}
   
Vector UniDefTemp3DOrtho::ScanningRHS()
{
   Vector rhs(2);
   Matrix S=Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[2][2];

   return rhs;
}

Vector UniDefTemp3DOrtho::ScanningDef()
{
   Vector def(2);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[2][2];

   return def;
}

void UniDefTemp3DOrtho::ScanningUpdate(const Vector &newval)
{
   Matrix U = Lattice_->DefGrad();

   U[0][0] = U[1][1] -= newval[0];
   U[2][2] -= newval[1];
   
   Lattice_->SetDefGrad(U);
}

Matrix UniDefTemp3DOrtho::ScanningStiffness()
{
   Matrix K(2,2),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];

   return K;
}
