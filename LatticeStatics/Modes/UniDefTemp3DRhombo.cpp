#include "UniDefTemp3DRhombo.h"

UniDefTemp3DRhombo::UniDefTemp3DRhombo(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3DRhombo::UniDefTemp3DRhombo(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3DRhombo::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(3);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[0][1];
   rhs[2] = DS*DS - (Diff[2]*Diff[2])/(Aspect*Aspect)
      - Diff[0]*Diff[0] - Diff[1]*Diff[1];

   return rhs;
}

Vector UniDefTemp3DRhombo::ArcLenDef()
{
   Vector def(3);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[0][1];
   def[2] = Lattice_->Temp();

   return def;
}

void UniDefTemp3DRhombo::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1]=U[2][2] -= newval[0];
   U[0][1]=U[1][0]=U[0][2]=U[2][0]=U[1][2]=U[2][1] -= newval[1];

   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[2]);
}

double UniDefTemp3DRhombo::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[2] /= Aspect;
   New[2] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3DRhombo::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(3,3);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];
   K[0][1] = 2.0*(Stiff[0][3] + Stiff[0][4] + Stiff[0][5]);
   K[0][2] = SDT[0][0];
   K[1][0] = Stiff[3][0] + Stiff[3][1] + Stiff[3][2];
   K[1][1] = 2.0*(Stiff[3][3] + Stiff[3][4] + Stiff[3][5]);
   K[1][2] = SDT[0][1];
   K[2][0] = -2.0*Diff[0];
   K[2][1] = -2.0*Diff[1];
   K[2][2] = -2.0*Diff[2]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp3DRhombo::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][1];
}

void UniDefTemp3DRhombo::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][1]=U[1][0]=U[0][2]=U[2][0]=U[1][2]=U[2][1] -= newval;

   Lattice_->SetDefGrad(U);
}

double UniDefTemp3DRhombo::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3DRhombo::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double UniDefTemp3DRhombo::ScanningStressParameter()
{
   return Lattice_->Stress()[0][1];
}
   
Vector UniDefTemp3DRhombo::ScanningRHS()
{
   Vector rhs(1);

   rhs[0] = Lattice_->Stress()[0][0];

   return rhs;
}

Vector UniDefTemp3DRhombo::ScanningDef()
{
   Vector def(1);

   def[0] = Lattice_->DefGrad()[0][0];

   return def;
}

void UniDefTemp3DRhombo::ScanningUpdate(const Vector &newval)
{
   Matrix U = Lattice_->DefGrad();

   U[0][0] = U[1][1] = U[2][2] -= newval[0];

   Lattice_->SetDefGrad(U);
}

Matrix UniDefTemp3DRhombo::ScanningStiffness()
{
   Matrix K(1,1),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];

   return K;
}
