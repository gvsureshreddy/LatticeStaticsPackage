#include "UniDefTemp3DTetrag.h"

UniDefTemp3DTetrag::UniDefTemp3DTetrag(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3DTetrag::UniDefTemp3DTetrag(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3DTetrag::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(3);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[2][2];
   rhs[2] = DS*DS - (Diff[2]*Diff[2])/(Aspect*Aspect)
      - Diff[0]*Diff[0] -Diff[1]*Diff[1];

   return rhs;
}

Vector UniDefTemp3DTetrag::ArcLenDef()
{
   Vector def(3);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[2][2];
   def[2] = Lattice_->Temp();

   return def;
}

void UniDefTemp3DTetrag::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1] -= newval[0];
   U[2][2] -= newval[1];
   U[0][1]=U[1][0]=U[0][2]=U[2][0]=U[1][2]=U[2][1] = 0.0;
   
   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[2]);
}

double UniDefTemp3DTetrag::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[2] /= Aspect;
   New[2] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3DTetrag::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(3,3);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[0][2] = SDT[0][0];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];
   K[1][2] = SDT[2][2];
   K[2][0] = -2.0*Diff[0];
   K[2][1] = -2.0*Diff[1];
   K[2][2] = -2.0*Diff[2]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp3DTetrag::ScanningDefParameter()
{
   return Lattice_->DefGrad()[2][2];
}

void UniDefTemp3DTetrag::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[2][2] -= newval;

   Lattice_->SetDefGrad(U);
}

double UniDefTemp3DTetrag::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3DTetrag::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double UniDefTemp3DTetrag::ScanningStressParameter()
{
   return Lattice_->Stress()[2][2];
}
   
Vector UniDefTemp3DTetrag::ScanningRHS()
{
   Vector rhs(1);
   Matrix S=Lattice_->Stress();

   rhs[0] = S[0][0];

   return rhs;
}

Vector UniDefTemp3DTetrag::ScanningDef()
{
   Vector def(1);

   def[0] = Lattice_->DefGrad()[0][0];

   return def;
}

void UniDefTemp3DTetrag::ScanningUpdate(const Vector &newval)
{
   Matrix U = Lattice_->DefGrad();

   U[0][0] = U[1][1] -= newval[0];
   
   Lattice_->SetDefGrad(U);
}

Matrix UniDefTemp3DTetrag::ScanningStiffness()
{
   Matrix K(1,1),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1];

   return K;
}
