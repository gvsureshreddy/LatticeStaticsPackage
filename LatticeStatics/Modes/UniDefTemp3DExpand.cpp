#include "UniDefTemp3DExpand.h"

UniDefTemp3DExpand::UniDefTemp3DExpand(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3DExpand::UniDefTemp3DExpand(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3DExpand::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - Diff[1]*Diff[1]/(Aspect*Aspect) - Diff[0]*Diff[0];

   return rhs;
}

Vector UniDefTemp3DExpand::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->Temp();

   return def;
}

void UniDefTemp3DExpand::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1]=U[2][2] -= newval[0];
   U[0][1]=U[1][0]=U[0][2]=U[2][0]=U[1][2]=U[2][1] = 0.0;

   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double UniDefTemp3DExpand::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3DExpand::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];
   K[0][1] = Lattice_->StressDT()[0][0];
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp3DExpand::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][0];
}

void UniDefTemp3DExpand::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1]=U[2][2] += newval;

   Lattice_->SetDefGrad(U);
}

double UniDefTemp3DExpand::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3DExpand::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() + newval);
}

double UniDefTemp3DExpand::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector UniDefTemp3DExpand::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector UniDefTemp3DExpand::ScanningDef()
{
   return Vector(1,0.0);
}

void UniDefTemp3DExpand::ScanningUpdate(const Vector &newval) {}

Matrix UniDefTemp3DExpand::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
