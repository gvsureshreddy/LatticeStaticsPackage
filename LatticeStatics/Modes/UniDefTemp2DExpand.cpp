#include "UniDefTemp2DExpand.h"

UniDefTemp2DExpand::UniDefTemp2DExpand(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp2DExpand::UniDefTemp2DExpand(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp2DExpand::ArcLenRHS(double DS,const Vector &Diff,
				   double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - Diff[1]*Diff[1]/(Aspect*Aspect) - Diff[0]*Diff[0];

   return rhs;
}

Vector UniDefTemp2DExpand::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->Temp();

   return def;
}

void UniDefTemp2DExpand::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1] -= newval[0];
   U[0][1]=U[1][0] = 0.0;

   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double UniDefTemp2DExpand::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp2DExpand::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Lattice_->StressDT()[0][0];
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double UniDefTemp2DExpand::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][0];
}

void UniDefTemp2DExpand::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1] -= newval;

   Lattice_->SetDefGrad(U);
}

double UniDefTemp2DExpand::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp2DExpand::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double UniDefTemp2DExpand::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector UniDefTemp2DExpand::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector UniDefTemp2DExpand::ScanningDef()
{
   return Vector(1,0.0);
}

void UniDefTemp2DExpand::ScanningUpdate(const Vector &newval) {}

Matrix UniDefTemp2DExpand::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
