#include "UniDefTemp3DNaCl.h"

UniDefTemp3DNaCl::UniDefTemp3DNaCl(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3DNaCl::UniDefTemp3DNaCl(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3DNaCl::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(2);

   rhs[0] = Lattice_->Stress()[0][0];
   rhs[1] = DS*DS - Diff[1]*Diff[1]/(Aspect*Aspect) - Diff[0]*Diff[0];

   return rhs;
}

Vector UniDefTemp3DNaCl::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->Temp();

   return def;
}

void UniDefTemp3DNaCl::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1]=U[2][2] -= newval[0];
   U[0][1]=U[1][0]=U[0][2]=U[2][0]=U[1][2]=U[2][1] = U[0][0]/4.0;

   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double UniDefTemp3DNaCl::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3DNaCl::ArcLenStiffness(const Vector &Diff,double Aspect)
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

double UniDefTemp3DNaCl::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][0];
}

void UniDefTemp3DNaCl::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0]=U[1][1]=U[2][2] -= newval;
   U[0][1]=U[1][0]=U[0][2]=U[2][0]=U[1][2]=U[2][1] = U[0][0]/4.0;

   Lattice_->SetDefGrad(U);
}

double UniDefTemp3DNaCl::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3DNaCl::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double UniDefTemp3DNaCl::ScanningStressParameter()
{
   return Lattice_->Stress()[0][0];
}
   
Vector UniDefTemp3DNaCl::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector UniDefTemp3DNaCl::ScanningDef()
{
   return Vector(1,0.0);
}

void UniDefTemp3DNaCl::ScanningUpdate(const Vector &newval) {}

Matrix UniDefTemp3DNaCl::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
