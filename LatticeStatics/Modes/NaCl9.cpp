#include "NaCl9.h"

#define DOFS 9

NaCl9::NaCl9(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector NaCl9::ArcLenRHS(double DS,const Vector &Diff,
			  double Aspect)
{
   Vector rhs(2);
   Matrix stress = Lattice_->Stress();

   rhs[0] = stress[0][0] + stress[0][1] + stress[0][2]
      + (stress[0][3] + stress[0][4] + stress[0][5])/4.0;
   rhs[1] = DS*DS - Diff[1]*Diff[1]/(Aspect*Aspect) - Diff[0]*Diff[0];

   return rhs;
}

Vector NaCl9::ArcLenDef()
{
   Vector def(2);

   def[0] = Lattice_->DOF()[0];
   def[1] = Lattice_->Temp();

   return def;
}

void NaCl9::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval[0];
   DOF[3]=DOF[4]=DOF[5] = DOF[0]/4.0;
   for (int i=6;i<DOFS;i++)
      DOF[i] = 0.0;
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[1]);
}

double NaCl9::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[1] /= Aspect;
   New[1] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector NaCl9::DrDt(const Vector &Diff)
{
   Vector ddt(DOFS,0.0);

   ddt[0] = Diff[0]/Diff[1];
   ddt[1] = ddt[2] = ddt[0];

   return ddt;
}

Matrix NaCl9::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(2,2);
   Matrix Stiff=Lattice_->Stiffness();
   Matrix stressdt = Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2]
      + (Stiff[0][3] + Stiff[0][4] + Stiff[0][5])/16.0
      + Stiff[1][0] + Stiff[1][1] + Stiff[1][2]
      + (Stiff[1][3] + Stiff[1][4] + Stiff[1][5])/16.0
      + Stiff[2][0] + Stiff[2][1] + Stiff[2][2]
      + (Stiff[2][3] + Stiff[2][4] + Stiff[2][5])/16.0
      + Stiff[3][0] + Stiff[3][1] + Stiff[3][2]
      + (Stiff[3][3] + Stiff[3][4] + Stiff[3][5])/16.0
      + Stiff[4][0] + Stiff[4][1] + Stiff[4][2]
      + (Stiff[4][3] + Stiff[4][4] + Stiff[4][5])/16.0
      + Stiff[5][0] + Stiff[5][1] + Stiff[5][2]
      + (Stiff[5][3] + Stiff[5][4] + Stiff[5][5])/16.0;
   K[0][1] = stressdt[0][0] + stressdt[0][1] + stressdt[0][2]
      + (2.0/4.0)*(stressdt[0][3] + stressdt[0][4] + stressdt[0][5]);
   K[1][0] = -2.0*Diff[0];
   K[1][1] = -2.0*Diff[1]/(Aspect*Aspect);
   
   return K;
}

double NaCl9::ScanningDefParameter()
{
   return Lattice_->DOF()[0];
}

void NaCl9::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval;
   DOF[3]=DOF[4]=DOF[5] = DOF[0]/4.0;

   Lattice_->SetDOF(DOF);
}

double NaCl9::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void NaCl9::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double NaCl9::ScanningStressParameter()
{
   Matrix stress = Lattice_->Stress();
   return stress[0][0] + stress[0][1] + stress[0][2]
      + (2.0/4.0)*(stress[0][3] + stress[0][4] + stress[0][5]);
}
   
Vector NaCl9::ScanningRHS()
{
   return Vector(1,0.0);
}

Vector NaCl9::ScanningDef()
{
   return Vector(1,0.0);
}

void NaCl9::ScanningUpdate(const Vector &newval) {}

Matrix NaCl9::ScanningStiffness()
{
   return Matrix(1,1,1.0);
}
