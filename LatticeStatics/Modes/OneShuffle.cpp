#include "OneShuffle.h"

OneShuffle::OneShuffle(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector OneShuffle::ArcLenRHS(double DS,const Vector &Diff,
				double Aspect)
{
   Vector rhs(8);
   Matrix Stress = Lattice_->Stress();

   rhs[0] = Stress[0][0];
   rhs[1] = Stress[0][1];
   rhs[2] = Stress[0][2];
   rhs[3] = Stress[0][3];
   rhs[4] = Stress[0][4];
   rhs[5] = Stress[0][5];
   rhs[6] = Stress[0][6];
   rhs[7] = DS*DS - Diff[7]*Diff[7]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1]
      - Diff[2]*Diff[2]
      - Diff[3]*Diff[3]
      - Diff[4]*Diff[4]
      - Diff[5]*Diff[5]
      - Diff[6]*Diff[6];

   return rhs;
}

Vector OneShuffle::ArcLenDef()
{
   Vector def(8),
      DOF = Lattice_->DOF();

   def[0] = DOF[0];
   def[1] = DOF[1];
   def[2] = DOF[2];
   def[3] = DOF[3];
   def[4] = DOF[4];
   def[5] = DOF[5];
   def[6] = DOF[6];
   def[7] = Lattice_->Temp();

   return def;
}

void OneShuffle::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] -= newval[1];
   DOF[2] -= newval[2];
   DOF[3] -= newval[3];
   DOF[4] -= newval[4];
   DOF[5] -= newval[5];
   DOF[6] -= newval[6];
   DOF[7] = DOF[8] = 0.0;
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[7]);
}

double OneShuffle::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[7] /= Aspect;
   New[7] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix OneShuffle::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(8,8);
   Matrix Stiff=Lattice_->Stiffness(),
      StressDT = Lattice_->StressDT();

   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1];
   K[0][2] = Stiff[0][2];
   K[0][3] = 2.0*Stiff[0][3];
   K[0][4] = 2.0*Stiff[0][4];
   K[0][5] = 2.0*Stiff[0][5];
   K[0][6] = Stiff[0][6];
   K[0][7] = StressDT[0][0];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = 2.0*Stiff[1][3];
   K[1][4] = 2.0*Stiff[1][4];
   K[1][5] = 2.0*Stiff[1][5];
   K[1][6] = Stiff[1][6];
   K[1][7] = StressDT[0][1];   
   K[2][0] = Stiff[2][0];
   K[2][1] = Stiff[2][1];
   K[2][2] = Stiff[2][2];
   K[2][3] = 2.0*Stiff[2][3];
   K[2][4] = 2.0*Stiff[2][4];
   K[2][5] = 2.0*Stiff[2][5];
   K[2][6] = Stiff[2][6];
   K[2][7] = StressDT[0][2];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1];
   K[3][2] = Stiff[3][2];
   K[3][3] = 2.0*Stiff[3][3];
   K[3][4] = 2.0*Stiff[3][4];
   K[3][5] = 2.0*Stiff[3][5];
   K[3][6] = Stiff[3][6];
   K[3][7] = StressDT[0][3];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = 2.0*Stiff[4][3];
   K[4][4] = 2.0*Stiff[4][4];
   K[4][5] = 2.0*Stiff[4][5];
   K[4][6] = Stiff[4][6];
   K[4][7] = StressDT[0][4];
   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = 2.0*Stiff[5][3];
   K[5][4] = 2.0*Stiff[5][4];
   K[5][5] = 2.0*Stiff[5][5];
   K[5][6] = Stiff[5][6];
   K[5][7] = StressDT[0][5];
   K[6][0] = Stiff[6][0];
   K[6][1] = Stiff[6][1];
   K[6][2] = Stiff[6][2];
   K[6][3] = 2.0*Stiff[6][3];
   K[6][4] = 2.0*Stiff[6][4];
   K[6][5] = 2.0*Stiff[6][5];
   K[6][6] = Stiff[6][6];
   K[6][7] = StressDT[0][6];
   K[7][0] = -2.0*Diff[0];
   K[7][1] = -2.0*Diff[1];
   K[7][2] = -2.0*Diff[2];
   K[7][3] = -2.0*Diff[3];
   K[7][4] = -2.0*Diff[4];
   K[7][5] = -2.0*Diff[5];
   K[7][6] = -2.0*Diff[6];
   K[7][7] = -2.0*Diff[7]/(Aspect*Aspect);
   
   return K;
}

double OneShuffle::ScanningDefParameter()
{
   return Lattice_->DOF()[6];
}

void OneShuffle::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[6] -= newval;

   Lattice_->SetDOF(DOF);
}

double OneShuffle::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void OneShuffle::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double OneShuffle::ScanningStressParameter()
{
   return Lattice_->Stress()[0][6];
}
   
Vector OneShuffle::ScanningRHS()
{
   Matrix Stress=Lattice_->Stress();
   Vector RHS(6);

   RHS[0] = Stress[0][0];
   RHS[1] = Stress[0][1];
   RHS[2] = Stress[0][2];
   RHS[3] = Stress[0][3];
   RHS[4] = Stress[0][4];
   RHS[5] = Stress[0][5];

   return RHS;
}

Vector OneShuffle::ScanningDef()
{
   Vector DOF=Lattice_->DOF();
   Vector Def(6);

   Def[0] = DOF[0];
   Def[1] = DOF[1];
   Def[2] = DOF[2];
   Def[3] = DOF[3];
   Def[4] = DOF[4];
   Def[5] = DOF[5];

   return Def;
}

void OneShuffle::ScanningUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] -= newval[1];
   DOF[2] -= newval[2];
   DOF[3] -= newval[3];
   DOF[4] -= newval[4];
   DOF[5] -= newval[5];
   
   Lattice_->SetDOF(DOF);
}

Matrix OneShuffle::ScanningStiffness()
{
   Matrix K(6,6);
   Matrix Stiff=Lattice_->Stiffness();
   
   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1];
   K[0][2] = Stiff[0][2];
   K[0][3] = 2.0*Stiff[0][3];
   K[0][4] = 2.0*Stiff[0][4];
   K[0][5] = 2.0*Stiff[0][5];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = 2.0*Stiff[1][3];
   K[1][4] = 2.0*Stiff[1][4];
   K[1][5] = 2.0*Stiff[1][5];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1];
   K[3][2] = Stiff[3][2];
   K[3][3] = 2.0*Stiff[3][3];
   K[3][4] = 2.0*Stiff[3][4];
   K[3][5] = 2.0*Stiff[3][5];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = 2.0*Stiff[4][3];
   K[4][4] = 2.0*Stiff[4][4];
   K[4][5] = 2.0*Stiff[4][5];
   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = 2.0*Stiff[5][3];
   K[5][4] = 2.0*Stiff[5][4];
   K[5][5] = 2.0*Stiff[5][5];

   return K;
}
