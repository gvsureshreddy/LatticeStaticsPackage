#include "FullShuffle3.h"

FullShuffle3::FullShuffle3(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector FullShuffle3::ArcLenRHS(double DS,const Vector &Diff,
				double Aspect)
{
   Vector rhs(10);
   Matrix Stress = Lattice_->Stress();

   rhs[0] = Stress[0][0];
   rhs[1] = Stress[0][1];
   rhs[2] = Stress[0][2];
   rhs[3] = Stress[0][3];
   rhs[4] = Stress[0][4];
   rhs[5] = Stress[0][5];
   rhs[6] = Stress[0][6];
   rhs[7] = Stress[0][7];
   rhs[8] = Stress[0][8];
   rhs[9] = DS*DS - Diff[9]*Diff[9]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1]
      - Diff[2]*Diff[2]
      - Diff[3]*Diff[3]
      - Diff[4]*Diff[4]
      - Diff[5]*Diff[5]
      - Diff[6]*Diff[6]
      - Diff[7]*Diff[7]
      - Diff[8]*Diff[8];

   return rhs;
}

Vector FullShuffle3::ArcLenDef()
{
   Vector def(10),
      DOF = Lattice_->DOF();

   def[0] = DOF[0];
   def[1] = DOF[1];
   def[2] = DOF[2];
   def[3] = DOF[3];
   def[4] = DOF[4];
   def[5] = DOF[5];
   def[6] = DOF[6];
   def[7] = DOF[7];
   def[8] = DOF[8];
   def[9] = Lattice_->Temp();

   return def;
}

void FullShuffle3::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] -= newval[1];
   DOF[2] -= newval[2];
   DOF[3] -= newval[3];
   DOF[4] -= newval[4];
   DOF[5] -= newval[5];
   DOF[6] -= newval[6];
   DOF[7] -= newval[7];
   DOF[8] -= newval[8];
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[9]);
}

double FullShuffle3::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[9] /= Aspect;
   New[9] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector FullShuffle3::DrDt(const Vector &Diff)
{
   int size=Diff.Dim()-1;
   Vector ddt(size);

   for (int i=0;i<size;i++)
   {
      ddt[i] = Diff[i]/Diff[size];
   }

   return ddt;
}

Matrix FullShuffle3::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(10,10);
   Matrix Stiff=Lattice_->Stiffness(),
      StressDT = Lattice_->StressDT();

   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1];
   K[0][2] = Stiff[0][2];
   K[0][3] = Stiff[0][3];
   K[0][4] = Stiff[0][4];
   K[0][5] = Stiff[0][5];
   K[0][6] = Stiff[0][6];
   K[0][7] = Stiff[0][7];
   K[0][8] = Stiff[0][8];
   K[0][9] = StressDT[0][0];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = Stiff[1][3];
   K[1][4] = Stiff[1][4];
   K[1][5] = Stiff[1][5];
   K[1][6] = Stiff[1][6];
   K[1][7] = Stiff[1][7];
   K[1][8] = Stiff[1][8];
   K[1][9] = StressDT[0][1];   
   K[2][0] = Stiff[2][0];
   K[2][1] = Stiff[2][1];
   K[2][2] = Stiff[2][2];
   K[2][3] = Stiff[2][3];
   K[2][4] = Stiff[2][4];
   K[2][5] = Stiff[2][5];
   K[2][6] = Stiff[2][6];
   K[2][7] = Stiff[2][7];
   K[2][8] = Stiff[2][8];
   K[2][9] = StressDT[0][2];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1];
   K[3][2] = Stiff[3][2];
   K[3][3] = Stiff[3][3];
   K[3][4] = Stiff[3][4];
   K[3][5] = Stiff[3][5];
   K[3][6] = Stiff[3][6];
   K[3][7] = Stiff[3][7];
   K[3][8] = Stiff[3][8];
   K[3][9] = StressDT[0][3];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = Stiff[4][3];
   K[4][4] = Stiff[4][4];
   K[4][5] = Stiff[4][5];
   K[4][6] = Stiff[4][6];
   K[4][7] = Stiff[4][7];
   K[4][8] = Stiff[4][8];
   K[4][9] = StressDT[0][4];
   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = Stiff[5][3];
   K[5][4] = Stiff[5][4];
   K[5][5] = Stiff[5][5];
   K[5][6] = Stiff[5][6];
   K[5][7] = Stiff[5][7];
   K[5][8] = Stiff[5][8];
   K[5][9] = StressDT[0][5];
   K[6][0] = Stiff[6][0];
   K[6][1] = Stiff[6][1];
   K[6][2] = Stiff[6][2];
   K[6][3] = Stiff[6][3];
   K[6][4] = Stiff[6][4];
   K[6][5] = Stiff[6][5];
   K[6][6] = Stiff[6][6];
   K[6][7] = Stiff[6][7];
   K[6][8] = Stiff[6][8];
   K[6][9] = StressDT[0][6];
   K[7][0] = Stiff[7][0];
   K[7][1] = Stiff[7][1];
   K[7][2] = Stiff[7][2];
   K[7][3] = Stiff[7][3];
   K[7][4] = Stiff[7][4];
   K[7][5] = Stiff[7][5];
   K[7][6] = Stiff[7][6];
   K[7][7] = Stiff[7][7];
   K[7][8] = Stiff[7][8];
   K[7][9] = StressDT[0][7];
   K[8][0] = Stiff[8][0];
   K[8][1] = Stiff[8][1];
   K[8][2] = Stiff[8][2];
   K[8][3] = Stiff[8][3];
   K[8][4] = Stiff[8][4];
   K[8][5] = Stiff[8][5];
   K[8][6] = Stiff[8][6];
   K[8][7] = Stiff[8][7];
   K[8][8] = Stiff[8][8];
   K[8][9] = StressDT[0][8];
   K[9][0] = -2.0*Diff[0];
   K[9][1] = -2.0*Diff[1];
   K[9][2] = -2.0*Diff[2];
   K[9][3] = -2.0*Diff[3];
   K[9][4] = -2.0*Diff[4];
   K[9][5] = -2.0*Diff[5];
   K[9][6] = -2.0*Diff[6];
   K[9][7] = -2.0*Diff[7];
   K[9][8] = -2.0*Diff[8];
   K[9][9] = -2.0*Diff[9]/(Aspect*Aspect);
   
   return K;
}

double FullShuffle3::ScanningDefParameter()
{
   return Lattice_->DOF()[6];
}

void FullShuffle3::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[6] -= newval;

   Lattice_->SetDOF(DOF);
}

double FullShuffle3::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void FullShuffle3::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double FullShuffle3::ScanningStressParameter()
{
   return Lattice_->Stress()[0][6];
}
   
Vector FullShuffle3::ScanningRHS()
{
   Matrix Stress=Lattice_->Stress();
   Vector RHS(8);

   RHS[0] = Stress[0][0];
   RHS[1] = Stress[0][1];
   RHS[2] = Stress[0][2];
   RHS[3] = Stress[0][3];
   RHS[4] = Stress[0][4];
   RHS[5] = Stress[0][5];
   RHS[6] = Stress[0][7];
   RHS[7] = Stress[0][8];

   return RHS;
}

Vector FullShuffle3::ScanningDef()
{
   Vector DOF=Lattice_->DOF();
   Vector Def(8);

   Def[0] = DOF[0];
   Def[1] = DOF[1];
   Def[2] = DOF[2];
   Def[3] = DOF[3];
   Def[4] = DOF[4];
   Def[5] = DOF[5];
   Def[6] = DOF[7];
   Def[7] = DOF[8];

   return Def;
}

void FullShuffle3::ScanningUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] -= newval[1];
   DOF[2] -= newval[2];
   DOF[3] -= newval[3];
   DOF[4] -= newval[4];
   DOF[5] -= newval[5];
   DOF[7] -= newval[6];
   DOF[8] -= newval[7];
   
   Lattice_->SetDOF(DOF);
}

Matrix FullShuffle3::ScanningStiffness()
{
   Matrix K(8,8);
   Matrix Stiff=Lattice_->Stiffness();
   
   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1];
   K[0][2] = Stiff[0][2];
   K[0][3] = Stiff[0][3];
   K[0][4] = Stiff[0][4];
   K[0][5] = Stiff[0][5];
   K[0][6] = Stiff[0][7];
   K[0][7] = Stiff[0][8];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = Stiff[1][3];
   K[1][4] = Stiff[1][4];
   K[1][5] = Stiff[1][5];
   K[1][6] = Stiff[1][7];
   K[1][7] = Stiff[1][8];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1];
   K[3][2] = Stiff[3][2];
   K[3][3] = Stiff[3][3];
   K[3][4] = Stiff[3][4];
   K[3][5] = Stiff[3][5];
   K[3][6] = Stiff[3][7];
   K[3][7] = Stiff[3][8];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = Stiff[4][3];
   K[4][4] = Stiff[4][4];
   K[4][5] = Stiff[4][5];
   K[4][6] = Stiff[4][7];
   K[4][7] = Stiff[4][8];
   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = Stiff[5][3];
   K[5][4] = Stiff[5][4];
   K[5][5] = Stiff[5][5];
   K[5][6] = Stiff[5][7];
   K[5][7] = Stiff[5][8];
   K[6][0] = Stiff[7][0];
   K[6][1] = Stiff[7][1];
   K[6][2] = Stiff[7][2];
   K[6][3] = Stiff[7][3];
   K[6][4] = Stiff[7][4];
   K[6][5] = Stiff[7][5];
   K[6][6] = Stiff[7][7];
   K[6][7] = Stiff[7][8];
   K[7][0] = Stiff[8][0];
   K[7][1] = Stiff[8][1];
   K[7][2] = Stiff[8][2];
   K[7][3] = Stiff[8][3];
   K[7][4] = Stiff[8][4];
   K[7][5] = Stiff[8][5];
   K[7][6] = Stiff[8][7];
   K[7][7] = Stiff[8][8];

   return K;
}
