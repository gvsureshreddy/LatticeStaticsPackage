#include "Full15.h"

Full15::Full15(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector Full15::ArcLenRHS(double DS,const Vector &Diff,
				double Aspect)
{
   Vector rhs(16);
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
   rhs[9] = Stress[0][9];
   rhs[10] = Stress[0][10];
   rhs[11] = Stress[0][11];
   rhs[12] = Stress[0][12];
   rhs[13] = Stress[0][13];
   rhs[14] = Stress[0][14];
   rhs[15] = DS*DS - Diff[15]*Diff[15]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1]
      - Diff[2]*Diff[2]
      - Diff[3]*Diff[3]
      - Diff[4]*Diff[4]
      - Diff[5]*Diff[5]
      - Diff[6]*Diff[6]
      - Diff[7]*Diff[7]
      - Diff[8]*Diff[8]
      - Diff[9]*Diff[9]
      - Diff[10]*Diff[10]
      - Diff[11]*Diff[11]
      - Diff[12]*Diff[12]
      - Diff[13]*Diff[13]
      - Diff[14]*Diff[14];

   return rhs;
}

Vector Full15::ArcLenDef()
{
   Vector def(16),
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
   def[9] = DOF[9];
   def[10] = DOF[10];
   def[11] = DOF[11];
   def[12] = DOF[12];
   def[13] = DOF[13];
   def[14] = DOF[14];
   def[15] = Lattice_->Temp();

   return def;
}

void Full15::ArcLenUpdate(const Vector &newval)
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
   DOF[9] -= newval[9];
   DOF[10] -= newval[10];
   DOF[11] -= newval[11];
   DOF[12] -= newval[12];
   DOF[13] -= newval[13];
   DOF[14] -= newval[14];
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[15]);
}

double Full15::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[15] /= Aspect;
   New[15] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix Full15::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(16,16);
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
   K[0][9] = Stiff[0][9];
   K[0][10] = Stiff[0][10];
   K[0][11] = Stiff[0][11];
   K[0][12] = Stiff[0][12];
   K[0][13] = Stiff[0][13];
   K[0][14] = Stiff[0][14];
   K[0][15] = StressDT[0][0];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = Stiff[1][3];
   K[1][4] = Stiff[1][4];
   K[1][5] = Stiff[1][5];
   K[1][6] = Stiff[1][6];
   K[1][7] = Stiff[1][7];
   K[1][8] = Stiff[1][8];
   K[1][9] = Stiff[1][9];
   K[1][10] = Stiff[1][10];
   K[1][11] = Stiff[1][11];
   K[1][12] = Stiff[1][12];
   K[1][13] = Stiff[1][13];
   K[1][14] = Stiff[1][14];
   K[1][15] = StressDT[0][1];
   K[2][0] = Stiff[2][0];
   K[2][1] = Stiff[2][1];
   K[2][2] = Stiff[2][2];
   K[2][3] = Stiff[2][3];
   K[2][4] = Stiff[2][4];
   K[2][5] = Stiff[2][5];
   K[2][6] = Stiff[2][6];
   K[2][7] = Stiff[2][7];
   K[2][8] = Stiff[2][8];
   K[2][9] = Stiff[2][9];
   K[2][10] = Stiff[2][10];
   K[2][11] = Stiff[2][11];
   K[2][12] = Stiff[2][12];
   K[2][13] = Stiff[2][13];
   K[2][14] = Stiff[2][14];
   K[2][15] = StressDT[0][2];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1];
   K[3][2] = Stiff[3][2];
   K[3][3] = Stiff[3][3];
   K[3][4] = Stiff[3][4];
   K[3][5] = Stiff[3][5];
   K[3][6] = Stiff[3][6];
   K[3][7] = Stiff[3][7];
   K[3][8] = Stiff[3][8];
   K[3][9] = Stiff[3][9];
   K[3][10] = Stiff[3][10];
   K[3][11] = Stiff[3][11];
   K[3][12] = Stiff[3][12];
   K[3][13] = Stiff[3][13];
   K[3][14] = Stiff[3][14];
   K[3][15] = StressDT[0][3];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = Stiff[4][3];
   K[4][4] = Stiff[4][4];
   K[4][5] = Stiff[4][5];
   K[4][6] = Stiff[4][6];
   K[4][7] = Stiff[4][7];
   K[4][8] = Stiff[4][8];
   K[4][9] = Stiff[4][9];
   K[4][10] = Stiff[4][10];
   K[4][11] = Stiff[4][11];
   K[4][12] = Stiff[4][12];
   K[4][13] = Stiff[4][13];
   K[4][14] = Stiff[4][14];
   K[4][15] = StressDT[0][4];
   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = Stiff[5][3];
   K[5][4] = Stiff[5][4];
   K[5][5] = Stiff[5][5];
   K[5][6] = Stiff[5][6];
   K[5][7] = Stiff[5][7];
   K[5][8] = Stiff[5][8];
   K[5][9] = Stiff[5][9];
   K[5][10] = Stiff[5][10];
   K[5][11] = Stiff[5][11];
   K[5][12] = Stiff[5][12];
   K[5][13] = Stiff[5][13];
   K[5][14] = Stiff[5][14];
   K[5][15] = StressDT[0][5];
   K[6][0] = Stiff[6][0];
   K[6][1] = Stiff[6][1];
   K[6][2] = Stiff[6][2];
   K[6][3] = Stiff[6][3];
   K[6][4] = Stiff[6][4];
   K[6][5] = Stiff[6][5];
   K[6][6] = Stiff[6][6];
   K[6][7] = Stiff[6][7];
   K[6][8] = Stiff[6][8];
   K[6][9] = Stiff[6][9];
   K[6][10] = Stiff[6][10];
   K[6][11] = Stiff[6][11];
   K[6][12] = Stiff[6][12];
   K[6][13] = Stiff[6][13];
   K[6][14] = Stiff[6][14];
   K[6][15] = StressDT[0][6];
   K[7][0] = Stiff[7][0];
   K[7][1] = Stiff[7][1];
   K[7][2] = Stiff[7][2];
   K[7][3] = Stiff[7][3];
   K[7][4] = Stiff[7][4];
   K[7][5] = Stiff[7][5];
   K[7][6] = Stiff[7][6];
   K[7][7] = Stiff[7][7];
   K[7][8] = Stiff[7][8];
   K[7][9] = Stiff[7][9];
   K[7][10] = Stiff[7][10];
   K[7][11] = Stiff[7][11];
   K[7][12] = Stiff[7][12];
   K[7][13] = Stiff[7][13];
   K[7][14] = Stiff[7][14];
   K[7][15] = StressDT[0][7];
   K[8][0] = Stiff[8][0];
   K[8][1] = Stiff[8][1];
   K[8][2] = Stiff[8][2];
   K[8][3] = Stiff[8][3];
   K[8][4] = Stiff[8][4];
   K[8][5] = Stiff[8][5];
   K[8][6] = Stiff[8][6];
   K[8][7] = Stiff[8][7];
   K[8][8] = Stiff[8][8];
   K[8][9] = Stiff[8][9];
   K[8][10] = Stiff[8][10];
   K[8][11] = Stiff[8][11];
   K[8][12] = Stiff[8][12];
   K[8][13] = Stiff[8][13];
   K[8][14] = Stiff[8][14];
   K[8][15] = StressDT[0][8];
   K[9][0] = Stiff[9][0];
   K[9][1] = Stiff[9][1];
   K[9][2] = Stiff[9][2];
   K[9][3] = Stiff[9][3];
   K[9][4] = Stiff[9][4];
   K[9][5] = Stiff[9][5];
   K[9][6] = Stiff[9][6];
   K[9][7] = Stiff[9][7];
   K[9][8] = Stiff[9][8];
   K[9][9] = Stiff[9][9];
   K[9][10] = Stiff[9][10];
   K[9][11] = Stiff[9][11];
   K[9][12] = Stiff[9][12];
   K[9][13] = Stiff[9][13];
   K[9][14] = Stiff[9][14];
   K[9][15] = StressDT[0][9];
   K[10][0] = Stiff[10][0];
   K[10][1] = Stiff[10][1];
   K[10][2] = Stiff[10][2];
   K[10][3] = Stiff[10][3];
   K[10][4] = Stiff[10][4];
   K[10][5] = Stiff[10][5];
   K[10][6] = Stiff[10][6];
   K[10][7] = Stiff[10][7];
   K[10][8] = Stiff[10][8];
   K[10][9] = Stiff[10][9];
   K[10][10] = Stiff[10][10];
   K[10][11] = Stiff[10][11];
   K[10][12] = Stiff[10][12];
   K[10][13] = Stiff[10][13];
   K[10][14] = Stiff[10][14];
   K[10][15] = StressDT[0][10];
   K[11][0] = Stiff[11][0];
   K[11][1] = Stiff[11][1];
   K[11][2] = Stiff[11][2];
   K[11][3] = Stiff[11][3];
   K[11][4] = Stiff[11][4];
   K[11][5] = Stiff[11][5];
   K[11][6] = Stiff[11][6];
   K[11][7] = Stiff[11][7];
   K[11][8] = Stiff[11][8];
   K[11][9] = Stiff[11][9];
   K[11][10] = Stiff[11][10];
   K[11][11] = Stiff[11][11];
   K[11][12] = Stiff[11][12];
   K[11][13] = Stiff[11][13];
   K[11][14] = Stiff[11][14];
   K[11][15] = StressDT[0][11];
   K[12][0] = Stiff[12][0];
   K[12][1] = Stiff[12][1];
   K[12][2] = Stiff[12][2];
   K[12][3] = Stiff[12][3];
   K[12][4] = Stiff[12][4];
   K[12][5] = Stiff[12][5];
   K[12][6] = Stiff[12][6];
   K[12][7] = Stiff[12][7];
   K[12][8] = Stiff[12][8];
   K[12][9] = Stiff[12][9];
   K[12][10] = Stiff[12][10];
   K[12][11] = Stiff[12][11];
   K[12][12] = Stiff[12][12];
   K[12][13] = Stiff[12][13];
   K[12][14] = Stiff[12][14];
   K[12][15] = StressDT[0][12];
   K[13][0] = Stiff[13][0];
   K[13][1] = Stiff[13][1];
   K[13][2] = Stiff[13][2];
   K[13][3] = Stiff[13][3];
   K[13][4] = Stiff[13][4];
   K[13][5] = Stiff[13][5];
   K[13][6] = Stiff[13][6];
   K[13][7] = Stiff[13][7];
   K[13][8] = Stiff[13][8];
   K[13][9] = Stiff[13][9];
   K[13][10] = Stiff[13][10];
   K[13][11] = Stiff[13][11];
   K[13][12] = Stiff[13][12];
   K[13][13] = Stiff[13][13];
   K[13][14] = Stiff[13][14];
   K[13][15] = StressDT[0][13];
   K[14][0] = Stiff[14][0];
   K[14][1] = Stiff[14][1];
   K[14][2] = Stiff[14][2];
   K[14][3] = Stiff[14][3];
   K[14][4] = Stiff[14][4];
   K[14][5] = Stiff[14][5];
   K[14][6] = Stiff[14][6];
   K[14][7] = Stiff[14][7];
   K[14][8] = Stiff[14][8];
   K[14][9] = Stiff[14][9];
   K[14][10] = Stiff[14][10];
   K[14][11] = Stiff[14][11];
   K[14][12] = Stiff[14][12];
   K[14][13] = Stiff[14][13];
   K[14][14] = Stiff[14][14];
   K[14][15] = StressDT[0][14];
   K[15][0] = -2.0*Diff[0];
   K[15][1] = -2.0*Diff[1];
   K[15][2] = -2.0*Diff[2];
   K[15][3] = -2.0*Diff[3];
   K[15][4] = -2.0*Diff[4];
   K[15][5] = -2.0*Diff[5];
   K[15][6] = -2.0*Diff[6];
   K[15][7] = -2.0*Diff[7];
   K[15][8] = -2.0*Diff[8];
   K[15][9] = -2.0*Diff[9];
   K[15][10] = -2.0*Diff[10];
   K[15][11] = -2.0*Diff[11];
   K[15][12] = -2.0*Diff[12];
   K[15][13] = -2.0*Diff[13];
   K[15][14] = -2.0*Diff[14];
   K[15][15] = -2.0*Diff[15]/(Aspect*Aspect);
   
   return K;
}

double Full15::ScanningDefParameter()
{
   return Lattice_->DOF()[3];
}

void Full15::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[3] -= newval;

   Lattice_->SetDOF(DOF);
}

double Full15::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void Full15::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double Full15::ScanningStressParameter()
{
   return Lattice_->Stress()[0][3];
}
   
Vector Full15::ScanningRHS()
{
   Matrix Stress=Lattice_->Stress();
   Vector RHS(14);

   RHS[0] = Stress[0][0];
   RHS[1] = Stress[0][1];
   RHS[2] = Stress[0][2];
   RHS[3] = Stress[0][4];
   RHS[4] = Stress[0][5];
   RHS[5] = Stress[0][6];
   RHS[6] = Stress[0][7];
   RHS[7] = Stress[0][8];
   RHS[8] = Stress[0][9];
   RHS[9] = Stress[0][10];
   RHS[10] = Stress[0][11];
   RHS[11] = Stress[0][12];
   RHS[12] = Stress[0][13];
   RHS[13] = Stress[0][14];

   return RHS;
}

Vector Full15::ScanningDef()
{
   Vector DOF=Lattice_->DOF();
   Vector Def(14);

   Def[0] = DOF[0];
   Def[1] = DOF[1];
   Def[2] = DOF[2];
   Def[3] = DOF[4];
   Def[4] = DOF[5];
   Def[5] = DOF[6];
   Def[6] = DOF[7];
   Def[7] = DOF[8];
   Def[8] = DOF[9];
   Def[9] = DOF[10];
   Def[10] = DOF[11];
   Def[11] = DOF[12];
   Def[12] = DOF[13];
   Def[13] = DOF[14];

   return Def;
}

void Full15::ScanningUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] -= newval[1];
   DOF[2] -= newval[2];
   DOF[4] -= newval[3];
   DOF[5] -= newval[4];
   DOF[6] -= newval[5];
   DOF[7] -= newval[6];
   DOF[8] -= newval[7];
   DOF[9] -= newval[8];
   DOF[10] -= newval[9];
   DOF[11] -= newval[10];
   DOF[12] -= newval[11];
   DOF[13] -= newval[12];
   DOF[14] -= newval[13];
   
   Lattice_->SetDOF(DOF);
}

Matrix Full15::ScanningStiffness()
{
   Matrix K(14,14);
   Matrix Stiff=Lattice_->Stiffness();
   
   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1];
   K[0][2] = Stiff[0][2];
   K[0][3] = Stiff[0][4];
   K[0][4] = Stiff[0][5];
   K[0][5] = Stiff[0][6];
   K[0][6] = Stiff[0][7];
   K[0][7] = Stiff[0][8];
   K[0][8] = Stiff[0][9];
   K[0][9] = Stiff[0][10];
   K[0][10] = Stiff[0][11];
   K[0][11] = Stiff[0][12];
   K[0][12] = Stiff[0][13];
   K[0][13] = Stiff[0][14];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = Stiff[1][4];
   K[1][4] = Stiff[1][5];
   K[1][5] = Stiff[1][6];
   K[1][6] = Stiff[1][7];
   K[1][7] = Stiff[1][8];
   K[1][8] = Stiff[1][9];
   K[1][9] = Stiff[1][10];
   K[1][10] = Stiff[1][11];
   K[1][11] = Stiff[1][12];
   K[1][12] = Stiff[1][13];
   K[1][13] = Stiff[1][14];
   K[2][0] = Stiff[2][0];
   K[2][1] = Stiff[2][1];
   K[2][2] = Stiff[2][2];
   K[2][3] = Stiff[2][4];
   K[2][4] = Stiff[2][5];
   K[2][5] = Stiff[2][6];
   K[2][6] = Stiff[2][7];
   K[2][7] = Stiff[2][8];
   K[2][8] = Stiff[2][9];
   K[2][9] = Stiff[2][10];
   K[2][10] = Stiff[2][11];
   K[2][11] = Stiff[2][12];
   K[2][12] = Stiff[2][13];
   K[2][13] = Stiff[2][14];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][1];
   K[3][2] = Stiff[3][2];
   K[3][3] = Stiff[3][4];
   K[3][4] = Stiff[3][5];
   K[3][5] = Stiff[3][6];
   K[3][6] = Stiff[3][7];
   K[3][7] = Stiff[3][8];
   K[3][8] = Stiff[3][9];
   K[3][9] = Stiff[3][10];
   K[3][10] = Stiff[3][11];
   K[3][11] = Stiff[3][12];
   K[3][12] = Stiff[3][13];
   K[3][13] = Stiff[3][14];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = Stiff[4][4];
   K[4][4] = Stiff[4][5];
   K[4][5] = Stiff[4][6];
   K[4][6] = Stiff[4][7];
   K[4][7] = Stiff[4][8];
   K[4][8] = Stiff[4][9];
   K[4][9] = Stiff[4][10];
   K[4][10] = Stiff[4][11];
   K[4][11] = Stiff[4][12];
   K[4][12] = Stiff[4][13];
   K[4][13] = Stiff[4][14];
   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = Stiff[5][4];
   K[5][4] = Stiff[5][5];
   K[5][5] = Stiff[5][6];
   K[5][6] = Stiff[5][7];
   K[5][7] = Stiff[5][8];
   K[5][8] = Stiff[5][9];
   K[5][9] = Stiff[5][10];
   K[5][10] = Stiff[5][11];
   K[5][11] = Stiff[5][12];
   K[5][12] = Stiff[5][13];
   K[5][13] = Stiff[5][14];
   K[6][0] = Stiff[6][0];
   K[6][1] = Stiff[6][1];
   K[6][2] = Stiff[6][2];
   K[6][3] = Stiff[6][4];
   K[6][4] = Stiff[6][5];
   K[6][5] = Stiff[6][6];
   K[6][6] = Stiff[6][7];
   K[6][7] = Stiff[6][8];
   K[6][8] = Stiff[6][9];
   K[6][9] = Stiff[6][10];
   K[6][10] = Stiff[6][11];
   K[6][11] = Stiff[6][12];
   K[6][12] = Stiff[6][13];
   K[6][13] = Stiff[6][14];
   K[7][0] = Stiff[7][0];
   K[7][1] = Stiff[7][1];
   K[7][2] = Stiff[7][2];
   K[7][3] = Stiff[7][4];
   K[7][4] = Stiff[7][5];
   K[7][5] = Stiff[7][6];
   K[7][6] = Stiff[7][7];
   K[7][7] = Stiff[7][8];
   K[7][8] = Stiff[7][9];
   K[7][9] = Stiff[7][10];
   K[7][10] = Stiff[7][11];
   K[7][11] = Stiff[7][12];
   K[7][12] = Stiff[7][13];
   K[7][13] = Stiff[7][14];
   K[8][0] = Stiff[8][0];
   K[8][1] = Stiff[8][1];
   K[8][2] = Stiff[8][2];
   K[8][3] = Stiff[8][4];
   K[8][4] = Stiff[8][5];
   K[8][5] = Stiff[8][6];
   K[8][6] = Stiff[8][7];
   K[8][7] = Stiff[8][8];
   K[8][8] = Stiff[8][9];
   K[8][9] = Stiff[8][10];
   K[8][10] = Stiff[8][11];
   K[8][11] = Stiff[8][12];
   K[8][12] = Stiff[8][13];
   K[8][13] = Stiff[8][14];
   K[9][0] = Stiff[9][0];
   K[9][1] = Stiff[9][1];
   K[9][2] = Stiff[9][2];
   K[9][3] = Stiff[9][4];
   K[9][4] = Stiff[9][5];
   K[9][5] = Stiff[9][6];
   K[9][6] = Stiff[9][7];
   K[9][7] = Stiff[9][8];
   K[9][8] = Stiff[9][9];
   K[9][9] = Stiff[9][10];
   K[9][10] = Stiff[9][11];
   K[9][11] = Stiff[9][12];
   K[9][12] = Stiff[9][13];
   K[9][13] = Stiff[9][14];
   K[10][0] = Stiff[10][0];
   K[10][1] = Stiff[10][1];
   K[10][2] = Stiff[10][2];
   K[10][3] = Stiff[10][4];
   K[10][4] = Stiff[10][5];
   K[10][5] = Stiff[10][6];
   K[10][6] = Stiff[10][7];
   K[10][7] = Stiff[10][8];
   K[10][8] = Stiff[10][9];
   K[10][9] = Stiff[10][10];
   K[10][10] = Stiff[10][11];
   K[10][11] = Stiff[10][12];
   K[10][12] = Stiff[10][13];
   K[10][13] = Stiff[10][14];
   K[11][0] = Stiff[11][0];
   K[11][1] = Stiff[11][1];
   K[11][2] = Stiff[11][2];
   K[11][3] = Stiff[11][4];
   K[11][4] = Stiff[11][5];
   K[11][5] = Stiff[11][6];
   K[11][6] = Stiff[11][7];
   K[11][7] = Stiff[11][8];
   K[11][8] = Stiff[11][9];
   K[11][9] = Stiff[11][10];
   K[11][10] = Stiff[11][11];
   K[11][11] = Stiff[11][12];
   K[11][12] = Stiff[11][13];
   K[11][13] = Stiff[11][14];
   K[12][0] = Stiff[12][0];
   K[12][1] = Stiff[12][1];
   K[12][2] = Stiff[12][2];
   K[12][3] = Stiff[12][4];
   K[12][4] = Stiff[12][5];
   K[12][5] = Stiff[12][6];
   K[12][6] = Stiff[12][7];
   K[12][7] = Stiff[12][8];
   K[12][8] = Stiff[12][9];
   K[12][9] = Stiff[12][10];
   K[12][10] = Stiff[12][11];
   K[12][11] = Stiff[12][12];
   K[12][12] = Stiff[12][13];
   K[12][13] = Stiff[12][14];
   K[13][0] = Stiff[13][0];
   K[13][1] = Stiff[13][1];
   K[13][2] = Stiff[13][2];
   K[13][3] = Stiff[13][4];
   K[13][4] = Stiff[13][5];
   K[13][5] = Stiff[13][6];
   K[13][6] = Stiff[13][7];
   K[13][7] = Stiff[13][8];
   K[13][8] = Stiff[13][9];
   K[13][9] = Stiff[13][10];
   K[13][10] = Stiff[13][11];
   K[13][11] = Stiff[13][12];
   K[13][12] = Stiff[13][13];
   K[13][13] = Stiff[13][14];
   K[14][0] = Stiff[14][0];
   K[14][1] = Stiff[14][1];
   K[14][2] = Stiff[14][2];
   K[14][3] = Stiff[14][4];
   K[14][4] = Stiff[14][5];
   K[14][5] = Stiff[14][6];
   K[14][6] = Stiff[14][7];
   K[14][7] = Stiff[14][8];
   K[14][8] = Stiff[14][9];
   K[14][9] = Stiff[14][10];
   K[14][10] = Stiff[14][11];
   K[14][11] = Stiff[14][12];
   K[14][12] = Stiff[14][13];
   K[14][13] = Stiff[14][14];
   
   return K;
}
