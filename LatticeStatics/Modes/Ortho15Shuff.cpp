#include "Ortho15Shuff.h"

Ortho15Shuff::Ortho15Shuff(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector Ortho15Shuff::ArcLenRHS(double DS,const Vector &Diff,
				double Aspect)
{
   Vector rhs(7);
   Matrix Stress = Lattice_->Stress();

   rhs[0] = Stress[0][0];
   rhs[1] = Stress[0][2];
   rhs[2] = Stress[0][3];
   rhs[3] = Stress[0][6];
   rhs[4] = Stress[0][9];
   rhs[5] = Stress[0][12];
   rhs[6] = DS*DS - Diff[6]*Diff[6]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1]
      - Diff[2]*Diff[2]
      - Diff[3]*Diff[3]
      - Diff[4]*Diff[4]
      - Diff[5]*Diff[5];

   return rhs;
}

Vector Ortho15Shuff::ArcLenDef()
{
   Vector def(7),
      DOF = Lattice_->DOF();

   def[0] = DOF[0];
   def[1] = DOF[2];
   def[2] = DOF[3];
   def[3] = DOF[6];
   def[4] = DOF[9];
   def[5] = DOF[12];
   def[6] = Lattice_->Temp();

   return def;
}

void Ortho15Shuff::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[1] = DOF[0] -= newval[0];
   DOF[2] -= newval[1];
   DOF[3] -= newval[2];
   DOF[6] -= newval[3];
   DOF[9] -= newval[4];
   DOF[12] -= newval[5];
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[6]);
}

double Ortho15Shuff::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[6] /= Aspect;
   New[6] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector Ortho15Shuff::DrDt(const Vector &Diff)
{
   int size=Diff.Dim()-1;
   int len=Lattice_->DOF().Dim();
   Vector ddt(len,0.0);

   ddt[0]=ddt[1] = Diff[0]/Diff[size];
   ddt[2] = Diff[1]/Diff[size];
   ddt[3] = Diff[2]/Diff[size];
   ddt[6] = Diff[3]/Diff[size];
   ddt[9] = Diff[4]/Diff[size];
   ddt[12] = Diff[5]/Diff[size];
   
   return ddt;
}

Matrix Ortho15Shuff::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(7,7);
   Matrix Stiff=Lattice_->Stiffness(),
      StressDT = Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1];
   K[0][1] = Stiff[0][2];
   K[0][2] = Stiff[0][3];
   K[0][3] = Stiff[0][6];
   K[0][4] = Stiff[0][9];
   K[0][5] = Stiff[0][12];
   K[0][6] = StressDT[0][0];
   K[1][0] = Stiff[2][0] + Stiff[2][1];
   K[1][1] = Stiff[2][2];
   K[1][2] = Stiff[2][3];
   K[1][3] = Stiff[2][6];
   K[1][4] = Stiff[2][9];
   K[1][5] = Stiff[2][12];
   K[1][6] = StressDT[0][2];
   K[2][0] = Stiff[3][0] + Stiff[3][1];
   K[2][1] = Stiff[3][2];
   K[2][2] = Stiff[3][3];
   K[2][3] = Stiff[3][6];
   K[2][4] = Stiff[3][9];
   K[2][5] = Stiff[3][12];
   K[2][6] = StressDT[0][3];
   K[3][0] = Stiff[6][0] + Stiff[6][1];
   K[3][1] = Stiff[6][2];
   K[3][2] = Stiff[6][3];
   K[3][3] = Stiff[6][6];
   K[3][4] = Stiff[6][9];
   K[3][5] = Stiff[6][12];
   K[3][6] = StressDT[0][6];
   K[4][0] = Stiff[9][0] + Stiff[9][1];
   K[4][1] = Stiff[9][2];
   K[4][2] = Stiff[9][3];
   K[4][3] = Stiff[9][6];
   K[4][4] = Stiff[9][9];
   K[4][5] = Stiff[9][12];
   K[4][6] = StressDT[0][9];
   K[5][0] = Stiff[12][0] + Stiff[12][1];
   K[5][1] = Stiff[12][2];
   K[5][2] = Stiff[12][3];
   K[5][3] = Stiff[12][6];
   K[5][4] = Stiff[12][9];
   K[5][5] = Stiff[12][12];
   K[5][6] = StressDT[0][12];
   K[6][0] = -2.0*Diff[0];
   K[6][1] = -2.0*Diff[1];
   K[6][2] = -2.0*Diff[2];
   K[6][3] = -2.0*Diff[3];
   K[6][4] = -2.0*Diff[4];
   K[6][5] = -2.0*Diff[5];
   K[6][6] = -2.0*Diff[6]/(Aspect*Aspect);

   return K;
}

double Ortho15Shuff::ScanningDefParameter()
{
   return Lattice_->DOF()[6];
}

void Ortho15Shuff::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[6] -= newval;

   Lattice_->SetDOF(DOF);
}

double Ortho15Shuff::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void Ortho15Shuff::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double Ortho15Shuff::ScanningStressParameter()
{
   return Lattice_->Stress()[0][6];
}
   
Vector Ortho15Shuff::ScanningRHS()
{
   Matrix Stress=Lattice_->Stress();
   Vector RHS(5);

   RHS[0] = Stress[0][0];
   RHS[1] = Stress[0][2];
   RHS[2] = Stress[0][3];
   RHS[3] = Stress[0][9];
   RHS[4] = Stress[0][12];

   return RHS;
}

Vector Ortho15Shuff::ScanningDef()
{
   Vector DOF=Lattice_->DOF();
   Vector Def(5);

   Def[0] = DOF[0];
   Def[1] = DOF[2];
   Def[2] = DOF[3];
   Def[3] = DOF[9];
   Def[4] = DOF[12];

   return Def;
}

void Ortho15Shuff::ScanningUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[1]=DOF[0] -= newval[0];
   DOF[2] -= newval[1];
   DOF[3] -= newval[2];
   DOF[9] -= newval[3];
   DOF[12] -= newval[4];
   
   Lattice_->SetDOF(DOF);
}

Matrix Ortho15Shuff::ScanningStiffness()
{
   Matrix K(5,5);
   Matrix Stiff=Lattice_->Stiffness();
   
   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][2];
   K[0][2] = Stiff[0][3];
   K[0][3] = Stiff[0][9];
   K[0][4] = Stiff[0][12];
   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][2];
   K[1][2] = Stiff[1][3];
   K[1][3] = Stiff[1][9];
   K[1][4] = Stiff[1][12];
   K[2][0] = Stiff[2][0];
   K[2][1] = Stiff[2][2];
   K[2][2] = Stiff[2][3];
   K[2][3] = Stiff[2][9];
   K[2][4] = Stiff[2][12];
   K[3][0] = Stiff[3][0];
   K[3][1] = Stiff[3][2];
   K[3][2] = Stiff[3][3];
   K[3][3] = Stiff[3][9];
   K[3][4] = Stiff[3][12];
   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][2];
   K[4][2] = Stiff[4][3];
   K[4][3] = Stiff[4][9];
   K[4][4] = Stiff[4][12];
   
   return K;
}
