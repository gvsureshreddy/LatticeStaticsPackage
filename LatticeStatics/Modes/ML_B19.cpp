#include "ML_B19.h"

ML_B19::ML_B19(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector ML_B19::ArcLenRHS(double DS,const Vector &Diff,
			  double Aspect)
{
   Vector rhs(7);
   Matrix stress = Lattice_->Stress();

   rhs[0] = stress[0][0];
   rhs[1] = stress[0][3];
   rhs[2] = stress[0][2];
   rhs[3] = stress[0][6];
   rhs[4] = stress[0][9];
   rhs[5] = stress[0][12];
   rhs[6] = DS*DS - Diff[6]*Diff[6]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1]
      - Diff[2]*Diff[2]
      - Diff[3]*Diff[3]
      - Diff[4]*Diff[4]
      - Diff[5]*Diff[5];

   return rhs;
}

Vector ML_B19::ArcLenDef()
{
   Vector def(7);
   Vector dof = Lattice_->DOF();

   def[0] = dof[0];
   def[1] = dof[3];
   def[2] = dof[2];
   def[3] = dof[6];
   def[4] = dof[9];
   def[5] = dof[12];
   def[6] = Lattice_->Temp();

   return def;
}

void ML_B19::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0] -= newval[0];
   DOF[1] = DOF[0];
   DOF[3] -= newval[1];
   DOF[2] -= newval[2];
   DOF[6] -= newval[3];
   DOF[9] -= newval[4];
   DOF[12] -= newval[5];
   
   for (int i=4;i<6;i++)
      DOF[i] = 0.0;
   for (int i=7;i<9;i++)
      DOF[i] = 0.0;
   for (int i=10;i<12;i++)
      DOF[i] = 0.0;
   for (int i=13;i<15;i++)
      DOF[i] = 0.0;
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[6]);
}

double ML_B19::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[6] /= Aspect;
   New[6] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector ML_B19::DrDt(const Vector &Diff)
{
   Vector ddt((Lattice_->DOF()).Dim(),0.0);

   ddt[0] = Diff[0]/Diff[6];
   ddt[1] = ddt[0];
   ddt[3] = Diff[1]/Diff[6];
   ddt[2] = Diff[2]/Diff[6];
   ddt[6] = Diff[3]/Diff[6];
   ddt[9] = Diff[4]/Diff[6];
   ddt[12] = Diff[5]/Diff[6];
   
   return ddt;
}

Matrix ML_B19::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(7,7);
   Matrix Stiff=Lattice_->Stiffness();
   Matrix stressdt = Lattice_->StressDT();


   K[0][0] = Stiff[0][0]+Stiff[0][1];
   K[0][1] = Stiff[0][3];
   K[0][2] = Stiff[0][2];
   K[0][3] = Stiff[0][6];
   K[0][4] = Stiff[0][9];
   K[0][5] = Stiff[0][12];
   K[0][6] = stressdt[0][0];

   K[1][0] = Stiff[3][0]+Stiff[3][1];
   K[1][1] = Stiff[3][3];
   K[1][2] = Stiff[3][2];
   K[1][3] = Stiff[3][6];
   K[1][4] = Stiff[3][9];
   K[1][5] = Stiff[3][12];
   K[1][6] = stressdt[0][3];

   K[2][0] = Stiff[2][0]+Stiff[2][1];
   K[2][1] = Stiff[2][3];
   K[2][2] = Stiff[2][2];
   K[2][3] = Stiff[2][6];
   K[2][4] = Stiff[2][9];
   K[2][5] = Stiff[2][12];
   K[2][6] = stressdt[0][2];

   K[3][0] = Stiff[6][0]+Stiff[6][1];
   K[3][1] = Stiff[6][3];
   K[3][2] = Stiff[6][2];
   K[3][3] = Stiff[6][6];
   K[3][4] = Stiff[6][9];
   K[3][5] = Stiff[6][12];
   K[3][6] = stressdt[0][6];
   
   K[4][0] = Stiff[9][0]+Stiff[9][1];
   K[4][1] = Stiff[9][3];
   K[4][2] = Stiff[9][2];
   K[4][3] = Stiff[9][6];
   K[4][4] = Stiff[9][9];
   K[4][5] = Stiff[9][12];
   K[4][6] = stressdt[0][9];

   K[5][0] = Stiff[12][0]+Stiff[12][1];
   K[5][1] = Stiff[12][3];
   K[5][2] = Stiff[12][2];
   K[5][3] = Stiff[12][6];
   K[5][4] = Stiff[12][9];
   K[5][5] = Stiff[12][12];
   K[5][6] = stressdt[0][12];
   
   K[6][0] = -2.0*Diff[0];
   K[6][1] = -2.0*Diff[1];
   K[6][2] = -2.0*Diff[2];
   K[6][3] = -2.0*Diff[3];
   K[6][4] = -2.0*Diff[4];
   K[6][5] = -2.0*Diff[5];
   K[6][6] = -2.0*Diff[6]/(Aspect*Aspect);
   
   return K;
}

double ML_B19::ScanningDefParameter()
{
   return Lattice_->DOF()[9];
}

void ML_B19::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[9] -= newval;

   Lattice_->SetDOF(DOF);
}

double ML_B19::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void ML_B19::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double ML_B19::ScanningStressParameter()
{
   Matrix stress = Lattice_->Stress();
   return stress[0][9];
}
   
Vector ML_B19::ScanningRHS()
{
   Matrix stress = Lattice_->Stress();
   Vector RHS(5,0.0);

   RHS[0] = stress[0][0];
   RHS[1] = stress[0][3];
   RHS[2] = stress[0][2];
   RHS[3] = stress[0][6];
   RHS[4] = stress[0][12];
   
   return RHS;
}

Vector ML_B19::ScanningDef()
{
   Vector DOF(5),
      dof = Lattice_->DOF();

   DOF[0] = dof[0];
   DOF[1] = dof[3];
   DOF[2] = dof[2];
   DOF[3] = dof[6];
   DOF[4] = dof[12];

   return DOF;
}

void ML_B19::ScanningUpdate(const Vector &newval)
{
   Vector dof=Lattice_->DOF();

   dof[0] -= newval[0];
   dof[1] = dof[0];
   dof[3] -= newval[1];
   dof[2] -= newval[2];
   dof[6] -= newval[3];
   dof[12] -= newval[4];

   Lattice_->SetDOF(dof);
}

Matrix ML_B19::ScanningStiffness()
{
   Matrix stiff = Lattice_->Stiffness(),
      K(5,5);

   K[0][0] = stiff[0][0]+stiff[0][1];
   K[0][1] = stiff[0][3];
   K[0][2] = stiff[0][2];
   K[0][3] = stiff[0][6];
   K[0][4] = stiff[0][12];

   K[1][0] = stiff[3][0]+stiff[3][1];
   K[1][1] = stiff[3][3];
   K[1][2] = stiff[3][2];
   K[1][3] = stiff[3][6];
   K[1][4] = stiff[3][12];

   K[2][0] = stiff[2][0]+stiff[2][1];
   K[2][1] = stiff[2][3];
   K[2][2] = stiff[2][2];
   K[2][3] = stiff[2][6];
   K[2][4] = stiff[2][12];

   K[3][0] = stiff[6][0]+stiff[6][1];
   K[3][1] = stiff[6][3];
   K[3][2] = stiff[6][2];
   K[3][3] = stiff[6][6];
   K[3][4] = stiff[6][12];

   K[4][0] = stiff[12][0]+stiff[12][1];
   K[4][1] = stiff[12][3];
   K[4][2] = stiff[12][2];
   K[4][3] = stiff[12][6];
   K[4][4] = stiff[12][12];

   return K;
}
