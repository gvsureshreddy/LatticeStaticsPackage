#include "ML_Rhombo.h"

ML_Rhombo::ML_Rhombo(Lattice *M)
{
   Lattice_ = (GenericLat *) M;
}

// Functions required by LatticeMode
Vector ML_Rhombo::ArcLenRHS(double DS,const Vector &Diff,
			  double Aspect)
{
   Vector rhs(3);
   Matrix stress = Lattice_->Stress();

   rhs[0] = stress[0][0] + stress[0][1] + stress[0][2];
   rhs[1] = stress[0][3] + stress[0][4] + stress[0][5];
   rhs[2] = DS*DS - Diff[2]*Diff[2]/(Aspect*Aspect)
      - Diff[0]*Diff[0]
      - Diff[1]*Diff[1];

   return rhs;
}

Vector ML_Rhombo::ArcLenDef()
{
   Vector def(3);
   Vector dof = Lattice_->DOF();

   def[0] = dof[0];
   def[1] = dof[3];
   def[2] = Lattice_->Temp();

   return def;
}

void ML_Rhombo::ArcLenUpdate(const Vector &newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[0]=DOF[1]=DOF[2] -= newval[0];
   DOF[3]=DOF[4]=DOF[5] -= newval[1];
   for (int i=6;i<DOF.Dim();i++)
      DOF[i] = 0.0;
   
   Lattice_->SetDOF(DOF);
   Lattice_->SetTemp(Lattice_->Temp() - newval[2]);
}

double ML_Rhombo::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[2] /= Aspect;
   New[2] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Vector ML_Rhombo::DrDt(const Vector &Diff)
{
   Vector ddt((Lattice_->DOF()).Dim(),0.0);

   ddt[0] = Diff[0]/Diff[2];
   ddt[1]=ddt[2]=ddt[0];
   ddt[3] = Diff[1]/Diff[2];
   ddt[4]=ddt[5]=ddt[3];
   
   return ddt;
}

Matrix ML_Rhombo::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(3,3);
   Matrix Stiff=Lattice_->Stiffness();
   Matrix stressdt = Lattice_->StressDT();

   K[0][0] = Stiff[0][0] + Stiff[0][1] + Stiff[0][2];
   K[0][1] = Stiff[0][3] + Stiff[0][4] + Stiff[0][5];
   K[0][2] = stressdt[0][0] + stressdt[0][1] + stressdt[0][2];
   K[1][0] = Stiff[3][0] + Stiff[3][1] + Stiff[3][2];
   K[1][1] = Stiff[3][3] + Stiff[3][4] + Stiff[3][5];
   K[1][2] = stressdt[0][3] + stressdt[0][4] + stressdt[0][5];
   K[2][0] = -2.0*Diff[0];
   K[2][1] = -2.0*Diff[1];
   K[2][2] = -2.0*Diff[2]/(Aspect*Aspect);
   
   return K;
}

double ML_Rhombo::ScanningDefParameter()
{
   return Lattice_->DOF()[3];
}

void ML_Rhombo::ScanningDefParamUpdate(const double newval)
{
   Vector DOF=Lattice_->DOF();

   DOF[3]=DOF[4]=DOF[5] -= newval;

   Lattice_->SetDOF(DOF);
}

double ML_Rhombo::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void ML_Rhombo::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double ML_Rhombo::ScanningStressParameter()
{
   Matrix stress = Lattice_->Stress();
   return stress[0][3] + stress[0][4] + stress[0][5];
}
   
Vector ML_Rhombo::ScanningRHS()
{
   Matrix stress = Lattice_->Stress();
   
   return Vector(1,stress[0][1]+stress[0][2]+stress[0][3]);
}

Vector ML_Rhombo::ScanningDef()
{
   return Vector(1,Lattice_->DOF()[0]);
}

void ML_Rhombo::ScanningUpdate(const Vector &newval) {}

Matrix ML_Rhombo::ScanningStiffness()
{
   Matrix stiff = Lattice_->Stiffness();
   
   return Matrix(1,1,stiff[0][0]+stiff[0][1]+stiff[0][2]);
}
