#include "UniDefTemp3DFull.h"

UniDefTemp3DFull::UniDefTemp3DFull(UniDefTempLat *M)
{
   Lattice_ = M;
}

UniDefTemp3DFull::UniDefTemp3DFull(Lattice *M)
{
   Lattice_ = (UniDefTempLat*) M;
}

// Functions required by LatticeMode
Vector UniDefTemp3DFull::ArcLenRHS(double DS,const Vector &Diff,
				     double Aspect)
{
   Vector rhs(7);
   Matrix S = Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[1][1];
   rhs[2] = S[2][2];
   rhs[3] = S[0][1];
   rhs[4] = S[0][2];
   rhs[5] = S[1][2];
   rhs[6] = DS*DS - (Diff[6]*Diff[6])/(Aspect*Aspect)
      - Diff[0]*Diff[0] - Diff[1]*Diff[1] - Diff[2]*Diff[2] - Diff[3]*Diff[3]
      - Diff[4]*Diff[4] - Diff[5]*Diff[5];

   return rhs;
}

Vector UniDefTemp3DFull::ArcLenDef()
{
   Vector def(7);
   Matrix Def=Lattice_->DefGrad();

   def[0] = Def[0][0];
   def[1] = Def[1][1];
   def[2] = Def[2][2];
   def[3] = Def[0][1];
   def[4] = Def[0][2];
   def[5] = Def[1][2];
   def[6] = Lattice_->Temp();

   return def;
}

void UniDefTemp3DFull::ArcLenUpdate(const Vector &newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][0] -= newval[0];
   U[1][1] -= newval[1];
   U[2][2] -= newval[2];
   U[0][1]=U[1][0] -= newval[3];
   U[0][2]=U[2][0] -= newval[4];
   U[1][2]=U[2][1] -= newval[5];
   
   Lattice_->SetDefGrad(U);
   Lattice_->SetTemp(Lattice_->Temp() - newval[6]);
}

double UniDefTemp3DFull::ArcLenAngle(Vector Old,Vector New,double Aspect)
{
   Old[6] /= Aspect;
   New[6] /= Aspect;

   return fabs(acos( (Old*New)/(Old.Norm()*New.Norm()) ));
}

Matrix UniDefTemp3DFull::ArcLenStiffness(const Vector &Diff,double Aspect)
{
   Matrix K(7,7);
   Matrix Stiff=Lattice_->Stiffness(),
      SDT=Lattice_->StressDT();

   for (int i=0;i<6;i++)
      for (int j=0;j<3;j++)
      {
	 K[i][j] = Stiff[i][j];
      }

   for (int i=0;i<6;i++)
      for (int j=3;j<6;j++)
      {
	 K[i][j] = 2.0*Stiff[i][j];
      }

   K[0][6] = SDT[0][0];
   K[1][6] = SDT[1][1];
   K[2][6] = SDT[2][2];
   K[3][6] = SDT[0][1];
   K[4][6] = SDT[0][2];
   K[5][6] = SDT[1][2];

   for (int i=0;i<7;i++)
      K[6][i] = -2.0*Diff[i];

   return K;
}

double UniDefTemp3DFull::ScanningDefParameter()
{
   return Lattice_->DefGrad()[0][1];
}

void UniDefTemp3DFull::ScanningDefParamUpdate(const double newval)
{
   Matrix U=Lattice_->DefGrad();

   U[0][1]=U[1][0] -= newval;
      
   Lattice_->SetDefGrad(U);
}

double UniDefTemp3DFull::ScanningLoadParameter()
{
   return Lattice_->Temp();
}

void UniDefTemp3DFull::ScanningLoadParamUpdate(const double newval)
{
   Lattice_->SetTemp(Lattice_->Temp() - newval);
}

double UniDefTemp3DFull::ScanningStressParameter()
{
   return Lattice_->Stress()[0][1];
}
   
Vector UniDefTemp3DFull::ScanningRHS()
{
   Vector rhs(5);
   Matrix S=Lattice_->Stress();

   rhs[0] = S[0][0];
   rhs[1] = S[1][1];
   rhs[2] = S[2][2];
   rhs[3] = S[0][2];
   rhs[4] = S[1][2];

   return rhs;
}

Vector UniDefTemp3DFull::ScanningDef()
{
   Vector def(5);

   def[0] = Lattice_->DefGrad()[0][0];
   def[1] = Lattice_->DefGrad()[1][1];
   def[2] = Lattice_->DefGrad()[2][2];
   def[3] = Lattice_->DefGrad()[0][2];
   def[4] = Lattice_->DefGrad()[1][2];

   return def;
}

void UniDefTemp3DFull::ScanningUpdate(const Vector &newval)
{
   Matrix U = Lattice_->DefGrad();

   U[0][0] -= newval[0];
   U[1][1] -= newval[1];
   U[2][2] -= newval[2];
   U[2][0]=U[0][2] -= newval[3];
   U[2][1]=U[1][2] -= newval[4];
   
   Lattice_->SetDefGrad(U);
}

Matrix UniDefTemp3DFull::ScanningStiffness()
{
   Matrix K(5,5),
      Stiff = Lattice_->Stiffness();

   K[0][0] = Stiff[0][0];
   K[0][1] = Stiff[0][1];
   K[0][2] = Stiff[0][2];
   K[0][3] = 2.0*Stiff[0][4];
   K[0][4] = 2.0*Stiff[0][5];

   K[1][0] = Stiff[1][0];
   K[1][1] = Stiff[1][1];
   K[1][2] = Stiff[1][2];
   K[1][3] = 2.0*Stiff[1][4];
   K[1][4] = 2.0*Stiff[1][5];

   K[2][0] = Stiff[2][0];
   K[2][1] = Stiff[2][1];
   K[2][2] = Stiff[2][2];
   K[2][3] = 2.0*Stiff[2][4];
   K[2][4] = 2.0*Stiff[2][5];

   K[4][0] = Stiff[4][0];
   K[4][1] = Stiff[4][1];
   K[4][2] = Stiff[4][2];
   K[4][3] = 2.0*Stiff[4][4];
   K[4][4] = 2.0*Stiff[4][5];

   K[5][0] = Stiff[5][0];
   K[5][1] = Stiff[5][1];
   K[5][2] = Stiff[5][2];
   K[5][3] = 2.0*Stiff[5][4];
   K[5][4] = 2.0*Stiff[5][5];

   return K;
}
