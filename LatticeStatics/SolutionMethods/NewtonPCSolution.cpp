#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "UtilityFunctions.h"
#include <cmath>

using namespace std;


NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   const Vector &one,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   // Set Lattice to solution "one"
   Mode_->SetModeDOF(one);
}

int NewtonPCSolution::AllSolutionsFound()
{
   if (CurrentSolution_ < 60)
   { 
      return 0;
   }
   else
   {
      return 1;
   }
   
}

NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   char *startfile,fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{

   const char *StartType[] = {"Bifurcation","Continuation","ConsistencyCheck"};
   switch (GetStringParameter(prefix,"StartType",startfile,StartType,3))
   {
      case -1:
      {
         cerr << "Unknown StartType!" << endl;
         exit(-1);
         break;
      }
      case 0:
      {

         break;
      }
      case 1:
      {
         // Get solution1
         Vector one(Mode_->ArcLenDef().Dim());
         cout << "yes";
         if(!GetVectorParameter(prefix,"Solution1",startfile,&one)) exit(-1);
         Mode_->SetModeDOF(one);

         break;
      }
      case 2:
      {
         break;
      }
   }
}


double NewtonPCSolution::FindNextSolution(int &good)
{
   //Finds the next solution 
   double norm=0,f=0;
   double delta=0,kappa=0,alpha=0,delta0=1.1,kappa0=1.1,alpha0=0.8;
   double d,k,a;
   //double uncertainity;
   //int itr = 0;
   if (CurrentSolution_ == 0)
   {
      h = 0.005;
   }
   int dim = Mode_->ModeDOF().Dim(),count =0;
   Vector Dxu(dim),Dxv(dim),Dxw(dim),Dx(dim),D1(dim),D2(dim);

   //Predictor Step
   Dxu = Mode_->ModeDOF();
   do
   {//Loop to find the solution again if the step size decreases by more than half
      Dxv = Dxu + h*tang(Dxu);
      //Corrector Step
      do
      {
	 Mode_->SetModeDOF(Dxv);   // Set mode DOF to Dxv
	 Dxw = Dxv - MPI(Dxv)*(Mode_->ModeForce());
	 Dx  = Dxw - Dxv;
	 norm = Norm(Dx); // remember Norm accepts only vectors not expressions in arguments
	 Dxv = Dxw;
	 if (count == 0)
	 {
	    count++;
	    D1 = MPI(Dxu)*Force(Dxw); D2 = MPI(Dxu)*Force(Dxv);
	    kappa = Norm(D1)/Norm(D2);
	    k = sqrt(kappa/kappa0);
	    D1 = MPI(Dxv)*Force(Dxv);
	    delta = Norm(D1);
	    d = sqrt(delta/delta0);
	    alpha = acos(tang(Dxu)*tang(Dxv));
	    a = alpha/alpha0;
	    f = max(k,max(d,a));
	    f = max(min(f,2.0),0.5);
	    h = h/f;
	 } 
      }
      while(norm>0.005);
    
   }
   while(f>2);
   CurrentSolution_++; 
   good = 1;
   Mode_->SetModeDOF(Dxv);    // Set mode DOF to Dxv
   return norm;
}

Vector NewtonPCSolution::Force(Vector &D)
{
   //Finds the force vector
   int N;
   N = Mode_->ModeDOF().Dim()-1;
   Vector D1(N+1),D2(N);
   D1 = Mode_->ModeDOF();
   Mode_->SetModeDOF(D);
   D2 = Mode_->ModeForce();
   Mode_->SetModeDOF(D1);
   return D2;
}

double NewtonPCSolution::Norm(Vector &D)
{
   //finds norm of a given vector
   int i,dim;
   double norm=0;
   dim = D.Dim();
   for (i=0;i<dim;i++)
   {
      norm = norm + D[i]*D[i];
   }
   return sqrt(norm);
}
	
Vector NewtonPCSolution::tang(Vector &D)
{
   //finds the tangent to the curve
   int j,k,N;
   double sum,norm=1;
   N=D.Dim()-1;
   Vector t(N+1);
   Matrix Q(N,N),R(N,N+1),A(N+1,N+1),Z(N,N+1);
   t[N]=1;
   Z = Mode_->ModeStiffness();
   QR(Z,Q,R);
   for (j=N;j>0;j--)
   {
      sum = 0;
      for (k=N+1;k>j;k--)
      {
	 sum = sum - R[j-1][k-1]*t[k-1];
      }
      t[j-1] = sum/R[j-1][j-1];
      norm = norm + t[j-1]*t[j-1];
   }
   norm = sqrt(norm);
   t = t/norm;  //normalizing tangent vector 

   //Adjusting the direction of tangent vector
   for (j=0;j<N;j++)
   {
      for (k=0;k<N+1;k++)
      {
	 A[j][k]= Z[j][k];
	 if (j == 0)
	 {
	    A[N][k]=t[k];
	 }
      }
   }
   if (A.Det()>0)
   {
      if (CurrentSolution_ > 45)
      {
	 t = -1*t;
      }
      else
      {
	 t = -1*t;
      }
   }
              
   return t;
}

Matrix NewtonPCSolution::MPI(Vector &D)
{
   //finds moore-penrose inverse of a stiffness matrix
   int N=D.Dim()-1;
   Matrix A(N,N+1),mpi(N+1,N);
   Vector D1(N+1);

   D1 = Mode_->ModeDOF();
   Mode_->SetModeDOF(D);
   //Finding Moore-Penrose Inverse
   A=Mode_->ModeStiffness();
   mpi = A.Transpose()*((A*A.Transpose()).Inverse());
   Mode_->SetModeDOF(D1);
   return mpi;
}

int NewtonPCSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,int Width,
				  fstream &out)
{
   return 1;
}
