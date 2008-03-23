#include "Matrix.h"
#include "CMatrix.h"
#include "Vector.h"
#include "CVector.h"
#include "Vector3D.h"
#include "SparseMatrix.h"
#include <cstdlib>
#include <iomanip>
#include <ctime>

void sparsematrixtest();
Matrix project(int N);

int main()
{
   Matrix A(3,3),
      B(3,3),
      C(3,6),
      D(6,3),
      E(6,6),
      F(6,6);

   CMatrix CA(3,3),
      CB(3,3),
      CC(3,6),
      CD(6,3),
      CE(6,6),
      CF(6,6);   

   Vector a(3),
      b(3),
      c(6),
      d(6);

   CVector ca(3),
      cb(3),
      cc(6),
      cd(6);

   Vector3D e,f;
   
//   srand(time(NULL));
   srand(1);

//   Matrix::MathematicaPrintFlag=1;
   
   for (unsigned i=0;i<3;i++)
      for (unsigned j=0;j<3;j++)
      {
	 A[i][j] = double(rand()%10);
	 CA[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }

   for (unsigned i=0;i<3;i++)
      for (unsigned j=0;j<3;j++)
      {
	 B[i][j] = double(rand()%10);
	 CB[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }

   for (unsigned i=0;i<3;i++)
      for (unsigned j=0;j<6;j++)
      {
	 C[i][j] = double(rand()%10);
	 CC[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }
   
   for (unsigned i=0;i<6;i++)
      for (unsigned j=0;j<3;j++)
      {
	 D[i][j] = double(rand()%10);
	 CD[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }
   
   for (unsigned i=0;i<6;i++)
      for (unsigned j=0;j<6;j++)
      {
	 E[i][j] = double(rand()%10);
	 CE[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }
   
   for (unsigned i=0;i<3;i++)
   {
      a[i] = double(rand()%10);
      ca[i] = MyComplexDouble(rand()%10,rand()%10);
      b[i] = double(rand()%10);
      cb[i] = MyComplexDouble(rand()%10,rand()%10);
   }

   for (unsigned i=0;i<6;i++)
   {
      c[i] = double(rand()%10);
      cc[i] = MyComplexDouble(rand()%10,rand()%10);
      d[i] = double(rand()%10);
      cd[i] = MyComplexDouble(rand()%10,rand()%10);
   }

   for (unsigned i=0;i<3;i++)
   {
      e[i] = double(rand()%10);
      f[i] = double(rand()%10);
   }


   cout << setiosflags(ios::fixed) << setprecision(12);
//   cout << setiosflags(ios::scientific) << setprecision(12);

   cout << "A" << "\n"
	<< setw(20) << A << "\n"
	<< "CA" << "\n"
	<< setw(20) << CA << "\n"
	<< "B" << "\n"
	<< setw(20) << B << "\n"
	<< "CB" << "\n"
	<< setw(20) << CB << "\n"
	<< "C" << "\n"
	<< setw(20) << C << "\n"
	<< "CC" << "\n"
	<< setw(20) << CC << "\n"
	<< "D" << "\n"
	<< setw(20) << D << "\n"
	<< "CD" << "\n"
	<< setw(20) << CD << "\n"
	<< "a" << "\n"
	<< setw(20) << a << "\n"
	<< "ca" << "\n"
	<< setw(20) << ca << "\n"      
	<< "b" << "\n"
	<< setw(20) << b << "\n"
	<< "cb" << "\n"
	<< setw(20) << cb << "\n"
	<< "c" << "\n"
	<< setw(20) << c << "\n"
	<< "cc" << "\n"
	<< setw(20) << cc << "\n"
	<< "d" << "\n"
	<< setw(20) << d << "\n"
	<< "cd" << "\n"
	<< setw(20) << cd << "\n"
	<< "e" << "\n"
	<< setw(20) << e << "\n"
	<< "f" << "\n"
	<< setw(20) << f << "\n";

   Matrix q(A);

   cout << "q=" << "\n" << setw(20) << q <<"\n";

   q=B;

   cout << "q=B" << "\n" << setw(20) << q << "\n";

   cout << "A+B" << "\n" << setw(20) << A+B << "\n";

   cout << "CA+CB" << "\n" << setw(20) << CA+CB << "\n";

   cout << "A*B" << "\n" << setw(20) << A*B << "\n";

   cout << "CA*CB" << "\n" << setw(20) << CA*CB << "\n";

   cout << "D*C" << "\n" << setw(20) << D*C << "\n";

   cout << "D^T" << "\n" << setw(20) << D.Transpose() << "\n";

   cout << "CD^T" << "\n" << setw(20) << CD.Transpose() << "\n";

   cout << "Det(A)" << setw(20) << A.Det() << "\n";

   cout << "Det(CA)" << setw(20) << CA.Det() << "\n";

   cout << "A^-1" << setw(20) << A.Inverse() << "\n";

   cout << "A*A^-1" << setw(20) << A*A.Inverse() << "\n";

   cout << "CA^-1" << setw(20) << CA.Inverse() << "\n";

   cout << "CA*CA^-1" << setw(20) << CA*CA.Inverse() << "\n";

   PLU(A,B,C,D);
   PLU(CA,CB,CC,CD);

   cout << "A" << "\n" << setw(20) << A << "\n";

   cout << "P" << "\n" << setw(20) << B << "\n";

   cout << "L" << "\n" << setw(20) << C << "\n";

   cout << "U" << "\n" << setw(20) << D << "\n";

   cout << "PLU" << "\n" << setw(20) << B.Inverse()*C*D << "\n";

   cout << "CA" << "\n" << setw(20) << CA << "\n";

   cout << "CP" << "\n" << setw(20) << CB << "\n";

   cout << "CL" << "\n" << setw(20) << CC << "\n";

   cout << "CU" << "\n" << setw(20) << CD << "\n";

   cout << "CPLU" << "\n" << setw(20) << CB.Inverse()*CC*CD << "\n";

   cout << "E" << "\n" << setw(20) << E << "\n";

   cout << "c" << "\n" << setw(20) << c << "\n";

   cout << "x=E^-1*c" << setw(20) << SolvePLU(E,c) << "\n";
   cout << "E*x" << setw(20) << E*SolvePLU(E,c) << "\n";

   cout << "CE" << "\n" << setw(20) << CE << "\n";

   cout << "cc" << "\n" << setw(20) << cc << "\n";

   cout << "cx=CE^-1*cc" << setw(20) << SolvePLU(CE,cc) << "\n";
   cout << "CE*cx" << setw(20) << CE*SolvePLU(CE,cc) << "\n";

   SVD(A,B,C,D,MAXCONDITION,1);

   cout << "A" << "\n" << setw(20) << A << "\n";

   cout << "U" << "\n" << setw(20) << B << "\n";

   cout << "W" << "\n" << setw(20) << C << "\n";

   cout << "V" << "\n" << setw(20) << D << "\n";

   cout << "SVD" << "\n" << setw(20) << B*C*D.Transpose() << "\n";

   // check cholesky
   Matrix DD,UU;
   CMatrix CDD,CUU;
   // Make E symmetric
   E=E*E.Transpose();
   // Make CE Hermitian
   CE=CE*CE.ConjTrans();
   Cholesky(E,UU,DD);
   Cholesky(CE,CUU,CDD);
   cout << "Cholesky" << "\n" << setw(20) << E << "\n";
   cout << "result" << "\n" << setw(20) << UU.Transpose()*DD*UU << "\n";

   cout << "Complex Cholesky" << "\n" << setw(20) << CE << "\n";
   cout << "Complex result" << "\n" << setw(20)
	<< CUU.ConjTrans()*CDD*CUU << "\n";

   // check SymEigVal
   cout << "SymEigVal" << "\n" << setw(20) << SymEigVal(E) << "\n";

   // check SymEigVal eigenvector calculation
   d=SymEigVal(E,&F);
   cout << "SymEigVal(E,F)" << "\n" << setw(20) << d << "\n";
   cout << "F" << "\n" << setw(20) << F << "\n";
   cout << "F*F.Transpose()" << "\n" << setw(20) << F*F.Transpose() << "\n";

   for (unsigned i=0;i<6;i++)
   {
      for (unsigned j=0;j<6;j++)
      {
	 c[j] = F[j][i];
      }
      
      cout << setw(20) << E*c << "\n";
      cout << setw(20) << (E*c)/d[i] << "\n" << "\n";
   }

   // check HermiteigVal eigenvector calculation
   d=HermiteEigVal(CE,&CF);
   cout << "HermiteEigVal(CE,&CF)" << "\n" << setw(20) << d << "\n";
   cout << "CF" << "\n" << setw(20) << CF << "\n";
   cout << "CF.ConjTrans()*CF" << "\n" << setw(20) << CF.ConjTrans()*CF << "\n";

   for (unsigned i=0;i<6;i++)
   {
      for (unsigned j=0;j<6;j++)
      {
	 cc[j] = CF[j][i];
      }
      
      cout << setw(20) << CE*cc << "\n";
      cout << setw(20) << (CE*cc)/d[i] << "\n" << "\n";
   }

   E.Resize(3,3);
   E=A;


   // Create an ill conditioned matrix !!!!!
   A[0][0] = 1.0;
   A[0][1] = 0.0;
   A[0][2] = 2.0;
   A[1][0] = 1.0 + 1.0e-17;
   A[1][1] = 0.0;
   A[1][2] = 2.0;
   A[2][0] = 5.0;
   A[2][1] = 6.0;
   A[2][2] = 4.0;

   SVD(A,B,C,D,MAXCONDITION,1);

   cout << "A" << "\n" << setw(20) << A << "\n";

   cout << "U" << "\n" << setw(20) << B << "\n";

   cout << "W" << "\n" << setw(20) << C << "\n";

   cout << "V" << "\n" << setw(20) << D << "\n";

   cout << "SVD" << "\n" << setw(20) << B*C*D.Transpose() << "\n";

   for (unsigned i=0;i<C.Rows();i++)
      if (C[i][i]) C[i][i] = 1.0/C[i][i];
   
   cout << "A^-1 via SVD" << "\n" << setw(20) << D*C*B.Transpose() << "\n"
	<< setw(20) << A*(D*C*B.Transpose()) << "\n";

   cout << "A^-1 via PLU" << "\n" << setw(20) << A.Inverse() << "\n"
	<< setw(20) << A*(A.Inverse()) << "\n";

   SVD(A,B,C,D,10e17,1);

   cout << "W" << "\n" << setw(20) << C << "\n";

   for (unsigned i=0;i<C.Rows();i++)
      if (C[i][i]) C[i][i] = 1.0/C[i][i];
   
   cout << "A^-1 via SVD" << "\n" << setw(20) << D*C*B.Transpose() << "\n"
	<< setw(20) << A*(D*C*B.Transpose()) << "\n";

   A=E;

   cout << "a" << setw(20) << a << "\n";

   cout << "b" << setw(20) << b << "\n";

   cout << "c" << setw(20) << c << "\n";

   cout << "d" << setw(20) << d << "\n";

   cout << "a+b" << setw(20) << a+b << "\n";

   cout << "axb" << setw(20) << a%b << "\n";

   cout << "a.b" << setw(20) << a*b << "\n";

   cout << "2a" << setw(20) << 2.0*a << "\n";

   cout << "a/5" << setw(20) << a/5.0 << "\n";

   a=b;

   cout << "a=b" << setw(20) << a << "\n";

   cout << "A*b" << setw(20) << A*b << "\n";


   cout << "ca" << setw(20) << ca << "\n";

   cout << "cb" << setw(20) << cb << "\n";

   cout << "cc" << setw(20) << cc << "\n";

   cout << "cd" << setw(20) << cd << "\n";

   cout << "ca+cb" << setw(20) << ca+cb << "\n";

   cout << "ca.cb" << setw(20) << ca*cb << "\n";

   cout << "2ca" << setw(20) << MyComplexDouble(2.0)*ca << "\n";

   cout << "ca/5" << setw(20) << ca/5.0 << "\n";

   ca=cb;

   cout << "ca=cb" << setw(20) << ca << "\n";

   cout << "CA*cb" << setw(20) << CA*cb << "\n";

   f=Vector3D(b);
   
   cout << "A*f" << setw(20) << A*f << "\n";

   Matrix zz(1,1);
   zz[0][0] = 4.0;
   cout << "zz.Inverse()" << setw(20) << zz.Inverse() << "\n";

   Matrix AA(3,4),
   Q(3,3),R(3,4),Z(3,4),I(3,3);
   Matrix BB(4,3),BQ(4,4),BR(4,3),BZ(4,3),BI(4,4);

   AA[0][0]=12.0;AA[0][1]=-51.0;AA[0][2]=4.0;AA[0][3]=0.0;
   AA[1][0]=6.0;AA[1][1]=167.0;AA[1][2]=-68.0;AA[1][3]=1.0;
   AA[2][0]=-4.0;AA[2][1]=24.0;AA[2][2]=-41.0;AA[2][3]=2.0;

   BB=AA.Transpose();

   QR(AA,Q,R);

   cout << "Test QR" << "\n";

   cout << "A" << "\n" << setw(20) << AA << "\n";

   cout << "Q" << "\n" << setw(20) << Q  << "\n";

   cout << "R" << "\n" << setw(20) << R << "\n";

   Z = Q*R;

   cout << "A=Q*R" << "\n" << setw(20) << Z << "\n";
  
   I=Q*Q.Transpose();
  
   cout << "I=Q*Q^T" << "\n" << setw(20) << I << "\n";

   QR(AA,Q,R,1);

   cout << "Test QR with A transpose" << "\n";

   cout << "A" << "\n" << setw(20) << AA << "\n";

   cout << "Q" << "\n" << setw(20) << Q  << "\n";

   cout << "R" << "\n" << setw(20) << R << "\n";

   Z.Resize(4,3);
   Z = Q*R;

   cout << "A=Q*R" << "\n" << setw(20) << Z << "\n";

   I.Resize(4,4);
   I=Q*Q.Transpose();
  
   cout << "I=Q*Q^T" << "\n" << setw(20) << I << "\n";

   QR(BB,BQ,BR);

   cout << "BB" << "\n" << setw(20) << BB << "\n";

   cout << "BQ" << "\n" << setw(20) << BQ  << "\n";

   cout << "BR" << "\n" << setw(20) << BR << "\n";
   
   BZ = BQ*BR;

   cout << "BA=BQ*BR" << "\n" << setw(20) << BZ << "\n";
  
   BI=BQ*BQ.Transpose();
  
   cout << "BI=BQ*BQ^T" << "\n" << setw(20) << BI << "\n";


   // sparse matrix stuff
   sparsematrixtest();

   return 1;
}


void sparsematrixtest()
{
   int N = 4;
   
   Matrix u(N+1,N,0);	
   u=project(N);
   
   SparseMatrix A(u);	
   Matrix C(N+1,N,1);
   Matrix G(N,N,1);
   Matrix H(N,N+1,1);
   SparseMatrix I(H);
   SparseMatrix D(C);
   SparseMatrix B(G);
   SparseMatrix F(10, A.Rows(), A.Cols());
   Vector V1(N,1);	
   Vector V2(N+1,1);
   Vector3D Z(1);
   Matrix M(3,3,1);
   SparseMatrix Y(M);
	
   cout << "u = " << "\n" << setw(20) << u << "\n";
	
   //ADDITION AND SUBTRACTION TEST
   cout << "+A = u" << setw(20) << +A << "\n";
   cout << "-A = -u" << setw(20) << -A << "\n";
	
   cout << "C = " <<"\n" << setw(20) << C << "\n";
   cout << "D = C" << "\n" << setw(20) << D <<"\n";
	
   // SparseMatrix A + Matrix C + Matrix C
   cout << "A+C+C = " << "\n" << setw(20) << A+C+C << "\n";
	
   // SparseMatrix D + SparseMatrix A + SparseMatrix D
   cout << "D+A+D = " << "\n" << setw(20) << D+A+D << "\n";
	
   // SparseMatrix D + SparseMatrix A - Matrix C 
   cout << "D+A-C = " << "\n" << setw(20) << D+A-C << "\n";
	
   // SparseMatrix D + SparseMatrix A - Matrix C 
   cout << "-C+A+D = " << "\n" << setw(20) << -C+A+D << "\n";
	
   // -SparseMatrix A + SparseMatrix D
   cout << "-A+D  = " << "\n" << setw(20) << -A+D << "\n";
	
   // -SparseMatrix + Matrix C
   cout << "-A+C  = " << "\n" << setw(20) << -A+C << "\n";
	
   // -SparseMatrix A - SparseMatrix D
   cout << "-A-D  = " << "\n" << setw(20) << -A-D << "\n";
	
   // -SparseMatrix - Matrix C
   cout << "-A-C  = " << "\n" << setw(20) << -A-C << "\n";
	

	
   //MATRIX MULTIPLICATION TEST
	
   //SparseMatrix A, Matrix H, SparseMatrix I
   cout << "N = " << N << "\n"<< "\n";
   cout << "A = " << "\n" << setw(20) << A <<"\n";
   cout << "A = u" << "\n" << setw(20) << ReverseSparse(A) <<"\n";
   cout << "I = H " << "\n" << setw(20) <<I<< "\n";
   cout << "H = I" << "\n" << setw (20) <<H << "\n";

   //SparseMatrix A * Double N
   cout << "N*A = " << "\n" << setw(20) << N*A <<"\n";
   cout << "A*N = " << "\n" << setw(20) << A*N <<"\n";
	
   //SparseMatrix A * SparseMatrix I
   cout << "A*I = " << "\n" <<setw(20) << A*I << "\n";
   //SparseMatrix A * Matrix H 
   cout << "A*H = " << "\n" <<setw(20) << A*H << "\n";
   //SparseMatrix I * SparseMatrix A
   cout << "I*A = " << "\n" <<setw(20) << I*A << "\n";
   //Matrix H * SparseMatrix A
   cout << "H*A = " << "\n" <<setw(20) << H*A << "\n";
	
   //SparseMatrix B * SparseMatrix B
   cout << "B*B = " << "\n" <<setw(20) << B*B << "\n";
   cout << "B*G = " << "\n" << setw(20) << B*G <<"\n";
   cout << "G*B = " << "\n" << setw(20) << G*B <<"\n";
	
	
   //VECTOR MULTIPLICATION TEST
   cout << "A = u " << "\n" <<setw(20) << u << "\n";
   cout << "V1 = " << "\n" << setw(20) << V1 << "\n"<<"\n"<<"\n"<<"\n";
   cout << "V2 = " << "\n" << setw(20) << V2 << "\n"<<"\n"<<"\n"<<"\n";	
   cout << "A * V1 = " <<"\n" << setw(20) << A*V1 << "\n"<<"\n"<<"\n"<<"\n";
   cout << "V2 * A = " <<"\n" << setw(20) << V2*A << "\n"<<"\n"<<"\n"<<"\n";
   cout << "Y = " << "\n" << setw(20) << ReverseSparse(Y) << "\n"<<"\n"<<"\n"<<"\n";
   cout << "Z = " << "\n" << setw(20) << Z << "\n"<<"\n"<<"\n"<<"\n";
   cout << "Y*Z = " <<"\n" << setw(20) << Y*Z << "\n"<<"\n"<<"\n"<<"\n";
   cout << "Z*Y = " <<"\n" << setw(20) << Z*Y << "\n"<<"\n"<<"\n"<<"\n";	
	
   //TRANSPOSE AND REVERSESPARSE TEST
   cout << "A = " << "\n" <<setw(20) << A << "\n"<< "\n"
	<< "A.Transpose = " << "\n" << setw(20) << A.Transpose() <<"\n";
   cout << "A = " << "\n" <<setw(20) << ReverseSparse(A) << "\n" << "\n"
	<< "A.Transpose = " << "\n" << setw(20) << ReverseSparse(A.Transpose())<<"\n";
	
   //Identity SparseMatrix
   cout << "Identity SparseMatrix= " <<"\n" << setw(20)<<SparseIdentity(6)<< "\n";
   cout << "Identity Matrix= " << "\n" << setw(20) << ReverseSparse(SparseIdentity(6)) <<"\n";
}

//Returns matrix of size Rows_ x Cols_
//Function does not take arguments. Size of array determined by Rows_ and Cols_
Matrix project(int N)
{
   int k=0;
   int h,i, j, p, q, g, t;
   double m, sum, y, R;
   
   Matrix u(N+1,N, 0);
   
   // The following generates first column vector [1 0 0 0...0]
   u[0][k]=1;
   for (h=1; h < N+1; h++)              
   {
      u[h][k]=0;
   }
   k=k+1;                               //Finishes first vector
   
	
   // The following generates floor(log2(N)) vectors with floor(double(N)/pow(2,(i+1)))
   // entry pairs
	
   // Generates floor(log2(N)) entry pairs
   for (i=0; i < floor((log10(double(N))/log10(2.0)));i++)
   {
      m =  floor(double(N)/pow(2.0,(i+1)));
      q=0;

      for (j=0; j<m;j++)  //generates floor(double(N)/pow(2,(i+1))) column vectors
      {
	 if (j==0)
	 {
	    sum=0;
	    u[0][k]=0;
	    q=2*j+1;
	    for (p=0;p<pow(2.0,i);p++)
	    {
	       u[q][k]=1;
	       sum = sum +pow(u[q][k],2.0);
	       q=q+1;
	    }
	    for (p=0;p<pow(2.0,i);p++)
	    {
	       u[q][k]=-1;
	       sum = sum +pow(u[q][k],2.0);
	       q=q+1;
	    }				
	    sum = sqrt(sum);
	    for(h=0;h<N+1;h++)
	    {
	       u[h][k]=u[h][k]/sum;
	    }

	    k=k+1;
	 }
	 else
	 {
	    sum =0; 
	    for (p=0; p<pow(2.0,i);p++)
	    {
	       u[q][k]=1;
	       sum = sum +pow(u[q][k],2.0);
	       q=q+1;
	    }
	    for (p=0; p<pow(2.0,i);p++)
	    {
	       u[q][k]=-1;
	       sum = sum +pow(u[q][k],2.0);
	       q=q+1;
	    }
	    sum = sqrt(sum);
	    for(h=0;h<N+1;h++)
	    {
	       u[h][k]=u[h][k]/sum;
	    }		
	    k=k+1;
	 }
	 
	 //Generates remainder column vectors
	 R= N - ((floor(double(N)/(pow(2.0,(i+1)))))*pow(2.0,(i+1)));
	 if (j==m-1 & R >= pow(2.0,i))
	 {
	    u[0][k]=0;
	    sum = 0;
	    for (g=0;g<(((pow(2.0,i))*2*(m)));g++)
	    {
	       u[g+1][k]=1;
	       sum = sum +pow(u[g+1][k],2.0);
	    }
	    y=-(g)/pow(2.0,i);
	    for (t=0;t<(pow(2.0,(i)));t++)
	    {
	       u[g+1][k]=y;
	       sum = sum +pow(u[g+1][k],2.0);
	       g=g+1;
	    }
	    sum = sqrt(sum);
	    for(h=0;h<N+1;h++)
	    {
	       u[h][k]=u[h][k]/sum;
	    }
	    k=k+1;
	 }
      } // closes column vector (j) index
   } //closes entry pair (i) index

   return u;
}     
