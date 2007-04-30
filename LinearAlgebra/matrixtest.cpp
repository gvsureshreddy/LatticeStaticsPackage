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
   
   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
      {
	 A[i][j] = double(rand()%10);
	 CA[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }

   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
      {
	 B[i][j] = double(rand()%10);
	 CB[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }

   for (int i=0;i<3;i++)
      for (int j=0;j<6;j++)
      {
	 C[i][j] = double(rand()%10);
	 CC[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }
   
   for (int i=0;i<6;i++)
      for (int j=0;j<3;j++)
      {
	 D[i][j] = double(rand()%10);
	 CD[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }
   
   for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
      {
	 E[i][j] = double(rand()%10);
	 CE[i][j] = MyComplexDouble(rand()%10,rand()%10);
      }
   
   for (int i=0;i<3;i++)
   {
      a[i] = double(rand()%10);
      ca[i] = MyComplexDouble(rand()%10,rand()%10);
      b[i] = double(rand()%10);
      cb[i] = MyComplexDouble(rand()%10,rand()%10);
   }

   for (int i=0;i<6;i++)
   {
      c[i] = double(rand()%10);
      cc[i] = MyComplexDouble(rand()%10,rand()%10);
      d[i] = double(rand()%10);
      cd[i] = MyComplexDouble(rand()%10,rand()%10);
   }

   for (int i=0;i<3;i++)
   {
      e[i] = double(rand()%10);
      f[i] = double(rand()%10);
   }


   cout << setiosflags(ios::fixed) << setprecision(12);
//   cout << setiosflags(ios::scientific) << setprecision(12);

   cout << "A" << endl
	<< setw(20) << A << endl
	<< "CA" << endl
	<< setw(20) << CA << endl
	<< "B" << endl
	<< setw(20) << B << endl
	<< "CB" << endl
	<< setw(20) << CB << endl
	<< "C" << endl
	<< setw(20) << C << endl
	<< "CC" << endl
	<< setw(20) << CC << endl
	<< "D" << endl
	<< setw(20) << D << endl
	<< "CD" << endl
	<< setw(20) << CD << endl
	<< "a" << endl
	<< setw(20) << a << endl
	<< "ca" << endl
	<< setw(20) << ca << endl      
	<< "b" << endl
	<< setw(20) << b << endl
	<< "cb" << endl
	<< setw(20) << cb << endl
	<< "c" << endl
	<< setw(20) << c << endl
	<< "cc" << endl
	<< setw(20) << cc << endl
	<< "d" << endl
	<< setw(20) << d << endl
	<< "cd" << endl
	<< setw(20) << cd << endl
	<< "e" << endl
	<< setw(20) << e << endl
	<< "f" << endl
	<< setw(20) << f << endl;

   Matrix q(A);

   cout << "q=" << endl << setw(20) << q <<endl;

   q=B;

   cout << "q=B" << endl << setw(20) << q << endl;

   cout << "A+B" << endl << setw(20) << A+B << endl;

   cout << "CA+CB" << endl << setw(20) << CA+CB << endl;

   cout << "A*B" << endl << setw(20) << A*B << endl;

   cout << "CA*CB" << endl << setw(20) << CA*CB << endl;

   cout << "D*C" << endl << setw(20) << D*C << endl;

   cout << "D^T" << endl << setw(20) << D.Transpose() << endl;

   cout << "CD^T" << endl << setw(20) << CD.Transpose() << endl;

   cout << "Det(A)" << setw(20) << A.Det() << endl;

   cout << "Det(CA)" << setw(20) << CA.Det() << endl;

   cout << "A^-1" << setw(20) << A.Inverse() << endl;

   cout << "A*A^-1" << setw(20) << A*A.Inverse() << endl;

   cout << "CA^-1" << setw(20) << CA.Inverse() << endl;

   cout << "CA*CA^-1" << setw(20) << CA*CA.Inverse() << endl;

   PLU(A,B,C,D);
   PLU(CA,CB,CC,CD);

   cout << "A" << endl << setw(20) << A << endl;

   cout << "P" << endl << setw(20) << B << endl;

   cout << "L" << endl << setw(20) << C << endl;

   cout << "U" << endl << setw(20) << D << endl;

   cout << "PLU" << endl << setw(20) << B.Inverse()*C*D << endl;

   cout << "CA" << endl << setw(20) << CA << endl;

   cout << "CP" << endl << setw(20) << CB << endl;

   cout << "CL" << endl << setw(20) << CC << endl;

   cout << "CU" << endl << setw(20) << CD << endl;

   cout << "CPLU" << endl << setw(20) << CB.Inverse()*CC*CD << endl;

   cout << "E" << endl << setw(20) << E << endl;

   cout << "c" << endl << setw(20) << c << endl;

   cout << "x=E^-1*c" << setw(20) << SolvePLU(E,c) << endl;
   cout << "E*x" << setw(20) << E*SolvePLU(E,c) << endl;

   cout << "CE" << endl << setw(20) << CE << endl;

   cout << "cc" << endl << setw(20) << cc << endl;

   cout << "cx=CE^-1*cc" << setw(20) << SolvePLU(CE,cc) << endl;
   cout << "CE*cx" << setw(20) << CE*SolvePLU(CE,cc) << endl;

   SVD(A,B,C,D,MAXCONDITION,1);

   cout << "A" << endl << setw(20) << A << endl;

   cout << "U" << endl << setw(20) << B << endl;

   cout << "W" << endl << setw(20) << C << endl;

   cout << "V" << endl << setw(20) << D << endl;

   cout << "SVD" << endl << setw(20) << B*C*D.Transpose() << endl;

   // check cholesky
   Matrix DD,UU;
   CMatrix CDD,CUU;
   // Make E symmetric
   E=E*E.Transpose();
   // Make CE Hermitian
   CE=CE*CE.ConjTrans();
   Cholesky(E,UU,DD);
   Cholesky(CE,CUU,CDD);
   cout << "Cholesky" << endl << setw(20) << E << endl;
   cout << "result" << endl << setw(20) << UU.Transpose()*DD*UU << endl;

   cout << "Complex Cholesky" << endl << setw(20) << CE << endl;
   cout << "Complex result" << endl << setw(20)
	<< CUU.ConjTrans()*CDD*CUU << endl;

   // check SymEigVal
   cout << "SymEigVal" << endl << setw(20) << SymEigVal(E) << endl;

   // check SymEigVal eigenvector calculation
   d=SymEigVal(E,&F);
   cout << "SymEigVal(E,F)" << endl << setw(20) << d << endl;
   cout << "F" << endl << setw(20) << F << endl;
   cout << "F*F.Transpose()" << endl << setw(20) << F*F.Transpose() << endl;

   for (int i=0;i<6;i++)
   {
      for (int j=0;j<6;j++)
      {
	 c[j] = F[j][i];
      }
      
      cout << setw(20) << E*c << endl;
      cout << setw(20) << (E*c)/d[i] << endl << endl;
   }

   // check HermiteigVal eigenvector calculation
   d=HermiteEigVal(CE,&CF);
   cout << "HermiteEigVal(CE,&CF)" << endl << setw(20) << d << endl;
   cout << "CF" << endl << setw(20) << CF << endl;
   cout << "CF.ConjTrans()*CF" << endl << setw(20) << CF.ConjTrans()*CF << endl;

   for (int i=0;i<6;i++)
   {
      for (int j=0;j<6;j++)
      {
	 cc[j] = CF[j][i];
      }
      
      cout << setw(20) << CE*cc << endl;
      cout << setw(20) << (CE*cc)/d[i] << endl << endl;
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

   cout << "A" << endl << setw(20) << A << endl;

   cout << "U" << endl << setw(20) << B << endl;

   cout << "W" << endl << setw(20) << C << endl;

   cout << "V" << endl << setw(20) << D << endl;

   cout << "SVD" << endl << setw(20) << B*C*D.Transpose() << endl;

   for (int i=0;i<C.Rows();i++)
      if (C[i][i]) C[i][i] = 1.0/C[i][i];
   
   cout << "A^-1 via SVD" << endl << setw(20) << D*C*B.Transpose() << endl
	<< setw(20) << A*(D*C*B.Transpose()) << endl;

   cout << "A^-1 via PLU" << endl << setw(20) << A.Inverse() << endl
	<< setw(20) << A*(A.Inverse()) << endl;

   SVD(A,B,C,D,10e17,1);

   cout << "W" << endl << setw(20) << C << endl;

   for (int i=0;i<C.Rows();i++)
      if (C[i][i]) C[i][i] = 1.0/C[i][i];
   
   cout << "A^-1 via SVD" << endl << setw(20) << D*C*B.Transpose() << endl
	<< setw(20) << A*(D*C*B.Transpose()) << endl;

   A=E;

   cout << "a" << setw(20) << a << endl;

   cout << "b" << setw(20) << b << endl;

   cout << "c" << setw(20) << c << endl;

   cout << "d" << setw(20) << d << endl;

   cout << "a+b" << setw(20) << a+b << endl;

   cout << "axb" << setw(20) << a%b << endl;

   cout << "a.b" << setw(20) << a*b << endl;

   cout << "2a" << setw(20) << 2.0*a << endl;

   cout << "a/5" << setw(20) << a/5.0 << endl;

   a=b;

   cout << "a=b" << setw(20) << a << endl;

   cout << "A*b" << setw(20) << A*b << endl;


   cout << "ca" << setw(20) << ca << endl;

   cout << "cb" << setw(20) << cb << endl;

   cout << "cc" << setw(20) << cc << endl;

   cout << "cd" << setw(20) << cd << endl;

   cout << "ca+cb" << setw(20) << ca+cb << endl;

   cout << "ca.cb" << setw(20) << ca*cb << endl;

   cout << "2ca" << setw(20) << MyComplexDouble(2.0)*ca << endl;

   cout << "ca/5" << setw(20) << ca/5.0 << endl;

   ca=cb;

   cout << "ca=cb" << setw(20) << ca << endl;

   cout << "CA*cb" << setw(20) << CA*cb << endl;

   f=Vector3D(b);
   
   cout << "A*f" << setw(20) << A*f << endl;

   Matrix zz(1,1);
   zz[0][0] = 4.0;
   cout << "zz.Inverse()" << setw(20) << zz.Inverse() << endl;

   Matrix AA(3,4),
   Q(3,3),R(3,4),Z(3,4),I(3,3);
   Matrix BB(4,3),BQ(4,4),BR(4,3),BZ(4,3),BI(4,4);

   AA[0][0]=12.0;AA[0][1]=-51.0;AA[0][2]=4.0;AA[0][3]=0.0;
   AA[1][0]=6.0;AA[1][1]=167.0;AA[1][2]=-68.0;AA[1][3]=1.0;
   AA[2][0]=-4.0;AA[2][1]=24.0;AA[2][2]=-41.0;AA[2][3]=2.0;

   BB=AA.Transpose();

   QR(AA,Q,R);

   cout << "A" << endl << setw(20) << AA << endl;

   cout << "Q" << endl << setw(20) << Q  << endl;

   cout << "R" << endl << setw(20) << R << endl;
   
   Z = Q*R;

   cout << "A=Q*R" << endl << setw(20) << Z << endl;
  
   I=Q*Q.Transpose();
  
   cout << "I=Q*Q^T" << endl << setw(20) << I << endl;

   QR(BB,BQ,BR);

   cout << "BB" << endl << setw(20) << BB << endl;

   cout << "BQ" << endl << setw(20) << BQ  << endl;

   cout << "BR" << endl << setw(20) << BR << endl;
   
   BZ = BQ*BR;

   cout << "BA=BQ*BR" << endl << setw(20) << BZ << endl;
  
   BI=BQ*BQ.Transpose();
  
   cout << "BI=BQ*BQ^T" << endl << setw(20) << BI << endl;


   // sparse matrix stuff
   sparsematrixtest();

   return 1;
}


void sparsematrixtest()
{
   int count=0;
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
	
   cout << "u = " << endl << setw(20) << u << endl;
	
   //ADDITION AND SUBTRACTION TEST
   cout << "+A = u" << setw(20) << +A << endl;
   cout << "-A = -u" << setw(20) << -A << endl;
	
   cout << "C = " <<endl << setw(20) << C << endl;
   cout << "D = C" << endl << setw(20) << D <<endl;
	
   // SparseMatrix A + Matrix C + Matrix C
   cout << "A+C+C = " << endl << setw(20) << A+C+C << endl;
	
   // SparseMatrix D + SparseMatrix A + SparseMatrix D
   cout << "D+A+D = " << endl << setw(20) << D+A+D << endl;
	
   // SparseMatrix D + SparseMatrix A - Matrix C 
   cout << "D+A-C = " << endl << setw(20) << D+A-C << endl;
	
   // SparseMatrix D + SparseMatrix A - Matrix C 
   cout << "-C+A+D = " << endl << setw(20) << -C+A+D << endl;
	
   // -SparseMatrix A + SparseMatrix D
   cout << "-A+D  = " << endl << setw(20) << -A+D << endl;
	
   // -SparseMatrix + Matrix C
   cout << "-A+C  = " << endl << setw(20) << -A+C << endl;
	
   // -SparseMatrix A - SparseMatrix D
   cout << "-A-D  = " << endl << setw(20) << -A-D << endl;
	
   // -SparseMatrix - Matrix C
   cout << "-A-C  = " << endl << setw(20) << -A-C << endl;
	

	
   //MATRIX MULTIPLICATION TEST
	
   //SparseMatrix A, Matrix H, SparseMatrix I
   cout << "N = " << N << endl<< endl;
   cout << "A = " << endl << setw(20) << A <<endl;
   cout << "A = u" << endl << setw(20) << ReverseSparse(A) <<endl;
   cout << "I = H " << endl << setw(20) <<I<< endl;
   cout << "H = I" << endl << setw (20) <<H << endl;

   //SparseMatrix A * Double N
   cout << "N*A = " << endl << setw(20) << N*A <<endl;
   cout << "A*N = " << endl << setw(20) << A*N <<endl;
	
   //SparseMatrix A * SparseMatrix I
   cout << "A*I = " << endl <<setw(20) << A*I << endl;
   //SparseMatrix A * Matrix H 
   cout << "A*H = " << endl <<setw(20) << A*H << endl;
   //SparseMatrix I * SparseMatrix A
   cout << "I*A = " << endl <<setw(20) << I*A << endl;
   //Matrix H * SparseMatrix A
   cout << "H*A = " << endl <<setw(20) << H*A << endl;
	
   //SparseMatrix B * SparseMatrix B
   cout << "B*B = " << endl <<setw(20) << B*B << endl;
   cout << "B*G = " << endl << setw(20) << B*G <<endl;
   cout << "G*B = " << endl << setw(20) << G*B <<endl;
	
	
   //VECTOR MULTIPLICATION TEST
   cout << "A = u " << endl <<setw(20) << u << endl;
   cout << "V1 = " << endl << setw(20) << V1 << endl<<endl<<endl<<endl;
   cout << "V2 = " << endl << setw(20) << V2 << endl<<endl<<endl<<endl;	
   cout << "A * V1 = " <<endl << setw(20) << A*V1 << endl<<endl<<endl<<endl;
   cout << "V2 * A = " <<endl << setw(20) << V2*A << endl<<endl<<endl<<endl;
   cout << "Y = " << endl << setw(20) << ReverseSparse(Y) << endl<<endl<<endl<<endl;
   cout << "Z = " << endl << setw(20) << Z << endl<<endl<<endl<<endl;
   cout << "Y*Z = " <<endl << setw(20) << Y*Z << endl<<endl<<endl<<endl;
   cout << "Z*Y = " <<endl << setw(20) << Z*Y << endl<<endl<<endl<<endl;	
	
   //TRANSPOSE AND REVERSESPARSE TEST
   cout << "A = " << endl <<setw(20) << A << endl<< endl
	<< "A.Transpose = " << endl << setw(20) << A.Transpose() <<endl;
   cout << "A = " << endl <<setw(20) << ReverseSparse(A) << endl << endl
	<< "A.Transpose = " << endl << setw(20) << ReverseSparse(A.Transpose())<<endl;
	
   //Identity SparseMatrix
   cout << "Identity SparseMatrix= " <<endl << setw(20)<<SparseIdentity(6)<< endl;
   cout << "Identity Matrix= " << endl << setw(20) << ReverseSparse(SparseIdentity(6)) <<endl;
}

int pow(int a,int b)
{
   int c=1;
   for (int i=0;i<b;i++)
   {
      c *= a;
   }
   return c;
}

//Returns matrix of size Rows_ x Cols_
//Function does not take arguments. Size of array determined by Rows_ and Cols_
Matrix project(int N)
{
   int k=0;
   int l=0;
   int h,i, j, p, q, g, t;
   double m, sum, x, y, R;
   int count=0;	
   
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
      m =  floor(double(N)/pow(2,(i+1)));
      q=0;

      for (j=0; j<m;j++)  //generates floor(double(N)/pow(2,(i+1))) column vectors
      {
	 if (j==0)
	 {
	    sum=0;
	    u[0][k]=0;
	    q=2*j+1;
	    for (p=0;p<pow(2,i);p++)
	    {
	       u[q][k]=1;
	       sum = sum +pow(u[q][k],2);
	       q=q+1;
	    }
	    for (p=0;p<pow(2,i);p++)
	    {
	       u[q][k]=-1;
	       sum = sum +pow(u[q][k],2);
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
	    for (p=0; p<pow(2,i);p++)
	    {
	       u[q][k]=1;
	       sum = sum +pow(u[q][k],2);
	       q=q+1;
	    }
	    for (p=0; p<pow(2,i);p++)
	    {
	       u[q][k]=-1;
	       sum = sum +pow(u[q][k],2);
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
	 R= N - ((floor(double(N)/(pow(2,(i+1)))))*pow(2,(i+1)));
	 if (j==m-1 & R >= pow(2,i))
	 {
	    u[0][k]=0;
	    sum = 0;
	    for (g=0;g<(((pow(2,i))*2*(m)));g++)
	    {
	       u[g+1][k]=1;
	       sum = sum +pow(u[g+1][k],2);
	    }
	    y=-(g)/pow(2,i);
	    for (t=0;t<(pow(2,(i)));t++)
	    {
	       u[g+1][k]=y;
	       sum = sum +pow(u[g+1][k],2);
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
