#include "Matrix.h"
#include "CMatrix.h"
#include "Vector.h"
#include "Vector3D.h"
#include <stdlib.h>
#include <iomanip.h>
#include <time.h>

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

   Vector3D e,f;
   
   srand(time(NULL));

//   Matrix::MathematicaPrintFlag=1;
   
   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
      {
	 A[i][j] = double(rand()%10);
	 CA[i][j] = complex<double>(rand()%10,rand()%10);
      }

   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
      {
	 B[i][j] = double(rand()%10);
	 CB[i][j] = complex<double>(rand()%10,rand()%10);
      }

   for (int i=0;i<3;i++)
      for (int j=0;j<6;j++)
      {
	 C[i][j] = double(rand()%10);
	 CC[i][j] = complex<double>(rand()%10,rand()%10);
      }
   
   for (int i=0;i<6;i++)
      for (int j=0;j<3;j++)
      {
	 D[i][j] = double(rand()%10);
	 CD[i][j] = complex<double>(rand()%10,rand()%10);
      }
   
   for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
      {
	 E[i][j] = double(rand()%10);
	 CE[i][j] = complex<double>(rand()%10,rand()%10);
      }
   
   for (int i=0;i<3;i++)
   {
      a[i] = double(rand()%10);
      b[i] = double(rand()%10);
   }

   for (int i=0;i<6;i++)
   {
      c[i] = double(rand()%10);
      d[i] = double(rand()%10);
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
	<< "b" << endl
	<< setw(20) << b << endl
	<< "c" << endl
	<< setw(20) << c << endl
	<< "d" << endl
	<< setw(20) << d << endl
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

   PLU(A,B,C,D);

   cout << "A" << endl << setw(20) << A << endl;

   cout << "P" << endl << setw(20) << B << endl;

   cout << "L" << endl << setw(20) << C << endl;

   cout << "U" << endl << setw(20) << D << endl;

   cout << "PLU" << endl << setw(20) << B.Inverse()*C*D << endl;

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
   // Make CE symmetric
   CE=CE*CE.Transpose();
   Cholesky(E,UU,DD);
   Cholesky(CE,CUU,CDD);
   cout << "Cholesky" << endl << setw(20) << E << endl;
   cout << "result" << endl << setw(20) << UU.Transpose()*DD*UU << endl;

   cout << "Complex Cholesky" << endl << setw(20) << CE << endl;
   cout << "Complex result" << endl << setw(20) << CUU.Transpose()*CDD*CUU << endl;

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

   f=Vector3D(b);
   
   cout << "A*f" << setw(20) << A*f << endl;

   return 1;
}

   
