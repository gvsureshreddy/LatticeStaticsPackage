#include <Matrix.h>
#include <Vector.h>
#include <stdlib.h>
#include <iomanip.h>
#include <time.h>

int main()
{
   Matrix A(3,3),
      B(3,3),
      C(3,6),
      D(6,3),
      E(6,6);

   Vector a(3),
      b(3),
      c(6),
      d(6);

   srand(time(NULL));
   
   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	 A[i][j] = double(rand()%10);

   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
	 B[i][j] = double(rand()%10);
   

   for (int i=0;i<3;i++)
      for (int j=0;j<6;j++)
	 C[i][j] = double(rand()%10);

   for (int i=0;i<6;i++)
      for (int j=0;j<3;j++)
	 D[i][j] = double(rand()%10);

   for (int i=0;i<6;i++)
      for (int j=0;j<6;j++)
	 E[i][j] = double(rand()%10);

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


   cout << setiosflags(ios::fixed) << setprecision(6);

   cout << setw(20) << A << endl
	<< setw(20) << B << endl
	<< setw(20) << C << endl
	<< setw(20) << D << endl
	<< setw(20) << a << endl
	<< setw(20) << b << endl
	<< setw(20) << c << endl
	<< setw(20) << d << endl;

   Matrix q(A);

   cout << "q=" << endl << setw(20) << q <<endl;

   q=B;

   cout << "q=B" << endl << setw(20) << q << endl;

   cout << "A+B" << endl << setw(20) << A+B << endl;

   cout << "A*B" << endl << setw(20) << A*B << endl;

   cout << "D*C" << endl << setw(20) << D*C << endl;

   cout << "D^T" << endl << setw(20) << D.Transpose() << endl;

   cout << "Det(A)" << setw(20) << A.Det() << endl;

   cout << "A^-1" << setw(20) << A.Inverse() << endl;

   cout << "A*A^-1" << setw(20) << A*A.Inverse() << endl;

   PLU(A,B,C,D);

   cout << "A" << endl << setw(20) << A << endl;

   cout << "P" << endl << setw(20) << B << endl;

   cout << "L" << endl << setw(20) << C << endl;

   cout << "U" << endl << setw(20) << D << endl;

   cout << "PLU" << endl << setw(20) << B.Inverse()*C*D << endl;

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

   return 1;
}

   
