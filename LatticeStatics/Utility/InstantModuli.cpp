#include "UtilityFunctions.h"
#include <iomanip.h>


int IND3D(int i,int j);

int main(int argc, char *argv[])
{
   Matrix F(3,3),L(6,6),L0(6,6);

   cout << "CATION:  This assumes L does not have multiple of 2 and 4 in it!!!" << endl;

   cin >> F;
   cout << endl << endl;
   cin >> L;

   double Jinv=1.0/F.Det();

   for (int i=0;i<3;i++)
      for (int j=i;j<3;j++)
	 for (int k=0;k<3;k++)
	    for (int l=k;l<3;l++)
	    {
	       L0[IND3D(i,j)][IND3D(k,l)] = 0;

	       for (int m=0;m<3;m++)
		  for (int n=0;n<3;n++)
		  {
		     L0[IND3D(i,j)][IND3D(k,l)] += Jinv*F[i][m]*F[k][n]*L[IND3D(m,j)][IND3D(n,l)];
		  }
	    }

   cout << setiosflags(ios::fixed) << setprecision(12);

   cout << setw(20) << L0 << endl;

   return 0;
}
