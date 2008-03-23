#include "UtilityFunctions.h"
#include "CBKinematics.h"
#include "SymLagrangeCB.h"
#include "LagrangeCB.h"
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
   if (argc != 3)
   {
      cerr << "Usage: " << argv[0] << " <Dx (deg)> <U or F>" << endl;
      exit(1);
   }
   
   double pi = 4.0*atan(1.0);
   double Dx = strtod(argv[1],NULL);

   Vector pos(3,0.0);
   Matrix Id(3,3);
   Matrix A;

   Id.SetIdentity();
   CBKinematics *CBK;
   
   if (!strcmp("U",argv[2]))
   {
      A.Resize(6,6);
      CBK = new SymLagrangeCB(1,Id,&pos);
   }
   else
   {
      A.Resize(9,9);
      CBK = new LagrangeCB(1,Id,&pos);
   }
   
   cin >> A;

   int rk;
   
//   rk = Rank1Convex3D(CBK,A,Dx*pi/180);
   rk = FullScanRank1Convex3D(CBK,A,Dx*pi/180);
   
   cout << "Rank1Convex:" << setw(20) << rk << endl;
   
   return 0;
}
