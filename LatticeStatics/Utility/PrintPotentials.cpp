#include "KnownPairPotentials.h"


int main(int argc, char *argv[])
{
   if (argc != 8)
   {
      cerr << "Usage: " << argv[0] << " inputfile Aatom Batom Temp Rstart Rend points" <<endl;
      exit(-1);
   }

   PairPotentials *pot;

   int
      i=atoi(argv[2]),
      j=atoi(argv[3]);
   int tmp;
   if (i>j)
   {
      tmp=i; i=j; j=tmp;
   }
   double
      temp=atof(argv[4]),
      rstart=atof(argv[5]),
      rend=atof(argv[6]),
      pts=atof(argv[7]);

   pot = InitializePairPotential(argv[1],"^",i,j);

   cout << setiosflags(ios::fixed) << setprecision(10);
   double inc=(rend-rstart)/pts;

   cout << "#" << setw(19) << "Temperature"
	<< setw(20) << "Radius"
	<< setw(20) << "Potential"
	<< setw(20) << "Force"
	<< setw(20) << "Stiffness" << endl;
   
   for (double q=rstart;q<=rend;q+=inc)
   {
      cout << setw(20) << temp
	   << setw(20) << q
	   << setw(20) << pot->PairPotential(temp,q*q)
	   << setw(20) << pot->PairPotential(temp,q*q,PairPotentials::DY)*(2.0*q)
	   << setw(20) << pot->PairPotential(temp,q*q,PairPotentials::D2Y)*(4.0*q*q)
	 + pot->PairPotential(temp,q*q,PairPotentials::DY)*(2.0)
	   << endl;
   }

   delete pot;
}
