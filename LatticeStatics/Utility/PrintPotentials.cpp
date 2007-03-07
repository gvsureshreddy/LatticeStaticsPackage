#include "KnownPairPotentials.h"


int main(int argc, char *argv[])
{
   if (argc < 10)
   {
      cerr << "Usage: " << argv[0] << " inputfile Aatom Batom Tempstart Tempend temppts Rstart Rend points [print numerical derivatives?(1/0)]" <<endl;
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
      tempst=atof(argv[4]),
      tempend=atof(argv[5]),
      temppts=atof(argv[6]),
      rstart=atof(argv[7]),
      rend=atof(argv[8]),
      pts=atof(argv[9]);
   int deriv;
   if (argc == 11)
      deriv=atoi(argv[10]);
   else
      deriv = 0;

   pot = InitializePairPotential(argv[1],"^",i,j);

   cout << setiosflags(ios::fixed) << setprecision(10);
   double inc=(rend-rstart)/pts,
      tempinc=(tempend-tempst)/temppts;

   cout << "#" << setw(19) << "Temperature"
	<< setw(20) << "Radius2"
	<< setw(20) << "Potential"
	<< setw(20) << "Force";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "Stiffness";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D3Y";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D4Y";
   if (deriv) cout << setw(20) << "numerical";
   cout	<< setw(20) << "DT";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D2T";
   if (deriv) cout << setw(20) << "numerical";
   cout	<< setw(20) << "DYDT";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D2YDT";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D2YD2T";
   if (deriv) cout << setw(20) << "numerical";
   cout	<< setw(20) << "D3YDT";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D3YD2T";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D4YDT";
   if (deriv) cout << setw(20) << "numerical";
   cout << setw(20) << "D4YD2T";
   if (deriv) cout << setw(20) << "numerical";
   cout << endl;
   
   
   for (double t=tempst;t<=tempend;t+=tempinc)
   {
      for (double q=rstart;q<=rend;q+=inc)
      {
	 cout << setw(20) << t
	      << setw(20) << q
	      << setw(20) << pot->PairPotential(t,q)
	      << setw(20) << pot->PairPotential(t,q,PairPotentials::DY);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t,q+inc,PairPotentials::Y0,PairPotentials::T0)
		     -pot->PairPotential(t,q-inc,PairPotentials::Y0,PairPotentials::T0))
	       /(2.0*inc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D2Y);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t,q+inc,PairPotentials::DY,PairPotentials::T0)
		     -pot->PairPotential(t,q-inc,PairPotentials::DY,PairPotentials::T0))
	       /(2.0*inc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D3Y);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t,q+inc,PairPotentials::D2Y,PairPotentials::T0)
		     -pot->PairPotential(t,q-inc,PairPotentials::D2Y,PairPotentials::T0))
	       /(2.0*inc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D4Y);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t,q+inc,PairPotentials::D3Y,PairPotentials::T0)
		     -pot->PairPotential(t,q-inc,PairPotentials::D3Y,PairPotentials::T0))
	       /(2.0*inc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::Y0,PairPotentials::DT);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::Y0,PairPotentials::T0)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::Y0,PairPotentials::T0))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::Y0,PairPotentials::D2T);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::Y0,PairPotentials::DT)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::Y0,PairPotentials::DT))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::DY,PairPotentials::DT);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::DY,PairPotentials::T0)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::DY,PairPotentials::T0))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D2Y,PairPotentials::DT);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::D2Y,PairPotentials::T0)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::D2Y,PairPotentials::T0))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D2Y,PairPotentials::D2T);
	 if (deriv)
	    cout << setw(20)
		<< (pot->PairPotential(t+tempinc,q,PairPotentials::D2Y,PairPotentials::DT)
		    -pot->PairPotential(t-tempinc,q,PairPotentials::D2Y,PairPotentials::DT))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D3Y,PairPotentials::DT);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::D3Y,PairPotentials::T0)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::D3Y,PairPotentials::T0))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D3Y,PairPotentials::D2T);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::D3Y,PairPotentials::DT)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::D3Y,PairPotentials::DT))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D4Y,PairPotentials::DT);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::D4Y,PairPotentials::T0)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::D4Y,PairPotentials::T0))
	       /(2.0*tempinc);
	 cout << setw(20) << pot->PairPotential(t,q,PairPotentials::D4Y,PairPotentials::D2T);
	 if (deriv)
	    cout << setw(20)
		 << (pot->PairPotential(t+tempinc,q,PairPotentials::D4Y,PairPotentials::DT)
		     -pot->PairPotential(t-tempinc,q,PairPotentials::D4Y,PairPotentials::DT))
	       /(2.0*tempinc);
	 cout << endl;
      }
      cout << endl;
   }
   delete pot;
}
