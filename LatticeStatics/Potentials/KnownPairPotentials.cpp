#include "KnownPairPotentials.h"

#include "UtilityFunctions.h"

PairPotentials* InitializePairPotential(char *datafile,int i,int j)
{
   const int NoPotentials = NOPOTENTIALS;
   const char *Potentials[]={"RadiiMorse","TempMorse"};

   char tmp[LINELENGTH];
   double Tref,A0,B0,Alpha,Rref,Tmelt;
   double Rtheta;
   
   sprintf(tmp,"^PotentialType_%u_%u",i,j);
   switch (GetStringParameter(tmp,datafile,Potentials,NoPotentials))
   {
      case 0:
      {
	 sprintf(tmp,"^Tref",i,j);
	 GetParameter(tmp,datafile,"%lf",&Tref);
	 sprintf(tmp,"^A0_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&A0);
	 sprintf(tmp,"^B0_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&B0);
	 sprintf(tmp,"^Alpha_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Alpha);
	 sprintf(tmp,"^Rref_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Rref);
	 sprintf(tmp,"^Rtheta_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Rtheta);
	 sprintf(tmp,"^Tmelt_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Tmelt);
	 
	 return new RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref,Tmelt);
      }
      break;
      case 1:
      {
	 sprintf(tmp,"^Tref",i,j);
	 GetParameter(tmp,datafile,"%lf",&Tref);
	 sprintf(tmp,"^A0_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&A0);
	 sprintf(tmp,"^B0_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&B0);
	 sprintf(tmp,"^Alpha_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Alpha);
	 sprintf(tmp,"^Rref_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Rref);
	 sprintf(tmp,"^Tmelt_%u_%u",i,j);
	 GetParameter(tmp,datafile,"%lf",&Tmelt);
	 
	 return new TempMorse(A0,B0,Alpha,Rref,Tref,Tmelt);
      }
      break;
      case -1:
      {
	 cerr << "Unknown Potential Type " << endl;
	 exit(-1);
      }
   }

   return NULL;
}
