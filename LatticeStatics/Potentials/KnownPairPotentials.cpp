#include "KnownPairPotentials.h"

#include "UtilityFunctions.h"

PairPotentials* InitializePairPotential(char *datafile,const char *prefix,int i,int j)
{
   const int NoPotentials = NOPOTENTIALS;
   const char *Potentials[]={"RadiiMorse","TempMorse"};

   char tmp[LINELENGTH];
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2;
   
   sprintf(tmp,"PotentialType_%u_%u",i,j);
   switch (GetStringParameter(prefix,tmp,datafile,Potentials,NoPotentials))
   {
      case 0:
      {
	 sprintf(tmp,"Tref",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Tref)) exit(-1);
	 sprintf(tmp,"A0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&A0)) exit(-1);
	 sprintf(tmp,"B0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&B0)) exit(-1);
	 sprintf(tmp,"Alpha_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Alpha)) exit(-1);
	 sprintf(tmp,"Rref1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rref1)) exit(-1);
	 sprintf(tmp,"Rtheta1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rtheta1)) exit(-1);
	 sprintf(tmp,"Rref2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rref2)) exit(-1);
	 sprintf(tmp,"Rtheta2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rtheta2)) exit(-1);
	 
	 return new RadiiMorse(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2,Tref);
      }
      case 1:
      {
	 sprintf(tmp,"Tref",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Tref)) exit(-1);
	 sprintf(tmp,"A0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&A0)) exit(-1);
	 sprintf(tmp,"B0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&B0)) exit(-1);
	 sprintf(tmp,"Alpha_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Alpha)) exit(-1);
	 sprintf(tmp,"Rref_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rref1)) exit(-1);
	 sprintf(tmp,"Tmelt_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Tmelt)) exit(-1);
	 
	 return new TempMorse(A0,B0,Alpha,Rref1,Tref,Tmelt);
      }
      case -1:
      {
	 cerr << "Unknown Potential Type " << endl;
	 exit(-1);
      }
   }

   return NULL;
}
