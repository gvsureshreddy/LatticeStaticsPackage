#include "KnownPairPotentials.h"

#include "UtilityFunctions.h"

PairPotentials* InitializePairPotential(char *datafile,const char *prefix,int i,int j)
{
   const int NoPotentials = NOPOTENTIALS;
   const char *Potentials[]={"RadiiMorse","TempMorse"};

   char tmp[LINELENGTH];
   double Tref,A0,B0,Alpha,Rref,Tmelt;
   double Rtheta;
   
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
	 sprintf(tmp,"Rref_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rref)) exit(-1);
	 sprintf(tmp,"Rtheta_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rtheta)) exit(-1);
	 
	 return new RadiiMorse(A0,B0,Alpha,Rref,Rtheta,Tref);
      }
      break;
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
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Rref)) exit(-1);
	 sprintf(tmp,"Tmelt_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,"%lf",&Tmelt)) exit(-1);
	 
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
