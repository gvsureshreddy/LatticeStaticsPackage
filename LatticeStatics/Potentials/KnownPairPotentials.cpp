#include "KnownPairPotentials.h"

#include "UtilityFunctions.h"

const char *POTENTIALNAMES[] = {"RadiiMorse","TempMorse"};

PairPotentials* InitializePairPotential(char *datafile,const char *prefix,int i,int j)
{
   char tmp[LINELENGTH];
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2;
   
   sprintf(tmp,"PotentialType_%u_%u",i,j);
   switch (GetStringParameter(prefix,tmp,datafile,POTENTIALNAMES,NOPOTENTIALS))
   {
      case 0:
      {
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
	 
	 return new RadiiMorse(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2);
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

void UpdatePairPotential(char *datafile,const char *prefix,int i,int j,
			 PairPotentials *Potential)
{
   char tmp[LINELENGTH];
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2;

   RadiiMorse *RM;
   RM = dynamic_cast<RadiiMorse *>(Potential);
   if (!strcmp(RM->Type(),"RadiiMorse"))
   {
      sprintf(tmp,"A0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&A0,0))
	 RM->SetA0(A0);
      
      sprintf(tmp,"B0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&B0,0))
	 RM->SetB0(B0);
      
      sprintf(tmp,"Alpha_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Alpha,0))
	 RM->SetAlpha(Alpha);
      
      sprintf(tmp,"Rref1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Rref1,0))
	 RM->SetRref1(Rref1);
      
      sprintf(tmp,"Rtheta1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Rtheta1,0))
	 RM->SetRtheta1(Rtheta1);
      
      sprintf(tmp,"Rref2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Rref2,0))
	 RM->SetRref2(Rref2);
      
      sprintf(tmp,"Rtheta2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Rtheta2,0))
	 RM->SetRtheta2(Rtheta2);
   }
   else if (!strcmp(Potential->Type(),"TempMorse"))
   {
      TempMorse *TM;
      TM = dynamic_cast<TempMorse *>(Potential);
      
      sprintf(tmp,"Tref",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Tref,0))
	 TM->SetTref(Tref);
      
      sprintf(tmp,"A0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&A0,0))
	 TM->SetA0(A0);
      
      sprintf(tmp,"B0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&B0,0))
	 TM->SetB0(B0);
      
      sprintf(tmp,"Alpha_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Alpha,0))
	 TM->SetAlpha(Alpha);
      
      sprintf(tmp,"Rref_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Rref1,0))
	 TM->SetRref(Rref1);
      
      sprintf(tmp,"Tmelt_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,"%lf",&Tmelt,0))
	 TM->SetTmelt(Tmelt);
   }
   else
   {
      cerr << "Unknown Potential Type " << endl;
      exit(-1);
   }
}
