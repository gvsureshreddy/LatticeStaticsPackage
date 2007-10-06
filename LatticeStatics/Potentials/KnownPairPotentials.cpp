#include "KnownPairPotentials.h"

#include "UtilityFunctions.h"

const char *POTENTIALNAMES[] = {"LJ","LJCutoff","RadiiMorse","RadiiMorse2","RadiiMorseCutoff",
				"TempMorse","Dobson"};

PairPotentials* InitializePairPotential(char *datafile,const char *prefix,int i,int j)
{
   char tmp[LINELENGTH];
   double Eps0,Eps1,Sigma0,Sigma1,rcut;
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2,Cutoff;
   
   sprintf(tmp,"PotentialType_%u_%u",i,j);
   switch (GetStringParameter(prefix,tmp,datafile,POTENTIALNAMES,NOPOTENTIALS))
   {
      case 0:
      {
	 sprintf(tmp,"Eps0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Eps0)) exit(-1);
	 sprintf(tmp,"Eps1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Eps1)) exit(-1);
	 sprintf(tmp,"Sigma0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Sigma0)) exit(-1);
	 sprintf(tmp,"Sigma1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Sigma1)) exit(-1);
	 
	 return new LJ(Eps0,Eps1,Sigma0,Sigma1);
      }
      case 1:
      {
	 sprintf(tmp,"Eps0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Eps0)) exit(-1);
	 sprintf(tmp,"Eps1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Eps1)) exit(-1);
	 sprintf(tmp,"Sigma0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Sigma0)) exit(-1);
	 sprintf(tmp,"Sigma1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Sigma1)) exit(-1);
	 sprintf(tmp,"Cutoff_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Cutoff)) exit(-1);

	 return new LJCutoff(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
      }
      case 2:
	 {
	 sprintf(tmp,"A0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&A0)) exit(-1);
	 sprintf(tmp,"B0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&B0)) exit(-1);
	 sprintf(tmp,"Alpha_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Alpha)) exit(-1);
	 sprintf(tmp,"Rref1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref1)) exit(-1);
	 sprintf(tmp,"Rtheta1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rtheta1)) exit(-1);
	 sprintf(tmp,"Rref2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref2)) exit(-1);
	 sprintf(tmp,"Rtheta2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rtheta2)) exit(-1);
	 
	 return new RadiiMorse(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2);
	 }
      case 3:
	 {
	 sprintf(tmp,"A0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&A0)) exit(-1);
	 sprintf(tmp,"B0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&B0)) exit(-1);
	 sprintf(tmp,"Alpha_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Alpha)) exit(-1);
	 sprintf(tmp,"Rref1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref1)) exit(-1);
	 sprintf(tmp,"Rtheta1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rtheta1)) exit(-1);
	 sprintf(tmp,"Rref2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref2)) exit(-1);
	 sprintf(tmp,"Rtheta2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rtheta2)) exit(-1);
	 
	 return new RadiiMorse2(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2);
      }
      case 4:
	 {
	 sprintf(tmp,"A0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&A0)) exit(-1);
	 sprintf(tmp,"B0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&B0)) exit(-1);
	 sprintf(tmp,"Alpha_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Alpha)) exit(-1);
	 sprintf(tmp,"Rref1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref1)) exit(-1);
	 sprintf(tmp,"Rtheta1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rtheta1)) exit(-1);
	 sprintf(tmp,"Rref2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref2)) exit(-1);
	 sprintf(tmp,"Rtheta2_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rtheta2)) exit(-1);
	 sprintf(tmp,"Cutoff_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Cutoff)) exit(-1);
	 
	 return new RadiiMorseCutoff(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2,Cutoff);
      }
      case 5:
      {
	 sprintf(tmp,"Tref");
	 if(!GetParameter(prefix,tmp,datafile,'l',&Tref)) exit(-1);
	 sprintf(tmp,"A0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&A0)) exit(-1);
	 sprintf(tmp,"B0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&B0)) exit(-1);
	 sprintf(tmp,"Alpha_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Alpha)) exit(-1);
	 sprintf(tmp,"Rref_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Rref1)) exit(-1);
	 sprintf(tmp,"Tmelt_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Tmelt)) exit(-1);
	 
	 return new TempMorse(A0,B0,Alpha,Rref1,Tref,Tmelt);
      }
      case 6:
      {
	 sprintf(tmp,"Eps0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Eps0)) exit(-1);
	 sprintf(tmp,"Eps1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Eps1)) exit(-1);
	 sprintf(tmp,"Sigma0_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Sigma0)) exit(-1);
	 sprintf(tmp,"Sigma1_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&Sigma1)) exit(-1);
	 sprintf(tmp,"rcut_%u_%u",i,j);
	 if(!GetParameter(prefix,tmp,datafile,'l',&rcut)) exit(-1);
	 
	 return new Dobson(Eps0,Eps1,Sigma0,Sigma1,rcut);
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
   double Eps0,Eps1,Sigma0,Sigma1,rcut;
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2,Cutoff;
   
   if (!strcmp(Potential->Type(),"LJ"))
   {
      LJ *LJp;
      LJp = dynamic_cast<LJ *>(Potential);
      
      sprintf(tmp,"Eps0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Eps0,0))
	 LJp->SetEps0(Eps0);
      
      sprintf(tmp,"Eps1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Eps1,0))
	 LJp->SetEps1(Eps1);
      
      sprintf(tmp,"Sigma0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Sigma0,0))
	 LJp->SetSigma0(Sigma0);
      
      sprintf(tmp,"Sigma1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Sigma1,0))
	 LJp->SetSigma1(Sigma1);
   }
   else if (!strcmp(Potential->Type(),"LJCutoff"))
   {
      LJCutoff *LJCutoffp;
      LJCutoffp = dynamic_cast<LJCutoff *>(Potential);
      
      sprintf(tmp,"Eps0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Eps0,0))
	 LJCutoffp->SetEps0(Eps0);
      
      sprintf(tmp,"Eps1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Eps1,0))
	 LJCutoffp->SetEps1(Eps1);
      
      sprintf(tmp,"Sigma0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Sigma0,0))
	 LJCutoffp->SetSigma0(Sigma0);
      
      sprintf(tmp,"Sigma1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Sigma1,0))
	 LJCutoffp->SetSigma1(Sigma1);
      
      sprintf(tmp,"Cutoff_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Cutoff,0))
	 LJCutoffp->SetCutoff(Cutoff);
   }
   else if (!strcmp(Potential->Type(),"RadiiMorse"))
   {
      RadiiMorse *RM;
      RM = dynamic_cast<RadiiMorse *>(Potential);
      
      sprintf(tmp,"A0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&A0,0))
	 RM->SetA0(A0);
      
      sprintf(tmp,"B0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&B0,0))
	 RM->SetB0(B0);
      
      sprintf(tmp,"Alpha_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Alpha,0))
	 RM->SetAlpha(Alpha);
      
      sprintf(tmp,"Rref1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref1,0))
	 RM->SetRref1(Rref1);
      
      sprintf(tmp,"Rtheta1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rtheta1,0))
	 RM->SetRtheta1(Rtheta1);
      
      sprintf(tmp,"Rref2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref2,0))
	 RM->SetRref2(Rref2);
      
      sprintf(tmp,"Rtheta2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rtheta2,0))
	 RM->SetRtheta2(Rtheta2);
   }
   else if (!strcmp(Potential->Type(),"RadiiMorse2"))
   {
      RadiiMorse2 *RM;
      RM = dynamic_cast<RadiiMorse2 *>(Potential);
      
      sprintf(tmp,"A0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&A0,0))
	 RM->SetA0(A0);
      
      sprintf(tmp,"B0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&B0,0))
	 RM->SetB0(B0);
      
      sprintf(tmp,"Alpha_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Alpha,0))
	 RM->SetAlpha(Alpha);
      
      sprintf(tmp,"Rref1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref1,0))
	 RM->SetRref1(Rref1);
      
      sprintf(tmp,"Rtheta1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rtheta1,0))
	 RM->SetRtheta1(Rtheta1);
      
      sprintf(tmp,"Rref2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref2,0))
	 RM->SetRref2(Rref2);
      
      sprintf(tmp,"Rtheta2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rtheta2,0))
	 RM->SetRtheta2(Rtheta2);
   }
   else if (!strcmp(Potential->Type(),"RadiiMorseCutoff"))
   {
      RadiiMorseCutoff *RM;
      RM = dynamic_cast<RadiiMorseCutoff *>(Potential);
      
      sprintf(tmp,"A0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&A0,0))
	 RM->SetA0(A0);
      
      sprintf(tmp,"B0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&B0,0))
	 RM->SetB0(B0);
      
      sprintf(tmp,"Alpha_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Alpha,0))
	 RM->SetAlpha(Alpha);
      
      sprintf(tmp,"Rref1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref1,0))
	 RM->SetRref1(Rref1);
      
      sprintf(tmp,"Rtheta1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rtheta1,0))
	 RM->SetRtheta1(Rtheta1);
      
      sprintf(tmp,"Rref2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref2,0))
	 RM->SetRref2(Rref2);
      
      sprintf(tmp,"Rtheta2_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rtheta2,0))
	 RM->SetRtheta2(Rtheta2);

      sprintf(tmp,"Cutoff_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Cutoff,0))
	 RM->SetCutoff(Cutoff);
   }
   else if (!strcmp(Potential->Type(),"TempMorse"))
   {
      TempMorse *TM;
      TM = dynamic_cast<TempMorse *>(Potential);
      
      sprintf(tmp,"Tref");
      if(GetParameter(prefix,tmp,datafile,'l',&Tref,0))
	 TM->SetTref(Tref);
      
      sprintf(tmp,"A0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&A0,0))
	 TM->SetA0(A0);
      
      sprintf(tmp,"B0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&B0,0))
	 TM->SetB0(B0);
      
      sprintf(tmp,"Alpha_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Alpha,0))
	 TM->SetAlpha(Alpha);
      
      sprintf(tmp,"Rref_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Rref1,0))
	 TM->SetRref(Rref1);
      
      sprintf(tmp,"Tmelt_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Tmelt,0))
	 TM->SetTmelt(Tmelt);
   }
   else if (!strcmp(Potential->Type(),"Dobson"))
   {
      Dobson *Dobsonp;
      Dobsonp = dynamic_cast<Dobson *>(Potential);
      
      sprintf(tmp,"Eps0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Eps0,0))
	 Dobsonp->SetEps0(Eps0);
      
      sprintf(tmp,"Eps1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Eps1,0))
	 Dobsonp->SetEps1(Eps1);

      sprintf(tmp,"Sigma0_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Sigma0,0))
	 Dobsonp->SetSigma0(Sigma0);

      sprintf(tmp,"Sigma1_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&Sigma1,0))
	 Dobsonp->SetSigma1(Sigma1);

      sprintf(tmp,"rcut_%u_%u",i,j);
      if(GetParameter(prefix,tmp,datafile,'l',&rcut,0))
	 Dobsonp->Setrcut(rcut);
   }
   else
   {
      cerr << "Unknown Potential Type " << endl;
      exit(-1);
   }
}
