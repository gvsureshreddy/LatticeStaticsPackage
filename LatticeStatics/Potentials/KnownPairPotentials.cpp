#include "KnownPairPotentials.h"

PairPotentials* InitializePairPotential(char *HashName,PerlInput &Input,int i,int j)
{
   return InitializePairPotential(Input.getHash(HashName),Input,i,j);
}

PairPotentials* InitializePairPotential(PerlInput::HashStruct ParentHash,PerlInput &Input,
                                        int i,int j)
{
   char tmp[LINELENGTH];
   double Eps0,Eps1,Sigma0,Sigma1,rcut;
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2,Cutoff;
   
   sprintf(tmp,"PotentialType_%u_%u",i,j);
   PerlInput::HashStruct Hash = Input.getHash(ParentHash,tmp);
   if (!strcmp("LJ",Input.getString(Hash,"Type")))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");

      return new LJ(Eps0,Eps1,Sigma0,Sigma1);
   }
   else if (!strcmp("LJCutoff",Input.getString(Hash,"Type")))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      Cutoff = Input.getDouble(Hash,"Cutoff");
      
      return new LJCutoff(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
   }
   else if (!strcmp("RadiiMorse",Input.getString(Hash,"Type")))
   {
      A0 = Input.getDouble(Hash,"A0");
      B0 = Input.getDouble(Hash,"B0");
      Alpha = Input.getDouble(Hash,"Alpha");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      
      return new RadiiMorse(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2);
   }
   else if (!strcmp("RadiiMorse2",Input.getString(Hash,"Type")))
   {
      A0 = Input.getDouble(Hash,"A0");
      B0 = Input.getDouble(Hash,"B0");
      Alpha = Input.getDouble(Hash,"Alpha");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      
      return new RadiiMorse2(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2);
   }
   else if (!strcmp("RadiiMorseCutoff",Input.getString(Hash,"Type")))
   {
      A0 = Input.getDouble(Hash,"A0");
      B0 = Input.getDouble(Hash,"B0");
      Alpha = Input.getDouble(Hash,"Alpha");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Cutoff = Input.getDouble(Hash,"Cutoff");

      return new RadiiMorseCutoff(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2,Cutoff);
   }
   else if (!strcmp("RadiiMorseCutoff2",Input.getString(Hash,"Type")))
   {
      A0 = Input.getDouble(Hash,"A0");
      B0 = Input.getDouble(Hash,"B0");
      Alpha = Input.getDouble(Hash,"Alpha");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Cutoff = Input.getDouble(Hash,"Cutoff");

      return new RadiiMorseCutoff2(A0,B0,Alpha,Rref1,Rref2,Rtheta1,Rtheta2,Cutoff);
   }
   else if (!strcmp("TempMorse",Input.getString(Hash,"Type")))
   {
      Tref = Input.getDouble(Hash,"Tref");
      A0 = Input.getDouble(Hash,"A0");
      B0 = Input.getDouble(Hash,"B0");
      Alpha = Input.getDouble(Hash,"Alpha");
      Rref1 = Input.getDouble(Hash,"Ref1");
      Tmelt = Input.getDouble(Hash,"Tmelt");
      
      return new TempMorse(A0,B0,Alpha,Rref1,Tref,Tmelt);
   }
   else if (!strcmp("Dobson",Input.getString(Hash,"Type")))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      rcut = Input.getDouble(Hash,"rcut");
      
      return new Dobson(Eps0,Eps1,Sigma0,Sigma1,rcut);
   }
   else
   {
      cerr << "Unknown Potential Type " << "\n";
      exit(-1);
   }
   
   return NULL;
}


void UpdatePairPotential(char *HashName,PerlInput &Input,int i,int j,
                         PairPotentials *Potential)
{return UpdatePairPotential(Input.getHash(HashName),Input,i,j,Potential);}

void UpdatePairPotential(PerlInput::HashStruct ParentHash,PerlInput &Input,int i,int j,
                         PairPotentials *Potential)
{
   char tmp[LINELENGTH];
   double Eps0,Eps1,Sigma0,Sigma1,rcut;
   double Tref,A0,B0,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2,Cutoff;

   sprintf(tmp,"PotentialType_%u_%u",i,j);
   PerlInput::HashStruct Hash = Input.getHash(ParentHash,tmp);
   if (!strcmp(Potential->Type(),"LJ"))
   {
      LJ *LJp;
      LJp = dynamic_cast<LJ *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJp->SetEps0(Eps0);
      }

      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJp->SetSigma0(Sigma0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJp->SetSigma1(Sigma1);
      }
   }
   else if (!strcmp(Potential->Type(),"LJCutoff"))
   {
      LJCutoff *LJCutoffp;
      LJCutoffp = dynamic_cast<LJCutoff *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJCutoffp->SetEps0(Eps0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJCutoffp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJCutoffp->SetSigma0(Sigma0);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJCutoffp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         LJCutoffp->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"RadiiMorse"))
   {
      RadiiMorse *RM;
      RM = dynamic_cast<RadiiMorse *>(Potential);
      
      if (Input.ParameterOK(Hash,"Update-A0"))
      {
         A0 = Input.getDouble(Hash,"Update-A0");
         RM->SetA0(A0);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-Alpha"))
      {
         Alpha = Input.getDouble(Hash,"Update-Alpha");
         RM->SetAlpha(Alpha);
      }

      if (Input.ParameterOK(Hash,"Update-Rref1"))
      {
         Rref1 = Input.getDouble(Hash,"Update-Rref1");
         RM->SetRref1(Rref1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta1"))
      {
         Rtheta1 = Input.getDouble(Hash,"Update-Rtheta1");
         RM->SetRtheta1(Rtheta1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rref2"))
      {
         Rref2 = Input.getDouble(Hash,"Update-Rref2");
         RM->SetRref2(Rref2);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta2"))
      {
         Rtheta2 = Input.getDouble(Hash,"Update-Rtheta2");
         RM->SetRtheta2(Rtheta2);
      }
   }
   else if (!strcmp(Potential->Type(),"RadiiMorse2"))
   {
      RadiiMorse2 *RM;
      RM = dynamic_cast<RadiiMorse2 *>(Potential);
      
      if (Input.ParameterOK(Hash,"Update-A0"))
      {
         A0 = Input.getDouble(Hash,"Update-A0");
         RM->SetA0(A0);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-Alpha"))
      {
         Alpha = Input.getDouble(Hash,"Update-Alpha");
         RM->SetAlpha(Alpha);
      }

      if (Input.ParameterOK(Hash,"Update-Rref1"))
      {
         Rref1 = Input.getDouble(Hash,"Update-Rref1");
         RM->SetRref1(Rref1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta1"))
      {
         Rtheta1 = Input.getDouble(Hash,"Update-Rtheta1");
         RM->SetRtheta1(Rtheta1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rref2"))
      {
         Rref2 = Input.getDouble(Hash,"Update-Rref2");
         RM->SetRref2(Rref2);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta2"))
      {
         Rtheta2 = Input.getDouble(Hash,"Update-Rtheta2");
         RM->SetRtheta2(Rtheta2);
      }
   }
   else if (!strcmp(Potential->Type(),"RadiiMorseCutoff"))
   {
      RadiiMorseCutoff *RM;
      RM = dynamic_cast<RadiiMorseCutoff *>(Potential);
      
      if (Input.ParameterOK(Hash,"Update-A0"))
      {
         A0 = Input.getDouble(Hash,"Update-A0");
         RM->SetA0(A0);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-Alpha"))
      {
         Alpha = Input.getDouble(Hash,"Update-Alpha");
         RM->SetAlpha(Alpha);
      }

      if (Input.ParameterOK(Hash,"Update-Rref1"))
      {
         Rref1 = Input.getDouble(Hash,"Update-Rref1");
         RM->SetRref1(Rref1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta1"))
      {
         Rtheta1 = Input.getDouble(Hash,"Update-Rtheta1");
         RM->SetRtheta1(Rtheta1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rref2"))
      {
         Rref2 = Input.getDouble(Hash,"Update-Rref2");
         RM->SetRref2(Rref2);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta2"))
      {
         Rtheta2 = Input.getDouble(Hash,"Update-Rtheta2");
         RM->SetRtheta2(Rtheta2);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         RM->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"RadiiMorseCutoff2"))
   {
      RadiiMorseCutoff2 *RM;
      RM = dynamic_cast<RadiiMorseCutoff2 *>(Potential);
      
      if (Input.ParameterOK(Hash,"Update-A0"))
      {
         A0 = Input.getDouble(Hash,"Update-A0");
         RM->SetA0(A0);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-Alpha"))
      {
         Alpha = Input.getDouble(Hash,"Update-Alpha");
         RM->SetAlpha(Alpha);
      }

      if (Input.ParameterOK(Hash,"Update-Rref1"))
      {
         Rref1 = Input.getDouble(Hash,"Update-Rref1");
         RM->SetRref1(Rref1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta1"))
      {
         Rtheta1 = Input.getDouble(Hash,"Update-Rtheta1");
         RM->SetRtheta1(Rtheta1);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rref2"))
      {
         Rref2 = Input.getDouble(Hash,"Update-Rref2");
         RM->SetRref2(Rref2);
      }
      
      if (Input.ParameterOK(Hash,"Update-Rtheta2"))
      {
         Rtheta2 = Input.getDouble(Hash,"Update-Rtheta2");
         RM->SetRtheta2(Rtheta2);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         RM->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"TempMorse"))
   {
      TempMorse *TM;
      TM = dynamic_cast<TempMorse *>(Potential);
      
      if (Input.ParameterOK(Hash,"Update-Tref"))
      {
         Tref = Input.getDouble(Hash,"Update-Tref");
         TM->SetTref(Tref);
      }

      if (Input.ParameterOK(Hash,"Update-A0"))
      {
         A0 = Input.getDouble(Hash,"Update-A0");
         TM->SetA0(A0);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         TM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-Alpha"))
      {
         Alpha = Input.getDouble(Hash,"Update-Alpha");
         TM->SetAlpha(Alpha);
      }

      if (Input.ParameterOK(Hash,"Update-Rref"))
      {
         Rref1 = Input.getDouble(Hash,"Update-Rref");
         TM->SetRref(Rref1);
      }

      if (Input.ParameterOK(Hash,"Update-Tmelt"))
      {
         Tmelt = Input.getDouble(Hash,"Update-Tmelt");
         TM->SetTmelt(Tmelt);
      }
   }
   else if (!strcmp(Potential->Type(),"Dobson"))
   {
      Dobson *Dobsonp;
      Dobsonp = dynamic_cast<Dobson *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         Dobsonp->SetEps0(Eps0);
      }

      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         Dobsonp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         Dobsonp->SetSigma0(Sigma0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         Dobsonp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-rcut"))
      {
         rcut = Input.getDouble(Hash,"Update-rcut");
         Dobsonp->Setrcut(rcut);
      }
   }
   else
   {
      cerr << "Unknown Potential Type " << "\n";
      exit(-1);
   }
}
