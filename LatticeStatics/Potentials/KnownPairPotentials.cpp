#include "KnownPairPotentials.h"
#include <sstream>

PairPotentials* InitializePairPotential(char const* const HashName,PerlInput const& Input,
                                        int const& i,int const& j)
{
   return InitializePairPotential(Input.getHash(HashName),Input,i,j);
}

PairPotentials* InitializePairPotential(PerlInput::HashStruct const& ParentHash,
                                        PerlInput const& Input,int const& i,int const& j)
{
   stringstream tmp;
   double Eps0,Eps1,Sigma0,Sigma1,rcut;
   double Tref,A0,AT,B0,BT,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta1Pow,Rtheta2,Rtheta2Pow,Cutoff;

   tmp.str("");
   tmp << "PotentialType_" << i << "_" << j;
   PerlInput::HashStruct Hash = Input.getHash(ParentHash,tmp.str().c_str());

   char const* const pptype = Input.getString(Hash,"Type");
   if (!strcmp("LJ",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");

      return new LJ(Eps0,Eps1,Sigma0,Sigma1);
   }
   else if (!strcmp("LJCutoff",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      Cutoff = Input.getDouble(Hash,"Cutoff");
      
      return new LJCutoff(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
   }
   else if (!strcmp("RadiiMorse",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rtheta1Pow = Input.getDouble(Hash,"Rtheta1Pow");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Rtheta2Pow = Input.getDouble(Hash,"Rtheta2Pow");
      
      return new RadiiMorse(A0,AT,B0,BT,Rref1,Rtheta1,Rtheta1Pow,Rref2,Rtheta2,Rtheta2Pow);
   }
   else if (!strcmp("RadiiMorseCutoff",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rtheta1Pow = Input.getDouble(Hash,"Rtheta1Pow");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Rtheta2Pow = Input.getDouble(Hash,"Rtheta2Pow");
      Cutoff = Input.getDouble(Hash,"Cutoff");

      return new RadiiMorseCutoff(A0,AT,B0,BT,Rref1,Rtheta1,Rtheta1Pow,Rref2,Rtheta2,
                                  Rtheta2Pow,Cutoff);
   }
   else if (!strcmp("RadiiMorseCutoff2",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rtheta1Pow = Input.getDouble(Hash,"Rtheta1Pow");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Rtheta2Pow = Input.getDouble(Hash,"Rtheta2Pow");
      Cutoff = Input.getDouble(Hash,"Cutoff");

      return new RadiiMorseCutoff2(A0,AT,B0,BT,Rref1,Rtheta1,Rtheta1Pow,Rref2,Rtheta2,
                                  Rtheta2Pow,Cutoff);
   }
   else if (!strcmp("TempMorse",pptype))
   {
      Tref = Input.getDouble(Hash,"Tref");
      A0 = Input.getDouble(Hash,"A0");
      B0 = Input.getDouble(Hash,"B0");
      Alpha = Input.getDouble(Hash,"Alpha");
      Rref1 = Input.getDouble(Hash,"Ref1");
      Tmelt = Input.getDouble(Hash,"Tmelt");
      
      return new TempMorse(A0,B0,Alpha,Rref1,Tref,Tmelt);
   }
   else if (!strcmp("Dobson",pptype))
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
   
   return 0;
}


void UpdatePairPotential(char const* const HashName,PerlInput const& Input,int const& i,
                         int const& j,PairPotentials* const Potential)
{return UpdatePairPotential(Input.getHash(HashName),Input,i,j,Potential);}

void UpdatePairPotential(PerlInput::HashStruct const& ParentHash,PerlInput const& Input,
                         int const& i,int const& j,PairPotentials* const Potential)
{
   stringstream tmp;
   double Eps0,Eps1,Sigma0,Sigma1,rcut;
   double Tref,A0,AT,B0,BT,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta2,Cutoff;

   tmp.str("");
   tmp << "PotentialType_" << i << "_" << j;
   PerlInput::HashStruct Hash = Input.getHash(ParentHash,tmp.str().c_str());
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

      if (Input.ParameterOK(Hash,"Update-AT"))
      {
         AT = Input.getDouble(Hash,"Update-AT");
         RM->SetAT(AT);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-BT"))
      {
         BT = Input.getDouble(Hash,"Update-BT");
         RM->SetBT(BT);
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

      if (Input.ParameterOK(Hash,"Update-AT"))
      {
         AT = Input.getDouble(Hash,"Update-AT");
         RM->SetAT(AT);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-BT"))
      {
         BT = Input.getDouble(Hash,"Update-BT");
         RM->SetBT(BT);
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

      if (Input.ParameterOK(Hash,"Update-AT"))
      {
         AT = Input.getDouble(Hash,"Update-AT");
         RM->SetAT(AT);
      }

      if (Input.ParameterOK(Hash,"Update-B0"))
      {
         B0 = Input.getDouble(Hash,"Update-B0");
         RM->SetB0(B0);
      }

      if (Input.ParameterOK(Hash,"Update-BT"))
      {
         BT = Input.getDouble(Hash,"Update-BT");
         RM->SetBT(BT);
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
