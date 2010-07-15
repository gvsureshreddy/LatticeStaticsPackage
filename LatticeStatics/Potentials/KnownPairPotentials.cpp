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
   double Tref,A0,AT,ATPow,B0,BT,BTPow,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta1Pow,Rtheta2,Rtheta2Pow,Cutoff,CutoffStart,CutoffEnd;

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
   else if (!strcmp("LJConstCutoff",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      Cutoff = Input.getDouble(Hash,"Cutoff");
      
      return new LJConstCutoff(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
   }
   else if (!strcmp("LJLinearCutoff",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      Cutoff = Input.getDouble(Hash,"Cutoff");
      
      return new LJLinearCutoff(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
   }
   else if (!strcmp("LJQuadraticCutoff",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      Cutoff = Input.getDouble(Hash,"Cutoff");
      
      return new LJQuadraticCutoff(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
   }
   else if (!strcmp("LJDobson",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      Cutoff = Input.getDouble(Hash,"Cutoff");
      
      return new LJDobson(Eps0,Eps1,Sigma0,Sigma1,Cutoff);
   }
   else if (!strcmp("LJSplineCutoff",pptype))
   {
      Eps0 = Input.getDouble(Hash,"Eps0");
      Eps1 = Input.getDouble(Hash,"Eps1");
      Sigma0 = Input.getDouble(Hash,"Sigma0");
      Sigma1 = Input.getDouble(Hash,"Sigma1");
      CutoffStart = Input.getDouble(Hash,"CutoffStart");
      CutoffEnd = Input.getDouble(Hash,"CutoffEnd");
      
      return new LJSplineCutoff(Eps0,Eps1,Sigma0,Sigma1,CutoffStart,CutoffEnd);
   }
   else if (!strcmp("RadiiMorse",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      
      return new RadiiMorse(A0,AT,B0,BT,Rref1,Rtheta1,Rref2,Rtheta2);
   }
   else if (!strcmp("RadiiMorseCutoff",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Cutoff = Input.getDouble(Hash,"Cutoff");

      return new RadiiMorseCutoff(A0,AT,B0,BT,Rref1,Rtheta1,Rref2,Rtheta2,Cutoff);
   }
   else if (!strcmp("RadiiMorseCutoff2",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Cutoff = Input.getDouble(Hash,"Cutoff");

      return new RadiiMorseCutoff2(A0,AT,B0,BT,Rref1,Rtheta1,Rref2,Rtheta2,Cutoff);
   }
   else if (!strcmp("RadiiMorseOriginal",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      
      return new RadiiMorseOriginal(A0,AT,B0,BT,Rref1,Rtheta1,Rref2,Rtheta2);
   }
   else if (!strcmp("GVMorse",pptype))
   {
      A0 = Input.getDouble(Hash,"A0");
      AT = Input.getDouble(Hash,"AT");
      ATPow = Input.getDouble(Hash,"ATPow");
      B0 = Input.getDouble(Hash,"B0");
      BT = Input.getDouble(Hash,"BT");
      BTPow = Input.getDouble(Hash,"BTPow");
      Rref1 = Input.getDouble(Hash,"Rref1");
      Rtheta1 = Input.getDouble(Hash,"Rtheta1");
      Rtheta1Pow = Input.getDouble(Hash,"Rtheta1Pow");
      Rref2 = Input.getDouble(Hash,"Rref2");
      Rtheta2 = Input.getDouble(Hash,"Rtheta2");
      Rtheta2Pow = Input.getDouble(Hash,"Rtheta2Pow");
      
      return new GVMorse(A0,AT,ATPow,B0,BT,BTPow,Rref1,Rtheta1,Rtheta1Pow,Rref2,Rtheta2,
                         Rtheta2Pow);
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
   double Tref,A0,AT,ATPow,B0,BT,BTPow,Alpha,Rref1,Rref2,Tmelt;
   double Rtheta1,Rtheta1Pow,Rtheta2,Rtheta2Pow,Cutoff,CutoffStart,CutoffEnd;

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
   else if (!strcmp(Potential->Type(),"LJConstCutoff"))
   {
      LJConstCutoff *LJConstCutoffp;
      LJConstCutoffp = dynamic_cast<LJConstCutoff *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJConstCutoffp->SetEps0(Eps0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJConstCutoffp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJConstCutoffp->SetSigma0(Sigma0);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJConstCutoffp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         LJConstCutoffp->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"LJLinearCutoff"))
   {
      LJLinearCutoff *LJLinearCutoffp;
      LJLinearCutoffp = dynamic_cast<LJLinearCutoff *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJLinearCutoffp->SetEps0(Eps0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJLinearCutoffp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJLinearCutoffp->SetSigma0(Sigma0);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJLinearCutoffp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         LJLinearCutoffp->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"LJDobson"))
   {
      LJDobson *LJDobsonp;
      LJDobsonp = dynamic_cast<LJDobson *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJDobsonp->SetEps0(Eps0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJDobsonp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJDobsonp->SetSigma0(Sigma0);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJDobsonp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         LJDobsonp->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"LJQuadraticCutoff"))
   {
      LJQuadraticCutoff *LJQuadraticCutoffp;
      LJQuadraticCutoffp = dynamic_cast<LJQuadraticCutoff *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJQuadraticCutoffp->SetEps0(Eps0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJQuadraticCutoffp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJQuadraticCutoffp->SetSigma0(Sigma0);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJQuadraticCutoffp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-Cutoff"))
      {
         Cutoff = Input.getDouble(Hash,"Update-Cutoff");
         LJQuadraticCutoffp->SetCutoff(Cutoff);
      }
   }
   else if (!strcmp(Potential->Type(),"LJSplineCutoff"))
   {
      LJSplineCutoff *LJSplineCutoffp;
      LJSplineCutoffp = dynamic_cast<LJSplineCutoff *>(Potential);

      if (Input.ParameterOK(Hash,"Update-Eps0"))
      {
         Eps0 = Input.getDouble(Hash,"Update-Eps0");
         LJSplineCutoffp->SetEps0(Eps0);
      }
      
      if (Input.ParameterOK(Hash,"Update-Eps1"))
      {
         Eps1 = Input.getDouble(Hash,"Update-Eps1");
         LJSplineCutoffp->SetEps1(Eps1);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma0"))
      {
         Sigma0 = Input.getDouble(Hash,"Update-Sigma0");
         LJSplineCutoffp->SetSigma0(Sigma0);
      }

      if (Input.ParameterOK(Hash,"Update-Sigma1"))
      {
         Sigma1 = Input.getDouble(Hash,"Update-Sigma1");
         LJSplineCutoffp->SetSigma1(Sigma1);
      }

      if (Input.ParameterOK(Hash,"Update-CutoffStart"))
      {
         CutoffStart = Input.getDouble(Hash,"Update-CutoffStart");
         LJSplineCutoffp->SetCutoffStart(CutoffStart);
      }

      if (Input.ParameterOK(Hash,"Update-CutoffEnd"))
      {
         CutoffEnd = Input.getDouble(Hash,"Update-CutoffEnd");
         LJSplineCutoffp->SetCutoffEnd(CutoffEnd);
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
   else if (!strcmp(Potential->Type(),"GVMorse"))
   {
      GVMorse *RM;
      RM = dynamic_cast<GVMorse *>(Potential);
      
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

      if (Input.ParameterOK(Hash,"Update-ATPow"))
      {
         ATPow = Input.getDouble(Hash,"Update-ATPow");
         RM->SetATPow(ATPow);
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

      if (Input.ParameterOK(Hash,"Update-BTPow"))
      {
         BTPow = Input.getDouble(Hash,"Update-BTPow");
         RM->SetBTPow(BTPow);
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

      if (Input.ParameterOK(Hash,"Update-Rtheta1Pow"))
      {
         Rtheta1Pow = Input.getDouble(Hash,"Update-Rtheta1Pow");
         RM->SetRtheta1Pow(Rtheta1Pow);
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

      if (Input.ParameterOK(Hash,"Update-Rtheta2Pow"))
      {
         Rtheta2Pow = Input.getDouble(Hash,"Update-Rtheta2Pow");
         RM->SetRtheta2Pow(Rtheta2Pow);
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
