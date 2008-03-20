#include "ChainSum.h"

int CHAINSUMcomp(const void *a,const void *b);
int CHAINSUMind(double i,double j);
#include <cstdlib>

using namespace std;

ChainSum::ChainSum(Vector *DOF,int LagrangeCB,int Translations,Matrix *RefLat,
                   int InternalAtoms,Vector *InternalPOS,PairPotentials ***PairPot,
                   double *InfluDist,double *Ntemp)
   : Recalc_(0),
     InfluanceDist_(InfluDist),
     DOF_(DOF),
     LagrangeCB_(LagrangeCB),
     RefLattice_(RefLat),
     InternalAtoms_(InternalAtoms),
     InternalPOS_(InternalPOS),
     Potential_(PairPot),
     Ntemp_(Ntemp),
     CurrentPOS_(0),Pairs_(0),
     Translations_(Translations),
     V_(InternalAtoms),
     RelPosDATA_(int(2*(*InfluDist)*InternalAtoms*InternalAtoms),CHAINSUMdatalen)
{
   Initialize();
}

void ChainSum::operator()(Vector *DOF,int LagrangeCB,int Translations,Matrix *RefLat,
                          int InternalAtoms,Vector *InternalPOS,PairPotentials ***PairPot,
                          double *InfluDist,double *Ntemp)
{
   DOF_ = DOF;
   LagrangeCB_ = LagrangeCB;
   Translations_ = Translations;
   RefLattice_ = RefLat;
   InternalAtoms_= InternalAtoms;
   InternalPOS_ = InternalPOS;
   InfluanceDist_ = InfluDist;
   V_.Resize(InternalAtoms);
   Recalc_ = 0;
   CurrentPOS_ = 0;
   Pairs_ = 0;
   Potential_ = PairPot;
   Ntemp_ = Ntemp;
   RelPosDATA_.Resize(int(2*(*InfluDist)*InternalAtoms*InternalAtoms),CHAINSUMdatalen);
   
   Initialize();
}

void ChainSum::Reset()
{
   if (Recalc_)
   {
      Initialize();
   }
   else
   {
      CurrentPOS_ = 0;
   }
}

void ChainSum::Initialize()
{
   static double X;
   static double Influancedist,tmp;
   static int p,q;
   static int Top,Bottom,CurrentInfluanceDist;
   
   F_ = (*DOF_)[0];
   
   if (Translations_)
   {
      for (q=0;q<InternalAtoms_;++q)
      {
         V_[q] = (*DOF_)[q+1];
      }
   }
   else
   {
      V_[0] = 0.0;
      for (q=1;q<InternalAtoms_;++q)
      {
         V_[q] = (*DOF_)[q];
      }
   }
   
   // Set to inverse eigenvalue
   tmp = 1.0/F_;
   Influancedist=tmp*(*InfluanceDist_);
   
   tmp = 1;
   // set influance distance based on cube size
   CurrentInfluanceDist = int(ceil(Influancedist));
   tmp *= 2.0*CurrentInfluanceDist;
   
   Top = CurrentInfluanceDist;
   Bottom = -CurrentInfluanceDist;
   
   // set tmp to the number of pairs in the cell to be scanned
   tmp *= InternalAtoms_*InternalAtoms_;
   // Vol of sphere of R=0.5 is 0.52
   // make sure there is enough memory to store a sphere (which fits inside
   // the cube) of pairs.
   if (RelPosDATA_.Rows() < tmp)
   {
      RelPosDATA_.Resize(int(tmp),CHAINSUMdatalen);
      cerr << "Resizing RELPOSDATA matrix in ChainSum object to " << tmp << endl;
   }
   
   Pairs_ = 0;
   for (p=0;p<InternalAtoms_;p++)
   {
      for (q=0;q<InternalAtoms_;q++)
      {
         for (X = Bottom;X <= Top;X++)
         {
            RelPosDATA_[Pairs_][CHAINSUMatomstart] = p;
            RelPosDATA_[Pairs_][CHAINSUMatomstart+1] = q;
            
            // reference position
            if (LagrangeCB_)
            {
               RelPosDATA_[Pairs_][CHAINSUMdXstart] =
                  (X + ((InternalPOS_[q][0] + V_[q])
                        - (InternalPOS_[p][0] + V_[p])))
                  *(*RefLattice_)[0][0];
            }
            else
            {
               RelPosDATA_[Pairs_][CHAINSUMdXstart] =
                  (X + InternalPOS_[q][0] - InternalPOS_[p][0])*(*RefLattice_)[0][0];
            }
            
            if (LagrangeCB_)
            {
               RelPosDATA_[Pairs_][CHAINSUMdxstart] = F_ * RelPosDATA_[Pairs_][CHAINSUMdXstart];
            }
            else
            {
               RelPosDATA_[Pairs_][CHAINSUMdxstart] = F_ * RelPosDATA_[Pairs_][CHAINSUMdXstart]
                  + V_[q] - V_[p];
            }
            RelPosDATA_[Pairs_][CHAINSUMr2start] =
               RelPosDATA_[Pairs_][CHAINSUMdxstart]*RelPosDATA_[Pairs_][CHAINSUMdxstart];
            
            // Only use Sphere of Influance (current)
            if ((RelPosDATA_[Pairs_][CHAINSUMr2start] != 0)
                && (RelPosDATA_[Pairs_][CHAINSUMr2start]
                    <= (*InfluanceDist_)*(*InfluanceDist_)))
            {
               // calculate phi1 and phi2
               RelPosDATA_[Pairs_][CHAINSUMphi1start] = Potential_[p][q]->PairPotential(
                  *Ntemp_,RelPosDATA_[Pairs_][CHAINSUMr2start],PairPotentials::DY);
               RelPosDATA_[Pairs_][CHAINSUMphi2start] = Potential_[p][q]->PairPotential(
                  *Ntemp_,RelPosDATA_[Pairs_][CHAINSUMr2start],PairPotentials::D2Y);
               
               ++Pairs_;
            }
         }
      }
   }
   
   Recalc_ = 0;
   CurrentPOS_ = 0;
}

Matrix ChainSum::NeighborDistances(int cutoff,double eps)
{
   Reset();
   Matrix NeighborInfo(Pairs_,3);
   
   for (unsigned i=0;i<Pairs_;++i)
   {
      NeighborInfo[i][0] = RelPosDATA_[i][CHAINSUMr2start];
      NeighborInfo[i][1] = RelPosDATA_[i][CHAINSUMatomstart];
      NeighborInfo[i][2] = RelPosDATA_[i][CHAINSUMatomstart+1];
   }
   qsort(&(NeighborInfo[0][0]),Pairs_,3*sizeof(double),&CHAINSUMcomp);
   
   Matrix NeighborDist(cutoff,(InternalAtoms_*(InternalAtoms_+1))/2 + 1,0.0);
   int i=0,j=0;
   double CurrentDist;
   while (i<cutoff)
   {
      CurrentDist = NeighborInfo[j][0];
      NeighborDist[i][0] = sqrt(CurrentDist);
      while (fabs(CurrentDist - NeighborInfo[j][0]) < eps)
      {
         ++NeighborDist[i][CHAINSUMind(NeighborInfo[j][1],NeighborInfo[j][2])];
         ++j;
      }
      ++i;
   }
   
   return NeighborDist;
}


int CHAINSUMcomp(const void *a,const void *b)
{
   double t;
   if( *((double*) a) == *((double*) b) ) return 0;
   else
   {
      t= *((double*) a) - *((double*) b);
      t/=fabs(t);
      return int(t);
   }
}

int CHAINSUMind(double i,double j)
{
   int I=int(i+1),J=int(j+1);
   
   // Make sure I<=J
   if (I>J)
   {
      int s;
      s=I;
      I=J;
      J=s;
   }
   
   return ((J-1)*J)/2 + I;
}
