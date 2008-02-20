#include "PPSum.h"

int PPSUMcomp(const void *a,const void *b);
int PPSUMind(double i,double j);
#include <cstdlib>

using namespace std;

PPSum::PPSum(CBKinematics *CBK,int InternalAtoms,PairPotentials ***PairPot,double *InfluDist,
             double *Ntemp)
   : CBK_(CBK),InternalAtoms_(InternalAtoms),Ntemp_(Ntemp),Potential_(PairPot),
     InfluenceDist_(InfluDist),Recalc_(0),CurrentPOS_(0),Pairs_(0),
     RelPosDATA_(int(pow(2*(*InfluDist),double(3))*pow(double(InternalAtoms),
                                                       double(2))),PPSUMdatalen)
{
   Initialize();
}

void PPSum::operator()(CBKinematics *CBK,int InternalAtoms,PairPotentials ***PairPot,
                       double *InfluDist,double *Ntemp)
{
   CBK_=CBK;
   InternalAtoms_= InternalAtoms;
   InfluenceDist_ = InfluDist;
   Recalc_ = 0;
   CurrentPOS_ = 0;
   Pairs_ = 0;
   Potential_ = PairPot;
   Ntemp_ = Ntemp;
   RelPosDATA_.Resize(
      int(pow(2*(*InfluDist),double(3))*pow(double(InternalAtoms),double(2))),
      PPSUMdatalen);
   
   Initialize();
}

void PPSum::Reset()
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

void PPSum::Initialize()
{
   static double X[3];
   static double Influencedist[3],tmp;
   static int p,q,i;
   static int Top[3],Bottom[3],CurrentInfluenceDist;
   
   CBK_->InfluenceRegion(Influencedist);
   for (i=0;i<3;i++)
   {
      Influencedist[i] *= (*InfluenceDist_);
   }
   
   tmp = 1;
   for (p=0;p<3;p++)
   {
      // set influance distance based on cube size
      //
      // also setup to be large enough to encompass Eulerian sphere
      CurrentInfluenceDist = int(ceil(Influencedist[p]));
      tmp *= 2.0*CurrentInfluenceDist;
      
      Top[p] = CurrentInfluenceDist;
      Bottom[p] = -CurrentInfluenceDist;
   }
   
   // set tmp to the number of pairs in the cube to be scanned
   tmp *= InternalAtoms_*InternalAtoms_;
   // Vol of sphere of R=0.5 is 0.52
   // make sure there is enough memory to store a sphere (which fits inside
   // the cube) of pairs.
   if (RelPosDATA_.Rows() < 0.55*tmp)
   {
      RelPosDATA_.Resize(int(tmp),PPSUMdatalen);
      cerr << "Resizing RELPOSDATA matrix in PPSum object to " << tmp << endl;
   }
   
   Pairs_ = 0;
   for (p=0;p<InternalAtoms_;p++)
   {
      for (q=0;q<InternalAtoms_;q++)
      {
         for (X[0] = Bottom[0];X[0] <= Top[0];X[0]++)
         {
            for (X[1] = Bottom[1];X[1] <= Top[1];X[1]++)
            {
               for (X[2] = Bottom[2];X[2] <= Top[2];X[2]++)
               {
                  RelPosDATA_[Pairs_][PPSUMatomstart] = p;
                  RelPosDATA_[Pairs_][PPSUMatomstart+1] = q;
                  
                  for (i=0;i<3;i++)
                  {
                     RelPosDATA_[Pairs_][PPSUMdXstart+i] = CBK_->DX(X,p,q,i);
                  }
                  RelPosDATA_[Pairs_][PPSUMr2start] = 0.0;
                  for (i=0;i<3;i++)
                  {
                     RelPosDATA_[Pairs_][PPSUMdxstart+i] = CBK_->Dx(X,p,q,i);
                     RelPosDATA_[Pairs_][PPSUMr2start]
                        += RelPosDATA_[Pairs_][PPSUMdxstart+i]*RelPosDATA_[Pairs_][PPSUMdxstart+i];
                  }
                  // Only use Sphere of Influence (current)
                  if ((RelPosDATA_[Pairs_][PPSUMr2start] != 0) &&
                      (RelPosDATA_[Pairs_][PPSUMr2start]
                       <= (*InfluenceDist_)*(*InfluenceDist_)))
                  {
                     // calculate phi1 and phi2
                     RelPosDATA_[Pairs_][PPSUMphi1start] = Potential_[p][q]->PairPotential(
                        *Ntemp_,RelPosDATA_[Pairs_][PPSUMr2start],PairPotentials::DY);
                     RelPosDATA_[Pairs_][PPSUMphi2start] = Potential_[p][q]->PairPotential(
                        *Ntemp_,RelPosDATA_[Pairs_][PPSUMr2start],PairPotentials::D2Y);
                     
                     ++Pairs_;
                  }
               }
            }
         }
      }
   }
   
   Recalc_ = 0;
   CurrentPOS_ = 0;
}

Matrix PPSum::NeighborDistances(int cutoff,double eps)
{
   Reset();
   Matrix NeighborInfo(Pairs_,3);
   
   for (int i=0;i<Pairs_;++i)
   {
      NeighborInfo[i][0] = RelPosDATA_[i][PPSUMr2start];
      NeighborInfo[i][1] = RelPosDATA_[i][PPSUMatomstart];
      NeighborInfo[i][2] = RelPosDATA_[i][PPSUMatomstart+1];
   }
   qsort(&(NeighborInfo[0][0]),Pairs_,3*sizeof(double),&PPSUMcomp);
   
   Matrix NeighborDist(cutoff,(InternalAtoms_*(InternalAtoms_+1))/2 + 1,0.0);
   int i=0,j=0;
   double CurrentDist;
   while (i<cutoff)
   {
      CurrentDist = NeighborInfo[j][0];
      NeighborDist[i][0] = sqrt(CurrentDist);
      while (fabs(CurrentDist - NeighborInfo[j][0]) < eps)
      {
         ++NeighborDist[i][PPSUMind(NeighborInfo[j][1],NeighborInfo[j][2])];
         ++j;
      }
      ++i;
   }
   
   return NeighborDist;
}


int PPSUMcomp(const void *a,const void *b)
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

int PPSUMind(double i,double j)
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
