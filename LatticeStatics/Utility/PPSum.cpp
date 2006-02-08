#include "PPSum.h"

int PPSUMcomp(const void *a,const void *b);
int PPSUMind(double i,double j);
#include <cstdlib>

using namespace std;

PPSum::PPSum(Vector *DOF,Matrix *RefLat,int InternalAtoms,Vector *InternalPOS,
	     PairPotentials ***PairPot,unsigned *InfluDist,double *Ntemp)
   : DOF_(DOF),RefLattice_(RefLat),InternalAtoms_(InternalAtoms),Ntemp_(Ntemp),
     InternalPOS_(InternalPOS),Potential_(PairPot),InfluanceDist_(InfluDist),
     U_(3,3),V_(InternalAtoms,3),Recalc_(0),CurrentPOS_(0),Pairs_(0),
     RelPosDATA_(int(pow(double(2*(*InfluDist)),double(3))*pow(double(InternalAtoms),
							       double(2))),PPSUMdatalen)
{
   Initialize();
}

void PPSum::operator()(Vector *DOF,Matrix *RefLat,int InternalAtoms,
		       Vector *InternalPOS,PairPotentials ***PairPot,
		       unsigned *InfluDist,double *Ntemp)
{
   DOF_ = DOF;
   RefLattice_ = RefLat;
   InternalAtoms_= InternalAtoms;
   InternalPOS_ = InternalPOS;
   InfluanceDist_ = InfluDist;
   U_.Resize(3,3);
   V_.Resize(InternalAtoms,3);
   Recalc_ = 0;
   CurrentPOS_ = 0;
   Pairs_ = 0;
   Potential_ = PairPot;
   Ntemp_ = Ntemp;
   RelPosDATA_.Resize(
      int(pow(double(2*(*InfluDist)),double(3))*pow(double(InternalAtoms),double(2))),
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
   static Matrix Eigvals(1,3);
   static double X[3];
   static double Influancedist[3],tmp;
   static int p,q,i,j;
   static int Top[3],Bottom[3],CurrentInfluanceDist;

   U_[0][0] = (*DOF_)[0];
   U_[1][1] = (*DOF_)[1];
   U_[2][2] = (*DOF_)[2];
   U_[0][1] = U_[1][0] = (*DOF_)[3];
   U_[0][2] = U_[2][0] = (*DOF_)[4];
   U_[1][2] = U_[2][1] = (*DOF_)[5];

   V_[0][0] = 0.0;
   V_[0][1] = 0.0;
   V_[0][2] = 0.0;
   i=6;
   for (q=1;q<InternalAtoms_;++q)
   {
      for (p=0;p<3;p++)
      {
	 V_[q][p] = (*DOF_)[i++];
      }
   }

   // find largest eigenvalue of the inverse transformation
   // (i.e. from current to ref) and use influence cube of
   // that size...
   //
   // Use the fact that eigs of Uinv = 1/ eigs of U.
   //
   // Use U_*RefLattice_ as def grad.  This takes an
   // orthonormal lattice to the current config.
   // Thus, allowing non-square unit cells....
   //
   // Use F*F^T and take sqrt of eigvecs.
   Eigvals = SymEigVal(U_*(*RefLattice_)*((U_*(*RefLattice_)).Transpose()));
   tmp = sqrt(Eigvals[0][0]);
   for (i=0;i<3;i++)
      if (sqrt(Eigvals[0][i]) < tmp) tmp = sqrt(Eigvals[0][i]);
   
   // Set to inverse eigenvalue
   tmp = 1.0/tmp;
   for (i=0;i<3;i++)
   {
      Influancedist[i]=tmp*(*InfluanceDist_);
   }

   tmp = 1;
   for (p=0;p<3;p++)
   {
      // set influance distance based on cube size
      //
      // also setup to be large enough to encompass Eulerian sphere
      CurrentInfluanceDist = int(ceil(Influancedist[p]));
      tmp *= 2.0*CurrentInfluanceDist;

      Top[p] = CurrentInfluanceDist;
      Bottom[p] = -CurrentInfluanceDist;
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
		     RelPosDATA_[Pairs_][PPSUMdXstart+i] = 0.0;

		     // "SHIFTED reference position"
		     for (j=0;j<3;j++)
		     {
			RelPosDATA_[Pairs_][PPSUMdXstart+i] +=
			   (X[j] + ((InternalPOS_[q][j] + V_[q][j])
				    - (InternalPOS_[p][j] + V_[p][j])))
			   *(*RefLattice_)[j][i];

		     }
		  }
		  RelPosDATA_[Pairs_][PPSUMr2start] = 0.0;
		  for (i=0;i<3;i++)
		  {
		     RelPosDATA_[Pairs_][PPSUMdxstart+i] = 0.0;

		     for (j=0;j<3;j++)
		     {
			RelPosDATA_[Pairs_][PPSUMdxstart+i]
			   += U_[i][j] * RelPosDATA_[Pairs_][PPSUMdXstart+j];
		     }
		     RelPosDATA_[Pairs_][PPSUMr2start]
			+= RelPosDATA_[Pairs_][PPSUMdxstart+i]*RelPosDATA_[Pairs_][PPSUMdxstart+i];
		  }
		  // Only use Sphere of Influance (current)
		  if ((RelPosDATA_[Pairs_][PPSUMr2start] != 0) &&
		      (RelPosDATA_[Pairs_][PPSUMr2start]
		       <= (*InfluanceDist_)*(*InfluanceDist_)))
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
