#include "CMatrix.h"
#include <iomanip>
#include <string>
#include <cmath>

// Global IDString
char CMatrixID[]="$Id: CMatrix.cpp,v 1.15 2008/03/23 21:01:27 elliott Exp $";

// Private Methods...

// Returns matrix of size Rows_-1 x Cols_-1 with ith row and
//   jth column removed
CMatrix CMatrix::Minor(unsigned i,unsigned j) const
{
   CMatrix A(Rows_-1,Cols_-1);
   
   if (!IsNull() || !A.IsNull())
   {
      for (register unsigned a=0;a<A.Rows_;a++)
      {
	 for (register unsigned b=0;b<A.Cols_;b++)
	 {
	    if (a < i)
	    {
	       if (b < j)
	       {
		  A.Elements_[a][b]=Elements_[a][b];
	       }
	       else
	       {
		  A.Elements_[a][b]=Elements_[a][b+1];
	       }
	    }
	    else
	    {
	       if (b < j)
	       {
		  A.Elements_[a][b]=Elements_[a+1][b];
	       }
	       else
	       {
		  A.Elements_[a][b]=Elements_[a+1][b+1];
	       }
	    }
	 }
      }
   }
   
   return A;
}

// Public Methods...

int CMatrix::MathematicaPrintFlag = 0;

CMatrix::CMatrix(unsigned Rows,unsigned Cols,CMatrix::Elm InitVal):
   Rows_(Rows), Cols_(Cols)
{
   if (IsNull())
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new CMatrix::Elm*[Rows_];
      Elements_[0]=new CMatrix::Elm[Rows_*Cols_];

      for (register unsigned i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      if (InitVal!=SENTINAL)
      {
	 for (register unsigned i=0;i<Rows_;i++)
	 {
	    for (register unsigned j=0;j<Cols_;j++)
	       Elements_[i][j]=InitVal;
	 }
      }
   }
   return;
}

CMatrix::CMatrix(const CMatrix& A):
   Rows_(A.Rows_), Cols_(A.Cols_)
{
   if (IsNull())
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new CMatrix::Elm*[Rows_];
      Elements_[0]=new CMatrix::Elm[Rows_*Cols_];

      for (register unsigned i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      memmove(Elements_[0],A.Elements_[0],sizeof(CMatrix::Elm)*Rows_*Cols_);
   }

   return;
}

CMatrix::CMatrix(const Matrix& A):
   Rows_(A.Rows()), Cols_(A.Cols())
{
   if (IsNull())
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new CMatrix::Elm*[Rows_];
      Elements_[0]=new CMatrix::Elm[Rows_*Cols_];

      for (register unsigned i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      for (register unsigned i=0;i<Rows_;i++)
	 for (register unsigned j=0;j<Cols_;j++)
	 {
	    Elements_[i][j] = CMatrix::Elm(A[i][j],0.0);
	 }
   }
   return;
}

CMatrix::~CMatrix()
{
   if (!IsNull())
   {
      delete [] Elements_[0];
      delete [] Elements_;
   }
   
   return;
}

CMatrix operator+(const CMatrix& A,const CMatrix& B)
{
   if (A.Rows_!=B.Rows_ || A.Cols_!=B.Cols_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix Operator+() Diff Size Matrices or Null CMatrix!!!"
	   << "\n";
      exit(-1);
   }

   CMatrix C(A.Rows_,A.Cols_);
   
   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=A.Elements_[i][j]+B.Elements_[i][j];
      }
   }
   
   return C;
}

CMatrix operator-(const CMatrix& A)
{
   CMatrix B(A.Rows_,A.Cols_);

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for(unsigned j=0;j<A.Cols_;j++)
      {
	 B.Elements_[i][j]=-A.Elements_[i][j];
      }
   }

   return B;
}

CMatrix operator-(const CMatrix& A,const CMatrix& B)
{
   if (A.Rows_!=B.Rows_ || A.Cols_!=B.Cols_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix Operator-() Diff Size Matrices or Null CMatrix!!!"
	   << "\n";
      exit(-1);
   }
      
   CMatrix C(A.Rows_,A.Cols_);
   
   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=A.Elements_[i][j]-B.Elements_[i][j];
      }
   }
   
   return C;
}

CMatrix operator*(const CMatrix& A,const CMatrix& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In CMatrix Operator* : A.Cols!=B.Rows or Null CMatrix"
	   <<"\n";
      exit(-1);
   }

   CMatrix C(A.Rows_,B.Cols_,0.0);
   
   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<B.Cols_;j++)
      {
	 for (register unsigned k=0;k<A.Cols_;k++)
	 {
	    C.Elements_[i][j]+=A.Elements_[i][k]*B.Elements_[k][j];
	 }
      }
   }
   
   return C;
}

CMatrix operator*(const CMatrix::Elm& A,const CMatrix& B)
{
   CMatrix C(B.Rows_,B.Cols_);

   for (register unsigned i=0;i<B.Rows_;i++)
   {
      for (register unsigned j=0;j<B.Cols_;j++)
      {
	 C.Elements_[i][j]=A*B.Elements_[i][j];
      }
   }

   return C;
}

CMatrix operator*(const CMatrix& A,const CMatrix::Elm& B)
{
   CMatrix C(A.Rows_,A.Cols_);

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=B*A.Elements_[i][j];
      }
   }

   return C;
}

CMatrix operator/(const CMatrix& A,const CMatrix::Elm& B)
{
   if (B == CMatrix::Elm(0.0,0.0))
   {
      cerr << "Divide By Zero Error in CMatrix operator/()"
	   << "\n";
      exit(-1);
   }
   
   CMatrix C(A.Rows_,A.Cols_);

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=A.Elements_[i][j]/B;
      }
   }

   return C;
}

#ifdef CHECK_BOUNDS
CMatrix::Elm* CMatrix::operator[](unsigned i)
{
   if (i >= Rows_)
   {
      cerr << "CMatrix Index Overflow -- CMatrix::Elm* operator[]()" << "\n";
      exit(-1);
   }

   return Elements_[i];
}

CMatrix::Elm* CMatrix::operator[](unsigned i) const
{
   if (i >= Rows_)
   {
      cerr << "CMatrix Index Overflow -- CMatrix::Elm* operator[]()" << "\n";
      exit(-1);
   }

   return Elements_[i];
}
#endif
   
CMatrix& CMatrix::operator=(const CMatrix& B)
{
   if (Rows_!=B.Rows_ || Cols_!=B.Cols_ || IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix& operator=() : Matricies not same size "
           << "or Null CMatrix"
           << "\n";
      exit(-1);
   }

   memmove(Elements_[0],B.Elements_[0],sizeof(CMatrix::Elm)*Rows_*Cols_);

   return *this;
}

CMatrix& CMatrix::SetIdentity(unsigned Size)
{
   if (Rows_!=Size || Cols_!=Size)
      Resize(Size,Size);

   for (register unsigned i=0;i<Size;i++)
   {
      for (register unsigned j=0;j<Size;j++)
      {
	 if (i==j)
	    Elements_[i][i]=1.0;
	 else
	    Elements_[i][j]=0.0;
      } 
   }

   return *this;
}

CMatrix CMatrix::Transpose() const
{
   CMatrix A(Cols_,Rows_);

   for (register unsigned i=0;i<Rows_;i++)
   {
      for (register unsigned j=0;j<Cols_;j++)
      {
	 A.Elements_[j][i]=Elements_[i][j];
      }
   }

   return A;
}

CMatrix CMatrix::Conjugate() const
{
   CMatrix A(Cols_,Rows_);

   for (register unsigned i=0;i<Rows_;++i)
   {
      for (register unsigned j=0;j<Cols_;++j)
      {
	 A.Elements_[i][j] = Elements_[i][j].conj();
      }
   }

   return A;
}

CMatrix CMatrix::Inverse() const
{
   if (!IsSquare() || IsNull())
   {
      cerr << "Error in CMatrix::Inverse() : Non-Square or Null CMatrix" << "\n";
      exit(-1);
   }

   CMatrix B(Rows_,1,0),X(Rows_,1),C(Rows_,Cols_);

   B.Elements_[0][0]=1.0;
   for (register unsigned i=0;i<Cols_;i++)
   {
      X=SolvePLU(*this,B);

      for (register unsigned j=0;j<Rows_;j++)
	 C.Elements_[j][i]=X.Elements_[j][0];

      B.Elements_[i][0]=0;
      if (i!=Cols_-1) B.Elements_[i+1][0]=1.0;
   }

   return C;
}

void CMatrix::Resize(unsigned Rows,unsigned Cols,CMatrix::Elm InitVal)
{
   if (Rows!=Rows_ || Cols!=Cols_)
   {
      if (Elements_!=NULL)
      {
	 delete [] Elements_[0];
	 delete [] Elements_;
      }

      Rows_=Rows;
      Cols_=Cols;
   
      if (IsNull())
      {
	 Elements_=NULL;
      }
      else
      {
	 Elements_=new CMatrix::Elm*[Rows_];
	 Elements_[0]=new CMatrix::Elm[Rows_*Cols_];

	 for (register unsigned i=1;i<Rows_;i++)
	 {
	    Elements_[i]=Elements_[i-1]+Cols_;
	 }
      }
   }

   if (InitVal != SENTINAL)
      for (register unsigned i=0;i<Rows_;i++)
      {
	 for (register unsigned j=0;j<Cols_;j++)
	 {
	    Elements_[i][j]=InitVal;
	 }
      }

   return;
}


// Recursivly calculate determinent
CMatrix::Elm CMatrix::Det() const
{
   if (IsNull() || !IsSquare())
   {
      cerr << "Error in CMatrix::Det() : Null or Non-Square CMatrix" << "\n";
      exit(-1);
   }

   if (Rows_==1)
      return Elements_[0][0];
   else
   {
      Elm det(0.0,0.0);

      for (register unsigned i=0;i<Cols_;i++)
      {
	 det += CMatrix::Elm(1-2*(i%2))*Elements_[0][i]*(Minor(0,i).Det());
      }

      return det;
   }
}

// Decompose PA=LU using scaled partial pivoting.
void PLU(const CMatrix& A,CMatrix& P,CMatrix& L,CMatrix& U)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in PLU -- Non-Square or Null CMatrix to decompose..." << "\n";
      exit(-1);
   }

   P.Resize(A.Rows_,A.Cols_,0);
   L.SetIdentity(A.Rows_);
   U.Resize(A.Rows_,A.Cols_,0);

   CMatrix Temp=A,
      S(A.Rows_,1,0);
   int *Ipivot;
   Ipivot = new int[A.Rows_];

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      Ipivot[i]=i;
   }

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<A.Cols_;j++)
      {
	 if (abs(Temp.Elements_[i][j]) > abs(S.Elements_[i][0]))
	    S.Elements_[i][0]=Temp.Elements_[i][j];
      }
   }

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      CMatrix::Elm temp1;
      temp1=Temp.Elements_[i][i]/S.Elements_[i][0];

      unsigned k=i;
      for (register unsigned j=i;j<A.Rows_;j++)
      {
	 if (abs(Temp.Elements_[j][i]) > abs(temp1))
	 {
	    temp1=Temp.Elements_[j][i]/S.Elements_[j][0];
	    k=j;
	 }
      }

      if (k>i)
      {
	 CMatrix::Elm *Switch;
	 Switch = new CMatrix::Elm[A.Rows_];
	 for (register unsigned j=i;j<A.Rows_;j++)
	 {
	    Switch[j]=Temp.Elements_[i][j];
	    Temp.Elements_[i][j]=Temp.Elements_[k][j];
	    Temp.Elements_[k][j]=Switch[j];
	 }

	 for (register unsigned j=0;j<i;j++)
	 {
	    Switch[j]=L.Elements_[i][j];
	    L.Elements_[i][j]=L.Elements_[k][j];
	    L.Elements_[k][j]=Switch[j];
	 }

	 delete [] Switch;

	 temp1=S.Elements_[i][0];
	 S.Elements_[i][0]=S.Elements_[k][0];
	 S.Elements_[k][0]=temp1;

	 int tempi1=Ipivot[i];
	 Ipivot[i]=Ipivot[k];
	 Ipivot[k]=tempi1;
      }

      for (register unsigned j=i+1;j<A.Rows_;j++)
      {
	 L.Elements_[j][i]=Temp.Elements_[j][i]/Temp.Elements_[i][i];
      }

      for (register unsigned j=i+1;j<A.Rows_;j++)
      {
	 for (register unsigned k=i+1;k<A.Rows_;k++)
	 {
	    Temp.Elements_[j][k]=Temp.Elements_[j][k]-(L.Elements_[j][i]*Temp.Elements_[i][k]);
	 }
      }
   }

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=i;j<A.Rows_;j++)
      {
	 U.Elements_[i][j]=Temp.Elements_[i][j];
      }
   }

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      P.Elements_[i][Ipivot[i]]=1;
   }

   delete [] Ipivot;

   return;
}

Matrix HermiteEigVal(CMatrix A,CMatrix *B,const unsigned MaxItr,const double Tol)
{
   unsigned count=0,
      converged=0;
   Matrix EigVals(1,A.Cols_);
   MyComplexDouble tau,t1,t2,t,c,s,cc,ss,cs,ssbar,csbar,
      aij1,aii1,ajj1,aki1,akj1,tmp;

   if (B != NULL)
   {
      B->SetIdentity(A.Cols_);
   }
   
   while ((count < MaxItr) && (!converged))
   {
      for (unsigned i=0;i<A.Cols_;i++)
      {
	 for (unsigned j=i+1;j<A.Cols_;j++)
	 {
	    if (abs(A.Elements_[i][j]) < Tol )
	       continue;
	    else
	    {
	       tau = (A.Elements_[i][i] - A.Elements_[j][j])
		  /(2.0*(A.Elements_[i][j].conj()));

	       // take care to make the sqrt well conditioned! (not overflow)
	       if (abs(tau) > 1.0)
	       {
		  t1 = tau*(-1.0 - sqrt( 1.0 + (A.Elements_[i][j]
						/((A.Elements_[i][j].conj())*tau*tau))));
		  t2 = tau*(-1.0 + sqrt( 1.0 + (A.Elements_[i][j]
						/((A.Elements_[i][j].conj())*tau*tau))));
	       }
	       else
	       {
		  t1 = -tau - sqrt( A.Elements_[i][j]/(A.Elements_[i][j].conj())
				    + tau*tau);
		  t2 = -tau + sqrt( A.Elements_[i][j]/(A.Elements_[i][j].conj())
				    + tau*tau);
	       }

	       if (abs(t1) >= abs(t2))
		  t = t2;
	       else
		  t = t1;
	       
	       c = 1.0/sqrt(1.0 + (t*(t.conj())).real());
	       s = t * c;
	       cc = c*c;
	       ss = s*s;
	       ssbar =(s*(s.conj())).real();
	       cs = c*s;
	       csbar = c*(s.conj());
	       aij1 = A.Elements_[i][j];
	       aii1 = A.Elements_[i][i];
	       ajj1 = A.Elements_[j][j];
	       
	       A.Elements_[i][i] = (aii1*cc + cs*(aij1.conj())
				    + csbar*aij1 + ssbar*ajj1).real();
	       A.Elements_[j][j] = (aii1*ssbar - cs*(aij1.conj())
				    - csbar*aij1 + cc*ajj1).real();
	       A.Elements_[i][j] = A.Elements_[j][i] = 0.0;
	       
	       for (unsigned k=0;k<A.Cols_;k++)
	       {
		  if (B != NULL)
		  {
		     tmp = B->Elements_[k][i]*c + B->Elements_[k][j]*(s.conj());
		     B->Elements_[k][j] = -B->Elements_[k][i]*s + B->Elements_[k][j]*c;
		     B->Elements_[k][i] = tmp;
		  }
		  
		  
		  if ( k==i || k==j)
		  {
		     continue;
		  }
		     
		  aki1 = A.Elements_[k][i];
		  akj1 = A.Elements_[k][j];
		  
		  A.Elements_[k][i] = aki1*c + akj1*(s.conj());
		  A.Elements_[i][k] = (A.Elements_[k][i].conj());
		  A.Elements_[k][j] = -aki1*s + akj1*c;
		  A.Elements_[j][k] = (A.Elements_[k][j].conj());
	       }
	    }
	 }
      }
      count++;
      
      converged = 1;
      for (unsigned i=0;i<A.Cols_;i++)
	 for (unsigned j=i+1;j<A.Cols_;j++)
	 {
	    if (abs(A.Elements_[i][j]) > Tol)
	    {
	       converged = 0;
	    }
	 }
   }

   if (!converged)
   {
      cerr << "Error: HermiteEigVal(): Failed - No convergence!" << "\n";
      exit(-1);
   }

   for (unsigned i=0;i<A.Cols_;i++)
   {
      EigVals[0][i] = (A.Elements_[i][i]).real();
   }
   
   return EigVals;
}

void Cholesky(const CMatrix& A,CMatrix& U,CMatrix& D)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in Cholesky() -- Non-Square or Null CMatrix" << "\n";
      exit(-1);
   }

   U.SetIdentity(A.Rows_);
   D.Resize(A.Rows_,A.Cols_,0);
   

   for (register unsigned i=0;i<A.Rows_;i++)
   {
      D.Elements_[i][i]=A.Elements_[i][i];
      if (i > 0)
      {
	 for (register unsigned k=0;k<i;k++)
	 {

	    D.Elements_[i][i]-=
	       D.Elements_[k][k]*U.Elements_[k][i]*(U.Elements_[k][i].conj());
	 }
      }
      for (register unsigned j=i+1;j<A.Rows_;j++)
      {
	 U.Elements_[i][j]=A.Elements_[i][j];
	 if (i > 0)
	 {
	    for (register unsigned k=0;k<i;k++)
	    {
	       U.Elements_[i][j]-=
		  D.Elements_[k][k]*(U.Elements_[k][i].conj())*U.Elements_[k][j];
	    }
	 }
	 U.Elements_[i][j]/=D.Elements_[i][i];
      }
   }
   
   return;
}

CMatrix SolvePLU(const CMatrix& A,const CMatrix& B)
{
   if (!A.IsSquare() || A.IsNull() || A.Cols_!=B.Rows_)
   {
      cerr << "Error in Solve() - Non-Square CMatrix, Null CMatrix, or system of != dimension"
	   << "\n";
      exit(-1);
   }
   
   CMatrix
      P,L,U; //PLU() resizes P,L,& U thus do not waste time initializing them.

   PLU(A,P,L,U);

   CMatrix::Elm *Y;
   Y = new CMatrix::Elm[B.Rows_];
   CMatrix Temp=P*B;

   Y[0]=Temp.Elements_[0][0];
   for (register unsigned i=1;i<B.Rows_;i++)
   {
      Y[i]=Temp[i][0];
      for (register unsigned j=0;j<i;j++)
      {
	 Y[i]-=Y[j]*L.Elements_[i][j];
      }
   }

   CMatrix X(B.Rows_,1,0);

   X.Elements_[B.Rows_-1][0]=Y[B.Rows_-1]/U.Elements_[A.Rows_-1][A.Cols_-1];
   for (register int i=A.Rows_-2;i>=0;i--)
   {
      X.Elements_[i][0]=Y[i];
      for (register int j=A.Rows_-1;j>i;j--)
      {
	 X.Elements_[i][0]-=X.Elements_[j][0]*U.Elements_[i][j];
      }
      X.Elements_[i][0]=X.Elements_[i][0]/U.Elements_[i][i];
   }

   delete [] Y;
      
   return X;
}

ostream& operator<<(ostream& out,const CMatrix& A)
{
   int W=out.width();
   
   out << "\n";

   if (CMatrix::MathematicaPrintFlag) out << setw(0) << "{{";
   for (register unsigned i=0;i<A.Rows_;i++)
   {
      for (register unsigned j=0;j<A.Cols_;j++)
      {
	 out << setw(W) << A.Elements_[i][j];
	 if ((CMatrix::MathematicaPrintFlag) && (j!=(A.Cols_-1)))
	    out << ",";
      }
      
      if (CMatrix::MathematicaPrintFlag)
      {
	 if (i!=(A.Rows_-1))
	    out << "},\n {";
	 else
	    out << "}}";
      }
      else
	 out << "\n";
   }

   out << "\n";

   return out;
}

istream& operator>>(istream& in,CMatrix& A)
{
   for (register unsigned i=0;i<A.Rows_;i++)
      for (register unsigned j=0;j<A.Cols_;j++)
	 in >> A.Elements_[i][j];

   return in;
}

char* CMatrix::Revision()
{
   return CMatrixID;
}
