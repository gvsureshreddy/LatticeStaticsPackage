#include "CMatrix.h"
#include "Matrix.h"
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <math.h>

// Global IDString
char CMatrixID[]="$Header: /h/elliott/.cvsroot/ttttttttttt/LinearAlgebra/CMatrix.cpp,v 1.1 2002/10/15 20:30:20 elliottr Exp $";

// Private Methods...


// Returns matrix of size Rows_-1 x Cols_-1 with ith row and
//   jth column removed
CMatrix CMatrix::Minor(unsigned i,unsigned j) const
{
   CMatrix A(Rows_-1,Cols_-1);
   
   if (!IsNull() || !A.IsNull())
   {
      for (register int a=0;a<A.Rows_;a++)
      {
	 for (register int b=0;b<A.Cols_;b++)
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

      for (register int i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      if (InitVal!=SENTINAL)
      {
	 for (register int i=0;i<Rows_;i++)
	 {
	    for (register int j=0;j<Cols_;j++)
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

      for (register int i=1;i<Rows_;i++)
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

      for (register int i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      for (register int i=0;i<Rows_;i++)
	 for (register int j=0;j<Cols_;j++)
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
	   << endl;
      exit(-1);
   }
   else
   {
      CMatrix C(A.Rows_,A.Cols_);

      for (register int i=0;i<A.Rows_;i++)
      {
	 for (register int j=0;j<A.Cols_;j++)
	 {
	    C.Elements_[i][j]=A.Elements_[i][j]+B.Elements_[i][j];
	 }
      }

      return C;
   }
}

CMatrix operator-(const CMatrix& A)
{
   CMatrix B(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for(int j=0;j<A.Cols_;j++)
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
	   << endl;
      exit(-1);
   }
   else
   {
      CMatrix C(A.Rows_,A.Cols_);

      for (register int i=0;i<A.Rows_;i++)
      {
	 for (register int j=0;j<A.Cols_;j++)
	 {
	    C.Elements_[i][j]=A.Elements_[i][j]-B.Elements_[i][j];
	 }
      }

      return C;
   }
}

CMatrix operator*(const CMatrix& A,const CMatrix& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In CMatrix Operator* : A.Cols!=B.Rows or Null CMatrix"
	   <<endl;
      exit(-1);
   }
   else
   {
      CMatrix C(A.Rows_,B.Cols_,0);
      
      for (register int i=0;i<A.Rows_;i++)
      {
	 for (register int j=0;j<B.Cols_;j++)
	 {
	    for (register int k=0;k<A.Cols_;k++)
	    {
	       C.Elements_[i][j]+=A.Elements_[i][k]*B.Elements_[k][j];
	    }
	 }
      }

      return C;
   }
}

CMatrix operator*(const CMatrix::Elm& A,const CMatrix& B)
{
   CMatrix C(B.Rows_,B.Cols_);

   for (register int i=0;i<B.Rows_;i++)
   {
      for (register int j=0;j<B.Cols_;j++)
      {
	 C.Elements_[i][j]=A*B.Elements_[i][j];
      }
   }

   return C;
}

CMatrix operator*(const CMatrix& A,const CMatrix::Elm& B)
{
   CMatrix C(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
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
	   << endl;
      exit(-1);
   }
   
   CMatrix C(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=A.Elements_[i][j]/B;
      }
   }

   return C;
}

#ifdef CHECK_BOUNDS
CMatrix::Elm* CMatrix::operator[](unsigned i)
{
   if (i < Rows_)
      return Elements_[i];
   else
   {
      cerr << "CMatrix Index Overflow -- CMatrix::Elm* operator[]()" << endl;
      exit(-1);
   }
}

CMatrix::Elm* CMatrix::operator[](unsigned i) const
{
   if (i < Rows_)
      return Elements_[i];
   else
   {
      cerr << "CMatrix Index Overflow -- CMatrix::Elm* operator[]()" << endl;
      exit(-1);
   }
}
#endif
   
CMatrix& CMatrix::operator=(const CMatrix& B)
{
   if (Rows_!=B.Rows_ || Cols_!=B.Cols_ || IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix& operator=() : Matricies not same size "
           << "or Null CMatrix"
           << endl;
      exit(-1);
   }

   memmove(Elements_[0],B.Elements_[0],sizeof(CMatrix::Elm)*Rows_*Cols_);

   return *this;
}

CMatrix& CMatrix::SetIdentity(unsigned Size)
{
   if (Rows_!=Size || Cols_!=Size)
      Resize(Size,Size);

   for (register int i=0;i<Size;i++)
   {
      for (register int j=0;j<Size;j++)
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

   for (register int i=0;i<Rows_;i++)
   {
      for (register int j=0;j<Cols_;j++)
      {
	 A.Elements_[j][i]=Elements_[i][j];
      }
   }

   return A;
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

	 for (register int i=1;i<Rows_;i++)
	 {
	    Elements_[i]=Elements_[i-1]+Cols_;
	 }
      }

      if (InitVal!=SENTINAL)
      {
	 for (register int i=0;i<Rows_;i++)
	 {
	    for (register int j=0;j<Cols_;j++)
	    {
	       Elements_[i][j]=InitVal;
	    }
	 }
      }
   }
   else
   {
      if (InitVal != SENTINAL)
	 for (register int i=0;i<Rows_;i++)
	 {
	    for (register int j=0;j<Cols_;j++)
	    {
	       Elements_[i][j]=InitVal;
	    }
	 }
   }
   return;
}


// Recursivly calculate determinent
CMatrix::Elm CMatrix::Det() const
{
   if (IsNull() || !IsSquare())
   {
      cerr << "Error in CMatrix::Det() : Null or Non-Square CMatrix" << endl;
      exit(-1);
   }

   if (Rows_==1)
      return Elements_[0][0];
   else
   {
      Elm det(0.0,0.0);

      for (register int i=0;i<Cols_;i++)
      {
	 det += CMatrix::Elm(1-2*(i%2))*Elements_[0][i]*(Minor(0,i).Det());
      }

      return det;
   }
}

void Cholesky(const CMatrix& A,CMatrix& U,CMatrix& D)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in Cholesky() -- Non-Square or Null CMatrix" << endl;
      exit(-1);
   }

   U.SetIdentity(A.Rows_);
   D.Resize(A.Rows_,A.Cols_,0);
   

   for (register int i=0;i<A.Rows_;i++)
   {
      D.Elements_[i][i]=A.Elements_[i][i];
      if (i > 0)
      {
	 for (register int k=0;k<i;k++)
	 {
	    D.Elements_[i][i]-=D.Elements_[k][k]*U.Elements_[k][i]*U.Elements_[k][i];
	 }
      }
      for (register int j=i+1;j<A.Rows_;j++)
      {
	 U.Elements_[i][j]=A.Elements_[i][j];
	 if (i > 0)
	 {
	    for (register int k=0;k<i;k++)
	    {
	       U.Elements_[i][j]-=D.Elements_[k][k]*U.Elements_[k][i]*U.Elements_[k][j];
	    }
	 }
	 U.Elements_[i][j]/=D.Elements_[i][i];
      }
   }
   
   return;
}

ostream& operator<<(ostream& out,const CMatrix& A)
{
   int W=out.width();
   
   out << endl;

   if (CMatrix::MathematicaPrintFlag) out << setw(0) << "{{";
   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
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
	 out << endl;
   }

   out << endl;

   return out;
}

istream& operator>>(istream& in,CMatrix& A)
{
   for (register int i=0;i<A.Rows_;i++)
      for (register int j=0;j<A.Cols_;j++)
	 in >> A.Elements_[i][j];

   return in;
}

char* CMatrix::Revision()
{
   return CMatrixID;
}
