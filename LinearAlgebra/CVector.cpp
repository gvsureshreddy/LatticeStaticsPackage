#include "CVector.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Vector.h"
#include "Matrix.h"

// Global IDString
char CVectorID[]="$Id: CVector.cpp,v 1.4 2006/04/18 14:27:21 elliott Exp $";

// Private Functions...

// Public Functions...

CVector::CVector(const unsigned& Cols,const CVector::Elm& InitVal)
{
   Cols_=Cols;

   if (Cols_==0)
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new CVector::Elm[Cols_];

      if (InitVal!=SENTINAL)
      {
	 for (register int j=0;j<Cols_;j++)
	    Elements_[j]=InitVal;
      }
   }

   return;
}

CVector::CVector(const CVector& A)
{
   Cols_=A.Cols_;
   
   if (Cols_==0)
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new CVector::Elm[Cols_];
   }

   for (int i=0;i<Cols_;i++)
   {
      Elements_[i] = A.Elements_[i];
   }

   return;
}

CVector::CVector(const CMatrix& A)
{
   if (A.IsNull() || (A.Rows()!=1 && A.Cols()!=1))
   {
      cerr << "Error in CVector::CVector(CMatrix& A) -- Null CMatrix or Non-CVector"
	   << endl;
      exit(-1);
   }

   if (A.Rows()==1)
   {
      Cols_=A.Cols();
   
      if (Cols_==0)
      {
         Elements_=NULL;
      }
      else
      {
         Elements_=new CVector::Elm[Cols_];
      }

      memmove(Elements_,A[0],sizeof(CVector::Elm)*Cols_);
   }
   else
   {
      Cols_=A.Rows();
   
      if (Cols_==0)
      {
         Elements_=NULL;
      }
      else
      {
         Elements_=new CVector::Elm[Cols_];
      }

      memmove(Elements_,A[0],sizeof(CVector::Elm)*Cols_);
   }

   return;
}

CVector::CVector(const Vector& A)
{
   Cols_=A.Dim();
   
   if (Cols_==0)
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new CVector::Elm[Cols_];
   }

   for (int i=0;i<Cols_;i++)
   {
      Elements_[i] = CVector::Elm(A[i]);
   }

   return;
}

CVector::CVector(const Matrix& A)
{
   if (A.IsNull() || (A.Rows()!=1 && A.Cols()!=1))
   {
      cerr << "Error in CVector::CVector(Matrix& A) -- Null Matrix or Non-CVector"
	   << endl;
      exit(-1);
   }

   if (A.Rows()==1)
   {
      Cols_=A.Cols();
   
      if (Cols_==0)
      {
         Elements_=NULL;
      }
      else
      {
         Elements_=new CVector::Elm[Cols_];
      }

      for (int i=0;i<Cols_;i++)
      {
	 Elements_[i] = CVector::Elm(A[0][i]);
      }
   }
   else
   {
      Cols_=A.Rows();
   
      if (Cols_==0)
      {
         Elements_=NULL;
      }
      else
      {
         Elements_=new CVector::Elm[Cols_];
      }

      for (int i=0;i<Cols_;i++)
      {
	 Elements_[i] = CVector::Elm(A[0][i]);
      }
   }

   return;
}

CVector::~CVector()
{
   delete [] Elements_;
   
   return;   
}

CVector operator+(const CVector& A,const CVector& B)
{
   if (A.Cols_!=B.Cols_ || A.Cols_==0 || B.Cols_==0)
   {
      cerr << "Error in CVector Operator+() Diff Size CVectors or Null CVector!!!"
	   << endl;
      exit(-1);
   }

   CVector C(A.Cols_);
   
   for (register int j=0;j<A.Cols_;j++)
   {
      C.Elements_[j]=A.Elements_[j]+B.Elements_[j];
   }
   
   return C;
}

CVector operator-(const CVector& A,const CVector& B)
{
   if (A.Cols_!=B.Cols_ || A.Cols_==0 || B.Cols_==0)
   {
      cerr << "Error in CVector Operator-() Diff Size CVectors or Null CVector!!!"
	   << endl;
      exit(-1);
   }

   CVector C(A.Cols_);
   
   for (register int j=0;j<A.Cols_;j++)
   {
      C.Elements_[j]=A.Elements_[j]-B.Elements_[j];
   }
   
   return C;
}

CVector operator-(const CVector& A)
{
   CVector B(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
      B.Elements_[i]=-A.Elements_[i];

   return B;
}

// Dot Product
CVector::Elm operator*(const CVector& A,const CVector& B)
{
   if (A.Cols_==0 || B.Cols_==0 || A.Cols_!=B.Cols_)
   {
      cerr <<"Error in Dot Product -- Null CVector or different Dimensions."
	   << endl;
      exit(-1);
   }
   
   CVector::Elm sum=0;

   for (register int i=0;i<A.Cols_;i++)
   {
      sum+=A.Elements_[i]*(B.Elements_[i].conj());
   }

   return sum;
}

CVector operator*(const CMatrix& A,const CVector& B)
{
   if (A.Cols()!=B.Cols_ || A.IsNull() || B.Cols_==0)
   {
      cerr << "Error In CVector Operator* : A.Cols!=B.Cols or Null CMatrix or CVector"
	   <<endl;
      exit(-1);
   }

   CVector C(A.Rows(),0);
   
   for (register int i=0;i<A.Rows();i++)
   {
      for (register int k=0;k<B.Cols_;k++)
      {
	 C.Elements_[i]+=A[i][k]*B.Elements_[k];
      }
   }
   
   return C;
}

CVector operator*(const CVector& A,const CMatrix& B)
{
   if (B.Cols()!=A.Cols_ || B.IsNull() || A.Cols_==0)
   {
      cerr << "Error In CVector Operator* : A.Cols!=B.Cols or Null CMatrix or CVector"
	   <<endl;
      exit(-1);
   }

   CVector C(A.Cols_,0);
   
   for (register int i=0;i<B.Cols();i++)
   {
      for (register int k=0;k<A.Cols_;k++)
      {
	 C.Elements_[i]+=A.Elements_[k]*B[k][i];
      }
   }
   
   return C;
}

CVector operator*(const CVector::Elm& A,const CVector& B)
{
   CVector C(B.Cols_);

   for (register int i=0;i<B.Cols_;i++)
   {
      C.Elements_[i]=A*B.Elements_[i];
   }

   return C;
}

CVector operator*(const CVector& A,const CVector::Elm& B)
{
   CVector C(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
   {
      C.Elements_[i]=B*A.Elements_[i];
   }

   return C;
}

CVector operator/(const CVector& A,const CVector::Elm& B)
{
   CVector C(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
   {
      C.Elements_[i]=A.Elements_[i]/B;
   }

   return C;
}

#ifdef CHECK_BOUNDS
CVector::Elm& CVector::operator[](const unsigned& i)
{
   if (i>=Cols_)
   {
      cerr << "Error in CVector::operator[]() -- Index Overflow"
	   << endl;
      exit(-1);
   }

   return Elements_[i];
}

const CVector::Elm CVector::operator[](const unsigned& i) const
{
   if (i>=Cols_)
   {
      cerr << "Error in CVector::operator[]() -- Index Overflow"
	   << endl;
      exit(-1);
   }

   return Elements_[i];
}
#endif

CVector& CVector::operator=(const CVector& B)
{
   if (Cols_!=B.Cols_ || Cols_==0 || B.Cols_==0)
   {
      cerr << "Error in CVector& operator=() : CVectors not same size "
           << "or Null CVector"
           << endl;
      exit(-1);
   }

   memmove(Elements_,B.Elements_,sizeof(CVector::Elm)*Cols_);

   return *this;
}

void CVector::Resize(const unsigned& Cols,const CVector::Elm& InitVal)
{
   if (Cols!=Cols_)
   {
      delete [] Elements_;

      Cols_=Cols;
      
      if (Cols_==0)
      {
	 Elements_=NULL;
      }
      else
      {
	 Elements_=new CVector::Elm[Cols_];
      }

      if (InitVal!=SENTINAL)
      {
	 for (register int j=0;j<Cols_;j++)
	 {
	    Elements_[j]=InitVal;
	 }
      }
   }

   return;
}

CVector::Elm CVector::Norm()
{
   return sqrt(*this*(*this));
}

CVector SolvePLU(const CMatrix& A,const CVector& B)
{
   CMatrix C(B.Cols_,1);

   for(register int i=0;i<B.Cols_;i++)
   {
      C[i][0]=B.Elements_[i];
   }

      return SolvePLU(A,C);
}

ostream& operator<<(ostream& out,const CVector& A)
{
   int W=out.width();

   for (register int i=0;i<A.Cols_;i++)
   {
      out << setw(W) << A.Elements_[i];
   }

   return out;
}

istream& operator>>(istream& in,CVector& A)
{
   for (register int i=0;i<A.Cols_;i++)
      in >> A.Elements_[i];

   return in;
}

char* CVector::Revision()
{
   return CVectorID;
}
