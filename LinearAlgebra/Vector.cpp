#include "Vector.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>

// Global IDString
char VectorID[]="$Id: Vector.cpp,v 1.13 2008/03/23 02:26:56 elliott Exp $";

// Private Functions...

// Public Functions...

Vector::Vector(const unsigned& Cols,const Vector::Elm& InitVal)
{
   Cols_=Cols;

   if (Cols_==0)
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new Vector::Elm[Cols_];

      if (InitVal!=SENTINAL)
      {
	 for (register int j=0;j<Cols_;j++)
	    Elements_[j]=InitVal;
      }
   }

   return;
}

Vector::Vector(const Vector& A)
{
   Cols_=A.Cols_;
   
   if (Cols_==0)
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new Vector::Elm[Cols_];
   }

   memmove(Elements_,A.Elements_,sizeof(Vector::Elm)*Cols_);

   return;
}

Vector::Vector(const Vector3D& A)
{
   Cols_=3;
   
   if (Cols_==0)
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new Vector::Elm[Cols_];
   }

   for (int i=0;i<3;i++)
   {
      Elements_[i] = A[i];
   }

   return;
}

Vector::Vector(const Matrix& A)
{
   if (A.IsNull() || (A.Rows()!=1 && A.Cols()!=1))
   {
      cerr << "Error in Vector::Vector(Matrix& A) -- Null Matrix or Non-Vector"
	   << "\n";
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
         Elements_=new Vector::Elm[Cols_];
      }

      memmove(Elements_,A[0],sizeof(Vector::Elm)*Cols_);
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
         Elements_=new Vector::Elm[Cols_];
      }

      memmove(Elements_,A[0],sizeof(Vector::Elm)*Cols_);
   }

   return;
}

Vector::~Vector()
{
   delete [] Elements_;
   
   return;   
}

Vector operator+(const Vector& A,const Vector& B)
{
   if (A.Cols_!=B.Cols_ || A.Cols_==0 || B.Cols_==0)
   {
      cerr << "Error in Vector Operator+() Diff Size Vectors or Null Vector!!!"
	   << "\n";
      exit(-1);
   }
   else
   {
      Vector C(A.Cols_);

      for (register int j=0;j<A.Cols_;j++)
      {
	 C.Elements_[j]=A.Elements_[j]+B.Elements_[j];
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator-(const Vector& A,const Vector& B)
{
   if (A.Cols_!=B.Cols_ || A.Cols_==0 || B.Cols_==0)
   {
      cerr << "Error in Vector Operator-() Diff Size Vectors or Null Vector!!!"
	   << "\n";
      exit(-1);
   }
   else
   {
      Vector C(A.Cols_);

      for (register int j=0;j<A.Cols_;j++)
      {
	 C.Elements_[j]=A.Elements_[j]-B.Elements_[j];
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator-(const Vector& A)
{
   Vector B(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
      B.Elements_[i]=-A.Elements_[i];

   return B;
}

// Dot Product
Vector::Elm operator*(const Vector& A,const Vector& B)
{
   if (A.Cols_==0 || B.Cols_==0 || A.Cols_!=B.Cols_)
   {
      cerr <<"Error in Dot Product -- Null Vector or different Dimensions."
	   << "\n";
      exit(-1);
   }
   
   Vector::Elm sum=0;

   for (register int i=0;i<A.Cols_;i++)
   {
      sum+=A.Elements_[i]*B.Elements_[i];
   }

   return sum;
}

// Cross Product
Vector operator%(const Vector& A,const Vector& B)
{
   if (A.Cols_==0 || B.Cols_==0 || A.Cols_!=B.Cols_ || A.Cols_!=3)
   {
      cerr << "Error in Cross Product -- Null Vector of different Dimensions"
	   << " or Cols_!=3" << "\n";
      exit(-1);
   }
   
   Vector C(A.Cols_);

   C.Elements_[0]=A.Elements_[1]*B.Elements_[2]
                  -A.Elements_[2]*B.Elements_[1];
   C.Elements_[1]=-(A.Elements_[0]*B.Elements_[2]
                  -A.Elements_[2]*B.Elements_[0]);
   C.Elements_[2]=A.Elements_[0]*B.Elements_[1]
                  -A.Elements_[1]*B.Elements_[0];

   return C;
}

Vector operator*(const Matrix& A,const Vector& B)
{
   if (A.Cols()!=B.Cols_ || A.IsNull() || B.Cols_==0)
   {
      cerr << "Error In Vector Operator* : A.Cols!=B.Cols or Null Matrix or Vector"
	   <<"\n";
      exit(-1);
   }
   else
   {
      Vector C(A.Rows(),0);
      
      for (register int i=0;i<A.Rows();i++)
      {
	 for (register int k=0;k<B.Cols_;k++)
	 {
	    C.Elements_[i]+=A.Elements_[i][k]*B.Elements_[k];
	 }
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator*(const Vector& A,const Matrix& B)
{
   if (B.Cols()!=A.Cols_ || B.IsNull() || A.Cols_==0)
   {
      cerr << "Error In Vector Operator* : A.Cols!=B.Cols or Null Matrix or Vector"
	   <<"\n";
      exit(-1);
   }
   else
   {
      Vector C(A.Cols_,0);
      
      for (register int i=0;i<B.Cols();i++)
      {
	 for (register int k=0;k<A.Cols_;k++)
	 {
	    C.Elements_[i]+=A.Elements_[k]*B.Elements_[k][i];
	 }
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator*(const Vector::Elm& A,const Vector& B)
{
   Vector C(B.Cols_);

   for (register int i=0;i<B.Cols_;i++)
   {
      C.Elements_[i]=A*B.Elements_[i];
   }

   return C;
}

Vector operator*(const Vector& A,const Vector::Elm& B)
{
   Vector C(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
   {
      C.Elements_[i]=B*A.Elements_[i];
   }

   return C;
}

Vector operator/(const Vector& A,const Vector::Elm& B)
{
   Vector C(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
   {
      C.Elements_[i]=A.Elements_[i]/B;
   }

   return C;
}

#ifdef CHECK_BOUNDS
Vector::Elm& Vector::operator[](const unsigned& i)
{
   if (i>=Cols_)
   {
      cerr << "Error in Vector::operator[]() -- Index Overflow"
	   << "\n";
      exit(-1);
   }
   return Elements_[i];
}

Vector::Elm Vector::operator[](const unsigned& i) const
{
   if (i>=Cols_)
   {
      cerr << "Error in Vector::operator[]() -- Index Overflow"
	   << "\n";
      exit(-1);
   }
   return Elements_[i];
}
#endif

Vector& Vector::operator=(const Vector& B)
{
   if (Cols_!=B.Cols_ || Cols_==0 || B.Cols_==0)
   {
      cerr << "Error in Vector& operator=() : Vectors not same size "
           << "or Null Vector"
           << "\n";
      exit(-1);
   }

   memmove(Elements_,B.Elements_,sizeof(Vector::Elm)*Cols_);

   return *this;
}

void Vector::Resize(const unsigned& Cols,const Matrix::Elm& InitVal)
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
	 Elements_=new Vector::Elm[Cols_];
      }
   }

   if (InitVal!=SENTINAL)
   {
      for (register int j=0;j<Cols_;j++)
      {
	 Elements_[j]=InitVal;
      }
   }

   return;
}

Matrix::Elm Vector::Norm()
{
   return sqrt(*this*(*this));
}

Vector SolvePLU(const Matrix& A,const Vector& B)
{
   Matrix C(B.Cols_,1);

   for(register int i=0;i<B.Cols_;i++)
   {
      C[i][0]=B.Elements_[i];
   }

      return SolvePLU(A,C);
}

Vector SolveSVD(const Matrix& A,const Vector& B,const Vector::Elm MaxCond,
		const int PrintFlag)
{
   Matrix C(B.Cols_,1);

   for(register int i=0;i<B.Cols_;i++)
   {
      C[i][0]=B.Elements_[i];
   }

      return SolveSVD(A,C,MaxCond,PrintFlag);
}

ostream& operator<<(ostream& out,const Vector& A)
{
   int W=out.width();

   for (register int i=0;i<A.Cols_;i++)
   {
      out << setw(W) << A.Elements_[i];
   }

   return out;
}

istream& operator>>(istream& in,Vector& A)
{
   for (register int i=0;i<A.Cols_;i++)
      in >> A.Elements_[i];

   return in;
}

char* Vector::Revision()
{
   return VectorID;
}
