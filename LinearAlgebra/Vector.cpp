#include "Vector.h"
#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

// Global IDString
char VectorID[]="$Id: Vector.cpp,v 1.19 2009/09/09 18:19:49 elliott Exp $";

// Private Functions...

// Public Functions...

Vector::Vector(int const& Cols)
{
   Cols_=Cols;

   if (Cols_==0)
   {
      Elements_=0;
   }
   else
   {
      Elements_=new Vector::Elm[Cols_];
   }

   return;
}

Vector::Vector(int const& Cols,Vector::Elm const& InitVal)
{
   Cols_=Cols;

   if (Cols_==0)
   {
      Elements_=0;
   }
   else
   {
      Elements_=new Vector::Elm[Cols_];

      for (register int j=0;j<Cols_;j++)
         Elements_[j]=InitVal;
   }

   return;
}

Vector::Vector(Vector const& A)
{
   Cols_=A.Cols_;
   
   if (Cols_==0)
   {
      Elements_=0;
   }
   else
   {
      Elements_=new Vector::Elm[Cols_];
   }

   memmove(Elements_,A.Elements_,sizeof(Vector::Elm)*Cols_);

   return;
}

Vector::Vector(Vector3D const& A)
{
   Cols_=3;
   
   if (Cols_==0)
   {
      Elements_=0;
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

Vector::Vector(Matrix const& A)
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
         Elements_=0;
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
         Elements_=0;
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

Vector operator+(Vector const& A,Vector const& B)
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
	 C[j]=A[j]+B[j];
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator-(Vector const& A,Vector const& B)
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
	 C[j]=A[j]-B[j];
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator-(Vector const& A)
{
   Vector B(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
      B[i]=-A[i];

   return B;
}

// Dot Product
Vector::Elm const operator*(Vector const& A,Vector const& B)
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
      sum+=A[i]*B[i];
   }

   return sum;
}

// Cross Product
Vector operator%(Vector const& A,Vector const& B)
{
   if (A.Cols_==0 || B.Cols_==0 || A.Cols_!=B.Cols_ || A.Cols_!=3)
   {
      cerr << "Error in Cross Product -- Null Vector of different Dimensions"
	   << " or Cols_!=3" << "\n";
      exit(-1);
   }
   
   Vector C(A.Cols_);

   C[0]=A[1]*B[2]-A[2]*B[1];
   C[1]=-(A[0]*B[2]-A[2]*B[0]);
   C[2]=A[0]*B[1]-A[1]*B[0];

   return C;
}

Vector operator*(Matrix const& A,Vector const& B)
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
	    C[i]+=A[i][k]*B[k];
	 }
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator*(Vector const& A,Matrix const& B)
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
	    C[i]+=A[k]*B[k][i];
	 }
      }

      return C;
   }
   return Vector(); // dummy statement to avoid warning
}

Vector operator*(Vector::Elm const& A,Vector const& B)
{
   Vector C(B.Cols_);

   for (register int i=0;i<B.Cols_;i++)
   {
      C[i]=A*B[i];
   }

   return C;
}

Vector operator*(Vector const& A,Vector::Elm const& B)
{
   Vector C(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
   {
      C[i]=B*A[i];
   }

   return C;
}

Vector operator/(Vector const& A,Vector::Elm const& B)
{
   Vector C(A.Cols_);

   for (register int i=0;i<A.Cols_;i++)
   {
      C[i]=A[i]/B;
   }

   return C;
}

#ifdef CHECK_BOUNDS
Vector::Elm& Vector::operator[](int const& i)
{
   if (i>=Cols_)
   {
      cerr << "Error in Vector::operator[]() -- Index Overflow"
	   << "\n";
      exit(-1);
   }
   return Elements_[i];
}

Vector::Elm const& Vector::operator[](int const& i) const
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

Vector& Vector::operator=(Vector const& B)
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

void Vector::Resize(int const& Cols)
{
   if (Cols!=Cols_)
   {
      delete [] Elements_;

      Cols_=Cols;
      
      if (Cols_==0)
      {
	 Elements_=0;
      }
      else
      {
	 Elements_=new Vector::Elm[Cols_];
      }
   }

   return;
}

void Vector::Resize(int const& Cols,Matrix::Elm const& InitVal)
{
   if (Cols!=Cols_)
   {
      delete [] Elements_;

      Cols_=Cols;
      
      if (Cols_==0)
      {
	 Elements_=0;
      }
      else
      {
	 Elements_=new Vector::Elm[Cols_];
      }
   }
   
   for (register int j=0;j<Cols_;j++)
   {
      Elements_[j]=InitVal;
   }

   return;
}

Matrix::Elm Vector::Norm() const
{
   return sqrt(*this*(*this));
}

void SolveQR(Matrix const& Q,Matrix const& R,Vector& x,Vector const& B)
{
   Matrix C(B.Cols_,1);
   for (register int i=0;i<B.Cols_;++i)
   {
      C[i][0] = B[i];
   }
   Matrix X(x.Cols_,1);
   SolveQR(Q,R,X,C);
   for (register int i=0;i<x.Cols_;++i)
   {
      x[i] = X[i][0];
   }
}

Vector SolvePLU(Matrix const& A,Vector const& B)
{
   Matrix C(B.Cols_,1);

   for(register int i=0;i<B.Cols_;i++)
   {
      C[i][0]=B[i];
   }

      return SolvePLU(A,C);
}

Vector SolveSVD(Matrix const& A,Vector const& B,Vector::Elm const& MaxCond,
		int const& PrintFlag)
{
   Matrix C(B.Cols_,1);

   for(register int i=0;i<B.Cols_;i++)
   {
      C[i][0]=B[i];
   }

      return SolveSVD(A,C,MaxCond,PrintFlag);
}

void BroydenQRUpdate(Matrix& Q,Matrix& R,Vector const& y,Vector const& x)
{
   Matrix Y(y.Cols_,1);
   Matrix X(x.Cols_,1);
   for (int i=0;i<y.Cols_;++i)
   {
      Y[i][0] = y[i];
   }
   for (int i=0;i<x.Cols_;++i)
   {
      X[i][0] = x[i];
   }

   BroydenQRUpdate(Q,R,Y,X);
}

ostream& operator<<(ostream& out,Vector const& A)
{
   int W=out.width();

   for (register int i=0;i<A.Cols_;i++)
   {
      out << setw(W) << A[i];
   }

   return out;
}

istream& operator>>(istream& in,Vector& A)
{
   for (register int i=0;i<A.Cols_;i++)
      in >> A[i];

   return in;
}

char const* const Vector::Revision()
{
   return VectorID;
}
