#include "Matrix.h"
#include <iostream.h>
#include <iomanip.h>
#include <string.h>
#include <math.h>

// ***********************************************************************************
// $Log: Matrix.cpp,v $
// Revision 1.2  1999/08/16 14:49:26  elliottr
// Added Identification (Revision) Routine
//
// Revision 1.1  1999/07/26 17:36:21  elliottr
// Initial revision
//
// ***********************************************************************************


// Global IDString
char MatrixID[]="$Id: Matrix.cpp,v 1.2 1999/08/16 14:49:26 elliottr Exp $";

// Private Methods...


// Returns matrix of size Rows_-1 x Cols_-1 with ith row and
//   jth column removed
Matrix Matrix::Minor(unsigned i,unsigned j) const
{
   Matrix A(Rows_-1,Cols_-1);
   
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

Matrix::Matrix(unsigned Rows,unsigned Cols,Matrix::Elm InitVal)
{
   Rows_=Rows;
   Cols_=Cols;

   if (IsNull())
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new Matrix::Elm*[Rows_];
      Elements_[0]=new Matrix::Elm[Rows_*Cols_];

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

Matrix::Matrix(const Matrix& A)
{
   Rows_=A.Rows_;
   Cols_=A.Cols_;
   
   if (IsNull())
   {
      Elements_=NULL;
   }
   else
   {
      Elements_=new Matrix::Elm*[Rows_];
      Elements_[0]=new Matrix::Elm[Rows_*Cols_];

      for (register int i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      memmove(Elements_[0],A.Elements_[0],sizeof(Matrix::Elm[Rows_*Cols_]));
   }

   return;
}

Matrix::~Matrix()
{
   if (!IsNull())
   {
      delete [] Elements_[0];
      delete [] Elements_;
   }
   
   return;
}

Matrix operator+(const Matrix& A,const Matrix& B)
{
   if (A.Rows_!=B.Rows_ || A.Cols_!=B.Cols_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix Operator+() Diff Size Matrices or Null Matrix!!!"
	   << endl;
      exit(-1);
   }
   else
   {
      Matrix C(A.Rows_,A.Cols_);

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

Matrix operator-(const Matrix& A)
{
   Matrix B(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for(int j=0;j<A.Cols_;j++)
      {
	 B.Elements_[i][j]=-A.Elements_[i][j];
      }
   }

   return B;
}

Matrix operator-(const Matrix& A,const Matrix& B)
{
   if (A.Rows_!=B.Rows_ || A.Cols_!=B.Cols_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix Operator-() Diff Size Matrices or Null Matrix!!!"
	   << endl;
      exit(-1);
   }
   else
   {
      Matrix C(A.Rows_,A.Cols_);

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

Matrix operator*(const Matrix& A,const Matrix& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In Matrix Operator* : A.Cols!=B.Rows or Null Matrix"
	   <<endl;
      exit(-1);
   }
   else
   {
      Matrix C(A.Rows_,B.Cols_,0);
      
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

Matrix operator*(const Matrix::Elm& A,const Matrix& B)
{
   Matrix C(B.Rows_,B.Cols_);

   for (register int i=0;i<B.Rows_;i++)
   {
      for (register int j=0;j<B.Cols_;j++)
      {
	 C.Elements_[i][j]=A*B.Elements_[i][j];
      }
   }

   return C;
}

Matrix operator*(const Matrix& A,const Matrix::Elm& B)
{
   Matrix C(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=B*A.Elements_[i][j];
      }
   }

   return C;
}

Matrix operator/(const Matrix& A,const Matrix::Elm& B)
{
   if (B==0)
   {
      cerr << "Divide By Zero Error in Matrix operator/()"
	   << endl;
      exit(-1);
   }
   
   Matrix C(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C.Elements_[i][j]=A.Elements_[i][j]/B;
      }
   }

   return C;
}

Matrix::Elm* Matrix::operator[](unsigned i)
{
   if (i < Rows_)
      return Elements_[i];
   else
   {
      cerr << "Matrix Index Overflow -- Matrix::Elm* operator[]()" << endl;
      exit(-1);
   }
}

Matrix::Elm* Matrix::operator[](unsigned i) const
{
   if (i < Rows_)
      return Elements_[i];
   else
   {
      cerr << "Matrix Index Overflow -- Matrix::Elm* operator[]()" << endl;
      exit(-1);
   }
}
   
Matrix& Matrix::operator=(const Matrix& B)
{
   if (Rows_!=B.Rows_ || Cols_!=B.Cols_ || IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix& operator=() : Matricies not same size "
           << "or Null Matrix"
           << endl;
      exit(-1);
   }

   memmove(Elements_[0],B.Elements_[0],sizeof(Matrix::Elm[Rows_*Cols_]));

   return *this;
}

Matrix& Matrix::SetIdentity(unsigned Size)
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
	    Elements_[i][j]=0;
      } 
   }

   return *this;
}

Matrix Matrix::Transpose() const
{
   Matrix A(Cols_,Rows_);

   for (register int i=0;i<Rows_;i++)
   {
      for (register int j=0;j<Cols_;j++)
      {
	 A.Elements_[j][i]=Elements_[i][j];
      }
   }

   return A;
}

Matrix Matrix::Inverse() const
{
   if (!IsSquare() || IsNull())
   {
      cerr << "Error in Matrix::Inverse() : Non-Square or Null Matrix" << endl;
      exit(-1);
   }

   Matrix B(Rows_,1,0),X(Rows_,1),C(Rows_,Cols_);

   B.Elements_[0][0]=1.0;
   for (register int i=0;i<Cols_;i++)
   {
      X=Solve(*this,B);

      for (register int j=0;j<Rows_;j++)
	 C.Elements_[j][i]=X.Elements_[j][0];

      B.Elements_[i][0]=0;
      if (i!=Cols_-1) B.Elements_[i+1][0]=1.0;
   }

   return C;
}

void Matrix::Resize(unsigned Rows,unsigned Cols,Matrix::Elm InitVal)
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
	 Elements_=new Matrix::Elm*[Rows_];
	 Elements_[0]=new Matrix::Elm[Rows_*Cols_];

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
   return;
}


// Recursivly calculate determinent
Matrix::Elm Matrix::Det() const
{
   if (IsNull() || !IsSquare())
   {
      cerr << "Error in Matrix::Det() : Null or Non-Square Matrix" << endl;
      exit(-1);
   }

   if (Rows_==1)
      return Elements_[0][0];
   else
   {
      Matrix::Elm det=0;

      for (register int i=0;i<Cols_;i++)
      {
	 det+=(1-2*(i%2))*Elements_[0][i]*(Minor(0,i).Det());
      }

      return det;
   }
}


// Decompose PA=LU using scaled partial pivoting.
void PLU(const Matrix& A,Matrix& P,Matrix& L,Matrix& U)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in PLU -- Non-Square or Null Matrix to decompose..." << endl;
      exit(-1);
   }

   P.Resize(A.Rows_,A.Cols_,0);
   L.SetIdentity(A.Rows_);
   U.Resize(A.Rows_,A.Cols_,0);

   Matrix Temp=A,
      S(A.Rows_,1,0);
   int Ipivot[A.Rows_];

   for (register int i=0;i<A.Rows_;i++)
   {
      Ipivot[i]=i;
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 if (fabs(Temp.Elements_[i][j]) > S.Elements_[i][0])
	    S.Elements_[i][0]=fabs(Temp.Elements_[i][j]);
      }
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      Matrix::Elm temp1;
      temp1=fabs(Temp.Elements_[i][i]/S.Elements_[i][0]);

      int k=i;
      for (register int j=i;j<A.Rows_;j++)
      {
	 if (fabs(Temp.Elements_[j][i]) > temp1)
	 {
	    temp1=fabs(Temp.Elements_[j][i]/S.Elements_[j][0]);
	    k=j;
	 }
      }

      if (k>i)
      {
	 for (register int j=i;j<A.Rows_;j++)
	 {
	    Matrix::Elm Switch[A.Rows_];
	    Switch[j]=Temp.Elements_[i][j];
	    Temp.Elements_[i][j]=Temp.Elements_[k][j];
	    Temp.Elements_[k][j]=Switch[j];
	 }

	 for (register int j=0;j<i;j++)
	 {
	    Matrix::Elm Switch[A.Rows_];
	    Switch[j]=L.Elements_[i][j];
	    L.Elements_[i][j]=L.Elements_[k][j];
	    L.Elements_[k][j]=Switch[j];
	 }

	 temp1=S.Elements_[i][0];
	 S.Elements_[i][0]=S.Elements_[k][0];
	 S.Elements_[k][0]=temp1;

	 int tempi1=Ipivot[i];
	 Ipivot[i]=Ipivot[k];
	 Ipivot[k]=tempi1;
      }

      for (register int j=i+1;j<A.Rows_;j++)
      {
	 L.Elements_[j][i]=Temp.Elements_[j][i]/Temp.Elements_[i][i];
      }

      for (register int j=i+1;j<A.Rows_;j++)
      {
	 for (register int k=i+1;k<A.Rows_;k++)
	 {
	    Temp.Elements_[j][k]=Temp.Elements_[j][k]-(L.Elements_[j][i]*Temp.Elements_[i][k]);
	 }
      }
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=i;j<A.Rows_;j++)
      {
	 U.Elements_[i][j]=Temp.Elements_[i][j];
      }
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      P.Elements_[i][Ipivot[i]]=1;
   }

   return;
}

void Cholesky(const Matrix& A,Matrix& U,Matrix& D)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in Cholesky() -- Non-Square or Null Matrix" << endl;
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

Matrix Solve(const Matrix& A,const Matrix& B)
{
   if (!A.IsSquare() || A.IsNull() || A.Cols_!=B.Rows_)
   {
      cerr << "Error in Solve() - Non-Square Matrix, Null Matrix, or system of != dimension"
	   << endl;
      exit(-1);
   }
   
   Matrix
      P,L,U; //PLU() resizes P,L,& U thus do not waste time initializing them.

   PLU(A,P,L,U);

   Matrix::Elm Y[B.Rows_];
   Matrix Temp=P*B;

   Y[0]=Temp.Elements_[0][0];
   for (register int i=1;i<B.Rows_;i++)
   {
      Y[i]=Temp[i][0];
      for (register int j=0;j<i;j++)
      {
	 Y[i]-=Y[j]*L.Elements_[i][j];
      }
   }

   Matrix X(B.Rows_,1,0);

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

   return X;
}

ostream& operator<<(ostream& out,const Matrix& A)
{
   int W=out.width();
   
   out << endl;

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 out << setw(W) << A.Elements_[i][j];
      }
      out << endl;
   }

   out << endl;

   return out;
}

istream& operator>>(istream& in,Matrix& A)
{
   for (register int i=0;i<A.Rows_;i++)
      for (register int j=0;j<A.Cols_;j++)
	 in >> A.Elements_[i][j];

   return in;
}

char* Matrix::Revision()
{
   return MatrixID;
}
