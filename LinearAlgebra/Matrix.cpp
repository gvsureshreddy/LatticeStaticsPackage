#include "Matrix.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cstdlib>

// Global IDString
char MatrixID[]="$Id: Matrix.cpp,v 1.15 2005/03/17 19:33:25 elliott Exp $";

// Private Methods...


// Computes sqrt(a^2 + b^2) without destructive underflow or overflow.
Matrix::Elm pythag(Matrix::Elm a,Matrix::Elm b);
Matrix::Elm pythag(Matrix::Elm a,Matrix::Elm b)
{
   Matrix::Elm absa,absb;
   absa=fabs(a);
   absb=fabs(b);
   if (absa > absb)
      return absa*sqrt(1.0+((absb/absa)*(absb/absa)));
   else
      return (absb == 0.0 ?
	      0.0 :
	      absb*sqrt(1.0+((absa/absb)*(absa/absb))));
}


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

int Matrix::MathematicaPrintFlag = 0;

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

      memmove(Elements_[0],A.Elements_[0],sizeof(Matrix::Elm)*Rows_*Cols_);
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

Matrix operator*(const Matrix& A,const Matrix& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In Matrix Operator* : A.Cols!=B.Rows or Null Matrix"
	   <<endl;
      exit(-1);
   }

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

#ifdef CHECK_BOUNDS
Matrix::Elm* Matrix::operator[](unsigned i)
{
   if (i >= Rows_)
   {
      cerr << "Matrix Index Overflow -- Matrix::Elm* operator[]()" << endl;
      exit(-1);
   }

   return Elements_[i];
}

Matrix::Elm* Matrix::operator[](unsigned i) const
{
   if (i >= Rows_)
   {
      cerr << "Matrix Index Overflow -- Matrix::Elm* operator[]()" << endl;
      exit(-1);
   }

   return Elements_[i];
}
#endif
   
Matrix& Matrix::operator=(const Matrix& B)
{
   if (Rows_!=B.Rows_ || Cols_!=B.Cols_ || IsNull() || B.IsNull())
   {
      cerr << "Error in Matrix& operator=() : Matricies not same size "
           << "or Null Matrix"
           << endl;
      exit(-1);
   }

   memmove(Elements_[0],B.Elements_[0],sizeof(Matrix::Elm)*Rows_*Cols_);

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
      X=SolvePLU(*this,B);

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
   int *Ipivot;
   Ipivot = new int[A.Rows_];

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
	 Matrix::Elm *Switch;
	 Switch = new Matrix::Elm[A.Rows_];
	 for (register int j=i;j<A.Rows_;j++)
	 {
	    Switch[j]=Temp.Elements_[i][j];
	    Temp.Elements_[i][j]=Temp.Elements_[k][j];
	    Temp.Elements_[k][j]=Switch[j];
	 }

	 for (register int j=0;j<i;j++)
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

   delete [] Ipivot;

   return;
}

// Singular Value Decomposition of A -- Algorithm from Numerical Recipies
//
// return value - condition number of A
// A - mxn matrix to decompose
// U - mxn "column-orthogonal" matrix
// W - nxn diagonal matrix (singular values)
// V - nxn orthogonal matrix
//
// A = U*W*V.Transpose();
//
Matrix::Elm SVD(const Matrix& A,Matrix& U,Matrix& W,Matrix& V,
		const Matrix::Elm MaxCond,const int PrintFlag)
{
   // Initialize U = A
   U.Resize(A.Rows_,A.Cols_);
   U = A;

   // Initialize W and V
   W.Resize(A.Cols_,A.Cols_,0.0);
   V.Resize(A.Cols_,A.Cols_);

   // allocate temp storage space
   Matrix::Elm *temp;
   temp = new Matrix::Elm[A.Cols_];

   // define local variables
   int flag,
      l,
      nm;
   Matrix::Elm anorm,
      c,f,g,h,s,scale,x,y,z;

   g=scale=anorm=0.0;

   // Householder reduction to bidiagonal form.
   for (int i=0;i<A.Cols_;i++)
   {
      l = i+1;
      temp[i] = scale*g;
      g=s=scale=0.0;

      if (i < A.Rows_)
      {
	 for (int k=i;k<A.Rows_;k++) scale += fabs(U.Elements_[k][i]);
	 if (scale)
	 {
	    for (int k=i;k<A.Rows_;k++)
	    {
	       U.Elements_[k][i] /= scale;
	       s += U.Elements_[k][i]*U.Elements_[k][i];
	    }
	    f = U.Elements_[i][i];
	    g = - ( f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)) );
	    h = f*g - s;
	    U.Elements_[i][i] = f-g;
	    for (int j=l;j<A.Cols_;j++)
	    {
	       s = 0.0;
	       for (int k=i;k<A.Rows_;k++)
		  s += U.Elements_[k][i]*U.Elements_[k][j];
	       f = s/h;
	       for (int k=i;k<A.Rows_;k++)
		  U.Elements_[k][j] += f*U.Elements_[k][i];
	    }

	    for (int k=i;k<A.Rows_;k++) U.Elements_[k][i] *= scale;
	 }
      }
      W.Elements_[i][i] = scale*g;

      g=s=scale=0.0;
      if ((i < A.Rows_) && (i != A.Cols_-1))
      {
	 for (int k=l;k<A.Cols_;k++) scale += fabs(U.Elements_[i][k]);
	 if (scale)
	 {
	    for (int k=l;k<A.Cols_;k++)
	    {
	       U.Elements_[i][k] /= scale;
	       s += U.Elements_[i][k]*U.Elements_[i][k];
	    }
	    f = U.Elements_[i][l];
	    g = - ( f >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)) );
	    h = f*g - s;
	    U.Elements_[i][l] = f-g;
	    for (int k=l;k<A.Cols_;k++) temp[k] = U.Elements_[i][k]/h;
	    for (int j=l;j<A.Rows_;j++)
	    {
	       s = 0.0;
	       for (int k=l;k<A.Cols_;k++)
		  s += U.Elements_[j][k]*U.Elements_[i][k];
	       for (int k=l;k<A.Cols_;k++) U.Elements_[j][k] += s*temp[k];
	    }
	    for (int k=l;k<A.Cols_;k++) U.Elements_[i][k] *= scale;
	 }
      }
      anorm = ( anorm > (fabs(W.Elements_[i][i])+fabs(temp[i])) ?
		anorm :
		(fabs(W.Elements_[i][i])+fabs(temp[i])));
   }

   // Accumulation of right-hand transformations.
   for (int i=A.Cols_-1;i>-1;i--)
   {
      if (i < A.Cols_-1)
      {
	 if (g)
	 {
	    for (int j=l;j<A.Cols_;j++)
	       // Double division to avoid possible underflow
	       V.Elements_[j][i] = (U.Elements_[i][j]/U.Elements_[i][l])/g;
	    for (int j=l;j<A.Cols_;j++)
	    {
	       s = 0.0;
	       for (int k=l;k<A.Cols_;k++)
		  s += U.Elements_[i][k]*V.Elements_[k][j];
	       for (int k=l;k<A.Cols_;k++)
		  V.Elements_[k][j] += s*V.Elements_[k][i];
	    }
	 }
	 for (int j=l;j<A.Cols_;j++)
	    V.Elements_[i][j] = V.Elements_[j][i]=0.0;
      }
      V.Elements_[i][i]=1.0;
      g = temp[i];
      l = i;
   }

   // Accumulation of left-hand transformations.
   for (int i=(A.Cols_ < A.Rows_ ? A.Cols_ : A.Rows_)-1;i>-1;i--)
   {

      l = i+1;
      g = W.Elements_[i][i];
      for (int j=l;j<A.Cols_;j++) U.Elements_[i][j] = 0.0;
      if (g)
      {
	 g = 1.0/g;
	 for (int j=l;j<A.Cols_;j++)
	 {
	    s = 0.0;
	    for (int k=l;k<A.Rows_;k++)
	       s += U.Elements_[k][i]*U.Elements_[k][j];
	    f = (s/U.Elements_[i][i])*g;
	    for (int k=i;k<A.Rows_;k++)
	       U.Elements_[k][j] += f*U.Elements_[k][i];
	 }
	 for (int j=i;j<A.Rows_;j++) U.Elements_[j][i] *= g;
      }
      else
      {
	 for (int j=i;j<A.Rows_;j++) U.Elements_[j][i] = 0.0;
      }
      ++U.Elements_[i][i];
   }

   // Diagonalization of the bidiagonal form: Loop over singular values, and
   // -- over allowed iterations.
   for (int k=A.Cols_-1;k>-1;k--)
   {
      for (int its=0;its<30;its++)
      {
	 flag = 1;
	 // Test for splitting
	 for (l=k;l>-1;l--)
	 {
	    // Note that temp[0] is always zero
	    nm=l-1;
	    if (fabs(temp[l])+anorm == anorm)
	    {
	       flag = 0;
	       break;
	    }
	    if (fabs(W.Elements_[nm][nm])+anorm == anorm) break;
	 }
	 if (flag)
	 {
	    // Cancellation of temp[l], if l > 0
	    c = 0.0;
	    s = 1.0;
	    for (int i=l;i<=k;i++)
	    {
	       f = s*temp[i];
	       temp[i]=c*temp[i];
	       if (fabs(f)+anorm == anorm) break;
	       g = W.Elements_[i][i];
	       h = pythag(f,g);
	       W.Elements_[i][i] = h;
	       h = 1.0/h;
	       c = g*h;
	       s = -f*h;
	       for (int j=0;j<A.Rows_;j++)
	       {
		  y = U.Elements_[j][nm];
		  z = U.Elements_[j][i];
		  U.Elements_[j][nm] = y*c + z*s;
		  U.Elements_[j][i] = z*c - y*s;
	       }
	    }
	 }
	 z = W.Elements_[k][k];

	 // Convergence
	 if (l == k)
	 {
	    // Singular value is made nonnegative.
	    if (z < 0.0)
	    {
	       W.Elements_[k][k] = -z;
	       for (int j=0;j<A.Cols_;j++)
		  V.Elements_[j][k] = -V.Elements_[j][k];
	    }
	    break;
	 }
	 if (its == 30)
	 {
	    cerr << "no convergence in 30 SVD iterations" << endl;
	    exit(-1);
	 }
	 // Shift from bottom 2-by-2 minor
	 x = W.Elements_[l][l];
	 nm = k-1;
	 y = W.Elements_[nm][nm];
	 g = temp[nm];
	 h = temp[k];
	 f = ((y-z)*(y+z) + (g-h)*(g+h))/(2.0*h*y);
	 g = pythag(f,1.0);
	 f = ((x-z)*(x+z) + h*((y/(f + (f >=0.0 ? fabs(g) : -fabs(g)))) - h))/x;
	 // Next QR transformation:
	 c=s=1.0;
	 for (int j=l;j<=nm;j++)
	 {
	    int i;
	    i = j+1;
	    g = temp[i];
	    y = W.Elements_[i][i];
	    h = s*g;
	    g = c*g;
	    z = pythag(f,h);
	    temp[j] = z;
	    c = f/z;
	    s = h/z;
	    f = x*c + g*s;
	    g = g*c - x*s;
	    h = y*s;
	    y *= c;
	    for (int jj=0;jj<A.Cols_;jj++)
	    {
	       x = V.Elements_[jj][j];
	       z = V.Elements_[jj][i];
	       V.Elements_[jj][j] = x*c + z*s;
	       V.Elements_[jj][i] = z*c - x*s;
	    }
	    z = pythag(f,h);
	    // Rotation can be arbitrary if z = 0
	    W.Elements_[j][j] = z;
	    if (z)
	    {
	       z = 1.0/z;
	       c = f*z;
	       s = h*z;
	    }
	    f = c*g + s*y;
	    x = c*y - s*g;
	    for (int jj=0;jj<A.Rows_;jj++)
	    {
	       y = U.Elements_[jj][j];
	       z = U.Elements_[jj][i];
	       U.Elements_[jj][j] = y*c + z*s;
	       U.Elements_[jj][i] = z*c - y*s;
	    }
	 }
	 temp[l] = 0.0;
	 temp[k] = f;
	 W.Elements_[k][k] = x;
      }
   }

   delete [] temp;

   // Condition number stuff...
   // Remember the singular values are >= 0.0
   Matrix::Elm
      ConditionNumber,
      wmax=0.0,
      wmin;
   for (int j=0;j<A.Cols_;j++)
      if (W.Elements_[j][j] > wmax) wmax = W.Elements_[j][j];
   wmin = wmax;
   for (int j=0;j<A.Cols_;j++)
      if (W.Elements_[j][j] < wmin) wmin = W.Elements_[j][j];
   
   ConditionNumber = wmax/wmin;
   if (PrintFlag)
   {
      cerr << "SVD: Condition Number is : " << ConditionNumber << endl;
   }

   // Fix up any singular values that are "too small"
   for (int j=0;j<A.Cols_;j++)
      if (W.Elements_[j][j] < wmax/MaxCond)
      {
	 W.Elements_[j][j] = 0.0;
	 cerr << "SVD: Explicitly set Singular Value #" << j
	      << " to 0.0  !!!" << endl;
      }

   return ConditionNumber;
}

Matrix SymEigVal(Matrix A,Matrix *B,const int MaxItr,const double Tol)
{
   int count=0,
      converged=0;
   Matrix EigVals(1,A.Cols_);
   double theta,c,s,cc,ss,cs,aij1,aii1,ajj1,aki1,akj1,tmp;
   const double PIby4=0.25*acos(-1.0);

   if (B != NULL)
   {
      B->SetIdentity(A.Cols_);
   }
   
   while ((count < MaxItr) && (!converged))
   {
      for (int i=0;i<A.Cols_;i++)
      {
	 for (int j=i+1;j<A.Cols_;j++)
	 {
	    if (fabs(A.Elements_[i][j]) < Tol )
	       continue;
	    
	    if (fabs(A.Elements_[i][i] - A.Elements_[j][j]) > Tol)
	    {
	       theta = 0.5*atan(2.0*A.Elements_[i][j]/(A.Elements_[i][i] - A.Elements_[j][j]));
	    }
	    else
	    {
	       theta = PIby4*(A.Elements_[i][j]/fabs(A.Elements_[i][j]));
	    }

	    if (fabs(theta) > PIby4)
	    {
	       theta -= 2.0*PIby4*(theta/fabs(theta));
	    }

	    c = cos(theta);
	    s = sin(theta);
	    cc = c*c;
	    ss = s*s;
	    cs = c*s;
	    aij1 = A.Elements_[i][j];
	    aii1 = A.Elements_[i][i];
	    ajj1 = A.Elements_[j][j];


	    A.Elements_[i][i] = aii1*cc + 2.0*aij1*cs + ajj1*ss;
	    A.Elements_[j][j] = aii1*ss - 2.0*aij1*cs + ajj1*cc;
	    A.Elements_[i][j] = A.Elements_[j][i] = 0.0;

	    for (int k=0;k<A.Cols_;k++)
	    {
	       if (B != NULL)
	       {
		  tmp = B->Elements_[k][i]*c + B->Elements_[k][j]*s;
		  B->Elements_[k][j] = -B->Elements_[k][i]*s + B->Elements_[k][j]*c;
		  B->Elements_[k][i] = tmp;
	       }

	       
	       if ( k==i || k==j)
	       {
		  continue;
	       }

	       aki1 = A.Elements_[k][i];
	       akj1 = A.Elements_[k][j];

	       A.Elements_[i][k] = A.Elements_[k][i] = aki1*c + akj1*s;
	       A.Elements_[j][k] = A.Elements_[k][j] = -aki1*s + akj1*c;
	    }
	 }
      }
      count++;

      converged = 1;
      for (int i=0;i<A.Cols_;i++)
	 for (int j=0;j<A.Cols_;j++)
	 {
	    if ((i != j) && (fabs(A.Elements_[i][j]) > Tol))
	    {
	       converged = 0;
	    }
	 }
   }

   if (!converged)
   {
      cerr << "Error: SymEigVal(): Failed - No convergence!" << endl;
      exit(-1);
   }

   for (int i=0;i<A.Cols_;i++)
      EigVals.Elements_[0][i] = A.Elements_[i][i];
   
   return EigVals;
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

Matrix SolvePLU(const Matrix& A,const Matrix& B)
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

   Matrix::Elm *Y;
   Y = new Matrix::Elm[B.Rows_];
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

   delete [] Y;
      
   return X;
}

Matrix SolveSVD(const Matrix& A,const Matrix& B,
		const Matrix::Elm MaxCond,const int PrintFlag)
{
   // SVD resizes U,W,V so don't bother now
   Matrix
      U,W,V;
   Matrix x(B.Rows_,B.Cols_);

   SVD(A,U,W,V,MaxCond,PrintFlag);

   int jj,j,i;
   Matrix::Elm s,*tmp;

   // Allocate temp space
   tmp = new Matrix::Elm[A.Cols_];

   // Claculate U.Transpose()*B
   for (j=0;j<A.Cols_;j++)
   {
      s = 0.0;
      // Nonzero result only if W[j][j] is nonzero
      if (W.Elements_[j][j])
      {
	 for (i=0;i<A.Rows_;i++)
	    s += U.Elements_[i][j]*B.Elements_[i][0];
	 // This is the divide by W[j][j];
	 s /= W.Elements_[j][j];
      }
      tmp[j] = s;
   }

   // Matrix multiply by V to get answer
   for (j=0;j<A.Cols_;j++)
   {
      s = 0.0;
      for (jj=0;jj<A.Cols_;jj++)
	 s += V.Elements_[j][jj]*tmp[jj];
      x[j][0] = s;
   }

   // release temp space
   delete [] tmp;

   return x;
}

ostream& operator<<(ostream& out,const Matrix& A)
{
   int W=out.width();
   
   out << endl;

   if (Matrix::MathematicaPrintFlag) out << setw(0) << "{{";
   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 out << setw(W) << A.Elements_[i][j];
	 if ((Matrix::MathematicaPrintFlag) && (j!=(A.Cols_-1)))
	    out << ",";
      }
      
      if (Matrix::MathematicaPrintFlag)
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
