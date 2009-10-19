#include "CMatrix.h"
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>

// Global IDString
char CMatrixID[]="$Id: CMatrix.cpp,v 1.22 2009/10/19 18:47:19 elliott Exp $";

// Private Methods...

// Returns matrix of size Rows_-1 x Cols_-1 with ith row and
//   jth column removed
CMatrix CMatrix::Minor(int const& i,int const& j) const
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
		  A[a][b]=Elements_[a][b];
	       }
	       else
	       {
		  A[a][b]=Elements_[a][b+1];
	       }
	    }
	    else
	    {
	       if (b < j)
	       {
		  A[a][b]=Elements_[a+1][b];
	       }
	       else
	       {
		  A[a][b]=Elements_[a+1][b+1];
	       }
	    }
	 }
      }
   }
   
   return A;
}

// Public Methods...

int CMatrix::MathematicaPrintFlag = 0;

CMatrix::CMatrix(int const& Rows,int const& Cols):
   Rows_(Rows), Cols_(Cols)
{
   if (IsNull())
   {
      Elements_=0;
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
   return;
}

CMatrix::CMatrix(int const& Rows,int const& Cols,CMatrix::Elm const& InitVal):
   Rows_(Rows), Cols_(Cols)
{
   if (IsNull())
   {
      Elements_=0;
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
      {
         for (register int j=0;j<Cols_;j++)
            Elements_[i][j]=InitVal;
      }
   }
   return;
}

CMatrix::CMatrix(CMatrix const& A):
   Rows_(A.Rows_), Cols_(A.Cols_)
{
   if (IsNull())
   {
      Elements_=0;
   }
   else
   {
      Elements_=new CMatrix::Elm*[Rows_];
      Elements_[0]=new CMatrix::Elm[Rows_*Cols_];

      for (register int i=1;i<Rows_;i++)
      {
	 Elements_[i]=Elements_[i-1]+Cols_;
      }

      memmove(Elements_[0],A[0],sizeof(CMatrix::Elm)*Rows_*Cols_);
   }

   return;
}

CMatrix::CMatrix(Matrix const& A):
   Rows_(A.Rows()), Cols_(A.Cols())
{
   if (IsNull())
   {
      Elements_=0;
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

CMatrix operator+(CMatrix const& A,CMatrix const& B)
{
   if (A.Rows_!=B.Rows_ || A.Cols_!=B.Cols_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix Operator+() Diff Size Matrices or Null CMatrix!!!"
	   << "\n";
      exit(-1);
   }

   CMatrix C(A.Rows_,A.Cols_);
   
   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C[i][j]=A[i][j]+B[i][j];
      }
   }
   
   return C;
}

CMatrix operator-(CMatrix const& A)
{
   CMatrix B(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for(int j=0;j<A.Cols_;j++)
      {
	 B[i][j]=-A[i][j];
      }
   }

   return B;
}

CMatrix operator-(CMatrix const& A,CMatrix const& B)
{
   if (A.Rows_!=B.Rows_ || A.Cols_!=B.Cols_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix Operator-() Diff Size Matrices or Null CMatrix!!!"
	   << "\n";
      exit(-1);
   }
      
   CMatrix C(A.Rows_,A.Cols_);
   
   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C[i][j]=A[i][j]-B[i][j];
      }
   }
   
   return C;
}

CMatrix operator*(CMatrix const& A,CMatrix const& B)
{
   if (A.Cols_!=B.Rows_ || A.IsNull() || B.IsNull())
   {
      cerr << "Error In CMatrix Operator* : A.Cols!=B.Rows or Null CMatrix"
	   <<"\n";
      exit(-1);
   }

   CMatrix C(A.Rows_,B.Cols_,0.0);
   
   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<B.Cols_;j++)
      {
	 for (register int k=0;k<A.Cols_;k++)
	 {
	    C[i][j]+=A[i][k]*B[k][j];
	 }
      }
   }
   
   return C;
}

CMatrix operator*(CMatrix::Elm const& A,CMatrix const& B)
{
   CMatrix C(B.Rows_,B.Cols_);

   for (register int i=0;i<B.Rows_;i++)
   {
      for (register int j=0;j<B.Cols_;j++)
      {
	 C[i][j]=A*B[i][j];
      }
   }

   return C;
}

CMatrix operator*(CMatrix const& A,CMatrix::Elm const& B)
{
   CMatrix C(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C[i][j]=B*A[i][j];
      }
   }

   return C;
}

CMatrix operator/(CMatrix const& A,CMatrix::Elm const& B)
{
   if (B == CMatrix::Elm(0.0,0.0))
   {
      cerr << "Divide By Zero Error in CMatrix operator/()"
	   << "\n";
      exit(-1);
   }
   
   CMatrix C(A.Rows_,A.Cols_);

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 C[i][j]=A[i][j]/B;
      }
   }

   return C;
}

#ifdef CHECK_BOUNDS
CMatrix::Elm* const CMatrix::operator[](int const& i)
{
   if (i >= Rows_)
   {
      cerr << "CMatrix Index Overflow -- CMatrix::Elm* operator[]()" << "\n";
      exit(-1);
   }

   return Elements_[i];
}

CMatrix::Elm const* const CMatrix::operator[](int const& i) const
{
   if (i >= Rows_)
   {
      cerr << "CMatrix Index Overflow -- CMatrix::Elm* operator[]()" << "\n";
      exit(-1);
   }

   return Elements_[i];
}
#endif
   
CMatrix& CMatrix::operator=(CMatrix const& B)
{
   if (Rows_!=B.Rows_ || Cols_!=B.Cols_ || IsNull() || B.IsNull())
   {
      cerr << "Error in CMatrix& operator=() : Matricies not same size "
           << "or Null CMatrix"
           << "\n";
      exit(-1);
   }

   memmove(Elements_[0],B[0],sizeof(CMatrix::Elm)*Rows_*Cols_);

   return *this;
}

CMatrix& CMatrix::SetIdentity(int const& Size)
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
	 A[j][i]=Elements_[i][j];
      }
   }

   return A;
}

CMatrix CMatrix::Conjugate() const
{
   CMatrix A(Rows_,Cols_);

   for (register int i=0;i<Rows_;++i)
   {
      for (register int j=0;j<Cols_;++j)
      {
	 A[i][j] = Elements_[i][j].conj();
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

   B[0][0]=1.0;
   for (register int i=0;i<Cols_;i++)
   {
      X=SolvePLU(*this,B);

      for (register int j=0;j<Rows_;j++)
	 C[j][i]=X[j][0];

      B[i][0]=0;
      if (i!=Cols_-1) B[i+1][0]=1.0;
   }

   return C;
}

void CMatrix::Resize(int const& Rows,int const& Cols)
{
   if (Rows!=Rows_ || Cols!=Cols_)
   {
      if (Elements_!=0)
      {
	 delete [] Elements_[0];
	 delete [] Elements_;
      }

      Rows_=Rows;
      Cols_=Cols;
   
      if (IsNull())
      {
	 Elements_=0;
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
   }

   return;
}

void CMatrix::Resize(int const& Rows,int const& Cols,CMatrix::Elm const& InitVal)
{
   if (Rows!=Rows_ || Cols!=Cols_)
   {
      if (Elements_!=0)
      {
	 delete [] Elements_[0];
	 delete [] Elements_;
      }

      Rows_=Rows;
      Cols_=Cols;
   
      if (IsNull())
      {
	 Elements_=0;
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
   }

   for (register int i=0;i<Rows_;i++)
   {
      for (register int j=0;j<Cols_;j++)
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

      for (register int i=0;i<Cols_;i++)
      {
	 det += CMatrix::Elm(1-2*(i%2))*Elements_[0][i]*(Minor(0,i).Det());
      }

      return det;
   }
}

// Decompose PA=LU using scaled partial pivoting.
void PLU(CMatrix const& A,CMatrix& P,CMatrix& L,CMatrix& U)
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

   for (register int i=0;i<A.Rows_;i++)
   {
      Ipivot[i]=i;
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 if (abs(Temp[i][j]) > abs(S[i][0]))
	    S[i][0]=Temp[i][j];
      }
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      CMatrix::Elm temp1;
      temp1=Temp[i][i]/S[i][0];

      int k=i;
      for (register int j=i;j<A.Rows_;j++)
      {
	 if (abs(Temp[j][i]) > abs(temp1))
	 {
	    temp1=Temp[j][i]/S[j][0];
	    k=j;
	 }
      }

      if (k>i)
      {
	 CMatrix::Elm *Switch;
	 Switch = new CMatrix::Elm[A.Rows_];
	 for (register int j=i;j<A.Rows_;j++)
	 {
	    Switch[j]=Temp[i][j];
	    Temp[i][j]=Temp[k][j];
	    Temp[k][j]=Switch[j];
	 }

	 for (register int j=0;j<i;j++)
	 {
	    Switch[j]=L[i][j];
	    L[i][j]=L[k][j];
	    L[k][j]=Switch[j];
	 }

	 delete [] Switch;

	 temp1=S[i][0];
	 S[i][0]=S[k][0];
	 S[k][0]=temp1;

	 int tempi1=Ipivot[i];
	 Ipivot[i]=Ipivot[k];
	 Ipivot[k]=tempi1;
      }

      for (register int j=i+1;j<A.Rows_;j++)
      {
	 L[j][i]=Temp[j][i]/Temp[i][i];
      }

      for (register int j=i+1;j<A.Rows_;j++)
      {
	 for (register int k=i+1;k<A.Rows_;k++)
	 {
	    Temp[j][k]=Temp[j][k]-(L[j][i]*Temp[i][k]);
	 }
      }
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=i;j<A.Rows_;j++)
      {
	 U[i][j]=Temp[i][j];
      }
   }

   for (register int i=0;i<A.Rows_;i++)
   {
      P[i][Ipivot[i]]=1;
   }

   delete [] Ipivot;

   return;
}

Matrix HermiteEigVal(CMatrix A,CMatrix* const B,int const& MaxItr,double const& Tol)
{
   int count=0,
      converged=0;
   Matrix EigVals(1,A.Cols_);
   MyComplexDouble tau,t1,t2,t,c,s,cc,ss,cs,ssbar,csbar,
      aij1,aii1,ajj1,aki1,akj1,tmp;

   if (B != 0)
   {
      B->SetIdentity(A.Cols_);
   }
   
   while ((count < MaxItr) && (!converged))
   {
      for (int i=0;i<A.Cols_;i++)
      {
	 for (int j=i+1;j<A.Cols_;j++)
	 {
	    if (abs(A[i][j]) < Tol )
	       continue;
	    else
	    {
	       tau = (A[i][i] - A[j][j])
		  /(2.0*(A[i][j].conj()));

	       // take care to make the sqrt well conditioned! (not overflow)
	       if (abs(tau) > 1.0)
	       {
		  t1 = tau*(-1.0 - sqrt( 1.0 + (A[i][j]
						/((A[i][j].conj())*tau*tau))));
		  t2 = tau*(-1.0 + sqrt( 1.0 + (A[i][j]
						/((A[i][j].conj())*tau*tau))));
	       }
	       else
	       {
		  t1 = -tau - sqrt( A[i][j]/(A[i][j].conj())
				    + tau*tau);
		  t2 = -tau + sqrt( A[i][j]/(A[i][j].conj())
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
	       aij1 = A[i][j];
	       aii1 = A[i][i];
	       ajj1 = A[j][j];
	       
	       A[i][i] = (aii1*cc + cs*(aij1.conj())
                          + csbar*aij1 + ssbar*ajj1).real();
	       A[j][j] = (aii1*ssbar - cs*(aij1.conj())
                          - csbar*aij1 + cc*ajj1).real();
	       A[i][j] = A[j][i] = 0.0;
	       
	       for (int k=0;k<A.Cols_;k++)
	       {
		  if (B != 0)
		  {
		     tmp = (*B)[k][i]*c + (*B)[k][j]*(s.conj());
		     (*B)[k][j] = -(*B)[k][i]*s + (*B)[k][j]*c;
		     (*B)[k][i] = tmp;
		  }
		  
		  
		  if ( k==i || k==j)
		  {
		     continue;
		  }
		     
		  aki1 = A[k][i];
		  akj1 = A[k][j];
		  
		  A[k][i] = aki1*c + akj1*(s.conj());
		  A[i][k] = (A[k][i].conj());
		  A[k][j] = -aki1*s + akj1*c;
		  A[j][k] = (A[k][j].conj());
	       }
	    }
	 }
      }
      count++;
      
      converged = 1;
      for (int i=0;i<A.Cols_;i++)
	 for (int j=i+1;j<A.Cols_;j++)
	 {
	    if (abs(A[i][j]) > Tol)
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

   for (int i=0;i<A.Cols_;i++)
   {
      EigVals[0][i] = (A[i][i]).real();
   }
   
   return EigVals;
}

// find QR factorization of A or A.Transpose()
void QR(CMatrix const& A,CMatrix& Q,CMatrix& R,int const& CalcTranspose)
{
   int i,j,k,m,n;
   CMatrix::Elm c,s,r,A1,A2,tau;

   if (CalcTranspose)
   {
      m=A.Cols_;
      n=A.Rows_;
   }
   else
   {
      m=A.Rows_;
      n=A.Cols_;
   }

   //initialize R and Q
   R.Resize(m,n);
   for (i=0;i<m;++i)
   {
      for (j=0;j<n;++j)
      {
	 if (CalcTranspose)
	    R[i][j] = A[j][i];
	 else
	    R[i][j] = A[i][j];
      }
   }
   Q.SetIdentity(m);

   // Perform Givens rotations G
   for (j=0;j<m;++j)
   {
      for (i=m-1;i>j;--i)
      {
         // calculate G such that G^T*R sets R[i][j] = 0
         if (abs(R[i][j]) == 0.0)
         {
            c=1.0; s=0.0;
         }
         else if (abs(R[i-1][j]) == 0.0)
         {
            c=0.0; s=1.0;
         }
         else
         {
            if (abs(R[i][j]) > abs(R[i-1][j]))
            {
               tau = -R[i-1][j]/R[i][j];
               c = sqrt(tau*tau.conj()/(1.0+tau*tau.conj()));
               s = c/tau;
            }
            else
            {
               tau = -R[i][j]/R[i-1][j];
               c = 1.0/sqrt(1.0+tau*tau.conj());
               s = tau*c;
            }
         }

         for (k=j;k<n;++k)
         {
            // perform G^T*R
            A1 = R[i-1][k];
            A2 = R[i][k];
            R[i-1][k] = A1*c - A2*s.conj();
            R[i][k] = A1*s + A2*c;
         }

         for (k=0;k<m;++k)
         {
            // perform Q*G
            A1 = Q[k][i-1];
            A2 = Q[k][i];
            Q[k][i-1] = A1*c - A2*s;
            Q[k][i] = A1*s.conj() + A2*c;
         }
      }
   }
}

CMatrix RightEigVals(CMatrix const& A,int const& MaxItr,double const& Tol)
{
   if (!A.IsSquare())
   {
      cerr << "Error: RightEigVals --- A is not square." << "\n";
      exit(-51);
   }
   
   CMatrix T(A);
   CMatrix Q,R;
   int n=A.Rows();
   int converged=0;
   int iterations = 0;

   while ((!converged) && (iterations < MaxItr))
   {
      iterations++;
      QR(T,Q,R);
      T = R*Q;
      converged = 1;
      for (int i=n-1;i>0;--i)
      {
         for (int j=0;j<i;++j)
         {
            if (abs(T[i][j]) > Tol)
            {
               converged = 0;
            }
         }
      }
   }

   CMatrix REVs(1,n);
   for (int i=0;i<n;++i)
      REVs[0][i] = T[i][i];

   if (iterations >= MaxItr)
   {
      cerr << "Error: RightEigVals did not converge!" << "\n";
   }
   return REVs;
}

void Cholesky(CMatrix const& A,CMatrix& U,CMatrix& D)
{
   if (!A.IsSquare() || A.IsNull())
   {
      cerr << "Error in Cholesky() -- Non-Square or Null CMatrix" << "\n";
      exit(-1);
   }

   U.SetIdentity(A.Rows_);
   D.Resize(A.Rows_,A.Cols_,0);
   

   for (register int i=0;i<A.Rows_;i++)
   {
      D[i][i]=A[i][i];
      if (i > 0)
      {
	 for (register int k=0;k<i;k++)
	 {

	    D[i][i]-=
	       D[k][k]*U[k][i]*(U[k][i].conj());
	 }
      }
      for (register int j=i+1;j<A.Rows_;j++)
      {
	 U[i][j]=A[i][j];
	 if (i > 0)
	 {
	    for (register int k=0;k<i;k++)
	    {
	       U[i][j]-=
		  D[k][k]*(U[k][i].conj())*U[k][j];
	    }
	 }
	 U[i][j]/=D[i][i];
      }
   }
   
   return;
}

CMatrix SolvePLU(CMatrix const& A,CMatrix const& B)
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

   Y[0]=Temp[0][0];
   for (register int i=1;i<B.Rows_;i++)
   {
      Y[i]=Temp[i][0];
      for (register int j=0;j<i;j++)
      {
	 Y[i]-=Y[j]*L[i][j];
      }
   }

   CMatrix X(B.Rows_,1,0);

   X[B.Rows_-1][0]=Y[B.Rows_-1]/U[A.Rows_-1][A.Cols_-1];
   for (register int i=A.Rows_-2;i>=0;i--)
   {
      X[i][0]=Y[i];
      for (register int j=A.Rows_-1;j>i;j--)
      {
	 X[i][0]-=X[j][0]*U[i][j];
      }
      X[i][0]=X[i][0]/U[i][i];
   }

   delete [] Y;
      
   return X;
}

ostream& operator<<(ostream& out,CMatrix const& A)
{
   int W=out.width();
   
   out << "\n";

   if (CMatrix::MathematicaPrintFlag) out << setw(0) << "{{";
   for (register int i=0;i<A.Rows_;i++)
   {
      for (register int j=0;j<A.Cols_;j++)
      {
	 out << setw(W) << A[i][j];
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
   for (register int i=0;i<A.Rows_;i++)
      for (register int j=0;j<A.Cols_;j++)
	 in >> A[i][j];

   return in;
}

char* CMatrix::Revision()
{
   return CMatrixID;
}
