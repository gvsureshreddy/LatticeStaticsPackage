#ifndef __SparseMatrix
#define __SparseMatrix

#ifndef __LinearAlgebraBuildDate
#define __LinearAlgebraBuildDate
char *LinearAlgebraBuildDate();
#endif

#include <iostream>

using namespace std;

class Vector3D;
class Vector;
class Matrix;

class SparseMatrix
{
public:
   typedef double Elm;

//protected:
public:
   int NoNonZero_;
   int Rows_;
   int Cols_;
   int *Row_id_;
   int *Column_id_;
   Elm *Nonzero_entry_;

public:

   // Constructor...
   // Precond. Matrix object has been declared
   // Receive. Rows,Cols,Initial Value
   // Output.  None
   // Postcondition. Matrix of size RowsXCols
   //   allocated and each element set to Initial Value
   // Defaults: Rows=0,Cols=0,Initial Value= (Uninitialized)
   SparseMatrix() {NoNonZero_ = 0;}
   SparseMatrix(Matrix const& A);
   SparseMatrix(Matrix const& A,int const& NoEntries);
   SparseMatrix(SparseMatrix const& A);
   SparseMatrix(int const& NoNonZero,int const& Rows,int const& Cols);
	
	
   // Deconstructor...
   ~SparseMatrix();

   // Size Access...
   int const& Rows() const {return Rows_;}
   int const& Cols() const {return Cols_;}
   int const& NoNonZero() const {return NoNonZero_;}
   
   // Mathematical Operations...
   friend SparseMatrix& operator+(SparseMatrix& A){return A;}
   friend Matrix operator+(SparseMatrix const& A,SparseMatrix const& B);
   friend SparseMatrix operator-(SparseMatrix const& A);
   friend Matrix operator-(SparseMatrix const& A,SparseMatrix const& B);
   friend SparseMatrix operator*(double const& A,SparseMatrix const& B);
   friend SparseMatrix operator*(SparseMatrix const& B,double const& A);
   friend Matrix operator*(SparseMatrix const& A,SparseMatrix const& B);
   friend Matrix operator*(SparseMatrix const& A,Matrix const& B);
   friend Matrix operator*(Matrix const& A,SparseMatrix const& B);
   friend Vector operator*(SparseMatrix const& A,Vector const& B);
   friend Vector operator*(Vector const& A,SparseMatrix const& B);
   friend Vector3D operator*(SparseMatrix const& A,Vector3D const& B);
   friend Vector3D operator*(Vector3D const& A,SparseMatrix const& B);
	

   //Assignment Operators
   SparseMatrix& operator=(Matrix const& A);
   
   
   //Miscellaneous Matrix operations
   SparseMatrix Transpose() const;
   SparseMatrix& SetSparseIdentity(int const& Size);
   friend Matrix ReverseSparse(SparseMatrix const& A);
	
   int IsNull() const {return NoNonZero_==0;}

   friend ostream& operator<<(ostream& out,SparseMatrix const& A);

   static char const* const Revision();
};

#endif
