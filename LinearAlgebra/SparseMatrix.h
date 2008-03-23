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
   unsigned NoNonZero_;
   unsigned Rows_;
   unsigned Cols_;
   unsigned *Row_id_;
   unsigned *Column_id_;
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
   SparseMatrix(const Matrix& A);
   SparseMatrix(const Matrix& A,unsigned NoEntries);
   SparseMatrix(const SparseMatrix& A);
   SparseMatrix(unsigned NoNonZero, unsigned Rows, unsigned Cols);
	
	
   // Deconstructor...
   ~SparseMatrix();

   // Size Access...
   unsigned Rows() const {return Rows_;}
   unsigned Cols() const {return Cols_;}
   unsigned NoNonZero() const {return NoNonZero_;}
   
   // Mathematical Operations...
   friend SparseMatrix& operator+(SparseMatrix& A){return A;}
   friend Matrix operator+(const SparseMatrix& A,const SparseMatrix&B);
   friend SparseMatrix operator-(const SparseMatrix& A);
   friend Matrix operator-(const SparseMatrix& A,const SparseMatrix&B);
   friend SparseMatrix operator*(const double A, const SparseMatrix& B);
   friend SparseMatrix operator*(const SparseMatrix& B, const double A);
   friend Matrix operator*(const SparseMatrix& A, const SparseMatrix& B);
   friend Matrix operator*(const SparseMatrix& A, const Matrix& B);
   friend Matrix operator*(const Matrix& A, const SparseMatrix& B);
   friend Vector operator*(const SparseMatrix& A, const Vector& B);
   friend Vector operator*(const Vector& A, const SparseMatrix& B);
   friend Vector3D operator*(const SparseMatrix& A,const Vector3D& B);
   friend Vector3D operator*(const Vector3D& A,const SparseMatrix& B);
	

   //Assignment Operators
   SparseMatrix& operator=(const Matrix& A);
   
   
   //Miscellaneous Matrix operations
   SparseMatrix Transpose() const;
   friend SparseMatrix SparseIdentity(unsigned Size);
   friend Matrix ReverseSparse(const SparseMatrix& A);
	
   int IsNull() const {return NoNonZero_==0;}

   friend ostream& operator<<(ostream& out,const SparseMatrix& A);

   static char* Revision();
};

#endif
