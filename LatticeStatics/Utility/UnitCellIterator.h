#ifndef __UnitCellIterator
#define __UnitCellIterator

#include <Matrix.h>
#include <Vector.h>

class UnitCellIterator
{
private:
   // Cubic grid of size ((2^(GridSize_)) + 1)
   int GridSize_;
   int VectorsLen_;
   double **Vectors_;
   int CurrentPOS_;
   enum ScanType {Simple,Volume};
   
   void Initialize(int GridSize,ScanType ScnTyp,int DoHalfOnly,int SkipZero);

public:

   UnitCellIterator(): Vectors_(NULL) {}
   UnitCellIterator(int GridSize,ScanType ScnTyp=Simple,int DoHalfOnly=1,int SkipZero=1)
      :Vectors_(NULL) {Initialize(GridSize,ScnTyp,DoHalfOnly,SkipZero);}
   ~UnitCellIterator();

   void operator()(int GridSize,ScanType ScnTyp=Simple,int DoHalfOnly=1,int SkipZero=1)
   {Initialize(GridSize,ScnTyp,DoHalfOnly,SkipZero);}

   void Reset() {CurrentPOS_ = 0;}

   int Done() {return VectorsLen_ <= CurrentPOS_;}
   const double operator[](int i) {return Vectors_[CurrentPOS_][i];}
   void operator++() {++CurrentPOS_;}
};

#endif

