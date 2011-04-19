#ifndef RSE__LJSplineCutoff
#define RSE__LJSplineCutoff

#include "LJ.h"
#include <Vector.h>
#include <Matrix.h>

using namespace std;

class LJSplineCutoff : public LJ
{
private:
   double CutoffStart_;
   double CutoffEnd_;
   Matrix CoeffsMat_;
   mutable Vector CoeffsVec_;
   mutable Vector Coeffs_;

public:
   LJSplineCutoff()
   {
   }
   LJSplineCutoff(double const& Eps0, double const& Eps1, double const& Sigma0,
                  double const& Sigma1, double const& CutoffStart, double const& CutoffEnd);
   ~LJSplineCutoff()
   {
   }
   friend ostream& operator<<(ostream& out, LJSplineCutoff const& A);
   double PairPotential(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                        TDeriv const& dt = T0) const;
   virtual int GetNoParameters() const
   {
      return (2 + LJ::GetNoParameters());
   }
   virtual void SetParameters(double const* const Vals);
   virtual void Print(ostream& out) const;
   virtual char const* const Type() const
   {
      return "LJSplineCutoff";
   }

   double const& CutoffStart() const
   {
      return CutoffStart_;
   }
   double const& CutoffEnd() const
   {
      return CutoffEnd_;
   }

   void SetCutoffStart(double const& CutoffStart)
   {
      CutoffStart_ = CutoffStart;
   }
   void SetCutoffEnd(double const& CutoffEnd)
   {
      CutoffEnd_ = CutoffEnd;
   }
private:
   double CutoffFunction(double const& NTemp, double const& r2, YDeriv const& dy = Y0,
                         TDeriv const& dt = T0) const;
};

#endif

