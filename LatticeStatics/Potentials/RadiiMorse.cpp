#include "RadiiMorse.h"
#include <cstdlib>

RadiiMorse::RadiiMorse(double const& A0,double const& AT,double const& B0,double const& BT,
                       double const& Rref1,double const& Rtheta1,double const& Rref2,
                       double const& Rtheta2):
   A0_(A0),AT_(AT),B0_(B0),BT_(BT),Rref1_(Rref1),Rtheta1_(Rtheta1),Rref2_(Rref2),
   Rtheta2_(Rtheta2)
{
}

void RadiiMorse::SetParameters(double const* const Vals)
{
   SetA0(Vals[0]);
   SetAT(Vals[1]);
   SetB0(Vals[2]);
   SetBT(Vals[3]);
   SetRref1(Vals[4]);
   SetRtheta1(Vals[5]);
   SetRref2(Vals[6]);
   SetRtheta2(Vals[7]);
}

double RadiiMorse::A(double const& NTemp,TDeriv const& dt) const
{
   double retval;
   switch (dt)
   {
      case T0:
         retval = A0_ + AT_*(NTemp-1.0);
         break;
      case DT:
         retval = AT_;
         break;
      case D2T:
         retval = 0.0;
         break;
      default:
         cerr << "Error in RadiiMorse::A" << "\n";
         exit(-1);
   }
   
   return retval;
}

double RadiiMorse::Beta(double const& NTemp,TDeriv const& dt) const
{
   double retval;
   
   switch (dt)
   {
      case T0:
         retval = B0_ + BT_*(NTemp-1.0);
         break;
      case DT:
         retval = BT_;
         break;
      case D2T:
         retval = 0.0;
         break;
      default:
         cerr << "Error in RadiiMorse::Beta" << "\n";
         exit(-1);
   }
   
   return retval;
}

double RadiiMorse::Rhat(double const& NTemp,TDeriv const& dt) const
{
   double rhat;
   
   switch (dt)
   {
      case T0:
         rhat = (Rref1_ + Rtheta1_*(NTemp - 1.0))
            *(Rref2_ + Rtheta2_*(NTemp - 1.0));
         break;
      case DT:
         rhat = Rtheta1_*(Rref2_ + Rtheta2_*(NTemp - 1.0))
            +(Rref1_ + Rtheta1_*(NTemp - 1.0))*Rtheta2_;
         break;
      case D2T:
         rhat = 2.0*Rtheta1_*Rtheta2_;
         break;
      default:
         cerr << "Error in RadiiMorse::Rhat" << "\n";
         exit(-1);
   }
   
   return rhat;
}

double RadiiMorse::PairPotential(double const& NTemp,double const& r2,YDeriv const& dy,
                                 TDeriv const& dt) const
{
   for (int i=0;i<DTmax;++i)
   {
      Achk_[i]=
         Bchk_[i]=
         Rchk_[i]=
         Ichk_[i]=0;
      for (int j=0;j<DYmax;++j)
         Gchk_[j][i] = 0;
   }
   
   double val=0.0;
   switch (dy)
   {
      case Y0:
         switch (dt)
         {
            case T0:
               val = a(NTemp)*g(NTemp,r2,Y0,T0)*(g(NTemp,r2,Y0,T0)-2.0);
               break;
            case DT:
               val = a(NTemp,DT)*g(NTemp,r2,Y0,T0)*(g(NTemp,r2,Y0,T0)-2.0)
                  + a(NTemp)*g(NTemp,r2,Y0,DT)*(g(NTemp,r2,Y0,T0)-2.0)
                  + a(NTemp)*g(NTemp,r2,Y0,T0)*g(NTemp,r2,Y0,DT);
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,Y0,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,Y0,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,Y0,DT)*(g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,Y0,T0)*g(NTemp,r2,Y0,DT)
                      + a(NTemp,D2T)*g(NTemp,r2,Y0,T0)*(g(NTemp,r2,Y0,T0)-2.0));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- Y0,D3T" << "\n";
               exit(-1);
               break;
         }
         break;
      case DY:
         switch (dt)
         {
            case T0:
               val = a(NTemp)*g(NTemp,r2,DY,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0);
               break;
            case DT:
               val = (a(NTemp)*g(NTemp,r2,DY,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,DY,T0)*g(NTemp,r2,Y0,DT)
                      + a(NTemp,DT)*g(NTemp,r2,DY,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0));
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,DY,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,DY,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,DY,T0)*g(NTemp,r2,Y0,D2T)
                      + a(NTemp,D2T)*g(NTemp,r2,DY,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,DY,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,DY,T0)*g(NTemp,r2,Y0,DT));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- DY,D3T" << "\n";
               exit(-1);
               break;
         }
         break;
      case D2Y:
         switch (dt)
         {
            case T0:
               val = (a(NTemp)*g(NTemp,r2,D2Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,DY,T0)*g(NTemp,r2,DY,T0));
               break;
            case DT:
               val = (a(NTemp)*g(NTemp,r2,D2Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,Y0,DT)
                      + 4.0*a(NTemp)*g(NTemp,r2,DY,DT)*g(NTemp,r2,DY,T0)
                      + a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,DY,T0)*g(NTemp,r2,DY,T0));
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,D2Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 4.0*a(NTemp)*g(NTemp,r2,DY,D2T)*g(NTemp,r2,DY,T0)
                      + 4.0*a(NTemp)*g(NTemp,r2,DY,DT)*g(NTemp,r2,DY,DT)
                      + a(NTemp,D2T)*g(NTemp,r2,D2Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp,D2T)*g(NTemp,r2,DY,T0)*g(NTemp,r2,DY,T0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D2Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,Y0,DT)
                      + 8.0*a(NTemp,DT)*g(NTemp,r2,DY,DT)*g(NTemp,r2,DY,T0));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- D2Y,D3T" << "\n";
               exit(-1);
               break;
         }
         break;
      case D3Y:
         switch (dt)
         {
            case T0:
               val = (4.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp)*g(NTemp,r2,D3Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,DY,T0)*g(NTemp,r2,D2Y,T0));
               break;
            case DT:
               val = (a(NTemp)*g(NTemp,r2,D3Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 6.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,DY,T0)
                      + 6.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,Y0,DT)
                      + a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 6.0*a(NTemp,DT)*g(NTemp,r2,DY,T0)*g(NTemp,r2,D2Y,T0));
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,D3Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 6.0*a(NTemp)*g(NTemp,r2,D2Y,D2T)*g(NTemp,r2,DY,T0)
                      + 12.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,DY,DT)
                      + 6.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,D2T)
                      + a(NTemp,D2T)*g(NTemp,r2,D3Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 6.0*a(NTemp,D2T)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,T0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D3Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 12.0*a(NTemp,DT)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,DY,T0)
                      + 12.0*a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,DT)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,Y0,DT));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- D3Y,D3T" << "\n";
               exit(-1);
               break;
         }
         break;
      case D4Y:
         switch (dt)
         {
            case T0:
               val = (6.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp)*g(NTemp,r2,D4Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0));
               break;
            case DT:
               val = (a(NTemp)*g(NTemp,r2,D4Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,Y0,DT)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,DY,T0)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,DT)
                      + 12.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,D2Y,T0)
                      + a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 8.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,T0)
                      + 6.0*a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,D2Y,T0));
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,D4Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,D2T)*g(NTemp,r2,DY,T0)
                      + 16.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,DY,DT)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,D2T)
                      + 12.0*a(NTemp)*g(NTemp,r2,D2Y,D2T)*g(NTemp,r2,D2Y,T0)
                      + 12.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,D2Y,DT)
                      + a(NTemp,D2T)*g(NTemp,r2,D4Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 8.0*a(NTemp,D2T)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,T0)
                      + 6.0*a(NTemp,D2T)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D4Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,Y0,DT)
                      + 16.0*a(NTemp,DT)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,DY,T0)
                      + 16.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,DT)
                      + 24.0*a(NTemp,DT)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,D2Y,T0));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- D4Y,D3T" << "\n";
               exit(-1);
               break;
         }
         break;
      case D5Y:
         switch (dt)
         {
            case T0:
               val = (20.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 10.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp)*g(NTemp,r2,D5Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0));
               break;
            case DT:
               val = (a(NTemp)*g(NTemp,r2,D5Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,D5Y,DT)*g(NTemp,r2,Y0,T0)
                      + 10.0*a(NTemp)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,DY,T0)
                      + 10.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,DY,DT)
                      + 20.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,D2Y,T0)
                      + 20.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D2Y,DT)
                      + a(NTemp,DT)*g(NTemp,r2,D5Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 10.0*a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,DY,T0)
                      + 20.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D2Y,T0));
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,D5Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D5Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 10.0*a(NTemp)*g(NTemp,r2,D4Y,D2T)*g(NTemp,r2,DY,T0)
                      + 20.0*a(NTemp)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,DY,DT)
                      + 10.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,DY,D2T)
                      + 20.0*a(NTemp)*g(NTemp,r2,D3Y,D2T)*g(NTemp,r2,D2Y,T0)
                      + 40.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,D2Y,DT)
                      + 20.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D2Y,D2T)
                      + a(NTemp,D2T)*g(NTemp,r2,D5Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 10.0*a(NTemp,D2T)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,DY,T0)
                      + 20.0*a(NTemp,D2T)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D5Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,Y0,DT)
                      + 20.0*a(NTemp,DT)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,DY,T0)
                      + 20.0*a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,DY,DT)
                      + 40.0*a(NTemp,DT)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,D2Y,T0)
                      + 40.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D2Y,DT));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- D5Y,D2T" << "\n";
               exit(-1);
               break;
         }
         break;
      case D6Y:
         switch (dt)
         {
            case T0:
               val = (20.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D3Y,T0)
                      + 30.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 12.0*a(NTemp)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp)*g(NTemp,r2,D6Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0));
               break;
            case DT:
               val = (a(NTemp)*g(NTemp,r2,D6Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,D6Y,T0)*g(NTemp,r2,Y0,DT)
                      + 12.0*a(NTemp)*g(NTemp,r2,D5Y,DT)*g(NTemp,r2,DY,T0)
                      + 12.0*a(NTemp)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,DY,DT)
                      + 30.0*a(NTemp)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,D2Y,T0)
                      + 30.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,D2Y,DT)
                      + 40.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,D3Y,T0)
                      + a(NTemp,DT)*g(NTemp,r2,D6Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 12.0*a(NTemp,DT)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,DY,T0)
                      + 30.0*a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 20.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D3Y,T0));
               break;
            case D2T:
               val = (a(NTemp)*g(NTemp,r2,D6Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D6Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D6Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 12.0*a(NTemp)*g(NTemp,r2,D5Y,D2T)*g(NTemp,r2,DY,T0)
                      + 24.0*a(NTemp)*g(NTemp,r2,D5Y,DT)*g(NTemp,r2,DY,DT)
                      + 12.0*a(NTemp)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,DY,D2T)
                      + 30.0*a(NTemp)*g(NTemp,r2,D4Y,D2T)*g(NTemp,r2,D2Y,T0)
                      + 60.0*a(NTemp)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,D2Y,DT)
                      + 30.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,D2Y,D2T)
                      + 40.0*a(NTemp)*g(NTemp,r2,D3Y,D2T)*g(NTemp,r2,D3Y,T0)
                      + 40.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,D3Y,DT)
                      + a(NTemp,D2T)*g(NTemp,r2,D6Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 12.0*a(NTemp,D2T)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,DY,T0)
                      + 30.0*a(NTemp,D2T)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 20.0*a(NTemp,D2T)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,D3Y,T0)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D6Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,D6Y,T0)*g(NTemp,r2,Y0,DT)
                      + 24.0*a(NTemp,DT)*g(NTemp,r2,D5Y,DT)*g(NTemp,r2,DY,T0)
                      + 24.0*a(NTemp,DT)*g(NTemp,r2,D5Y,T0)*g(NTemp,r2,D5Y,DT)
                      + 60.0*a(NTemp,DT)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,D2Y,T0)
                      + 60.0*a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,D2Y,DT)
                      + 80.0*a(NTemp,DT)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,D3Y,T0));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- D6Y,DT" << "\n";
               exit(-1);
               break;
         }
         break;
      default:
         cerr << "Error in RadiiMorse::PairPotential -- D5Y" << "\n";
         exit(-1);
         break;
   }
   
   return val;
}

double RadiiMorse::I(double const& NTemp,TDeriv const& dt) const
{
   double val;
   
   switch (dt)
   {
      case T0:
         val = - b(NTemp)/(2.0*rhat(NTemp));
         break;
      case DT:
         val = -(rhat(NTemp)*b(NTemp,DT) - rhat(NTemp,DT)*b(NTemp))
            /(2.0*rhat(NTemp)*rhat(NTemp));
         break;
      case D2T:
         val = -((b(NTemp,D2T)/(2.0*rhat(NTemp)))
                 - (b(NTemp,DT)*rhat(NTemp,DT))/(rhat(NTemp)*rhat(NTemp))
                 - (b(NTemp)*rhat(NTemp,D2T))/(2.0*rhat(NTemp)*rhat(NTemp))
                 + ((b(NTemp)*rhat(NTemp,DT)*rhat(NTemp,DT))
                    /(rhat(NTemp)*rhat(NTemp)*rhat(NTemp)))
                );
         break;
      default:
         cerr << "Error in RadiiMorse::I -- D3T" << "\n";
         exit(-1);
         break;
   }
   
   return val;
}

double RadiiMorse::G(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt)
   const
{
   double val;
   double r=sqrt(r2);
   
   switch (dy)
   {
      case Y0:
         switch (dt)
         {
            case T0:
               val = exp(2.0*i(NTemp)*r + b(NTemp));
               break;
            case DT:
               val = (2.0*i(NTemp,DT)*r + b(NTemp,DT))*g(NTemp,r2,Y0,T0);
               break;
            case D2T:
               val = (2.0*i(NTemp,DT)*r + b(NTemp,DT))*g(NTemp,r2,Y0,DT)
                  + (2.0*i(NTemp,D2T)*r + b(NTemp,D2T))*g(NTemp,r2,Y0,T0);
               break;
            default:
               cerr << "Error in RadiiMorse::G -- Y0,D3T" << "\n";
               exit(-1);
         }
         break;
      case DY:
         switch (dt)
         {
            case T0:
               val = i(NTemp)*(g(NTemp,r2,Y0,T0)/r);
               break;
            case DT:
               val = (i(NTemp,DT)*g(NTemp,r2,Y0,T0) + i(NTemp)*g(NTemp,r2,Y0,DT))/r;
               break;
            case D2T:
               val = (i(NTemp,D2T)*g(NTemp,r2,Y0,T0) + 2.0*i(NTemp,DT)*g(NTemp,r2,Y0,DT)
                      + i(NTemp)*g(NTemp,r2,Y0,D2T))/r;
               break;
            default:
               cerr << "Error in RadiiMorse::G -- DY,D3T" << "\n";
               exit(-1);
         }
         break;
      case D2Y:
         switch (dt)
         {
            case T0:
               val = i(NTemp)*(g(NTemp,r2,DY,T0)/r - g(NTemp,r2,Y0,T0)/(2.0*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(g(NTemp,r2,DY,T0)/r - g(NTemp,r2,Y0,T0)/(2.0*r2*r))
                      + i(NTemp)*(g(NTemp,r2,DY,DT)/r - g(NTemp,r2,Y0,DT)/(2.0*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(g(NTemp,r2,DY,T0)/r - g(NTemp,r2,Y0,T0)/(2.0*r2*r))
                      + 2.0*i(NTemp,DT)*(g(NTemp,r2,DY,DT)/r - g(NTemp,r2,Y0,DT)/(2.0*r2*r))
                      + i(NTemp)*(g(NTemp,r2,DY,D2T)/r - g(NTemp,r2,Y0,D2T)/(2.0*r2*r)));
               break;
            default:
               cerr << "Error in RadiiMorse::G -- D2Y,D3T" << "\n";
               exit(-1);
         }
         break;
      case D3Y:
         switch (dt)
         {
            case T0:
               val = i(NTemp)*(g(NTemp,r2,D2Y,T0)/r - g(NTemp,r2,DY,T0)/(r2*r)
                               + 3.0*g(NTemp,r2,Y0,T0)/(4.0*r2*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(g(NTemp,r2,D2Y,T0)/r - g(NTemp,r2,DY,T0)/(r2*r)
                                   + 3.0*g(NTemp,r2,Y0,T0)/(4.0*r2*r2*r))
                      + i(NTemp)*(g(NTemp,r2,D2Y,DT)/r - g(NTemp,r2,DY,DT)/(r2*r)
                                  + 3.0*g(NTemp,r2,Y0,DT)/(4.0*r2*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(g(NTemp,r2,D2Y,T0)/r - g(NTemp,r2,DY,T0)/(r2*r)
                                    + 3.0*g(NTemp,r2,Y0,T0)/(4.0*r2*r2*r))
                      + 2.0*i(NTemp,DT)*(g(NTemp,r2,D2Y,DT)/r - g(NTemp,r2,DY,DT)/(r2*r)
                                         + 3.0*g(NTemp,r2,Y0,DT)/(4.0*r2*r2*r))
                      + i(NTemp)*(g(NTemp,r2,D2Y,D2T)/r - g(NTemp,r2,DY,D2T)/(r2*r)
                                  + 3.0*g(NTemp,r2,Y0,D2T)/(4.0*r2*r2*r)));
               break;
            default:
               cerr << "Error in RadiiMorse::G -- D3Y,D3T" << "\n";
               exit(-1);
         }
         break;
      case D4Y:
         switch (dt)
         {
            case T0:
               val = i(NTemp)*(g(NTemp,r2,D3Y,T0)/r - 3.0*g(NTemp,r2,D2Y,T0)/(2.0*r2*r)
                               + 9.0*g(NTemp,r2,DY,T0)/(4.0*r2*r2*r)
                               - 15.0*g(NTemp,r2,Y0,T0)/(8.0*r2*r2*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(g(NTemp,r2,D3Y,T0)/r - 3.0*g(NTemp,r2,D2Y,T0)/(2.0*r2*r)
                                   + 9.0*g(NTemp,r2,DY,T0)/(4.0*r2*r2*r)
                                   - 15.0*g(NTemp,r2,Y0,T0)/(8.0*r2*r2*r2*r))
                      + i(NTemp)*(g(NTemp,r2,D3Y,DT)/r - 3.0*g(NTemp,r2,D2Y,DT)/(2.0*r2*r)
                                  + 9.0*g(NTemp,r2,DY,DT)/(4.0*r2*r2*r)
                                  - 15.0*g(NTemp,r2,Y0,DT)/(8.0*r2*r2*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(g(NTemp,r2,D3Y,T0)/r - 3.0*g(NTemp,r2,D2Y,T0)/(2.0*r2*r)
                                    + 9.0*g(NTemp,r2,DY,T0)/(4.0*r2*r2*r)
                                    - 15.0*g(NTemp,r2,Y0,T0)/(8.0*r2*r2*r2*r))
                      + 2.0*i(NTemp,DT)*(g(NTemp,r2,D3Y,DT)/r
                                         - 3.0*g(NTemp,r2,D2Y,DT)/(2.0*r2*r)
                                         + 9.0*g(NTemp,r2,DY,DT)/(4.0*r2*r2*r)
                                         - 15.0*g(NTemp,r2,Y0,DT)/(8.0*r2*r2*r2*r))
                      + i(NTemp)*(g(NTemp,r2,D3Y,D2T)/r - 3.0*g(NTemp,r2,D2Y,D2T)/(2.0*r2*r)
                                  + 9.0*g(NTemp,r2,DY,D2T)/(4.0*r2*r2*r)
                                  - 15.0*g(NTemp,r2,Y0,D2T)/(8.0*r2*r2*r2*r)));
               break;
            default:
               cerr << "Error in RadiiMorse::G -- D4Y,D3T" << "\n";
               exit(-1);
         }
         break;
      case D5Y:
         switch (dt)
         {
            case T0:
               val = i(NTemp)*(G(NTemp,r2,D4Y,T0)/r - 2.0*G(NTemp,r2,D3Y,T0)/(r2*r)
                               + 9.0*G(NTemp,r2,D2Y,T0)/(2.0*r2*r2*r)
                               - 15.0*G(NTemp,r2,DY,T0)/(2.0*r2*r2*r2*r)
                               + 105.0*G(NTemp,r2,Y0,T0)/(16.0*r2*r2*r2*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(G(NTemp,r2,D4Y,T0)/r - 2.0*G(NTemp,r2,D3Y,T0)/(r2*r)
                                   + 9.0*G(NTemp,r2,D2Y,T0)/(2.0*r2*r2*r)
                                   - 15.0*G(NTemp,r2,DY,T0)/(2.0*r2*r2*r2*r)
                                   + 105.0*G(NTemp,r2,Y0,T0)/(16.0*r2*r2*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D4Y,DT)/r - 2.0*G(NTemp,r2,D3Y,DT)/(r2*r)
                                  + 9.0*G(NTemp,r2,D2Y,DT)/(2.0*r2*r2*r)
                                  - 15.0*G(NTemp,r2,DY,DT)/(2.0*r2*r2*r2*r)
                                  + 105.0*G(NTemp,r2,Y0,DT)/(16.0*r2*r2*r2*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(G(NTemp,r2,D4Y,T0)/r - 2.0*G(NTemp,r2,D3Y,T0)/(r2*r)
                                    + 9.0*G(NTemp,r2,D2Y,T0)/(2.0*r2*r2*r)
                                    - 15.0*G(NTemp,r2,DY,T0)/(2.0*r2*r2*r2*r)
                                    + 105.0*G(NTemp,r2,Y0,T0)/(16.0*r2*r2*r2*r2*r))
                      + 2.0*i(NTemp,DT)*(G(NTemp,r2,D4Y,DT)/r - 2.0*G(NTemp,r2,D3Y,DT)/(r2*r)
                                         + 9.0*G(NTemp,r2,D2Y,DT)/(2.0*r2*r2*r)
                                         - 15.0*G(NTemp,r2,DY,DT)/(2.0*r2*r2*r2*r)
                                         + 105.0*G(NTemp,r2,Y0,DT)/(16.0*r2*r2*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D4Y,D2T)/r - 2.0*G(NTemp,r2,D3Y,D2T)/(r2*r)
                                  + 9.0*G(NTemp,r2,D2Y,D2T)/(2.0*r2*r2*r)
                                  - 15.0*G(NTemp,r2,DY,D2T)/(2.0*r2*r2*r2*r)
                                  + 105.0*G(NTemp,r2,Y0,D2T)/(16.0*r2*r2*r2*r2*r)));
               break;
            default:
               cerr << "Error in RadiiMorse::G -- D5Y,D2T" << "\n";
               exit(-1);
         }
         break;
      case D6Y:
         switch (dt)
         {
            case T0:
               val = i(NTemp)*(G(NTemp,r2,D5Y,T0)/r - 5.0*G(NTemp,r2,D4Y,T0)/(2.0*r2*r)
                               + 15.0*G(NTemp,r2,D3Y,T0)/(2.0*r2*r2*r)
                               - 75.0*G(NTemp,r2,D2Y,T0)/(4.0*r2*r2*r2*r)
                               + 105.0*5.0*G(NTemp,r2,DY,T0)/(16.0*r2*r2*r2*r2*r)
                               - 105.0*9.0*G(NTemp,r2,Y0,T0)/(32.0*r2*r2*r2*r2*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(G(NTemp,r2,D5Y,T0)/r - 5.0*G(NTemp,r2,D4Y,T0)/(2.0*r2*r)
                                   + 15.0*G(NTemp,r2,D3Y,T0)/(2.0*r2*r2*r)
                                   - 75.0*G(NTemp,r2,D2Y,T0)/(4.0*r2*r2*r2*r)
                                   + 105.0*5.0*G(NTemp,r2,DY,T0)/(16.0*r2*r2*r2*r2*r)
                                   - 105.0*9.0*G(NTemp,r2,Y0,T0)/(32.0*r2*r2*r2*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D5Y,DT)/r - 5.0*G(NTemp,r2,D4Y,DT)/(2.0*r2*r)
                                  + 15.0*G(NTemp,r2,D3Y,DT)/(2.0*r2*r2*r)
                                  - 75.0*G(NTemp,r2,D2Y,DT)/(4.0*r2*r2*r2*r)
                                  + 105.0*5.0*G(NTemp,r2,DY,DT)/(16.0*r2*r2*r2*r2*r)
                                  - 105.0*9.0*G(NTemp,r2,Y0,DT)/(32.0*r2*r2*r2*r2*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(G(NTemp,r2,D5Y,T0)/r - 5.0*G(NTemp,r2,D4Y,T0)/(2.0*r2*r)
                                    + 15.0*G(NTemp,r2,D3Y,T0)/(2.0*r2*r2*r)
                                    - 75.0*G(NTemp,r2,D2Y,T0)/(4.0*r2*r2*r2*r)
                                    + 105.0*5.0*G(NTemp,r2,DY,T0)/(16.0*r2*r2*r2*r2*r)
                                    - 105.0*9.0*G(NTemp,r2,Y0,T0)/(32.0*r2*r2*r2*r2*r2*r))
                      + 2.0*i(NTemp,DT)*(G(NTemp,r2,D5Y,DT)/r - 5.0*G(NTemp,r2,D4Y,DT)/(2.0*r2*r)
                                         + 15.0*G(NTemp,r2,D3Y,DT)/(2.0*r2*r2*r)
                                         - 75.0*G(NTemp,r2,D2Y,DT)/(4.0*r2*r2*r2*r)
                                         + 105.0*5.0*G(NTemp,r2,DY,DT)/(16.0*r2*r2*r2*r2*r)
                                         - 105.0*9.0*G(NTemp,r2,Y0,DT)/(32.0*r2*r2*r2*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D5Y,D2T)/r - 5.0*G(NTemp,r2,D4Y,D2T)/(2.0*r2*r)
                                  + 15.0*G(NTemp,r2,D3Y,D2T)/(2.0*r2*r2*r)
                                  - 75.0*G(NTemp,r2,D2Y,D2T)/(4.0*r2*r2*r2*r)
                                  + 105.0*5.0*G(NTemp,r2,DY,D2T)/(16.0*r2*r2*r2*r2*r)
                                  - 105.0*9.0*G(NTemp,r2,Y0,D2T)/(32.0*r2*r2*r2*r2*r2*r)));
               break;
            default:
               cerr << "Error in RadiiMorse::G -- D6Y,DT" << "\n";
               exit(-1);
         }
         break;
      default:
         cerr << "Error in RadiiMorse::G -- D5Y" << "\n";
         exit(-1);
   }
   
   return val;
}

void RadiiMorse::Print(ostream& out) const
{
   int W=out.width();
   
   out.width(0);
   
   out << "A0=" << setw(W) << A0_
       << "; AT=" << setw(W) << AT_
       << "; B0=" << setw(W) << B0_
       << "; BT=" << setw(W) << BT_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_;
}

ostream& operator<<(ostream& out,RadiiMorse const& A)
{
   A.Print(out);
   return out;
}
