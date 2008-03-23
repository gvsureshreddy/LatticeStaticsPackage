#include "RadiiMorse.h"

RadiiMorse::RadiiMorse(double A0,double B0,double Alpha,double Rref1,double Rref2,
                       double Rtheta1,double Rtheta2):
   A0_(A0),B0_(B0),Alpha_(Alpha),Rref1_(Rref1),Rref2_(Rref2),
   Rtheta1_(Rtheta1),Rtheta2_(Rtheta2)
{
}

double RadiiMorse::A(double NTemp,TDeriv dt)
{
   double retval;
   switch (dt)
   {
      case T0:
         retval = A0_;
         break;
      case DT:
         retval = 0.0;
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

double RadiiMorse::Beta(double NTemp,TDeriv dt)
{
   double retval;
   
   switch (dt)
   {
      case T0:
         retval = B0_ + Alpha_*(NTemp-1.0);
         break;
      case DT:
         retval = Alpha_;
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

double RadiiMorse::Rhat(double NTemp,TDeriv dt)
{
   double rhat;
   
   switch (dt)
   {
      case T0:
         rhat = (Rref1_ + Rtheta1_*(NTemp - 1.0))*(Rref2_ + Rtheta2_*(NTemp - 1.0));
         break;
      case DT:
         rhat = Rtheta1_*(Rref2_ + Rtheta2_*(NTemp - 1.0))
            + Rtheta2_*(Rref1_ + Rtheta1_*(NTemp - 1.0));
         break;
      case D2T:
         rhat = Rtheta1_*Rtheta2_ + Rtheta2_*Rtheta1_;
         break;
      default:
         cerr << "Error in RadiiMorse::Rhat" << "\n";
         exit(-1);
   }
   
   return rhat;
}

double RadiiMorse::PairPotential(double NTemp,double r2,YDeriv dy,TDeriv dt)
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
                  +  a(NTemp)*g(NTemp,r2,Y0,DT)*(g(NTemp,r2,Y0,T0)-2.0)
                  +  a(NTemp)*g(NTemp,r2,Y0,T0)*g(NTemp,r2,Y0,DT);
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
               val = (3.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,DY,T0)
                      + 3.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,DT)
                      + a(NTemp)*g(NTemp,r2,D3Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,Y0,DT)
                      + 3.0*a(NTemp)*g(NTemp,r2,DY,DT)*g(NTemp,r2,D2Y,T0)
                      + 3.0*a(NTemp)*g(NTemp,r2,DY,T0)*g(NTemp,r2,D2Y,DT)
                      + 3.0*a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 3.0*a(NTemp,DT)*g(NTemp,r2,DY,T0)*g(NTemp,r2,D2Y,T0));
               break;
            case D2T:
               val = (6.0*a(NTemp)*g(NTemp,r2,D2Y,D2T)*g(NTemp,r2,DY,T0)
                      + 12.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,DY,DT)
                      + 6.0*a(NTemp)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,D2T)
                      + a(NTemp)*g(NTemp,r2,D3Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 6.0*a(NTemp,D2T)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp,D2T)*g(NTemp,r2,D3Y,T0)*(g(NTemp,r2,Y0,T0)-2.0)
                      + a(NTemp,D2T)*g(NTemp,r2,Y0,T0)*g(NTemp,r2,D3Y,T0)
                      + 12.0*a(NTemp,DT)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,DY,T0)
                      + 12.0*a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,DY,DT)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D3Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
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
               val = (12.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,D2Y,T0)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,DY,T0)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,DT)
                      + a(NTemp)*g(NTemp,r2,D4Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 2.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,Y0,DT)
                      + 6.0*a(NTemp,DT)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 8.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0));
               break;
            case D2T:
               val = (12.0*a(NTemp)*g(NTemp,r2,D2Y,D2T)*g(NTemp,r2,D2Y,T0)
                      + 12.0*a(NTemp)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,D2Y,DT)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,D2T)*g(NTemp,r2,DY,T0)
                      + 16.0*a(NTemp)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,DY,DT)
                      + 8.0*a(NTemp)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,D2T)
                      + a(NTemp)*g(NTemp,r2,D4Y,D2T)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp)*g(NTemp,r2,D4Y,DT)*g(NTemp,r2,Y0,DT)
                      + 2.0*a(NTemp)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,Y0,D2T)
                      + 6.0*a(NTemp,D2T)*g(NTemp,r2,D2Y,T0)*g(NTemp,r2,D2Y,T0)
                      + 8.0*a(NTemp,D2T)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,T0)
                      + a(NTemp,D2T)*g(NTemp,r2,D4Y,T0)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 24.0*a(NTemp,DT)*g(NTemp,r2,D2Y,DT)*g(NTemp,r2,D2Y,T0)
                      + 16.0*a(NTemp,DT)*g(NTemp,r2,D3Y,DT)*g(NTemp,r2,DY,T0)
                      + 16.0*a(NTemp,DT)*g(NTemp,r2,D3Y,T0)*g(NTemp,r2,DY,DT)
                      + 2.0*a(NTemp,DT)*g(NTemp,r2,D4Y,DT)*(2.0*g(NTemp,r2,Y0,T0)-2.0)
                      + 4.0*a(NTemp,DT)*g(NTemp,r2,D4Y,T0)*g(NTemp,r2,Y0,DT));
               break;
            default:
               cerr << "Error in RadiiMorse::PairPotential -- D4Y,D3T" << "\n";
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

double RadiiMorse::I(double NTemp,TDeriv dt)
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

double RadiiMorse::G(double NTemp,double r2,YDeriv dy,TDeriv dt)
{
   double val;
   double r=sqrt(r2);
   
   switch (dy)
   {
      case Y0:
         switch (dt)
         {
            case T0:
               val = exp(-b(NTemp)*(r/rhat(NTemp) - 1.0));
               break;
            case DT:
               val = ((b(NTemp)*rhat(NTemp,DT)*r)/(rhat(NTemp)*rhat(NTemp))
                      - (b(NTemp,DT)*(r/rhat(NTemp) - 1.0)))*G(NTemp,r2,Y0,T0);
               break;
            case D2T:
               val = ((b(NTemp)*rhat(NTemp,DT)*r)/(rhat(NTemp)*rhat(NTemp))
                      - (b(NTemp,DT)*(r/rhat(NTemp) - 1.0)))*G(NTemp,r2,Y0,DT)
                  + ((r/rhat(NTemp)*rhat(NTemp))*(2.0*b(NTemp,DT)*rhat(NTemp,DT)
                                                  + b(NTemp)*rhat(NTemp,D2T))
                     - ((2.0*b(NTemp)*r*rhat(NTemp,DT)*rhat(NTemp,DT))
                        /(rhat(NTemp)*rhat(NTemp)*rhat(NTemp)))
                     - b(NTemp,D2T)*(r/rhat(NTemp) - 1.0))*G(NTemp,r2,Y0,T0);
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
               val = i(NTemp)*(G(NTemp,r2,Y0,T0)/r);
               break;
            case DT:
               val = (i(NTemp,DT)*G(NTemp,r2,Y0,T0) + i(NTemp)*G(NTemp,r2,Y0,DT))/r;
               break;
            case D2T:
               val = (i(NTemp,D2T)*G(NTemp,r2,Y0,T0) + 2.0*i(NTemp,DT)*G(NTemp,r2,Y0,DT)
                      + i(NTemp)*G(NTemp,r2,Y0,D2T))/r;
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
               val = i(NTemp)*(G(NTemp,r2,DY,T0)/r - G(NTemp,r2,Y0,T0)/(2.0*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(G(NTemp,r2,DY,T0)/r - G(NTemp,r2,Y0,T0)/(2.0*r2*r))
                      + i(NTemp)*(G(NTemp,r2,DY,DT)/r - G(NTemp,r2,Y0,DT)/(2.0*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(G(NTemp,r2,DY,T0)/r - G(NTemp,r2,Y0,T0)/(2.0*r2*r))
                      + 2.0*i(NTemp,DT)*(G(NTemp,r2,DY,DT)/r - G(NTemp,r2,Y0,DT)/(2.0*r2*r))
                      + i(NTemp)*(G(NTemp,r2,DY,D2T)/r - G(NTemp,r2,Y0,D2T)/(2.0*r2*r)));
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
               val = i(NTemp)*(G(NTemp,r2,D2Y,T0)/r - G(NTemp,r2,DY,T0)/(r2*r)
                               + 3.0*G(NTemp,r2,Y0,T0)/(4.0*r2*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(G(NTemp,r2,D2Y,T0)/r - G(NTemp,r2,DY,T0)/(r2*r)
                                   + 3.0*G(NTemp,r2,Y0,T0)/(4.0*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D2Y,DT)/r - G(NTemp,r2,DY,DT)/(r2*r)
                                  + 3.0*G(NTemp,r2,Y0,DT)/(4.0*r2*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(G(NTemp,r2,D2Y,T0)/r - G(NTemp,r2,DY,T0)/(r2*r)
                                    + 3.0*G(NTemp,r2,Y0,T0)/(4.0*r2*r2*r))
                      + 2.0*i(NTemp,DT)*(G(NTemp,r2,D2Y,DT)/r - G(NTemp,r2,DY,DT)/(r2*r)
                                         + 3.0*G(NTemp,r2,Y0,DT)/(4.0*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D2Y,D2T)/r - G(NTemp,r2,DY,D2T)/(r2*r)
                                  + 3.0*G(NTemp,r2,Y0,D2T)/(4.0*r2*r2*r)));
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
               val = i(NTemp)*(G(NTemp,r2,D3Y,T0)/r - 3.0*G(NTemp,r2,D2Y,T0)/(2.0*r2*r)
                               + 9.0*G(NTemp,r2,DY,T0)/(4.0*r2*r2*r)
                               - 15.0*G(NTemp,r2,Y0,T0)/(8.0*r2*r2*r2*r));
               break;
            case DT:
               val = (i(NTemp,DT)*(G(NTemp,r2,D3Y,T0)/r - 3.0*G(NTemp,r2,D2Y,T0)/(2.0*r2*r)
                                   + 9.0*G(NTemp,r2,DY,T0)/(4.0*r2*r2*r)
                                   - 15.0*G(NTemp,r2,Y0,T0)/(8.0*r2*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D3Y,DT)/r - 3.0*G(NTemp,r2,D2Y,DT)/(2.0*r2*r)
                                  + 9.0*G(NTemp,r2,DY,DT)/(4.0*r2*r2*r)
                                  - 15.0*G(NTemp,r2,Y0,DT)/(8.0*r2*r2*r2*r)));
               break;
            case D2T:
               val = (i(NTemp,D2T)*(G(NTemp,r2,D3Y,T0)/r - 3.0*G(NTemp,r2,D2Y,T0)/(2.0*r2*r)
                                    + 9.0*G(NTemp,r2,DY,T0)/(4.0*r2*r2*r)
                                    - 15.0*G(NTemp,r2,Y0,T0)/(8.0*r2*r2*r2*r))
                      + 2.0*i(NTemp,DT)*(G(NTemp,r2,D3Y,DT)/r
                                         - 3.0*G(NTemp,r2,D2Y,DT)/(2.0*r2*r)
                                         + 9.0*G(NTemp,r2,DY,DT)/(4.0*r2*r2*r)
                                         - 15.0*G(NTemp,r2,Y0,DT)/(8.0*r2*r2*r2*r))
                      + i(NTemp)*(G(NTemp,r2,D3Y,D2T)/r - 3.0*G(NTemp,r2,D2Y,D2T)/(2.0*r2*r)
                                  + 9.0*G(NTemp,r2,DY,D2T)/(4.0*r2*r2*r)
                                  - 15.0*G(NTemp,r2,Y0,D2T)/(8.0*r2*r2*r2*r)));
               break;
            default:
               cerr << "Error in RadiiMorse::G -- D4Y,D3T" << "\n";
               exit(-1);
         }
         break;
      default:
         cerr << "Error in RadiiMorse::G -- D5Y" << "\n";
         exit(-1);
   }
   
   return val;
}

void RadiiMorse::Print(ostream &out)
{
   int W=out.width();
   
   out.width(0);
   
   out << "A0=" << setw(W) << A0_
       << "; B0=" << setw(W) << B0_
       << "; Alpha=" << setw(W) << Alpha_
       << "; Rref1=" << setw(W) << Rref1_
       << "; Rtheta1=" << setw(W) << Rtheta1_
       << "; Rref2=" << setw(W) << Rref2_
       << "; Rtheta2=" << setw(W) << Rtheta2_;
}

ostream &operator<<(ostream &out,RadiiMorse &A)
{
   A.Print(out);
   return out;
}
