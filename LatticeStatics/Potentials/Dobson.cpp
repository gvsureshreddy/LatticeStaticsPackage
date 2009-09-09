#include "Dobson.h"
#include <cstdlib>

Dobson::Dobson(double const& Eps0,double const& Eps1,double const& Sigma0,double const& Sigma1,
               double const& rcut):
   Eps0_(Eps0),Eps1_(Eps1),Sigma0_(Sigma0),Sigma1_(Sigma1),rcut_(rcut)
{
   if (rcut_ < D3_RCUTU)
   {
      cerr << "Error: rcut < D3_RCUTU in Dobson3.  exiting." << "\n";
      exit(-1);
   }
}

void Dobson::SetParameters(double const* const Vals)
{
   SetEps0(Vals[0]);
   SetEps1(Vals[1]);
   SetSigma0(Vals[2]);
   SetSigma1(Vals[3]);
   Setrcut(Vals[4]);
}

double Dobson::Eps(double const& NTemp,TDeriv const& dt) const
{
   double retval;
   switch (dt)
   {
      case T0:
         retval = Eps0_ + Eps1_*(NTemp - 1.0);
         break;
      case DT:
         retval = Eps1_;
         break;
      case D2T:
         retval = 0.0;
         break;
      default:
         cerr << "Error in Dobson::Eps" << "\n";
         exit(-1);
   }
   
   return retval;
}

double Dobson::Sigma(double const& NTemp,TDeriv const& dt) const
{
   double retval;
   
   switch (dt)
   {
      case T0:
         retval = Sigma0_ + Sigma1_*(NTemp-1.0);
         break;
      case DT:
         retval = Sigma1_;
         break;
      case D2T:
         retval = 0.0;
         break;
      default:
         cerr << "Error in Dobson::Sigma" << "\n";
         exit(-1);
   }
   
   return retval;
}

double Dobson::j(double const& NTemp,double const& r2,YDeriv const& dy,TDeriv const& dt) const
{
   double val=0.0;
   
   switch (dy)
   {
      case Y0:
         switch (dt)
         {
            case T0:
               if ((sqrt(r2) > D3_RCUTL) && (sqrt(r2) < D3_RCUTU))
                  val = D3_jFACT * Eps(NTemp) * pow(sqrt(r2) - D3_RCUTL * Sigma(NTemp),2) *
                     pow(sqrt(r2) - D3_RCUTU * Sigma(NTemp),2);
               else
                  val = 0.0;
               break;
            case DT:
               // fix me
               val = 0.0;
               break;
            case D2T:
               // fix me
               val = 0.0;
               break;
            default:
               cerr << "Error in Dobson::j - Dj/DT and above not programmed" << "\n";
               exit(-1);
               break;
         }
         break;
      case DY:
         switch (dt)
         {
            case T0:
               if ((sqrt(r2) > D3_RCUTL) && (sqrt(r2) < D3_RCUTU))
                  val = (D3_jFACT * Eps(NTemp)/sqrt(r2))*(
                     (sqrt(r2) - 1.5*Sigma(NTemp))*pow(sqrt(r2)-1.7*Sigma(NTemp),2)
                     +pow(sqrt(r2) - 1.5*Sigma(NTemp),2)*(sqrt(r2)-1.7*Sigma(NTemp)));
               else
                  val = 0.0;
               break;
            case DT:
               val = 0.0;
               break;
            default:
               cerr << "Error in Dobson::j - D2j/DYDT and above not programmed" << "\n";
               exit(-1);
               break;
         }
         break;
      case D2Y:
         switch (dt)
         {
            case T0:
               if ((sqrt(r2) > D3_RCUTL) && (sqrt(r2) < D3_RCUTU))
                  val = (-D3_jFACT * Eps(NTemp)/(2.0*pow(sqrt(r2),3)))*(
                     (sqrt(r2) - 1.5*Sigma(NTemp))*pow(sqrt(r2)-1.7*Sigma(NTemp),2)
                     +pow(sqrt(r2) - 1.5*Sigma(NTemp),2)*(sqrt(r2)-1.7*Sigma(NTemp)))
                     +(D3_jFACT * Eps(NTemp)/r2)*(0.5*pow(sqrt(r2)-1.5*Sigma(NTemp),2)
                                                  + 2.0*(sqrt(r2)-1.5*Sigma(NTemp))*
                                                  (sqrt(r2)-1.7*Sigma(NTemp)) +
                                                  0.5*pow(sqrt(r2)-1.7*Sigma(NTemp),2));
               else
                  val = 0.0;
               break;
            default:
               cerr << "Error in Dobson::j - D3j/D2YDT and above not programmed" << "\n";
               exit(-1);
               break;
         }
         break;
      case D3Y:
         switch (dt)
         {
            case T0:
               if ((sqrt(r2) > D3_RCUTL) && (sqrt(r2) < D3_RCUTU))
                  val = (0.75*D3_jFACT * Eps(NTemp)/pow(r2,2.5))*(
                     (sqrt(r2)-1.5*Sigma(NTemp))*pow(sqrt(r2)-1.7*Sigma(NTemp),2)
                     +pow(sqrt(r2)-1.5*Sigma(NTemp),2)*(sqrt(r2)-1.7*Sigma(NTemp)))
                     - (1.5*D3_jFACT*Eps(NTemp)/(r2*r2))
                     *(0.5*pow(sqrt(r2)-1.7*Sigma(NTemp),2)
                       + 2.0*(sqrt(r2)-1.5*Sigma(NTemp))*(sqrt(r2)-1.7*Sigma(NTemp))
                       + 0.5*pow(sqrt(r2)-1.5*Sigma(NTemp),2))
                     + (1.5*D3_jFACT*Eps(NTemp)/pow(r2,1.5))*
                     ((sqrt(r2)-1.7*Sigma(NTemp)) + (sqrt(r2)-1.5*Sigma(NTemp)));
               else
                  val = 0.0;
               break;
            default:
               cerr << "Error in Dobson::j - D4j/D3YDT and above not programmed" << "\n";
               exit(-1);
               break;
         }
         break;
      case D4Y:
         switch (dt)
         {
            case T0:
               if ((sqrt(r2) > D3_RCUTL) && (sqrt(r2) < D3_RCUTU))
                  val = (-15.0*D3_jFACT*Eps(NTemp)*pow(r2,-3.2)/8.0)*(
                     (sqrt(r2)-1.5*Sigma(NTemp))*pow(sqrt(r2)-1.7*Sigma(NTemp),2)
                     +pow(sqrt(r2)-1.5*Sigma(NTemp),2)*(sqrt(r2)-1.7*Sigma(NTemp)))
                     +(27.0*D3_jFACT*Eps(NTemp)*pow(r2,-3)/8.0)*
                     (0.5*pow(sqrt(r2)-1.7*Sigma(NTemp),2)
                      + 2.0*(sqrt(r2)-1.5*Sigma(NTemp))*(sqrt(r2)-1.7*Sigma(NTemp))
                      + 0.5*pow(sqrt(r2)-1.5*Sigma(NTemp),2))
                     -(4.5*D3_jFACT*Eps(NTemp)*pow(r2,-2.5))*
                     ((sqrt(r2)-1.7*Sigma(NTemp)) + (sqrt(r2)-1.5*Sigma(NTemp)))
                     + 1.5*D3_jFACT*Eps(NTemp)*pow(r2,-2);
               else
                  val = 0.0;
               break;
            default:
               cerr << "Error in Dobson::j - D5j/D3YDT and above not programmed" << "\n";
               exit(-1);
               break;
         }
         break;
      default:
         cerr << "Error in Dobson::j" << "\n";
         exit(-1);
   }
   
   return val;
}


double Dobson::A(double const& NTemp,TDeriv const& dt) const
{
   double val=0.0;
   
   switch (dt)
   {
      case T0:
         val = (12.0*Eps(NTemp)/(rcut_*rcut_))*(
            2.0*pow(Sigma(NTemp)/rcut_,12) - pow(Sigma(NTemp)/rcut_,6));
         break;
      case DT:
         val = (12.0*Eps(NTemp,DT)/(rcut_*rcut_))*(
            2.0*pow(Sigma(NTemp)/rcut_,12) - pow(Sigma(NTemp)/rcut_,6))
            +(72.0*Eps(NTemp)*Sigma(NTemp,DT)/(rcut_*rcut_*rcut_))*(
               4.0*pow(Sigma(NTemp)/rcut_,11) - pow(Sigma(NTemp)/rcut_,5));
         break;
      case D2T:
         val = (12.0*Eps(NTemp,D2T)/(rcut_*rcut_))*(
            2.0*pow(Sigma(NTemp)/rcut_,12) - pow(Sigma(NTemp)/rcut_,6))
            + (144.0*Eps(NTemp,DT)*Sigma(NTemp,DT)/(rcut_*rcut_*rcut_)
               + 72.0*Eps(NTemp)*Sigma(NTemp,D2T)/(rcut_*rcut_*rcut_))*(
                  4.0*pow(Sigma(NTemp)/rcut_,11) - pow(Sigma(NTemp)/rcut_,5))
            + (72.0*Eps(NTemp)*Sigma(NTemp,DT)*Sigma(NTemp,DT)/pow(rcut_,4))*(
               44.0*pow(Sigma(NTemp)/rcut_,10) - 5.0*pow(Sigma(NTemp)/rcut_,4));
         break;
      default:
         cerr << "Error in Dobson::A" << "\n";
         exit(-1);
         break;
   }
   
   return val;
}

double Dobson::B(double const& NTemp,TDeriv const& dt) const
{
   double val=0.0;
   
   switch (dt)
   {
      case T0:
         val = -Eps(NTemp)*(
            28.0*pow(Sigma(NTemp)/rcut_,12) - 16.0*pow(Sigma(NTemp)/rcut_,6));
         break;
      case DT:
         val = -Eps(NTemp,DT)*(
            28.0*pow(Sigma(NTemp)/rcut_,12) - 16.0*pow(Sigma(NTemp)/rcut_,6))
            - (Eps(NTemp)*Sigma(NTemp,DT)/rcut_)*(
               336.0*pow(Sigma(NTemp)/rcut_,11) - 96.0*pow(Sigma(NTemp)/rcut_,5));
         break;
      case D2T:
         val = -Eps(NTemp,D2T)*(
            28.0*pow(Sigma(NTemp)/rcut_,12) - 16.0*pow(Sigma(NTemp)/rcut_,6))
            - (2.0*Eps(NTemp,DT)*Sigma(NTemp,DT)/rcut_ + Eps(NTemp)*Sigma(NTemp,D2T)/rcut_)*(
               336.0*pow(Sigma(NTemp)/rcut_,11) - 96.0*pow(Sigma(NTemp)/rcut_,5))
            - (Eps(NTemp,DT)*Sigma(NTemp,DT)*Sigma(NTemp,DT)/(rcut_*rcut_))*(
               3696.0*pow(Sigma(NTemp)/rcut_,10) - 480.0*pow(Sigma(NTemp)/rcut_,4));
         break;
      default:
         cerr << "Error in Dobson::B" << "\n";
         exit(-1);
         break;
   }
   
   return val;
}



double Dobson::PairPotential(double const& NTemp,double const& r2,YDeriv const& dy,
                             TDeriv const& dt) const
{
   double S=Sigma(NTemp),
      E=Eps(NTemp),
      r = sqrt(r2);
   
   double val=0;
   
   switch (dy)
   {
      case Y0:
         switch (dt)
         {
            case T0:
               val = 4.0*E*(pow(S/r,12.0) - pow(S/r,6.0))
                  + j(NTemp,r2) + A(NTemp)*r2 + B(NTemp);
               break;
            case DT:
               val = 4.0*Eps(NTemp,DT)*(pow(S/r,12.0) - pow(S/r,6.0))
                  + (24.0*E*Sigma(NTemp,DT)/r)*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  + j(NTemp,r2,Y0,DT) + A(NTemp,DT)*r2 + B(NTemp,DT);
               break;
            case D2T:
               val = 4.0*Eps(NTemp,D2T)*(pow(S/r,12.0) - pow(S/r,6.0))
                  + (24.0/r)*(2.0*Eps(NTemp,DT)*Sigma(NTemp,DT)
                              + E*Sigma(NTemp,D2T))
                  *(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  + (24*E*pow(Sigma(NTemp,DT),2.0)/r2)
                  *(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
                  + j(NTemp,r2,Y0,D2T) + A(NTemp,D2T)*r2 + B(NTemp,D2T);
               break;
            default:
               cerr << "Error in Dobson::PairPotential" << "\n";
               exit(-1);}
         break;
      case DY:
         switch (dt)
         {
            case T0:
               val = -(12.0*E*S/pow(r2,1.5))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  + j(NTemp,r2,DY) + A(NTemp);
               break;
            case DT:
               val = -(12.0/pow(r2,1.5))*(Eps(NTemp,DT)*S + E*Sigma(NTemp,DT))
                  *(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  - (12.0*E*S*Sigma(NTemp,DT)/pow(r2,2.0))
                  *(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
                  + j(NTemp,r2,DY,DT) + A(NTemp,DT);
               break;
            default:
               cerr << "Error in Dobson::PairPotential -- DY,D2T not programmed" << "\n";
               exit(-1);
         }
         break;
      case D2Y:
         switch (dt)
         {
            case T0:
               val = (18.0*E*S/pow(r2,2.5))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  +(6.0*E*S*S/pow(r2,3.0))*(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
                  + j(NTemp,r2,D2Y);
               break;
            case DT:
               val = (18.0/pow(r2,2.5))*(Eps(NTemp,DT)*S + E*Sigma(NTemp,DT))
                  *(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  +(6.0/pow(r2,3.0))*(Eps(NTemp,DT)*S*S + 5.0*E*S*Sigma(NTemp,DT))
                  *(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
                  +(120.0*E*S*S*Sigma(NTemp,DT)/pow(r2,3.5))
                  *(11.0*pow(S/r,9.0) - pow(S/r,3.0))
                  + j(NTemp,r2,D2Y,DT);
               break;
            default:
               cerr << "Error in Dobson::PairPotential -- D2Y,D2T not programmed" << "\n";
               exit(-1);
         }
         break;
      case D3Y:
         switch (dt)
         {
            case T0:
               val = (-45.0*E*S/pow(r2,3.5))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  -(27.0*E*S*S/pow(r2,4.0))*(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
                  -(60.0*E*S*S*S/pow(r2,4.5))*(11.0*pow(S/r,9.0) - pow(S/r,3.0))
                  + j(NTemp,r2,D3Y);
               break;
            case DT:
               cerr << "D4phi/Dy3DT Not Coded... " << "\n";
               exit(-1);
               break;
            default:
               cerr << "Error in Dobson::PairPotential -- D3Y,D2T not programmed" << "\n";
               exit(-1);
         }
         break;
      case D4Y:
         switch (dt)
         {
            case T0:
               val = (315.0*E*S/(2.0*pow(r2,4.5)))*(2.0*pow(S/r,11.0) - pow(S/r,5.0))
                  +(261.0*E*S*S/(2.0*pow(r2,5.0)))*(22.0*pow(S/r,10.0) - 5.0*pow(S/r,4.0))
                  +(540.0*E*S*S*S/pow(r2,5.5))*(11.0*pow(S/r,9.0) - pow(S/r,3.0))
                  +(90.0*E*S*S*S*S/pow(r2,6.0))*(33.0*pow(S/r,8.0) - pow(S/r,2.0))
                  + j(NTemp,r2,D4Y);
               break;
            case DT:
               cerr << "D5phi/Dy4DT Not Coded... " << "\n";
               exit(-1);
               break;
            default:
               cerr << "Error in Dobson::PairPotential -- D4Y,D2T not programmed" << "\n";
               exit(-1);
         }
         break;
      case DYmax:
      default:
         cerr << "Error in Dobson::PairPotential -- DYmax out of range\n";
         exit(-1);
   }
   return val;
}

void Dobson::Print(ostream& out) const
{
   int W=out.width();
   
   out.width(0);
   
   out << "Eps0=" << setw(W) << Eps0_
       << "; Eps1=" << setw(W) << Eps1_
       << "; Sigma0=" << setw(W) << Sigma0_
       << "; Sigma1=" << setw(W) << Sigma1_
       << "; rcut=" << setw(W) << rcut_;
}

ostream& operator<<(ostream& out,Dobson const& A)
{
   A.Print(out);
   return out;
}
