#include "PerlInput.h"
#include "KnownPairPotentials.h"

char* builddate();

int main(int argc, char* argv[])
{
   if (argc < 10)
   {
      cerr << "Usage: "
           << argv[0] << " inputfile Aatom Batom Tempstart Tempend temppts ystart "
                    "yend points [print numerical derivatives?(1/0)]" << "\n";
      exit(-1);
   }

   int i = atoi(argv[2]);
   int j = atoi(argv[3]);
   int tmp;
   if (i > j)
   {
      tmp = i; i = j; j = tmp;
   }
   double
      tempst = atof(argv[4]),
      tempend = atof(argv[5]),
      temppts = atof(argv[6]),
      ystart = atof(argv[7]),
      yend = atof(argv[8]),
      pts = atof(argv[9]);
   int deriv;
   if (argc == 11)
   {
      deriv = atoi(argv[10]);
   }
   else
   {
      deriv = 0;
   }

   PerlInput Input(argv[1]);
   PerlInput::HashStruct Hash = Input.getHash("Lattice");
   Hash = Input.getHash(Hash, Input.getString(Hash, "Type"));
   PairPotentials* const pot = InitializePairPotential(Hash, Input, i, j);
   Input.EndofInputSection();

   int Width, Precision;
   Width = Input.getPosInt("Main", "FieldWidth");
   if (Input.ParameterOK("Main", "Precision"))
   {
      Precision = Input.getPosInt("Main", "Precision");
   }
   else
   {
      Precision = Input.usePosInt(10, "Main", "Precision"); // Default Value
   }
   Input.EndofInputSection();

   cout << setiosflags(ios::fixed) << setprecision(Precision);
   double inc = (yend - ystart) / pts,
          tempinc = (tempend - tempst) / temppts;

   cout << "#" << setw(Width) << "Temperature"
        << setw(Width) << "Y"
        << setw(Width) << "Potential"
        << setw(Width) << "DY";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D2Y";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D3Y";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D4Y";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "DT";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D2T";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "DYDT";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D2YDT";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D2YD2T";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D3YDT";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D3YD2T";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D4YDT";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << setw(Width) << "D4YD2T";
   if (deriv)
   {
      cout << setw(Width) << "numerical";
   }
   cout << "\n";


   for (double t = tempst; t <= tempend; t += tempinc)
   {
      for (double q = ystart; q <= yend; q += inc)
      {
         cout << setw(Width) << t
              << setw(Width) << q
              << setw(Width) << pot->PairPotential(t, q)
              << setw(Width) << pot->PairPotential(t, q, PairPotentials::DY);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t, q + inc, PairPotentials::Y0, PairPotentials::T0)
                - pot->PairPotential(t, q - inc, PairPotentials::Y0, PairPotentials::T0))
            / (2.0 * inc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D2Y);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t, q + inc, PairPotentials::DY, PairPotentials::T0)
                - pot->PairPotential(t, q - inc, PairPotentials::DY, PairPotentials::T0))
            / (2.0 * inc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D3Y);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t, q + inc, PairPotentials::D2Y, PairPotentials::T0)
                - pot->PairPotential(t, q - inc, PairPotentials::D2Y, PairPotentials::T0))
            / (2.0 * inc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D4Y);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t, q + inc, PairPotentials::D3Y, PairPotentials::T0)
                - pot->PairPotential(t, q - inc, PairPotentials::D3Y, PairPotentials::T0))
            / (2.0 * inc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::Y0, PairPotentials::DT);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::Y0, PairPotentials::T0)
                - pot->PairPotential(t - tempinc, q, PairPotentials::Y0, PairPotentials::T0))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::Y0, PairPotentials::D2T);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::Y0, PairPotentials::DT)
                - pot->PairPotential(t - tempinc, q, PairPotentials::Y0, PairPotentials::DT))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::DY, PairPotentials::DT);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::DY, PairPotentials::T0)
                - pot->PairPotential(t - tempinc, q, PairPotentials::DY, PairPotentials::T0))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D2Y, PairPotentials::DT);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::D2Y, PairPotentials::T0)
                - pot->PairPotential(t - tempinc, q, PairPotentials::D2Y, PairPotentials::T0))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D2Y, PairPotentials::D2T);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::D2Y, PairPotentials::DT)
                - pot->PairPotential(t - tempinc, q, PairPotentials::D2Y, PairPotentials::DT))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D3Y, PairPotentials::DT);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::D3Y, PairPotentials::T0)
                - pot->PairPotential(t - tempinc, q, PairPotentials::D3Y, PairPotentials::T0))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D3Y, PairPotentials::D2T);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::D3Y, PairPotentials::DT)
                - pot->PairPotential(t - tempinc, q, PairPotentials::D3Y, PairPotentials::DT))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D4Y, PairPotentials::DT);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::D4Y, PairPotentials::T0)
                - pot->PairPotential(t - tempinc, q, PairPotentials::D4Y, PairPotentials::T0))
            / (2.0 * tempinc);
         }
         cout << setw(Width) << pot->PairPotential(t, q, PairPotentials::D4Y, PairPotentials::D2T);
         if (deriv)
         {
            cout << setw(Width)
                 << (pot->PairPotential(t + tempinc, q, PairPotentials::D4Y, PairPotentials::DT)
                - pot->PairPotential(t - tempinc, q, PairPotentials::D4Y, PairPotentials::DT))
            / (2.0 * tempinc);
         }
         cout << "\n";
      }
      cout << "\n";
   }
   delete pot;
}

