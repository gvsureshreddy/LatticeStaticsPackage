#include "UtilityFunctions.h"

#include <fstream>

#define DIM3 3

using namespace std;

#ifdef UNIX_TERMINAL  // Unix/Linux terminal: use non-blocking to check for debug mode
// ------------------------------------------------------------------------------------
#include <termios.h>
#include <fcntl.h>
#include <unistd.h>
/*
   The fcntl() function accepts a file descriptor and
   an action (depending on the action, it may accept
   a third parameter), and returns a result depending
   on that action.

   If things break bad, it returns -1.
 */

int setblock(int const& file_desc, int const& block)
{
   int flags;

   /* retrieve the file descriptor's flags */
   if ((flags = fcntl(file_desc, F_GETFL)) == -1)
   {
      return false;    /* something went wrong! */
   }

   if (block)
   {
      flags &= ~O_NONBLOCK;      /* we want blocking input */
   }
   else
   {
      flags |= O_NONBLOCK;       /* we want non-blocking input */
   }

   /* set the flags (note the third parameter) */
   fcntl(file_desc, F_SETFL, flags);

   return 1;
}

int EnterDebugMode()
{
   int n;
   char dbg[256];
   setblock(fileno(stdin), 0);
   n = read(fileno(stdin), dbg, 255);
   setblock(fileno(stdin), 1);
   // remove newline
   if (n > 0)
   {
      dbg[n - 1] = 0;

      if (!strcmp(dbg, "debug"))
      {
         return 1;
      }
   }
   return 0;
}
// ------------------------------------------------------------------------------------
#else  // Default debugmode behavior: enter debug first time then ignore calls
// ------------------------------------------------------------------------------------
int EnterDebugMode()
{
   static int flag = 0;

   if (flag)
   {
      return 0;
   }
   else
   {
      return ++flag;
   }
}
// ------------------------------------------------------------------------------------
#endif
// ------------------------------------------------------------------------------------


char kbhitWait()
{
   char t;
   cin.get(t);
   return t;
}

// ======================================================================
int IND2D(int const& i, int const& j);

int IND2D(int const& i, int const& j)
{
   if (i == j)
   {
      return i;
   }
   else
   {
      return 1 + i + j;
   }
}

int FullScanRank1Convex3D(CBKinematics const* const CBK, Matrix const& K, double const& dx)
{
   Matrix A(3, 3);
   Matrix Eigvals(1, 3);
   double Pi = 4.0 * atan(1.0);
   double Piby2 = Pi / 2,
          Pi2 = 2 * Pi;
   double phi, theta;
   double n[3];

   for (phi = -Piby2; phi < Piby2; phi += dx)
   {
      for (theta = 0; theta < Pi2; theta += dx)
      {
         n[0] = cos(phi) * cos(theta);
         n[1] = cos(phi) * sin(theta);
         n[2] = sin(phi);

         // Acoustic tensor
         for (int i = 0; i < 3; i++)
         {
            for (int j = 0; j < 3; j++)
            {
               A[i][j] = 0.0;
               for (int k = 0; k < 3; k++)
               {
                  for (int l = 0; l < 3; l++)
                  {
                     A[i][j] += K[CBK->INDF(i, k)][CBK->INDF(j, l)] * n[k] * n[l];
                  }
               }
            }
         }

         Eigvals = SymEigVal(A);
         for (int i = 0; i < 3; i++)
         {
            if (Eigvals[0][i] <= 0.0)
            {
               return 0;
            }
         }
      }
   }

   return 1;
}

int FullScanRank1Convex2D(Matrix const& K, double const& dx)
{
   Matrix A(2, 2);
   Matrix Eigvals(1, 2);
   double Pi = 4.0 * atan(1.0);
   double Pi2 = 2 * Pi;
   double theta;
   double n[2];

   for (theta = 0; theta < Pi2; theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);

      // Acoustic tensor
      for (int i = 0; i < 2; i++)
      {
         for (int j = 0; j < 2; j++)
         {
            A[i][j] = 0.0;
            for (int k = 0; k < 2; k++)
            {
               for (int l = 0; l < 2; l++)
               {
                  A[i][j] += K[IND2D(i, k)][IND2D(j, l)] * n[k] * n[l];
               }
            }
         }
      }

      Eigvals = SymEigVal(A);
      for (int i = 0; i < 2; i++)
      {
         if (Eigvals[0][i] <= 0.0)
         {
            return 0;
         }
      }
   }

   return 1;
}

int Rank1Convex3D(CBKinematics const* const CBK, Matrix const& K, double const& dx)
{
   double Pi = 4.0 * atan(1.0);
   MyComplexDouble A[3][3][3];
   double n[2];

   for (double theta = 0; theta < 2.0 * Pi; theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);

      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            for (int k = 0; k < 3; k++)
            {
               A[i][j][k] = MyComplexDouble(0.0, 0.0);
            }
         }
      }


      // Calculate Aij Polynomials
      for (int i = 0; i < 3; i++)
      {
         for (int j = 0; j < 3; j++)
         {
            for (int k = 0; k < 2; k++)
            {
               for (int l = 0; l < 2; l++)
               {
                  A[i][j][0] += K[CBK->INDF(i, k)][CBK->INDF(j, l)] * n[k] * n[l];
               }
               A[i][j][1] += K[CBK->INDF(i, k)][CBK->INDF(j, 2)] * n[k]
                             + K[CBK->INDF(i, 2)][CBK->INDF(j, k)] * n[k];
            }
            A[i][j][2] = K[CBK->INDF(i, 2)][CBK->INDF(j, 2)];
         }
      }


      MyComplexDouble Roots[6];
      MyComplexDouble SolveMe[7], PA[5], PB[5];

      PolyRootsLaguerre(A[0][0], 2, Roots, 1);
      for (int i = 0; i < 2; i++)
      {
         if (Roots[i].imag() == 0.0)
         {
            return 0;
         }
      }

      // Only need to check the leading principal minors
      // PolyRootsLaguerre(A[1][1],2,Roots,1);
      // for (int i=0;i<2;i++) if (Roots[i].imag() == 0.0) return 0;

      // Only need to check the leading principal minors
      // PolyRootsLaguerre(A[2][2],2,Roots,1);
      // for (int i=0;i<2;i++) if (Roots[i].imag() == 0.0) return 0;

      PolyMult(A[0][0], 2, A[1][1], 2, PA);
      PolyMult(A[0][1], 2, A[1][0], 2, PB);
      for (int i = 0; i <= 4; i++)
      {
         SolveMe[i] = PA[i] - PB[i];
      }
      PolyRootsLaguerre(SolveMe, 4, Roots, 1);
      for (int i = 0; i < 4; i++)
      {
         if (Roots[i].imag() == 0.0)
         {
            return 0;
         }
      }

      // Only need to check the leading principal minors
      // PolyMult(A[0][0],2,A[2][2],2,PA);
      // PolyMult(A[0][2],2,A[2][0],2,PB);
      // for (int i=0;i<=4;i++) SolveMe[i] = PA[i] - PB[i];
      // PolyRootsLaguerre(SolveMe,4,Roots,1);
      // for (int i=0;i<4;i++) if (Roots[i].imag() == 0.0) return 0;

      // Only need to check the leading principal minors
      // PolyMult(A[1][1],2,A[2][2],2,PA);
      // PolyMult(A[1][2],2,A[2][1],2,PB);
      // for (int i=0;i<=4;i++) SolveMe[i] = PA[i] - PB[i];
      // PolyRootsLaguerre(SolveMe,4,Roots,1);
      // for (int i=0;i<4;i++) if (Roots[i].imag() == 0.0) return 0;

      MyComplexDouble DA[7], DB[7], DC[7];


      PolyMult(A[1][1], 2, A[2][2], 2, PA);
      PolyMult(A[1][2], 2, A[2][1], 2, PB);
      for (int i = 0; i <= 4; i++)
      {
         PA[i] = PA[i] - PB[i];
      }
      PolyMult(A[0][0], 2, PA, 4, DA);

      PolyMult(A[1][0], 2, A[2][2], 2, PA);
      PolyMult(A[1][2], 2, A[2][0], 2, PB);
      for (int i = 0; i <= 4; i++)
      {
         PA[i] = PA[i] - PB[i];
      }
      PolyMult(A[0][1], 2, PA, 4, DB);

      PolyMult(A[1][0], 2, A[2][1], 2, PA);
      PolyMult(A[1][1], 2, A[2][0], 2, PB);
      for (int i = 0; i <= 4; i++)
      {
         PA[i] = PA[i] - PB[i];
      }
      PolyMult(A[0][2], 2, PA, 4, DC);

      for (int i = 0; i <= 6; i++)
      {
         SolveMe[i] = DA[i] - DB[i] + DC[i];
      }

      PolyRootsLaguerre(SolveMe, 6, Roots, 1);
      for (int i = 0; i < 6; i++)
      {
         if (Roots[i].imag() == 0.0)
         {
            return 0;
         }
      }
   }

   return 1;
}

int Rank1Convex2D(Matrix const& K, double const& dx)
{
   double Pi = 4.0 * atan(1.0);
   double A[2][2];
   double n[2];

   for (double theta = 0; theta < 2.0 * Pi; theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);

      for (int i = 0; i < 2; i++)
      {
         for (int j = 0; j < 2; j++)
         {
            A[i][j] = 0.0;
         }
      }


      // Calculate Aij
      for (int i = 0; i < 2; i++)
      {
         for (int j = 0; j < 2; j++)
         {
            for (int k = 0; k < 2; k++)
            {
               for (int l = 0; l < 2; l++)
               {
                  A[i][j] += K[IND2D(i, k)][IND2D(j, l)] * n[k] * n[l];
               }
            }
         }
      }


      if (A[0][0] <= 0.0)
      {
         return 0;
      }
      // Only need to check the leading principal minors
      // if (A[1][1] <= 0.0) return 0;

      if ((A[0][0] * A[1][1] - A[0][1] * A[1][0]) <= 0.0)
      {
         return 0;
      }
   }

   return 1;
}

int pow(int const& a, int const& b);
int pow(int const& a, int const& b)
{
   int c = 1;
   for (int i = 0; i < b; i++)
   {
      c *= a;
   }
   return c;
}

Matrix TranslationProjection1D(int const& NoAtoms)
{
   int k = 0;
   int h, i, j, p, q, r, g, t, m, Rows, Columns, pow2_i;
   double sum, y, R, log2_NoAtoms;

   Rows = NoAtoms;
   Columns = NoAtoms + 1;

   Matrix P(Rows, Columns, 0.0);
   // The following generates first column vector [1 0 0 0...0]

   P[k][0] = 1.0; // the rest of P[k][i] are already zero
   k += 1; // Finishes first vector

   // The following generates floor(log2(NoAtoms)) vectors with
   // floor(double(NoAtoms)/pow(2,(i+1))) entry pairs
   log2_NoAtoms = log2(NoAtoms);

   for (i = 0; i < log2_NoAtoms; i++) // Generates floor(log2(NoAtoms)) entry pairs
   {
      m = int(floor(double(NoAtoms) / pow(2, (i + 1))));
      q = 0;

      pow2_i = pow(2, i);

      R = NoAtoms - (m * pow(2, (i + 1)));
      q = 1;

      for (j = 0; j < m; j++) // generates floor(double(NoAtoms)/pow(2,(i+1))) column vectors
      {
         sum = 0.0;

         for (p = 0; p < pow2_i; p++)
         {
            P[k][q] = 1.0;
            sum += 1.0;
            q += 1;
         }

         for (p = 0; p < pow2_i; p++)
         {
            P[k][q] = -1.0;
            sum += 1.0;
            q += 1;
         }

         sum = sqrt(sum);

         for (h = 0; h < Columns; h++)
         {
            P[k][h] /= sum;
         }
         k += 1;

         // Generates remainder column vectors
         if ((j == m - 1) && (R >= pow2_i))
         {
            sum = 0.0;

            r = (((pow2_i) * 2 * (m)));

            for (g = 0; g < r; g++)
            {
               P[k][g + 1] = 1.0;
               sum += 1.0;
            }

            y = -(g) / (pow2_i);

            for (t = 0; t < (pow2_i); t++)
            {
               P[k][g + 1] = y;
               sum += pow(P[k][g + 1], 2);
               g += 1;
            }

            sum = sqrt(sum);

            for (h = 0; h < Columns; h++)
            {
               P[k][h] /= sum;
            }
            k += 1;
         }
      } // closes column vector (j) index
   } // closes entry pair (i) index

   return P;
}

Matrix TranslationProjection3D(int const& Fsize, int const& NoAtoms)
{
   int i, j, g, h, m, t, p, r, qx, qy, qz, kx, ky, kz, wx, wy, wz, pow2_i;
   double sum, y, log2_NoAtoms, R;
   int Rows = (Fsize - DIM3) + DIM3 * NoAtoms;
   int Columns = Fsize + DIM3 * NoAtoms;

   Matrix P(Rows, Columns, 0.0);

   // The following takes care of the uniform deformation part
   for (i = 0; i < Fsize; i++)
   {
      P[i][i] = 1.0;
   }

   kx = Fsize;
   ky = kx + 1;
   kz = ky + 1;

   // The following generates the Rigid Body Translation Projection Operator
   log2_NoAtoms = log2(NoAtoms);

   for (i = 0; i < log2_NoAtoms; i++)
   {
      m = int(floor(double(NoAtoms) / pow(2, (i + 1))));

      pow2_i = pow(2, i);

      qx = Fsize;
      qy = qx + 1;
      qz = qy + 1;

      R = NoAtoms - (m * pow(2, (i + 1)));

      for (j = 0; j < m; j++)
      {
         sum = 0.0;

         for (p = 0; p < pow2_i; p++)
         {
            P[kx][qx] = 1.0;
            P[ky][qy] = 1.0;
            P[kz][qz] = 1.0;
            sum += 1.0;
            qx += DIM3;
            qy += DIM3;
            qz += DIM3;
         }

         for (p = 0; p < pow2_i; p++)
         {
            P[kx][qx] = -1.0;
            P[ky][qy] = -1.0;
            P[kz][qz] = -1.0;
            sum += 1.0;
            qx += DIM3;
            qy += DIM3;
            qz += DIM3;
         }

         sum = sqrt(sum);

         // normalize
         for (h = Fsize; h < Columns; h++)
         {
            P[kx][h] /= sum;
            P[ky][h] /= sum;
            P[kz][h] /= sum;
         }
         kx += DIM3;
         ky += DIM3;
         kz += DIM3;

         // remainder vectors
         if ((j == m - 1) && (R >= pow2_i))
         {
            sum = 0.0;

            wx = Fsize;
            wy = wx + 1;
            wz = wy + 1;

            r = (((pow2_i) * 2 * (m)));

            for (g = 0; g < r; g++)
            {
               P[kx][wx] = 1.0;
               P[ky][wy] = 1.0;
               P[kz][wz] = 1.0;
               sum += 1.0;
               wx += DIM3;
               wy += DIM3;
               wz += DIM3;
            }

            y = -(g) / (pow2_i);

            for (t = 0; t < (pow2_i); t++)
            {
               P[kx][wx] = y;
               P[ky][wy] = y;
               P[kz][wz] = y;
               sum += pow(P[kx][wx], 2);
               wx += DIM3;
               wy += DIM3;
               wz += DIM3;
            }

            sum = sqrt(sum);

            // normalize
            for (h = Fsize; h < Columns; h++)
            {
               P[kx][h] /= sum;
               P[ky][h] /= sum;
               P[kz][h] /= sum;
            }

            kx += DIM3;
            ky += DIM3;
            kz += DIM3;
         }
      } // closes j
   } // closes i

   return P;
}
