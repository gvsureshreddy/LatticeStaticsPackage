#include "mex.h"

#include "PerlInput.h"
#include <Matrix.h>
#include <Vector.h>
#include "KnownLattices.h"

#define param0 6
#define param1 6
#define param2 6
#define totparams 18
// RadiiMorse has 8 parameters
// [A0, AT, B0, BT, Rref1, Rtheta1, Rref2, Rtheta2]
// Lat0 -- FCC Ni; Lat1 -- FCC Ti; Lat2 -- B2 NiTi

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   int j;
   double* output;

   static int flag = 1;
   cerr << setiosflags(ios::fixed) << setprecision(15);

   static PerlInput Input;
   if (flag)
   {
      cerr << "setup 1\n";
      Input.Readfile("Input1");
   }
   static MultiLatticeTPP Lat0(Input);
   if (flag)
   {
      Input.ClearHash("Main");
      Input.ClearHash("Lattice");
      cerr << "setup 2\n";
      Input.Readfile("Input2");
   }
   static MultiLatticeTPP Lat1(Input);
   if (flag)
   {
      Input.ClearHash("Main");
      Input.ClearHash("Lattice");
      cerr << "setup 3\n";
      Input.Readfile("Input3");
   }
   static MultiLatticeTPP Lat2(Input);

   flag = 0;

   int m, n;
   m = mxGetM(prhs[0]);
   n = mxGetN(prhs[0]);
   if ((m != 1) && (n != totparams))
   {
      mexPrintf("Input must be of size 1x%i\n", totparams);
      mexErrMsgTxt("Input wrong size. exiting");
   }

   double* params;

   params = mxGetPr(prhs[0]);

   for (int i = 0; i < 15; ++i)
   {
      mexPrintf("%17.15lf\t", params[i]);
      cerr << params[i] << " \t";
   }
   cerr << "\n";

   double p[8];
   double q[24];

   for (int i = 0; i < 6; i++)
   {
      p[i] = params[i];
   }
   p[6] = 1.0;
   p[7] = 0.0;
   Lat0.SetParameters(p);

   for (int i = 6; i < 12; i++)
   {
      p[i - 6] = params[i];
   }
   p[6] = 1.0;
   p[7] = 0.0;
   Lat1.SetParameters(p);

   for (int i = 0; i < 6; i++)
   {
      q[i] = params[i];
   }
   q[6] = 1.0;
   q[7] = 0.0;
   for (int i = 12; i < 18; i++)
   {
      q[i - 4] = params[i];
   }
   q[14] = 1.0;
   q[15] = 0.0;
   for (int i = 6; i < 12; i++)
   {
      q[i + 10] = params[i];
   }
   q[22] = 1.0;
   q[23] = 0.0;
   Lat2.SetParameters(q);

   // Create the array
   plhs[0] = mxCreateDoubleMatrix(1, 15, mxREAL);
   output = mxGetPr(plhs[0]);

   double RefVolume_0, RefVolume_1, RefVolume_2;

   double a_Au_exp, a_Cd_exp, a_AuCd_exp, bulk_Au_exp, bulk_Cd_exp, bulk_AuCd_exp, Ec_Au_exp, Ec_Cd_exp, Ec_AuCd_exp;
   double Cp_Au_exp, Cp_Cd_exp, Cp_AuCd_exp, alpha_Au_exp, alpha_Cd_exp, alpha_AuCd_exp;
   double T_ref, conversion1, conversion2, N_A;

   a_Au_exp = 4.08122143; // Angstrom
   a_Cd_exp = 2.97514419;
   a_AuCd_exp = 3.3165;

   bulk_Au_exp = 172.210933; // GPa
   bulk_Cd_exp = 48.3793758;
   bulk_AuCd_exp = 95.7854406;

   Ec_Au_exp = -3.8142161; // eV/atom
   Ec_Cd_exp = -1.1570427;
   Ec_AuCd_exp = 1 * Ec_Au_exp + 1 * Ec_Cd_exp - 0.414580521;

   Cp_Au_exp = 25.4255929; // J/(mol K)
   Cp_Cd_exp = 26.2368431;
   Cp_AuCd_exp = 26.3181507;

   alpha_Au_exp = 14.2355 * pow(10.0, -6.0); // K^-1
   alpha_Cd_exp = 20.460619 * pow(10.0, -6.0);
   alpha_AuCd_exp = 20.9 * pow(10.0, -6.0);

   T_ref = 323; // K
   conversion1 = 0.00624150974; // 1 GPa = 0.00624150974 eV/(angstrom^3)
   // conversion2 = 6.24150974*pow(10.0,18.0); // 1 J = conversion2 eV
   // N_A = 6.023*pow(10.0,23.0); // Avagadro number
   conversion2 = (6.24150974 / 6.023) * pow(10.0, 18 - 23);

   RefVolume_0 = Lat0.RefLattice()[0][0] * (Lat0.RefLattice()[1][1] * Lat0.RefLattice()[2][2] - Lat0.RefLattice()[1][2] * Lat0.RefLattice()[2][1]) - Lat0.RefLattice()[0][1] * (Lat0.RefLattice()[1][0] * Lat0.RefLattice()[2][2] - Lat0.RefLattice()[1][2] * Lat0.RefLattice()[2][1]) + Lat0.RefLattice()[0][2] * (Lat0.RefLattice()[1][0] * Lat0.RefLattice()[2][1] - Lat0.RefLattice()[1][1] * Lat0.RefLattice()[2][1]);
   RefVolume_1 = Lat1.RefLattice()[0][0] * (Lat1.RefLattice()[1][1] * Lat1.RefLattice()[2][2] - Lat1.RefLattice()[1][2] * Lat1.RefLattice()[2][1]) - Lat1.RefLattice()[0][1] * (Lat1.RefLattice()[1][0] * Lat1.RefLattice()[2][2] - Lat1.RefLattice()[1][2] * Lat1.RefLattice()[2][1]) + Lat1.RefLattice()[0][2] * (Lat1.RefLattice()[1][0] * Lat1.RefLattice()[2][1] - Lat1.RefLattice()[1][1] * Lat1.RefLattice()[2][1]);
   RefVolume_2 = Lat2.RefLattice()[0][0] * (Lat2.RefLattice()[1][1] * Lat2.RefLattice()[2][2] - Lat2.RefLattice()[1][2] * Lat2.RefLattice()[2][1]) - Lat2.RefLattice()[0][1] * (Lat2.RefLattice()[1][0] * Lat2.RefLattice()[2][2] - Lat2.RefLattice()[1][2] * Lat2.RefLattice()[2][1]) + Lat2.RefLattice()[0][2] * (Lat2.RefLattice()[1][0] * Lat2.RefLattice()[2][1] - Lat2.RefLattice()[1][1] * Lat2.RefLattice()[2][1]);

   output[0] = Lat0.RefLattice()[0][0] - a_Au_exp;
   output[1] = Lat1.RefLattice()[0][0] - a_Cd_exp;
   output[2] = Lat2.RefLattice()[0][0] - a_AuCd_exp;

   output[3] = -bulk_Au_exp * conversion1 + (1 / RefVolume_0) * (Lat0.CondensedModuli()[0][0] + 2 * Lat0.CondensedModuli()[0][1]) / 3;
   output[4] = -bulk_Cd_exp * conversion1 + (1 / RefVolume_1) * ((Lat1.CondensedModuli()[0][0] + Lat1.CondensedModuli()[0][1]) * Lat1.CondensedModuli()[2][2] - 2 * Lat1.CondensedModuli()[0][2] * Lat1.CondensedModuli()[0][2]) / (Lat1.CondensedModuli()[0][0] + Lat1.CondensedModuli()[0][1] + 2 * Lat1.CondensedModuli()[2][2] - 4 * Lat1.CondensedModuli()[0][2]);
   output[5] = -bulk_AuCd_exp * conversion1 + (1 / RefVolume_2) * (Lat2.CondensedModuli()[0][0] + 2 * Lat2.CondensedModuli()[0][1]) / 3;

   output[6] = Lat0.E0() / 4 - Ec_Au_exp;
   output[7] = Lat1.E0() / 2 - Ec_Cd_exp;
   output[8] = Lat2.E0() / 4 - Ec_AuCd_exp;


   Matrix TE0 = Lat0.ThermalExpansion();
   Matrix TE1 = Lat1.ThermalExpansion();
   Matrix TE2 = Lat2.ThermalExpansion();

   output[9] = -alpha_Au_exp * T_ref + TE0[0][0];
   output[10] = -alpha_Cd_exp * T_ref + TE1[0][0];
   output[11] = -alpha_AuCd_exp * T_ref + TE2[0][0];

   Matrix StressDT0 = Lat0.StressDT();
   Matrix StressDT1 = Lat1.StressDT();
   Matrix StressDT2 = Lat2.StressDT();

   output[12] = -Cp_Au_exp * conversion2 * T_ref + (Lat0.HeatCapacity() + StressDT0[0][0] * TE0[0][0] + StressDT0[0][1] * TE0[0][1] + StressDT0[0][2] * TE0[0][2] + StressDT0[0][3] * TE0[0][3] + StressDT0[0][4] * TE0[0][4] + StressDT0[0][5] * TE0[0][5] + StressDT0[0][6] * TE0[0][6] + StressDT0[0][7] * TE0[0][7] + StressDT0[0][8] * TE0[0][8] + StressDT0[0][9] * TE0[0][9] + StressDT0[0][10] * TE0[0][10] + StressDT0[0][11] * TE0[0][11] + StressDT0[0][12] * TE0[0][12] + StressDT0[0][13] * TE0[0][13] + StressDT0[0][14] * TE0[0][14]) / 4;
   output[13] = -Cp_Cd_exp * conversion2 * T_ref + (Lat1.HeatCapacity() + StressDT1[0][0] * TE1[0][0] + StressDT1[0][1] * TE1[0][1] + StressDT1[0][2] * TE1[0][2] + StressDT1[0][3] * TE1[0][3] + StressDT1[0][4] * TE1[0][4] + StressDT1[0][5] * TE1[0][5] + StressDT1[0][6] * TE1[0][6] + StressDT1[0][7] * TE1[0][7] + StressDT1[0][8] * TE1[0][8]) / 2;
   output[14] = -Cp_AuCd_exp * conversion2 * T_ref + (Lat2.HeatCapacity() + StressDT2[0][0] * TE2[0][0] + StressDT2[0][1] * TE2[0][1] + StressDT2[0][2] * TE2[0][2] + StressDT2[0][3] * TE2[0][3] + StressDT2[0][4] * TE2[0][4] + StressDT2[0][5] * TE2[0][5] + StressDT2[0][6] * TE2[0][6] + StressDT2[0][7] * TE2[0][7] + StressDT2[0][8] * TE2[0][8] + StressDT2[0][9] * TE2[0][9] + StressDT2[0][10] * TE2[0][10] + StressDT2[0][11] * TE2[0][11] + StressDT2[0][12] * TE2[0][12] + StressDT2[0][13] * TE2[0][13] + StressDT2[0][14] * TE2[0][14]) / 4;

   for (int i = 0; i < 15; ++i)
   {
      mexPrintf("%17.15lf\t", output[i]);
   }
   mexPrintf("\n");
}

