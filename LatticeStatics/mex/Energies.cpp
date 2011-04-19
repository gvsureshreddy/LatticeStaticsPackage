#include "mex.h"

#include "PerlInput.h"
#include <Matrix.h>
#include <Vector.h>
#include "KnownLattices.h"

#define param0 8
#define param1 8
#define param2 24
#define totparams 40
// RadiiMorse has 8 parameters
// [A0, AT, B0, BT, Rref1, Rtheta1, Rref2, Rtheta2]
// Lat0 -- FCC Ni; Lat1 -- FCC Ti; Lat2 -- B2 NiTi

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
   int j;
   double* output;

   static int flag = 1;

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

   Lat0.SetParameters(&(params[0]));
   Lat1.SetParameters(&(params[param0]));
   Lat2.SetParameters(&(params[param0 + param1]));

   // Create the array
   plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
   output = mxGetPr(plhs[0]);


   output[0] = Lat0.E0();
   output[1] = Lat1.E0();
   output[2] = Lat2.E0();
}

