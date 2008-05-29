#include "mex.h"

#include "PerlInput.h"
#include <Matrix.h>
#include <Vector.h>
#include "KnownLattices.h"

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
   int j;
   double *output;
   double energy[3];

   static int flag=1;

   static PerlInput Input;
   if (flag)
   {
      cerr << "setup 1\n";
      Input.Readfile("Input1");
   }
   static MultiLatticeTPP Lat0(Input);
   if (flag)
   {
      cerr << "setup 2\n";
      Input.Readfile("Input2");
   }
   static MultiLatticeTPP Lat1(Input);
   if (flag)
   {
      cerr << "setup 3\n";
      Input.Readfile("Input3");
   }
   static MultiLatticeTPP Lat2(Input);

   static int dofs[3];
   if (flag)
   {
      dofs[0] = (Lat0.DOF()).Dim();
      dofs[1] = (Lat1.DOF()).Dim();
      dofs[2] = (Lat2.DOF()).Dim();
   }
   flag=0;

   energy[0] = Lat0.E0();
   energy[1] = Lat1.E0();
   energy[2] = Lat2.E0();

// Create the array
   plhs[0] = mxCreateDoubleMatrix(1,3,mxREAL);
   output = mxGetPr(plhs[0]);


   for (j = 0; j < 3; j++)
   {
      output[j] = energy[j];
   }
}
