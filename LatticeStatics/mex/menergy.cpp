#include "mex.h"

extern int dofs[3];
extern double *energy;
extern void start();
extern double GetEnergies();
extern void SetDOF(double *DOFs);
extern void quit();

void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
   int j;
   double *output;
   
   // Create the array
   plhs[0] = mxCreateDoubleMatrix(1,3,mxREAL);
   output = mxGetPr(plhs[0]);

   start();

   GetEnergies();

   for (j = 0; j < 3; j++)
   {
      output[j] = energy[j];
   }

   quit();   
}
