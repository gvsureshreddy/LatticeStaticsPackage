#include "PerlInput.h"
#include <Matrix.h>
#include <Vector.h>
#include "KnownLattices.h"

PerlInput Input("Input1");

MultiLatticeTPP *Lat[3];
int dofs[3];
double *energy = NULL;


void start()
{
   Lat[0] = new MultiLatticeTPP(Input);
   dofs[0] = (Lat[0]->DOF()).Dim();

   Input.Readfile("Input2");
   Lat[1] = new MultiLatticeTPP(Input);
   dofs[1] = (Lat[1]->DOF()).Dim();

   Input.Readfile("Input3");
   Lat[2] = new MultiLatticeTPP(Input);
   dofs[2] = (Lat[2]->DOF()).Dim();
}

void GetEnergies()
{
   if (energy == NULL) energy = new double[3];
   for (int i=0;i<3;++i)
      energy[i] = Lat[i]->E0();
}

void SetDOF(double *DOFs)
{
   static Vector D[3];
   for (int i=0;i<3;++i)
      D[i].Resize(dofs[i]);

   int start=0;
   for (int i=0;i<3;++i)
   {
      for (int j=0;j<dofs[i];++j)
         D[i][j] = DOFs[start + j];
      start = dofs[i];
      Lat[i]->SetDOF(D[i]);
   }
}

void quit()
{
   if (energy != NULL)
      delete [] energy;

   for (int i=0;i<3;++i)
   {
      delete Lat[i];
   }
}
