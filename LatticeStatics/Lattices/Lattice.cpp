#include "Lattice.h"

int Lattice::StiffnessNulity(double *Min)
{
   int NoNegEigVal = 0;
   int index = 0;

   Matrix EigenValues=SymEigVal(Stiffness());
   int dofs=EigenValues.Cols();

   if (Min != NULL) *Min = fabs(EigenValues[0][0]);
   for (int i=0;i<dofs;i++)
   {
      if (EigenValues[0][i] < 0.0) NoNegEigVal++;
      if ((Min != NULL)
	  && (fabs(EigenValues[0][i]) < *Min))
      {
	 *Min = fabs(EigenValues[0][i]);
	 index = i;
      }
   }

   if (Min != NULL) *Min = EigenValues[0][index];
   return NoNegEigVal;
}
