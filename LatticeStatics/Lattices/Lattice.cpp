#include "Lattice.h"
#include "UtilityFunctions.h"


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

void Lattice::CriticalPointInfo(const Vector &DrDt,double Tolerance,
				char *datafile,const char *prefix,
				int Width,ostream &out)
{

   const char *thirdorderchoices[] = {"No","Yes"};
   int thirdorder=GetStringParameter(prefix,"CriticalPoint_T2",datafile,
				     thirdorderchoices,2);
   if (thirdorder < 0) exit(-1);
   
   Matrix 
      D3=E3(),
      D2=Stiffness(),
      D2T=StiffnessDT(),
      EigVec,
      EigVal=SymEigVal(D2,&EigVec);

   int dofs;
   
   if (DOFMAX < (dofs=D2.Rows()))
   {
      cerr << "Error: DOFMAX < " << dofs << " in GenericLat.h" << endl;
      exit(-5);
   }
   
   Matrix Mode;
   double
      Eijk[DOFMAX][DOFMAX][DOFMAX],
      EijT[DOFMAX][DOFMAX];

   // Find the modes
   int count = 0,
      Ind[DOFMAX];

   for (int i=0;i<dofs;i++)
   {
      Ind[i] = 0;
      if (fabs(EigVal[0][i]) < Tolerance)
      {
	 Ind[count++]=i;
      }
   }
   
   Mode.Resize(count,dofs);

   for (int i=0;i<count;i++)
   {
      for (int j=0;j<dofs;j++)
      {
	 Mode[i][j] = EigVec[j][Ind[i]];
      }
   }

   // Print out the Eigenvectors
   cout << "EigenVectors" << endl << setw(Width) << EigVec;
   out << "EigenVectors" << endl << setw(Width) << EigVec;

   // Eijk
   for (int i=0;i<count;i++)
      for (int j=0;j<count;j++)
	 for (int k=0;k<count;k++)
	 {
	    Eijk[i][j][k] = 0.0;
	    for (int a=0;a<dofs;a++)
	       for (int b=0;b<dofs;b++)
		  for (int c=0;c<dofs;c++)
		  {
		     Eijk[i][j][k] += D3[a*dofs + b][c]*Mode[i][a]*Mode[j][b]*Mode[k][c];
		  }
	 }
   
   //EijT
   for (int i=0;i<count;i++)
      for (int j=0;j<count;j++)
      {
	 EijT[i][j] = 0.0;
	 for (int a=0;a<dofs;a++)
	    for (int b=0;b<dofs;b++)
	    {
	       for (int c=0;c<dofs;c++)
	       {
		  EijT[i][j] += D3[a*dofs + b][c]*DrDt[c]*Mode[i][a]*Mode[j][b];
	       }
	       EijT[i][j] += D2T[a][b]*Mode[i][a]*Mode[j][b];
	    }
      }

   // Print out results
   for (int i=0;i<70;i++)
   {
      cout << "-";
      out << "-";
   }
   cout << endl << endl << "2nd Order Bifurcation Equations:" << endl;
   out << endl << endl << "2nd Order Bifurcation Equations:" << endl;

   int prec = out.precision();

   out << setiosflags(ios::scientific) << setprecision(prec/2);
   cout << setiosflags(ios::scientific) << setprecision(prec/2);

   for (int i=0;i<count;i++)
   {
      for (int j=0;j<count;j++)
	 for (int k=j;k<count;k++)
	 {
	    cout << "(" << setw(Width)
		 << Eijk[i][j][k] + ( (j==k)?0.0:Eijk[i][k][j] )
		 << ")a_" << j
		 << "a_"  << k
		 << " + ";
	    out << "(" << setw(Width)
		<< Eijk[i][j][k] + ( (j==k)?0.0:Eijk[i][k][j] )
		<< ")a_" << j
		<< "a_"  << k
		<< " + ";
	 }

      cout << "2T_1( ";
      out << "2T_1( ";
      
      for (int j=0;j<count-1;j++)
      {
	 cout << "(" << setw(Width)
	      << EijT[i][j]
	      << ")a_" << j
	      << " + ";
	 out << "(" << setw(Width)
	     << EijT[i][j]
	     << ")a_" << j
	     << " + ";
      }
      cout << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
	   << " ) = 0" << endl;
      out << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
	  << " ) = 0" << endl;
   }

   if (thirdorder)
   {
      double Eijkl[DOFMAX][DOFMAX][DOFMAX][DOFMAX],
	 Vij[DOFMAX][DOFMAX][DOFMAX];
      Matrix D4=E4(),
	 S(dofs-count,dofs);

      // Create projection operator
      for (int t=0,j=0,i=0;i<dofs;i++)
      {
	 if (Ind[t] == i)
	 {
	    t++;
	    continue;
	 }
	 else
	 {
	    for (int k=0;k<dofs;k++)
	    {
	       S[j][k] = EigVec[k][i];
	    }
	    j++;
	 }
      }

      Matrix Ainv=(S*D2*S.Transpose()).Inverse();
      Matrix b(dofs-count,1);

      for (int i=0;i<count;i++)
	 for (int j=0;j<count;j++)
	 {
	    for (int n=0;n<dofs-count;n++)
	    {
	       b[n][0] = 0.0;
	       for (int k=0;k<dofs;k++)
		  for (int l=0;l<dofs;l++)
		     for (int m=0;m<dofs;m++)
		     {
			b[n][0] += -S[n][k]*D3[k*dofs + l][m]*Mode[i][l]*Mode[j][m];
		     }
	    }

	    cout << endl << "V[" << i << "][" << j << "]=";
	    out << endl << "V[" << i << "][" << j << "]=";

	    for (int k=0;k<dofs;k++)
	    {
	       Vij[i][j][k] = 0.0;
	       for (int l=0;l<dofs-count;l++)
		  for (int m=0;m<dofs-count;m++)
		  {
		     Vij[i][j][k] += S[l][k]*Ainv[l][m]*b[m][0];
		  }
	       cout << setw(Width) << Vij[i][j][k];
	       out << setw(Width) << Vij[i][j][k];
	    }
	 }
      cout << endl;
      out << endl;

      // Eijkl
      for (int i=0;i<count;i++)
	 for (int j=0;j<count;j++)
	    for (int k=0;k<count;k++)
	       for (int l=0;l<count;l++)
	       {
		  Eijkl[i][j][k][l] = 0.0;
		  for (int m=0;m<dofs;m++)
		     for (int n=0;n<dofs;n++)
			for (int p=0;p<dofs;p++)
			   for (int q=0;q<dofs;q++)
			   {
			      Eijkl[i][j][k][l] += D4[m*dofs + n][p*dofs + q]
				 *Mode[i][m]*Mode[j][n]*Mode[k][p]*Mode[l][q];
			   }
		  for (int m=0;m<dofs;m++)
		     for (int n=0;n<dofs;n++)
			for (int p=0;p<dofs;p++)
			{
			   Eijkl[i][j][k][l] += D3[m*dofs + n][p]
			      *Mode[i][m]*(Mode[j][n]*Vij[k][l][p]
					    + Mode[k][n]*Vij[j][l][p]
					    + Mode[l][n]*Vij[j][k][p]);
			}
	       }

      // Print out results
      cout << endl << "3rd Order Bifurcation Equations:" << endl;
      out <<  endl << "3rd Order Bifurcation Equations:" << endl;
      
      for (int i=0;i<count;i++)
      {
	 for (int j=0;j<count;j++)
	    for (int k=j;k<count;k++)
	       for (int l=k;l<count;l++)
	       {
		  cout << "(" << setw(Width)
		       << Eijkl[i][j][k][l] + ( (j==k)?0.0:Eijkl[i][k][j][l] )
		     + ( (j==l)?0.0:Eijkl[i][l][k][j] ) + ( (k==l)?0.0:Eijkl[i][j][l][k] )
		       << ")a_" << j
		       << "a_"  << k
		       << "a_"  << l
		       << " + ";
		  out << "(" << setw(Width)
		      << Eijkl[i][j][k][l] + ( (j==k)?0.0:Eijkl[i][k][j][l] )
		     + ( (j==l)?0.0:Eijkl[i][l][k][j] ) + ( (k==l)?0.0:Eijkl[i][j][l][k] )
		      << ")a_" << j
		      << "a_"  << k
		      << "a_"  << l
		      << " + ";
	       }
	 
	 cout << "3T_2( ";
	 out << "3T_2( ";
	 
	 for (int j=0;j<count-1;j++)
	 {
	    cout << "(" << setw(Width)
		 << EijT[i][j]
		 << ")a_" << j
		 << " + ";
	    out << "(" << setw(Width)
		<< EijT[i][j]
		<< ")a_" << j
		<< " + ";
	 }
	 cout << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
	      << " ) = 0" << endl;
	 out << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
	     << " ) = 0" << endl;
      }      
   }

   cout << endl;
   out << endl;
   for (int i=0;i<70;i++)
   {
      cout << "-";
      out << "-";
   }
   cout << endl;
   out << endl;

   out << setiosflags(ios::fixed) << setprecision(prec);
   cout << setiosflags(ios::fixed) << setprecision(prec);
   
   return;
}
