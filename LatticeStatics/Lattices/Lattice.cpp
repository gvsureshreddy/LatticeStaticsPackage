#include "Lattice.h"
#include "UtilityFunctions.h"

void Lattice::SetLoadParameter(const double &load)
{
   if (LoadParameter_==Temperature)
   {
      SetTemp(load);
   }
   else
   {
      SetLambda(load);
   }
}

int Lattice::TestFunctions(Vector &TF1, StateType State , Vector *TF2)
{
   static int size =DOF().Dim();
   static Matrix Stiffness_1=E2();
   static Matrix Stiffness_2=Stiffness_1;
   static Matrix Stiffness_3(size, size);
   static Matrix Stiffness_temp(size, size);
   static Matrix Stiffness_diagonalized(size, size);
   static Matrix EigVect(size, size);
   static Matrix EV1(1,size),EV2(1,size);
   
   double sum;
   int NoNegEigVal = 0;
   int temp1, temp2;
   int Diff_NoNegEigVal;
   
   
   if (State == LHS)
   {
      Stiffness_1 = Stiffness_2;
      Stiffness_2 = E2();
      
      EV1 = SymEigVal(Stiffness_2);
      
      for (int i=0;i<size;i++)
      {
	 if (EV1[0][i] < 0.0) NoNegEigVal++;
	 TF1[i]=EV1[0][i];
      }
      
      return NoNegEigVal;
   }
   if (State == RHS)
   {
      if(TF2 == NULL)
      {
	 cerr << "Error in Lattice::TestFunctions(): TF2 == NULL" << endl;
	 exit(-53);
      }
      
      EV1 = SymEigVal(Stiffness_1, &EigVect);
      
      // Stiffness_diagonalized = EigVect.Transpose() * (Stiffness_2 * EigVect)
      //Stiffness_temp = Stiffness_2 * EigVect
      for (int i =0; i<size; i++)
      {
	 for (int j=0; j<size; j++)
	 {
	    sum = 0.0;
	    for(int k=0; k<size; k++)
	    {
	       sum += Stiffness_2[i][k] * EigVect[k][j];
	    }
	    Stiffness_temp[i][j] = sum;
	 }
      }
      //stiffness_diagonalized = Eigvect.Transpose() * Stiffness_temp
      for(int i=0; i<size; i++)
      {
	 for(int j=0; j<size; j++)
	 {
	    sum = 0.0;
	    for (int k=0; k<size; k++)
	    {
	       sum += EigVect[k][i] * Stiffness_temp[k][j];
	    }
	    Stiffness_diagonalized[i][j] = sum;
	 }
      }
      
      EV2 = SymEigVal(Stiffness_diagonalized);
      Diff_NoNegEigVal = 0;
      for (int i=0;i<size;i++)
      {
	 if ((EV1[0][i] * EV2[0][i]) < 0.0)
	 {
	    Diff_NoNegEigVal += 1;
	 }
	 TF1[i]=EV1[0][i];
	 (*TF2)[i]=EV2[0][i];
      }
      
      return Diff_NoNegEigVal;
   }
   if (State == CRITPT)
   {
      Stiffness_3 = E2();
      
      for (int i=0; i<size; i++)
      {
	 for (int j=0; j<size; j++)
	 {
	    sum = 0.0;
	    for(int k=0; k<size; k++)
	    {
	       sum += Stiffness_3[i][k] * EigVect[k][j];
	    }
	    Stiffness_temp[i][j] = sum;
	 }
      }
      //stiffness_diagonalized = Eigvect.Transpose() * Stiffness_temp
      for(int i=0; i<size; i++)
      {
	 for(int j=0; j<size; j++)
	 {
	    sum = 0.0;
	    for (int k=0; k<size; k++)
	    {
	       sum += EigVect[k][i] * Stiffness_temp[k][j];
	    }
	    Stiffness_diagonalized[i][j] = sum;
	 }
      }
      
      //cout << "STIFFNESS_DIAGONALIZED = " << endl << setw(15)
      //<< Stiffness_diagonalized << endl << endl;
      
      EV1 = SymEigVal(Stiffness_diagonalized);
      //EV1 = SymEigVal(Stiffness_diagonalized,&EigVect);
      
      for (int i=0;i<size;i++)
      {
	 if (EV1[0][i] < 0.0) NoNegEigVal++;
	 TF1[i]=EV1[0][i];
      }
      
      return NoNegEigVal;
   }
}

void Lattice::CriticalPointInfo(const Vector &DrDt,int NumZeroEigenVals,
				double Tolerance,char *datafile,const char *prefix,
				int Width,ostream &out)
{
   
   const char *thirdorderchoices[] = {"No","Yes"};
   int thirdorder=GetStringParameter(prefix,"CriticalPoint_T2",datafile,
				     thirdorderchoices,2);
   if (thirdorder < 0) exit(-1);
   
   Matrix
      D3=E3(),
      D2=E2(),
      D2T(D2.Rows(),D2.Cols()),
      D1T(1,D2.Cols()),
      EigVec,
      EigVal=SymEigVal(D2,&EigVec);
   if (LoadParameter_ == Temperature)
   {
      D1T=StressDT();
      D2T=StiffnessDT();
   }
   else if (LoadParameter_ == Load)
   {
      D1T=StressDL();
      D2T=StiffnessDL();
   }
   
   int dofs;
   
   if (DOFMAX < (dofs=D2.Rows()))
   {
      cerr << "Error: DOFMAX < " << dofs << " in Lattice.h" << endl;
      exit(-5);
   }
   
   Matrix Mode;
   double
      Eijk[BIFMAX][BIFMAX][BIFMAX],
      EijT[BIFMAX][BIFMAX];
   
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
   
   // Check for incorrect number of modes
   if (count != NumZeroEigenVals)
   {
      int skp;
      for (int j=count;j<NumZeroEigenVals;++j)
      {
	 Ind[count] = 0;
	 int a=0;
	 while ((Ind[count] == Ind[a]) && (a < j))
	 {
	    (Ind[count])++;
	    a++;
	 }
	 
	 for (int i=1;i<dofs;++i)
	 {
	    skp=0;
	    for (int k=0;k<count;++k)
	    {
	       if (Ind[k] == i) skp=1;
	    }
	    
	    if (!skp)
	    {
	       if (fabs(EigVal[0][i]) < fabs(EigVal[0][Ind[count]]))
		  Ind[count] = i;
	    }
	 }
	 count++;
      }
      
      out << "NOTE: Incorrect number of zero eigenvalues found. "
	  << "Modes with smallest abs. value used." << endl;
      if (Echo_)
      {
	 cout << "NOTE: Incorrect number of zero eigenvalues found. "
	      << "Modes with smallest abs. value used." << endl;
      }
   }
   
   for (int i=0;i<count;++i)
      out << "Mode[" << i << "] DOF: " << Ind[i] << ",  ";
   out << endl;
   if (Echo_)
   {
      for (int i=0;i<count;++i)
	 cout << "Mode[" << i << "] DOF: " << Ind[i] << ",  ";
      cout << endl;
   }
   
   if (BIFMAX < count)
   {
      cerr << "Error: BIFMAX < " << count << " in Lattice.h" << endl;
      exit(-6);
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
   out << "EigenVectors" << endl << setw(Width) << EigVec;
   if (Echo_) cout << "EigenVectors" << endl << setw(Width) << EigVec;
   
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
      out << "-";
      if (Echo_) cout << "-";
   }
   out << endl; if (Echo_) cout << endl;
   
   // Print out the critical point character test (Limit-load/Bifurcation)
   for (int i=0;i<count;++i)
   {
      double z=0.0;
      for (int j=0;j<dofs;++j)
      {
	 z+= Mode[i][j]*D1T[0][j];
      }
      out << "StressDT*Mode[" << i << "] = " << setw(Width) << z << endl;
      if (Echo_) cout << "StressDT*Mode[" << i << "] = " << setw(Width) << z << endl;
   }
   
   
   out << endl << endl << "2nd Order Bifurcation Equations:" << endl;
   if (Echo_) cout << endl << endl << "2nd Order Bifurcation Equations:" << endl;
   
   int prec = out.precision();
   
   out.flags(ios::scientific); out << setprecision(prec/2);
   if (Echo_) cout.flags(ios::scientific); cout << setprecision(prec/2);
   
   for (int i=0;i<count;i++)
   {
      for (int j=0;j<count;j++)
	 for (int k=j;k<count;k++)
	 {
	    if (Echo_) cout << "(" << setw(Width)
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
      
      out << "2T_1( ";
      if (Echo_) cout << "2T_1( ";
      
      for (int j=0;j<count-1;j++)
      {
	 if (Echo_) cout << "(" << setw(Width)
			 << EijT[i][j]
			 << ")a_" << j
			 << " + ";
	 out << "(" << setw(Width)
	     << EijT[i][j]
	     << ")a_" << j
	     << " + ";
      }
      if (Echo_) cout << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
		      << " ) = 0" << endl;
      out << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
	  << " ) = 0" << endl;
   }
   
   //------- output coefficients ----------
   out << endl << "DrDt = " << setw(Width) << DrDt << endl;
   if (Echo_) cout << endl << "DrDt = " << setw(Width) << DrDt << endl;
   out << "Eijk = " << endl;
   if (Echo_) cout << "Eijk = " << endl;
   for (int i=0;i<count;++i)
      for (int j=0;j<count;++j)
      {
	 for (int k=0;k<count;++k)
	 {
	    out << setw(Width) << Eijk[i][j][k];
	    if (Echo_) cout << setw(Width) << Eijk[i][j][k];
	 }
	 out << endl;
	 if (Echo_) cout << endl;
      }
   out << endl;
   if (Echo_) cout << endl;
   
   out << "EijT = " << endl;
   if (Echo_) cout << "EijT = " << endl;
   for (int i=0;i<count;++i)
   {
      for (int j=0;j<count;++j)
      {
	 out << setw(Width) << EijT[i][j];
	 if (Echo_) cout << setw(Width) << EijT[i][j];
      }
      out << endl;
      if (Echo_) cout << endl;
   }
   out << endl;
   if (Echo_) cout << endl;
   // ----------------------------
   
   
   if (thirdorder)
   {
      double Eijkl[BIFMAX][BIFMAX][BIFMAX][BIFMAX],
	 Vij[BIFMAX][BIFMAX][DOFMAX];
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
	    
	    if (Echo_) cout << endl << "V[" << i << "][" << j << "]=";
	    out << endl << "V[" << i << "][" << j << "]=";
	    
	    for (int k=0;k<dofs;k++)
	    {
	       Vij[i][j][k] = 0.0;
	       for (int l=0;l<dofs-count;l++)
		  for (int m=0;m<dofs-count;m++)
		  {
		     Vij[i][j][k] += S[l][k]*Ainv[l][m]*b[m][0];
		  }
	       if (Echo_) cout << setw(Width) << Vij[i][j][k];
	       out << setw(Width) << Vij[i][j][k];
	    }
	 }
      if (Echo_) cout << endl;
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
      if (Echo_) cout << endl << "3rd Order Bifurcation Equations:" << endl;
      out <<  endl << "3rd Order Bifurcation Equations:" << endl;
      
      for (int i=0;i<count;i++)
      {
	 for (int j=0;j<count;j++)
	    for (int k=j;k<count;k++)
	       for (int l=k;l<count;l++)
	       {
		  if (Echo_) cout << "(" << setw(Width)
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
	 
	 if (Echo_) cout << "3T_2( ";
	 out << "3T_2( ";
	 
	 for (int j=0;j<count-1;j++)
	 {
	    if (Echo_) cout << "(" << setw(Width)
			    << EijT[i][j]
			    << ")a_" << j
			    << " + ";
	    out << "(" << setw(Width)
		<< EijT[i][j]
		<< ")a_" << j
		<< " + ";
	 }
	 if (Echo_) cout << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
			 << " ) = 0" << endl;
	 out << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
	     << " ) = 0" << endl;
      }
      
      //-------- output coefficients -----------
      out << endl << "Eijkl = " << endl;
      if (Echo_) cout << endl << "Eijkl = " << endl;
      for (int i=0;i<count;++i)
	 for (int j=0;j<count;++j)
	 {
	    for (int k=0;k<count;++k)
	       for (int l=0;l<count;++l)
	       {
		  out << setw(Width) << Eijkl[i][j][k][l];
		  if (Echo_) cout << setw(Width) << Eijkl[i][j][k][l];
	       }
	    out << endl;
	    if (Echo_) cout << endl;
	 }
      out << endl;
      if (Echo_) cout << endl;
      // ----------------------------
   }
   
   out.flags(ios::fixed); out << setprecision(prec);
   if (Echo_) cout.flags(ios::fixed); cout << setprecision(prec);
   
   if (Echo_) cout << endl;
   out << endl;
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "-";
      out << "-";
   }
   if (Echo_) cout << endl;
   out << endl;
   
   
   if (dbg_)
   {
      if (EnterDebugMode())
      {
	 cout << setw(Width);
	 DebugMode();
      }
   }
   return;
}

void Lattice::ConsistencyCheck(double ConsistencyEpsilon,int Width,ostream &out)
{
   double potential;
   int Dim=DOF().Dim();
   Matrix
      Stiff(Dim,Dim),
      PerturbedStiff(Dim,Dim),
      D3(Dim*Dim,Dim),
      PerturbedD3(Dim*Dim,Dim),
      D4(Dim*Dim,Dim*Dim),
      PerturbedD4(Dim*Dim,Dim*Dim),
      Force(1,Dim),
      stress1(1,Dim),
      stress2(1,Dim),
      stiff1(Dim,Dim),
      stiff2(Dim,Dim),
      d31(Dim*Dim,Dim),
      d32(Dim*Dim,Dim);
   Vector
      PerturbedForce(Dim),
      OriginalState(Dim),
      pert(Dim,0.0);
   
   // Get current state
   OriginalState = DOF();
   
   out.flags(ios::scientific);
   if (Echo_) cout.flags(ios::scientific);
   
   // Do Consistency check
   int Do2=1,Do3=0,Do4=0;
   
   if (Echo_)
   {
      for (int i=0;i<70;i++) cout << "="; cout << endl;
      cout << "Consistency Check." << endl;
   }
   for (int i=0;i<70;i++) out << "="; out << endl;
   out << "Consistency Check." << endl;
   
   Force = ConsistencyEpsilon*E1();
   
   cout << "For E2 (1-yes,0-no) :"; cin >> Do2;
   cout << "For E3 (1-yes,0-no) :"; cin >> Do3;
   cout << "For E4 (1-yes,0-no) :"; cin >> Do4;
   
   if (Do2) Stiff = ConsistencyEpsilon*E2();
   if (Do3) D3 = ConsistencyEpsilon*E3();
   if (Do4) D4 = ConsistencyEpsilon*E4();
   
   potential = E0();
   if (Do2) stress1 = E1();
   if (Do3) stiff1 = E2();
   if (Do4) d31 = E3();
   for (int i=0;i<Dim;++i)
   {
      // Perturb the lattice state
      pert.Resize(Dim,0.0);
      pert[i]=1.0;
      SetDOF(OriginalState+ConsistencyEpsilon*pert);
      // Get Check
      PerturbedForce[i] = E0() - potential;
      if (Do2) stress2 = E1();
      if (Do3) stiff2 = E2();
      if (Do4) d32 = E3();
      for (int j=0;j<Dim;j++)
      {
	 if (Do2) PerturbedStiff[j][i] = stress2[0][j] - stress1[0][j];
	 
	 if (Do3 || Do4)
	 {
	    for (int k=0;k<Dim;++k)
	    {
	       if (Do3) PerturbedD3[k*Dim + j][i] = stiff2[k][j] - stiff1[k][j];
	       
	       if (Do4)
	       {
		  for (int l=0;l<Dim;++l)
		  {
		     PerturbedD4[l*Dim + k][j*Dim + i]
			= d32[l*Dim + k][j] - d31[l*Dim + k][j];
		  }
	       }
	    }
	 }
      }
   }
   // Set Lattice back to originalstate
   SetDOF(OriginalState);
   
   // Print out the facts
   if (Echo_)
   {
      cout << "E1(U) * Epsilon" << endl;
      cout << setw(Width) << Force << endl;
      cout << "E1(U + Epsilon*Vj)" << endl;
      cout << setw(Width) << PerturbedForce << endl;
      
      if (Do2)
      {
	 cout << "E2(U)*Epsilon" << endl;
	 cout << setw(Width) << Stiff << endl;
	 cout << "E2(U + Epsilon*Vj)" << endl;
	 cout << setw(Width) << PerturbedStiff << endl;
      }
      
      if (Do3)
      {
	 cout << "E3(U)*Epsilon" << endl;
	 cout << setw(Width) << D3 << endl;
	 cout << "E3(U + Epsilon*Vj)" << endl;
	 cout << setw(Width) << PerturbedD3 << endl;
      }
      
      if (Do4)
      {
	 cout << "E4(U)*Epsilon" << endl;
	 cout << setw(Width) << D4 << endl;
	 cout << "E4(U + Epsilon*Vj)" << endl;
	 cout << setw(Width) << PerturbedD4 << endl;
      }
      
      cout << "Difference" << endl;
      cout << setw(Width) << Force - PerturbedForce << endl << endl;
      if (Do2) cout << setw(Width) << Stiff - PerturbedStiff << endl;
      if (Do3) cout << setw(Width) << D3 - PerturbedD3 << endl;
      if (Do4) cout << setw(Width) << D4 - PerturbedD4 << endl;
   }
   
   out << "E1(U) * Epsilon" << endl;
   out << setw(Width) << Force << endl;
   out << "E1(U + Epsilon*Vj)" << endl;
   out << setw(Width) << PerturbedForce << endl;
   
   if (Do2)
   {
      out << "E2(U)*Epsilon" << endl;
      out << setw(Width) << Stiff << endl;
      out << "E2(U + Epsilon*Vj)" << endl;
      out << setw(Width) << PerturbedStiff << endl;
   }
   
   if (Do3)
   {
      out << "E3(U)*Epsilon" << endl;
      out << setw(Width) << D3 << endl;
      out << "E3(U + Epsilon*Vj)" << endl;
      out << setw(Width) << PerturbedD3 << endl;
   }
   
   if (Do4)
   {
      out << "E4(U)*Epsilon" << endl;
      out << setw(Width) << D4 << endl;
      out << "E4(U + Epsilon*Vj)" << endl;
      out << setw(Width) << PerturbedD4 << endl;
   }
   
   out << "Difference" << endl;
   out << setw(Width) << Force - PerturbedForce << endl << endl;
   if (Do2) out << setw(Width) << Stiff - PerturbedStiff << endl;
   if (Do3) out << setw(Width) << D3 - PerturbedD3 << endl;
   if (Do4) out << setw(Width) << D4 - PerturbedD4 << endl;
   
   out.flags(ios::fixed);
   if (Echo_) cout.flags(ios::fixed);
   
   
   if (Echo_)
   {
      for (int i=0;i<70;i++) cout << "="; cout << endl;
   }
   for (int i=0;i<70;i++) out << "="; out << endl;
}
