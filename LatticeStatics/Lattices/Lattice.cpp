#include <fstream>
#include "Lattice.h"
#include "UtilityFunctions.h"

Lattice::Lattice(PerlInput const& Input,int const& Echo):
   Echo_(Echo),
   test_flag_static(0)
{
   if (Input.ParameterOK("Lattice","OrderedTFs"))
   {
      char const* const orderedtfs = Input.getString("Lattice","OrderedTFs");
      if (!strcmp("No",orderedtfs))
      {
         OrderedTFs_ = 0;
      }
      else if (!strcmp("Yes",orderedtfs))
      {
         OrderedTFs_ = 1;
      }
      else
      {
         cerr << "Error: Unknown value for Lattice{OrderedTFs}.\n";
         exit(-1);
      }
   }
   else
   {
      OrderedTFs_ = 0;
      Input.useString("No","Lattice","OrderedTFs"); // Default Value
   }
   
   if (Input.ParameterOK("Lattice","LSKAnalysis"))
   {
      char const* const lsk = Input.getString("Lattice","LSKAnalysis");
      if (!strcmp("ThirdOrder",lsk))
      {
         LSKAnalysis_ = 2;
      }
      else if (!strcmp("SecondOrder",lsk))
      {
         LSKAnalysis_ = 1;
      }
      else if (!strcmp("None",lsk))
      {
         LSKAnalysis_ = 0;
      }
      else
      {
         cerr << "Error: Unknown value for Lattice{LSKAnalysis}.\n";
         exit(-1);
      }
   }
   else
   {
      LSKAnalysis_ = 1;
      Input.useString("SecondOrder","Lattice","LSKAnalysis"); // Default Value
   }

   if (Input.ParameterOK("Lattice","FullPrint"))
   {
      char const* const fullprint = Input.getString("Lattice","FullPrint");
      if (!strcmp("Yes",fullprint))
      {
         FullPrint_ = 1;
      }
      else if (!strcmp("No",fullprint))
      {
         FullPrint_ = 0;
      }
      else
      {
         cerr << "Error: Unknown value for Lattice{FullPrint}.\n";
         exit(-1);
      }
   }
   else
   {
      FullPrint_ = 1;
      Input.useString("Yes","Lattice","FullPrint"); // Default Value
   }

   if (Input.ParameterOK("Lattice","GuessModes"))
   {
      char const* const guessmodes = Input.getString("Lattice","GuessModes");
      if (!strcmp("Yes",guessmodes))
      {
         GuessModes_ = 1;
      }
      else if (!strcmp("No",guessmodes))
      {
         GuessModes_ = 0;
      }
      else
      {
         cerr << "Error: Unknown value for Lattice{GuessModes}.\n";
         exit(-1);
      }
   }
   else
   {
      GuessModes_ = 0;
      Input.useString("No","Lattice","GuessModes"); // Default Value
   }

   if (Input.ParameterOK("Lattice","UseExtension"))
   {
      UseExtension_ = Input.getString("Lattice","UseExtension");
      UseExtension_ = "." + UseExtension_;
   }
   else
   {
      UseExtension_ = "";
   }
   
   Input.EndofInputSection();
}

void Lattice::SetLoadParameter(double const& load)
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

int Lattice::TestFunctions(Vector &TF1,StateType const& State,Vector* const TF2) const
{
   double sum;
   int NoNegEigVal = 0;
   int Diff_NoNegEigVal = 0;
   int retval = 0;
   
   int size=DOF().Dim();
   if(test_flag_static == 0)
   {
      Stiffness_1_static.Resize(size,size);
      Stiffness_1_static=E2();
      
      Stiffness_2_static.Resize(size,size);
      Stiffness_2_static=Stiffness_1_static;
      
      Stiffness_3_static.Resize(size, size);
      Stiffness_temp_static.Resize(size, size);
      Stiffness_diagonalized_static.Resize(size, size);
      EigVect_static.Resize(size, size);
      EigVectRHS_static.Resize(size, size);
      EigVectLHS_static.SetIdentity(size);
      EV1_static.Resize(1,size);
      EV2_static.Resize(1,size);
      
      test_flag_static = test_flag_static+1;
   }
   
   if(OrderedTFs_ == 1)
   {
      if (State == LHS)
      {
         Stiffness_2_static = E2();
         
         //Stiffness_diagonalized = EigVect.Transpose() * (Stiffness_2 * EigVect)
         //Stiffness_temp = Stiffness_2 * EigVect
         for (int i =0; i<size; i++)
         {
            for (int j=0; j<size; j++)
            {
               sum = 0.0;
               for(int k=0; k<size; k++)
               {
                  sum += Stiffness_2_static[i][k] * EigVectLHS_static[k][j];
               }
               Stiffness_temp_static[i][j] = sum;
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
                  sum += EigVectLHS_static[k][i] * Stiffness_temp_static[k][j];
               }
               Stiffness_diagonalized_static[i][j] = sum;
            }
         }
         
         for(int i=0;i<size;i++)
         {
            EV1_static[0][i] = EV2_static[0][i];
            for (int j =0;j<size;j++)
            {
               EigVectRHS_static[i][j] = EigVectLHS_static[i][j];
            }
            
         }
         
         EV2_static = SymEigVal(Stiffness_diagonalized_static, &EigVect_static);
         
         for(int i=0; i<size; i++)
         {
            for(int j=0; j<size; j++)
            {
               sum = 0.0;
               for (int k=0; k<size; k++)
               {
                  sum += EigVectRHS_static[i][k] * EigVect_static[k][j];
               }
               EigVectLHS_static[i][j] = sum;
            }
         }
         
         for (int i=0;i<size;i++)
         {
            if (EV2_static[0][i] < 0.0) NoNegEigVal++;
            TF1[i]=EV2_static[0][i];
         }
         
         retval = NoNegEigVal;
      }
      if (State == RHS)
      {
         if(TF2 == 0)
         {
            cerr << "Error in Lattice::TestFunctions(): TF2 == 0" << "\n";
            exit(-53);
         }
         
         Diff_NoNegEigVal = 0;
         NoNegEigVal = 0;
         for (int i=0;i<size;i++)
         {
            NoNegEigVal += ((EV1_static[0][i] < 0.0) ? 1 : 0);
            NoNegEigVal -= ((EV2_static[0][i] < 0.0) ? 1 : 0);
            if ((EV1_static[0][i] * EV2_static[0][i]) < 0.0)
            {
               Diff_NoNegEigVal += 1;
            }
            TF1[i]=EV1_static[0][i];
            (*TF2)[i]=EV2_static[0][i];
         }
         NoNegEigVal = (NoNegEigVal < 0) ? -NoNegEigVal : NoNegEigVal;
         
         if (Diff_NoNegEigVal != NoNegEigVal)
         {
            // return a negative number to indicate the error
            Diff_NoNegEigVal = -Diff_NoNegEigVal;
         }
         
         retval = Diff_NoNegEigVal;
      }
      if (State == CRITPT)
      {
         Stiffness_3_static = E2();
         
         for (int i=0; i<size; i++)
         {
            for (int j=0; j<size; j++)
            {
               sum = 0.0;
               for(int k=0; k<size; k++)
               {
                  sum += Stiffness_3_static[i][k] * EigVectLHS_static[k][j];
               }
               Stiffness_temp_static[i][j] = sum;
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
                  sum += EigVectLHS_static[k][i] * Stiffness_temp_static[k][j];
               }
               Stiffness_diagonalized_static[i][j] = sum;
            }
         }
         
         EV1_static = SymEigVal(Stiffness_diagonalized_static);
         //EV1 = SymEigVal(Stiffness_diagonalized,&EigVectLHS);
         
         for (int i=0;i<size;i++)
         {
            if (EV1_static[0][i] < 0.0) NoNegEigVal++;
            TF1[i]=EV1_static[0][i];
         }
         
         retval = NoNegEigVal;
      }
   }
   else
   {
      if (State == LHS)
      {
         Stiffness_1_static = Stiffness_2_static;
         Stiffness_2_static = E2();
         
         EV1_static = SymEigVal(Stiffness_2_static);
         
         for (int i=0;i<size;i++)
         {
            if (EV1_static[0][i] < 0.0) NoNegEigVal++;
            TF1[i]=EV1_static[0][i];
         }
         
         retval = NoNegEigVal;
      }
      if (State == RHS)
      {
         if(TF2 == 0)
         {
            cerr << "Error in Lattice::TestFunctions(): TF2 == 0" << "\n";
            exit(-53);
         }
         
         EV1_static = SymEigVal(Stiffness_1_static, &EigVect_static);
         
         // Stiffness_diagonalized = EigVect.Transpose() * (Stiffness_2 * EigVect)
         //Stiffness_temp = Stiffness_2 * EigVect
         for (int i =0; i<size; i++)
         {
            for (int j=0; j<size; j++)
            {
               sum = 0.0;
               for(int k=0; k<size; k++)
               {
                  sum += Stiffness_2_static[i][k] * EigVect_static[k][j];
               }
               Stiffness_temp_static[i][j] = sum;
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
                  sum += EigVect_static[k][i] * Stiffness_temp_static[k][j];
               }
               Stiffness_diagonalized_static[i][j] = sum;
            }
         }
         
         EV2_static = SymEigVal(Stiffness_diagonalized_static);
         Diff_NoNegEigVal = 0;
         NoNegEigVal = 0;
         for (int i=0;i<size;i++)
         {
            NoNegEigVal += ((EV1_static[0][i] < 0.0) ? 1 : 0);
            NoNegEigVal -= ((EV2_static[0][i] < 0.0) ? 1 : 0);
            if ((EV1_static[0][i] * EV2_static[0][i]) < 0.0)
            {
               Diff_NoNegEigVal += 1;
            }
            TF1[i]=EV1_static[0][i];
            (*TF2)[i]=EV2_static[0][i];
         }
         NoNegEigVal = (NoNegEigVal < 0) ? -NoNegEigVal : NoNegEigVal;
         
         if (Diff_NoNegEigVal != NoNegEigVal)
         {
            // return a negative number to indicate the error
            Diff_NoNegEigVal = -Diff_NoNegEigVal;
         }
         
         retval = Diff_NoNegEigVal;
      }
      if (State == CRITPT)
      {
         Stiffness_3_static = E2();
         
         for (int i=0; i<size; i++)
         {
            for (int j=0; j<size; j++)
            {
               sum = 0.0;
               for(int k=0; k<size; k++)
               {
                  sum += Stiffness_3_static[i][k] * EigVect_static[k][j];
               }
               Stiffness_temp_static[i][j] = sum;
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
                  sum += EigVect_static[k][i] * Stiffness_temp_static[k][j];
               }
               Stiffness_diagonalized_static[i][j] = sum;
            }
         }
         
         //cout << "STIFFNESS_DIAGONALIZED = " << "\n" << setw(15)
         //<< Stiffness_diagonalized << "\n" << "\n";
         
         EV1_static = SymEigVal(Stiffness_diagonalized_static);
         //EV1 = SymEigVal(Stiffness_diagonalized,&EigVect);
         
         for (int i=0;i<size;i++)
         {
            if (EV1_static[0][i] < 0.0) NoNegEigVal++;
            TF1[i]=EV1_static[0][i];
         }
         
         retval = NoNegEigVal;
      }
   }
   
   return retval;
}

int Lattice::CriticalPointInfo(int const& CPCrossingNum,Vector const& DrDt,int const& CPorBif,
                               int const& NumZeroEigenVals,double const& Tolerance,
                               int const& Width,PerlInput const& Input,ostream& out)
{
   Matrix
      D2=E2(),
      EigVec,
      EigVal=SymEigVal(D2,&EigVec);
   Vector D1T(D2.Cols());

   // default to the passed in information
   int Bif = CPorBif;
   
   if (LoadParameter_ == Temperature)
   {
      D1T=StressDT();
   }
   else if (LoadParameter_ == Load)
   {
      D1T=StressDL();
   }
   
   int dofs = D2.Rows();
   
   out << "Critical Point Crossing Number: " << CPCrossingNum << "\n";
   if (Echo_)
   {
      cout << "Critical Point Crossing Number: " << CPCrossingNum << "\n";
   }
   
   // Find the modes
   int count = 0;
   int* Ind = new int[dofs];
   
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
          << "Modes with smallest abs. value used." << "\n";
      if (Echo_)
      {
         cout << "NOTE: Incorrect number of zero eigenvalues found. "
              << "Modes with smallest abs. value used." << "\n";
      }
   }
   
   for (int i=0;i<count;++i)
      out << "Mode[" << i << "] DOF: " << Ind[i] << ",  ";
   out << "\n";
   if (Echo_)
   {
      for (int i=0;i<count;++i)
         cout << "Mode[" << i << "] DOF: " << Ind[i] << ",  ";
      cout << "\n";
   }

   Matrix Mode(count,dofs);
   
   for (int i=0;i<count;i++)
   {
      for (int j=0;j<dofs;j++)
      {
         Mode[i][j] = EigVec[j][Ind[i]];
      }
   }
   
   // Print out the Eigenvalues and Eigenvectors
   if (FullPrint_)
   {
      out << "EigenValues" << "\n" << setw(Width) << EigVal;
      if (Echo_) cout << "EigenValues" << "\n" << setw(Width) << EigVal;
      out << "EigenVectors" << "\n" << setw(Width) << EigVec;
      if (Echo_) cout << "EigenVectors" << "\n" << setw(Width) << EigVec;
   }
   else
   {
      out << "EigenValues: ";
      if (Echo_) cout << "EigenValues: ";
      for (int i=0;i<count;++i)
      {
         out << setw(Width) << EigVal[0][Ind[i]];
         if (Echo_) cout << setw(Width) << EigVal[0][Ind[i]];
      }
      out << "\n";
      if (Echo_) cout << "\n";
   }

   // Print out results
   for (int i=0;i<70;i++)
   {
      out << "-";
      if (Echo_) cout << "-";
   }
   out << "\n"; if (Echo_) cout << "\n";
   
   // Print out the critical point character test (Limit-load/Bifurcation)
   Bif = 1;
   for (int i=0;i<count;++i)
   {
      double z=0.0;
      for (int j=0;j<dofs;++j)
      {
         z+= Mode[i][j]*D1T[j];
      }
      if (fabs(z) > Tolerance) Bif = 0;
      out << "StressDT*Mode[" << i << "] = " << setw(Width) << z << "\n";
      if (Echo_) cout << "StressDT*Mode[" << i << "] = " << setw(Width) << z << "\n";
   }
   if (Bif != CPorBif) // information conflicts. use CPorBif
   {
      out << "NOTE: Conflict between critical point identification methods.\n"
          << "      Using characterization provided to CriticalPointInfo()." << "\n";
      if (Echo_)
      {
         cout << "NOTE: Conflict between critical point identification methods.\n"
              << "      Using characterization provided to CriticalPointInfo()." << "\n";
      }
      Bif = CPorBif;
   }
   
   if (LSKAnalysis_ > 0)
   {
      Matrix D3=E3();
      Matrix D2T(D2.Rows(),D2.Cols());
      if (LoadParameter_ == Temperature)
      {
         D2T=StiffnessDT();
      }
      else if (LoadParameter_ == Load)
      {
         D2T=StiffnessDL();
      }
   
      Matrix* Eijk = new Matrix[count];
      for (int i=0;i<count;++i) Eijk[i].Resize(count,count);
      Matrix EijT(count,count);
      
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
      out << "\n" << "\n" << "2nd Order Bifurcation Equations:" << "\n";
      if (Echo_) cout << "\n" << "\n" << "2nd Order Bifurcation Equations:" << "\n";
      
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
                         << " ) = 0" << "\n";
         out << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
             << " ) = 0" << "\n";
      }
      
      //------- output coefficients ----------
      if (FullPrint_)
      {
         out << "\n" << "DrDt = " << setw(Width) << DrDt << "\n";
         if (Echo_) cout << "\n" << "DrDt = " << setw(Width) << DrDt << "\n";
         out << "Eijk = " << "\n";
         if (Echo_) cout << "Eijk = " << "\n";
         for (int i=0;i<count;++i)
            for (int j=0;j<count;++j)
            {
               for (int k=0;k<count;++k)
               {
                  out << setw(Width) << Eijk[i][j][k];
                  if (Echo_) cout << setw(Width) << Eijk[i][j][k];
               }
               out << "\n";
               if (Echo_) cout << "\n";
            }
         out << "\n";
         if (Echo_) cout << "\n";
         
         out << "EijT = " << "\n";
         if (Echo_) cout << "EijT = " << "\n";
         for (int i=0;i<count;++i)
         {
            for (int j=0;j<count;++j)
            {
               out << setw(Width) << EijT[i][j];
               if (Echo_) cout << setw(Width) << EijT[i][j];
            }
            out << "\n";
            if (Echo_) cout << "\n";
         }
         out << "\n";
         if (Echo_) cout << "\n";
      }
      // ----------------------------
      
      if (LSKAnalysis_ > 1)
      {
         Matrix** Eijkl = new Matrix*[count];
         Eijkl[0] = new Matrix[count*count];
         for (int i=1;i<count;++i)
         {
            Eijkl[i] = Eijkl[i-1] + count;
         }
         for (int i=0;i<count;++i)
            for (int j=0;j<count;++j)
            {
               Eijkl[i][j].Resize(count,count);
            }
         Matrix* Vij;
         Vij = new Matrix[count];
         for (int i=0;i<count;++i) Vij[i].Resize(count,dofs);
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

               if (FullPrint_)
               {
                  if (Echo_) cout << "\n" << "V[" << i << "][" << j << "]=";
                  out << "\n" << "V[" << i << "][" << j << "]=";
               }
               
               for (int k=0;k<dofs;k++)
               {
                  Vij[i][j][k] = 0.0;
                  for (int l=0;l<dofs-count;l++)
                     for (int m=0;m<dofs-count;m++)
                     {
                        Vij[i][j][k] += S[l][k]*Ainv[l][m]*b[m][0];
                     }
                  if (FullPrint_)
                  {
                     if (Echo_) cout << setw(Width) << Vij[i][j][k];
                     out << setw(Width) << Vij[i][j][k];
                  }
               }
            }
         if (FullPrint_)
         {
            if (Echo_) cout << "\n";
            out << "\n";
         }
         
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
         if (Echo_) cout << "\n" << "3rd Order Bifurcation Equations:" << "\n";
         out <<  "\n" << "3rd Order Bifurcation Equations:" << "\n";
         
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
                            << " ) = 0" << "\n";
            out << "(" << setw(Width) << EijT[i][count-1] << ")a_" << count-1
                << " ) = 0" << "\n";
         }
         
         //-------- output coefficients -----------
         if (FullPrint_)
         {
            out << "\n" << "Eijkl = " << "\n";
            if (Echo_) cout << "\n" << "Eijkl = " << "\n";
            for (int i=0;i<count;++i)
               for (int j=0;j<count;++j)
               {
                  for (int k=0;k<count;++k)
                     for (int l=0;l<count;++l)
                     {
                        out << setw(Width) << Eijkl[i][j][k][l];
                        if (Echo_) cout << setw(Width) << Eijkl[i][j][k][l];
                     }
                  out << "\n";
                  if (Echo_) cout << "\n";
               }
            out << "\n";
            if (Echo_) cout << "\n";
         }
         // ----------------------------
         
         // release memory
         delete [] Eijkl[0];
         delete [] Eijkl;
         delete [] Vij;
      }

      // release memory
      delete [] Eijk;

      out.flags(ios::fixed); out << setprecision(prec);
      if (Echo_) cout.flags(ios::fixed); cout << setprecision(prec);
   }

   // release memory
   delete [] Ind;
   
   for (int i=0;i<70;i++)
   {
      if (Echo_) cout << "-";
      out << "-";
   }
   if (Echo_) cout << endl;
   out << endl;

   // output a new input file to help restart at this critical point
   ostringstream cpfilename;
   fstream cpfile;
   cpfilename << Input.LastInputFileName();
   if ("" != UseExtension_)
   {
      unsigned pos = (cpfilename.str().rfind(UseExtension_,cpfilename.str().length()-1));
      if (string::npos != pos)
      {
         string a = cpfilename.str().substr(0,pos);
         cpfilename.str("");
         cpfilename << a;
      }
   }
   if (1 == Bif)
      cpfilename << ".BP.";
   else
      cpfilename << ".TP.";

   cpfilename << setw(2) << setfill('0') << CPCrossingNum << UseExtension_;
   cpfile.open(cpfilename.str().c_str(),ios::out);

   cpfile << setprecision(out.precision()) << scientific;
   Vector T(dofs+1);
   Vector M(dofs+1);
   Vector dof=DOF();
   for (int i=0;i<dofs;++i)
   {
      T[i] = dof[i];
   }
   T[dofs] = ((LoadParameter_ == Temperature) ? Temp() : Lambda());
   if (Bif > 0)
   {
      cpfile << Input.ReconstructedInput();
      cpfile << "\n\n";

      if (GuessModes_)
      {
         int** colms;
         colms = new int*[dofs];
         colms[0] = new int[dofs*dofs];
         for (int i=1;i<dofs;++i) colms[i] = colms[i-1] + dofs;
         int** colmssgn;
         colmssgn = new int*[dofs];
         colmssgn[0] = new int[dofs*dofs];
         for (int i=1;i<dofs;++i) colmssgn[i] = colmssgn[i-1] + dofs;
         int* colmscnt = new int[dofs];
         
         for (int i=0;i<count;++i)
         {
            int cnt=0;
            int foundflg;
            // make a guess at the new Restriction DOFS
            for (int j=0;j<dofs;++j)
            {
               colmscnt[j] = 0;
               for (int k=0;k<dofs;++k)
                  colms[k][j] = colmssgn[k][j] = 0;
            }
            for (int j=0;j<dofs;++j)
            {
               if ((Mode[i][j] > 20*Tolerance) || (Mode[i][j] < -20*Tolerance)) // M[j] != 0.0
               {
                  foundflg = 0;
                  for (int k=0;k<cnt;++k)
                  {
                     if ((fabs(Mode[i][j]) > fabs(Mode[i][colms[k][0]]) - 20*Tolerance) &&
                         (fabs(Mode[i][j]) < fabs(Mode[i][colms[k][0]]) + 20*Tolerance))
                     {
                        colms[k][colmscnt[k]] = j;
                        colmssgn[k][colmscnt[k]]
                           = int(round(Mode[i][j]/fabs(Mode[i][colms[k][0]])));
                        ++colmscnt[k];
                        foundflg = 1;
                     }
                  }
                  
                  if (!foundflg)
                  {
                     colms[cnt][0] = j;
                     colmssgn[cnt][0] = int(round(Mode[i][j]/fabs(Mode[i][j])));
                     ++colmscnt[cnt];
                     ++cnt;
                  }
               }
            }
            
            for (int j=0;j<cnt;++j)
            {
               Input.writePosIntVector(cpfile,&(colms[j][0]),colmscnt[j],
                                       Input.useHash("Restriction",
                                                     "RestrictToTranslatedSubSpace"),
                                       "DOF_?",0);
               Vector z(colmscnt[j]);
               for (int k=0;k<colmscnt[j];++k) z[k] = colmssgn[j][k];
               Input.writeVector(cpfile,z,
                                 Input.useHash("Restriction","RestrictToTranslatedSubSpace"),
                                 "DOF_?",1);
            }
         }
         cpfile << "\n";

         delete [] colms[0];
         delete [] colms;
         delete [] colmssgn[0];
         delete [] colmssgn;
         delete [] colmscnt;
      }
      
      Input.writeString(cpfile,"Bifurcation","StartType","Type");
      for (int i=0;i<count;++i)
      {
         for (int j=0;j<dofs;++j)
         {
            M[j] = Mode[i][j];
         }
         M[dofs] = 0.0;
         Input.writeVector(cpfile,M,"StartType","Tangent");
      }
      Input.writeVector(cpfile,T,"StartType","BifurcationPoint");
   }
   else
   {
      Input.writeString(cpfile,"Continuation","StartType","Type");
      for (int i=0;i<count;++i) // be safe, cover a miss-identified bifurcation point
      {
         for (int j=0;j<dofs;++j)
         {
            M[j] = Mode[i][j];
         }
         M[dofs] = 0.0;
         Input.writeVector(cpfile,M,"StartType","Tangent");
      }
      Input.writeVector(cpfile,T,"StartType","Solution");
   }

   cpfile.close();
   
   if (dbg_)
   {
      if (EnterDebugMode())
      {
         cout << setw(Width);
         DebugMode();
      }
   }
   
   return Bif;
}

void Lattice::ConsistencyCheck(double const& ConsistencyEpsilon,int const& Width,ostream& out)
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
      stiff1(Dim,Dim),
      stiff2(Dim,Dim),
      d31(Dim*Dim,Dim),
      d32(Dim*Dim,Dim);
   Vector
      Force(Dim),
      stress1(Dim),
      stress2(Dim),
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
      for (int i=0;i<70;i++) cout << "="; cout << "\n";
      cout << "Consistency Check." << "\n";
   }
   for (int i=0;i<70;i++) out << "="; out << "\n";
   out << "Consistency Check." << "\n";
   
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
         if (Do2) PerturbedStiff[j][i] = stress2[j] - stress1[j];
         
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
      cout << "E1(U) * Epsilon" << "\n";
      cout << setw(Width) << Force << "\n";
      cout << "E1(U + Epsilon*Vj)" << "\n";
      cout << setw(Width) << PerturbedForce << "\n";
      
      if (Do2)
      {
         cout << "E2(U)*Epsilon" << "\n";
         cout << setw(Width) << Stiff << "\n";
         cout << "E2(U + Epsilon*Vj)" << "\n";
         cout << setw(Width) << PerturbedStiff << "\n";
      }
      
      if (Do3)
      {
         cout << "E3(U)*Epsilon" << "\n";
         cout << setw(Width) << D3 << "\n";
         cout << "E3(U + Epsilon*Vj)" << "\n";
         cout << setw(Width) << PerturbedD3 << "\n";
      }
      
      if (Do4)
      {
         cout << "E4(U)*Epsilon" << "\n";
         cout << setw(Width) << D4 << "\n";
         cout << "E4(U + Epsilon*Vj)" << "\n";
         cout << setw(Width) << PerturbedD4 << "\n";
      }
      
      cout << "Difference" << "\n";
      cout << setw(Width) << Force - PerturbedForce << "\n" << "\n";
      if (Do2) cout << setw(Width) << Stiff - PerturbedStiff << "\n";
      if (Do3) cout << setw(Width) << D3 - PerturbedD3 << "\n";
      if (Do4) cout << setw(Width) << D4 - PerturbedD4 << "\n";
   }
   
   out << "E1(U) * Epsilon" << "\n";
   out << setw(Width) << Force << "\n";
   out << "E1(U + Epsilon*Vj)" << "\n";
   out << setw(Width) << PerturbedForce << "\n";
   
   if (Do2)
   {
      out << "E2(U)*Epsilon" << "\n";
      out << setw(Width) << Stiff << "\n";
      out << "E2(U + Epsilon*Vj)" << "\n";
      out << setw(Width) << PerturbedStiff << "\n";
   }
   
   if (Do3)
   {
      out << "E3(U)*Epsilon" << "\n";
      out << setw(Width) << D3 << "\n";
      out << "E3(U + Epsilon*Vj)" << "\n";
      out << setw(Width) << PerturbedD3 << "\n";
   }
   
   if (Do4)
   {
      out << "E4(U)*Epsilon" << "\n";
      out << setw(Width) << D4 << "\n";
      out << "E4(U + Epsilon*Vj)" << "\n";
      out << setw(Width) << PerturbedD4 << "\n";
   }
   
   out << "Difference" << "\n";
   out << setw(Width) << Force - PerturbedForce << "\n" << "\n";
   if (Do2) out << setw(Width) << Stiff - PerturbedStiff << "\n";
   if (Do3) out << setw(Width) << D3 - PerturbedD3 << "\n";
   if (Do4) out << setw(Width) << D4 - PerturbedD4 << "\n";
   
   out.flags(ios::fixed);
   if (Echo_) cout.flags(ios::fixed);
   
   
   if (Echo_)
   {
      for (int i=0;i<70;i++) cout << "="; cout << endl;
   }
   for (int i=0;i<70;i++) out << "="; out << endl;
}
