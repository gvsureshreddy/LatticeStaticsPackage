#include "UtilityFunctions.h"


void GetParameter(const char *tag,const char *datafile,const char *scanffmt,
		  void *parameter)
{
   char command[LINELENGTH];
   FILE *pipe;

   SetPerlCommand(command,datafile,tag);
   pipe = OpenPipe(command,"r");
   fscanf(pipe,scanffmt,parameter);
   if (pclose(pipe)) Errfun(tag);
}

int GetStringParameter(const char *tag,const char *datafile,
		       const char *choices[],const unsigned numb)
{
   int i;
   char strng[LINELENGTH];
   GetParameter(tag,datafile,"%s",strng);
   
   for (i=numb-1;i>=0;i--)
   {
      if (!strcasecmp(strng,choices[i]))
      {
	 return i;
      }
   }

   return i;
}

void SetPerlCommand(char *string,const char *datafile,const char *tag)
{
   char format[]=
     {"perl -e '$R=findref($ARGV[1],$ARGV[0]); print $R;"\
      "sub findref {my($tag,$df) = @_; my($fnd); $fnd=1; "\
      "open(R,$df); while (<R>) {if (/$tag/) {$fnd=0; "\
      "$_=deref($_,$df); split(\"=\",$_); return eval($_[1]);}} "\
      "close(R); if ($fnd == 1) {exit $fnd;}} sub deref "\
      "{my($fld,$df)=@_; my($t); while ($fld =~ m/<([^>]+)>/g) "\
      "{$t=$1; $v=findref(\"^$t\",$df); "\
      "$fld =~ s/<$t>/$v/} return $fld;}' %s \\%s"};
      sprintf(string,format,datafile,tag);
}

FILE *OpenPipe(const char *command,const char *mode)
{
   FILE *pipe;

   pipe=popen(command,mode);
   if (!pipe)
   {
      cerr << "popen failed! -- " << endl;
      exit(-1);
   }

   return pipe;
}

void Errfun(const char *string)
{
   cerr << "Error -- Unable to find : "
	<< string << endl;
   exit(-1);
}

//======================================================================
int IND3D(int i,int j);
int IND2D(int i,int j);

int IND3D(int i,int j)
{
   if (i==j)
      return i;
   else
      return 2+i+j;
}

int IND2D(int i,int j)
{
   if (i==j)
      return i;
   else
      return 1+i+j;
}

unsigned Rank1Convex3D(Matrix K,double dx)
{
   double Pi=4.0*atan(1.0);
   double_complex A[3][3][3];
   double n[2];
   
   for (double theta=0;theta<2.0*Pi;theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);

      for (int i=0;i<3;i++)
	 for (int j=0;j<3;j++)
	    for (int k=0;k<3;k++)
	       A[i][j][k] = double_complex(0.0,0.0);


      // Calculate Aij Polynomials
      for (int i=0;i<3;i++)
      {
	 for (int j=0;j<3;j++)
	 {
	    for (int k=0;k<2;k++)
	    {
	       for (int l=0;l<2;l++)
	       {
		  A[i][j][0] += K[IND3D(i,k)][IND3D(j,l)]*n[k]*n[l];
	       }
	       A[i][j][1] += K[IND3D(i,k)][IND3D(j,2)]*n[k]
		  + K[IND3D(i,2)][IND3D(j,k)]*n[k];
	    }
	    A[i][j][2] = K[IND3D(i,2)][IND3D(j,2)];
	 }
      }


      double_complex Roots[6];
      double_complex SolveMe[7], PA[5], PB[5];

      PolyRootsLaguerre(A[0][0],2,Roots,1);
      for (int i=0;i<2;i++) if (Roots[i].imag() == 0.0) return 0;
      
      // Only need to check the leading principal minors
      //PolyRootsLaguerre(A[1][1],2,Roots,1);
      //for (int i=0;i<2;i++) if (Roots[i].imag() == 0.0) return 0;

      // Only need to check the leading principal minors
      //PolyRootsLaguerre(A[2][2],2,Roots,1);
      //for (int i=0;i<2;i++) if (Roots[i].imag() == 0.0) return 0;

      PolyMult(A[0][0],2,A[1][1],2,PA);
      PolyMult(A[0][1],2,A[1][0],2,PB);
      for (int i=0;i<=4;i++) SolveMe[i] = PA[i] - PB[i];
      PolyRootsLaguerre(SolveMe,4,Roots,1);
      for (int i=0;i<4;i++) if (Roots[i].imag() == 0.0) return 0;

      // Only need to check the leading principal minors
      //PolyMult(A[0][0],2,A[2][2],2,PA);
      //PolyMult(A[0][2],2,A[2][0],2,PB);
      //for (int i=0;i<=4;i++) SolveMe[i] = PA[i] - PB[i];
      //PolyRootsLaguerre(SolveMe,4,Roots,1);
      //for (int i=0;i<4;i++) if (Roots[i].imag() == 0.0) return 0;
      
      // Only need to check the leading principal minors
      //PolyMult(A[1][1],2,A[2][2],2,PA);
      //PolyMult(A[1][2],2,A[2][1],2,PB);
      //for (int i=0;i<=4;i++) SolveMe[i] = PA[i] - PB[i];
      //PolyRootsLaguerre(SolveMe,4,Roots,1);
      //for (int i=0;i<4;i++) if (Roots[i].imag() == 0.0) return 0;

      double_complex DA[7],DB[7],DC[7];

      
      PolyMult(A[1][1],2,A[2][2],2,PA);
      PolyMult(A[1][2],2,A[2][1],2,PB);
      for (int i=0;i<=4;i++) PA[i] = PA[i] - PB[i];
      PolyMult(A[0][0],2,PA,4,DA);

      PolyMult(A[1][0],2,A[2][2],2,PA);
      PolyMult(A[1][2],2,A[2][0],2,PB);
      for (int i=0;i<=4;i++) PA[i] = PA[i] - PB[i];
      PolyMult(A[0][1],2,PA,4,DB);

      PolyMult(A[1][0],2,A[2][1],2,PA);
      PolyMult(A[1][1],2,A[2][0],2,PB);
      for (int i=0;i<=4;i++) PA[i] = PA[i] - PB[i];
      PolyMult(A[0][2],2,PA,4,DC);

      for (int i=0;i<=6;i++) SolveMe[i] = DA[i] - DB[i] + DC[i];

      PolyRootsLaguerre(SolveMe,6,Roots,1);
      for (int i=0;i<6;i++) if (Roots[i].imag() == 0.0) return 0;
   }

   return 1;
}

unsigned Rank1Convex2D(Matrix K,double dx)
{
   double Pi=4.0*atan(1.0);
   double A[2][2];
   double n[2];
   
   for (double theta=0;theta<2.0*Pi;theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);

      for (int i=0;i<2;i++)
	 for (int j=0;j<2;j++)
	    A[i][j] = 0.0;


      // Calculate Aij
      for (int i=0;i<2;i++)
      	 for (int j=0;j<2;j++)
	    for (int k=0;k<2;k++)
	       for (int l=0;l<2;l++)
	       {
		  A[i][j] += K[IND2D(i,k)][IND2D(j,l)]*n[k]*n[l];
	       }


      if (A[0][0] <= 0.0) return 0;
      // Only need to check the leading principal minors
      //if (A[1][1] <= 0.0) return 0;
      
      if ((A[0][0]*A[1][1] - A[0][1]*A[1][0]) <= 0.0) return 0;
   }
   
   return 1;
}
