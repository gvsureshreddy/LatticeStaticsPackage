#include "UtilityFunctions.h"

#include <termios.h>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

/*
  The fcntl() function accepts a file descriptor and 
  an action (depending on the action, it may accept 
  a third parameter), and returns a result depending 
  on that action.  

  If things break bad, it returns -1.
*/

int setblock  (int file_desc, int block)
{
 int flags;

   /* retrieve the file descriptor's flags */
   if ( (flags = fcntl (file_desc, F_GETFL)) == -1 )
   {
     return false;    /* something went wrong! */
   }

   if (block)
   {
    flags &= ~O_NONBLOCK;      /* we want blocking input */
   }
   else
   {
    flags |= O_NONBLOCK;       /* we want non-blocking input */
   }

   /* set the flags (note the third parameter) */
   fcntl (file_desc, F_SETFL, flags);

   return 1;
}

int kbhitWait(void)
{
  struct termios old, newstate;
  int ch;

  tcgetattr(0, &old);

  newstate = old;
  newstate.c_lflag &= ~(ICANON);
  newstate.c_lflag &= ~(ECHO);
  newstate.c_cc[VTIME] = 1;
  newstate.c_cc[VMIN] = 1;

  tcsetattr(0, TCSANOW, &newstate);

  ch = getchar();

  tcsetattr(0, TCSANOW, &old);

  if(ch == EOF)
    return 0;
  return ch;
}

int EnterDebugMode()
{
   int n;
   char dbg[256];
   setblock(fileno(stdin),0);
   n=read(fileno(stdin),dbg,255);
   setblock(fileno(stdin),1);
   // remove newline
   if (n>0)
   {
      dbg[n-1]=0;
      
      if (!strcmp(dbg,"debug"))
	 return 1;
   }
   return 0;
}

int GetParameter(const char *prefix,const char *tag,const char *datafile,
		 const char *scanffmt,void *parameter,int DispErr)
{
   char command[LINELENGTH];
   FILE *pipe;

   SetPerlCommand(command,datafile,prefix,tag);
   pipe = OpenPipe(command,"r");
   fscanf(pipe,scanffmt,parameter);
   if (pclose(pipe))
   {
      if (DispErr) Errfun(tag);
      return 0;
   }
   else
   {
      return 1;
   }
}

int GetVectorParameter(const char *prefix,const char *tag,const char *datafile,Vector *V,
		       int DispErr)
{
   char command[LINELENGTH];
   FILE *pipe;

   SetPerlCommand(command,datafile,prefix,tag);
   pipe = OpenPipe(command,"r");
   for (int i=0;i<V->Dim();++i)
   {
      fscanf(pipe,"%lf",&((*V)[i]));
   }
   if (pclose(pipe))
   {
      if (DispErr) Errfun(tag);
      return 0;
   }
   else
   {
      return 1;
   }
}

int GetIntVectorParameter(const char *prefix,const char *tag,
			  const char *datafile,int N,int *Vec,int DispErr)
{
   char command[LINELENGTH];
   FILE *pipe;

   SetPerlCommand(command,datafile,prefix,tag);
   pipe = OpenPipe(command,"r");
   for (int i=0;i<N;++i)
   {
      fscanf(pipe,"%u",&(Vec[i]));
   }
   if (pclose(pipe))
   {
      if (DispErr) Errfun(tag);
      return 0;
   }
   else
   {
      return 1;
   }
}

int GetMatrixParameter(const char *prefix,const char *tag,const char *datafile,Matrix *M,
		       int DispErr)
{
   char command[LINELENGTH];
   FILE *pipe;

   SetPerlCommand(command,datafile,prefix,tag);
   pipe = OpenPipe(command,"r");
   for (int i=0;i<M->Rows();++i)
   {
      for (int j=0;M->Cols();++j)
      {
	 fscanf(pipe,"%lf",&((*M)[i][j]));
      }
   }
   if (pclose(pipe))
   {
      if (DispErr) Errfun(tag);
      return 0;
   }
   else
   {
      return 1;
   }
}

int GetStringParameter(const char *prefix,const char *tag,const char *datafile,
		       const char *choices[],const unsigned numb,int DispErr)
{
   int i;
   char strng[LINELENGTH];
   if (!GetParameter(prefix,tag,datafile,"%s",strng,DispErr)) return -1;
   for (i=numb-1;i>=0;i--)
   {
      if (!strcasecmp(strng,choices[i]))
      {
	 return i;
      }
   }

   return i;
}

void SetPerlCommand(char *string,const char *datafile,const char *prefix,const char *tag)
{
   char format[]=
     {"perl -e 'use Math::Trig;"\
      "@__R=__findref($ARGV[1],$ARGV[2],$ARGV[0]);"\
      "for $__i(0..@__R-1){print ($__R[$__i],\"\n\");};"\
      "sub __findref {my($__prfx,$__tag,$__df) = @_; my($__fnd,$__reg); $__fnd=1;"\
      "$__reg = $__prfx . $__tag;"\
      "open(__R,$__df); while (<__R>) {if (/$__reg/) {$__fnd=0; "\
      "$_=__deref($__prfx,$_,$__df); $_=~s/$__tag\\s*=\\s*//; return eval($_);}} "\
      "close(__R); if ($__fnd == 1) {exit $__fnd;}} sub __deref "\
      "{my($__prfx,$__fld,$__df)=@_; my($__t); while ($__fld =~ m/<([^>]+)>/g) "\
      "{$__t=$1; $__v=__findref($__prfx,\"$__t\",$__df); "\
      "$__fld =~ s/<$__t>/$__v/} return $__fld;}' %s '%s' '%s'"};
      sprintf(string,format,datafile,prefix,tag);
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

unsigned FullScanRank1Convex3D(Matrix K, double dx)
{
   Matrix A(3,3);
   Matrix Eigvals(1,3);
   double Pi=4.0*atan(1.0);
   double Piby2 = Pi/2,
      Pi2 = 2*Pi;
   double phi,theta;
   double n[3];

   for (phi = -Piby2;phi < Piby2;phi += dx)
   {
      for (theta = 0;theta < Pi2;theta += dx)
      {
	 n[0] = cos(phi)*cos(theta);
	 n[1] = cos(phi)*sin(theta);
	 n[2] = sin(phi);

	 // Acoustic tensor
	 for (int i=0;i<3;i++)
	    for (int j=0;j<3;j++)
	    {
	       A[i][j] = 0.0;
	       for (int k=0;k<3;k++)
		  for (int l=0;l<3;l++)
		  {
		     A[i][j] += K[IND3D(k,i)][IND3D(j,l)] * n[k] * n[l];
		  }
	    }

	 Eigvals = SymEigVal(A);
	 for(int i=0;i<3;i++)
	    if (Eigvals[0][i] <= 0.0) return 0;
      }
   }

   return 1;
}

unsigned FullScanRank1Convex2D(Matrix K, double dx)
{
   Matrix A(2,2);
   Matrix Eigvals(1,2);
   double Pi=4.0*atan(1.0);
   double Pi2 = 2*Pi;
   double theta;
   double n[2];

   for (theta = 0;theta < Pi2;theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);
      
      // Acoustic tensor
      for (int i=0;i<2;i++)
	 for (int j=0;j<2;j++)
	 {
	    A[i][j] = 0.0;
	    for (int k=0;k<2;k++)
	       for (int l=0;l<2;l++)
	       {
		  A[i][j] += K[IND2D(k,i)][IND2D(j,l)] * n[k] * n[l];
	       }
	 }
      
      Eigvals = SymEigVal(A);
      for(int i=0;i<2;i++)
	 if (Eigvals[0][i] <= 0.0) return 0;
   }
   
   return 1;
}

unsigned Rank1Convex3D(Matrix K,double dx)
{
   double Pi=4.0*atan(1.0);
   MyComplexDouble A[3][3][3];
   double n[2];
   
   for (double theta=0;theta<2.0*Pi;theta += dx)
   {
      n[0] = cos(theta);
      n[1] = sin(theta);

      for (int i=0;i<3;i++)
	 for (int j=0;j<3;j++)
	    for (int k=0;k<3;k++)
	       A[i][j][k] = MyComplexDouble(0.0,0.0);


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


      MyComplexDouble Roots[6];
      MyComplexDouble SolveMe[7], PA[5], PB[5];

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

      MyComplexDouble DA[7],DB[7],DC[7];

      
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
