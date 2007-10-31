#include "NewtonPCSolution.h"
#include "Matrix.h"
#include "UtilityFunctions.h"
#include <cmath>

using namespace std;

#define CLOSEDDEFAULT 30


NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   const Vector &one,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{
   // get needed parameters
	if(!GetParameter(prefix,"PCNumSolutions",datafile,'u',&NumSolutions_)) exit(-1);
	if(!GetParameter(prefix,"PCStepLength",datafile,'l',&CurrentDS_)) exit(-1);
	if(!GetParameter(prefix,"PCContraction",datafile,'l',&cont_rate_nom_)) exit(-1);
	if(!GetParameter(prefix,"PCDistance",datafile,'l',&delta_nom_)) exit(-1);
   	if(!GetParameter(prefix,"PCAngle",datafile,'l',&alpha_nom_)) exit(-1);
	if(!GetParameter(prefix,"PCConvergeCriteria",datafile,'l',&Converge_)) exit(-1);
	if(!GetParameter(prefix,"PCClosedLoopStart",datafile,'u',&ClosedLoopStart_,0))
	
	{
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
	}

	
	FirstSolution_.Resize(one.Dim());
	FirstSolution_ = one;
	
	//cout << "size of one = " << one.Dim() << endl;
	
	// Set Lattice to solution "one"
	Mode_->SetModeDOF(one);	
	
	Stiffness_.Resize(Mode_->ModeStiffness().Rows(), Mode_->ModeStiffness().Cols());
	//cout << "size of stiffness rows= " << Stiffness_.Rows() << endl;
	//cout << "size of Mode_->ModeStiffness().Rows_= " << Mode_->ModeStiffness().Rows() << endl << endl;
	//cout << "size of stiffness columns= " << Stiffness_.Cols() << endl; 
	//cout << "size of Mode_->ModeStiffness().Cols_= " << Mode_->ModeStiffness().Cols() << endl;
   
	// set tangent. Force1_ might not be needed
	Stiffness_ = Mode_->ModeStiffness();
	int count = Mode_->ModeStiffness().Cols();
	int count_minus_one = Mode_->ModeStiffness().Rows();
	
	
	//QR Decomposition of Stiffness Matrix
	
	Matrix Q(count, count);
	Matrix R(count, count_minus_one);
			
			
	QR(Stiffness_, Q, R, 1);	//Performs QR decomposition using A^T = Q*R. Section 4.1 of ISBN 3-540-12760-7
	
	
	
	//Assigns Tangent Vector
	//Last column of Q is the tangent vector
	
	
	Tangent1_.Resize(count);
	Tangent2_.Resize(count);
		
	
	for(int i=0;i< count;i++)
	{
		Tangent1_[i] = Q[i][count_minus_one];
		Tangent2_[i] = Q[i][count_minus_one];
	}
	

	cout << "Stiffness = " << Stiffness_ << endl << endl;
	cout << "(Q*R)^T = " << (Q*R).Transpose() << endl << endl;
	cout << "Tangent1 = " << Tangent1_ << endl << endl;
	cout << "Tangent2 = " << Tangent2_ << endl << endl;	
	
}



NewtonPCSolution::NewtonPCSolution(LatticeMode *Mode,char *datafile,const char *prefix,
				   char *startfile,fstream &out,int Echo)
   : Mode_(Mode), CurrentSolution_(0), Echo_(Echo)
{

    // get needed parameters12
	if(!GetParameter(prefix,"PCNumSolutions",datafile,'u',&NumSolutions_)) exit(-1);
	if(!GetParameter(prefix,"PCStepLength",datafile,'l',&CurrentDS_)) exit(-1);
	if(!GetParameter(prefix,"PCContraction",datafile,'l',&cont_rate_nom_)) exit(-1);
	if(!GetParameter(prefix,"PCDistance",datafile,'l',&delta_nom_)) exit(-1);
   	if(!GetParameter(prefix,"PCAngle",datafile,'l',&alpha_nom_)) exit(-1);
	if(!GetParameter(prefix,"PCConvergeCriteria",datafile,'l',&Converge_)) exit(-1);
	if(!GetParameter(prefix,"PCClosedLoopStart",datafile,'u',&ClosedLoopStart_,0))
	{
      // Set default value
      ClosedLoopStart_ = CLOSEDDEFAULT;
	}


   const char *StartType[] = {"Bifurcation","Continuation","ConsistencyCheck"};
   switch (GetStringParameter(prefix,"StartType",startfile,StartType,3))
   {
      case -1:
      {
         cerr << "Unknown StartType!" << endl;
         exit(-1);
         break;
      }
      case 0:
      {
	 // Bifurcation

	 // read in bifurcation point and tangent then proceed as in case 1
         break;
      }
      case 1:
      {
	 // Continuation
	 
         // Get solution1
		Vector one((Mode_->ModeDOF()).Dim());
		Tangent1_.Resize((Mode_->ModeDOF()).Dim());
		Tangent2_.Resize((Mode_->ModeDOF()).Dim());
		
		//int count = one.Dim();
		//int count_plus_one = one.Dim()+1;
		//int count_minus_one = one.Dim()-1;
		
		if(!GetVectorParameter(prefix,"Solution1",startfile,&one)) exit(-1);
		if(!GetVectorParameter(prefix,"Tangent",startfile,&Tangent1_)) exit(-1);
		
		FirstSolution_.Resize(one.Dim());
		FirstSolution_ = one;
		
		Tangent2_ = Tangent1_;
		
		//cout << "size of one = " << one.Dim() << endl;
		
		Mode_->SetModeDOF(one);
		Stiffness_.Resize(Mode_->ModeStiffness().Rows(), Mode_->ModeStiffness().Cols());
		//cout << "size of stiffness rows= " << Stiffness_.Rows() << endl;
		//cout << "size of Mode_->ModeStiffness().Rows_= " << Mode_->ModeStiffness().Rows() << endl << endl;
		//cout << "size of stiffness columns= " << Stiffness_.Cols() << endl; 
		//cout << "size of Mode_->ModeStiffness().Cols_= " << Mode_->ModeStiffness().Cols() << endl;
		
		cout << "FirstSolution = "<< endl<< setw(15) << FirstSolution_ << endl << endl;
		cout << " Tangent1_ = " << endl << setw(15) << Tangent1_ << endl << endl;
		cout << "Tangent2_ = " << endl << setw(15) << Tangent2_ << endl << endl;
		
		cout << " NumSolutions_ = " << NumSolutions_ << endl << endl;
		cout << " CurrentDS_ = " << CurrentDS_ << endl << endl;
		cout << " cont_rate_nom_ = " << cont_rate_nom_ << endl << endl;
		cout << " delta_nom_ = " << delta_nom_ << endl << endl;
		cout << " alpha_nom_ = " << alpha_nom_ << endl << endl;
		cout << " Converge_ = "<< setprecision(20) <<Converge_ << endl << endl;
		cout << " ClosedLoopStart_ = " << setprecision(10) <<ClosedLoopStart_ << endl << endl;
		
		// set tangent. Force1_ Might not be needed.
		//Stiffness_ = Mode_->ModeStiffness();
		//int count = Mode_->ModeStiffness().Cols();
		//int count_minus_one = Mode_->ModeStiffness().Rows();
		
		//QR Decomposition of Stiffness Matrix
		//Matrix Q(count, count);
		//Matrix R(count, count_minus_one);

	
		//QR(Stiffness_, Q, R, 1);	//Performs QR decomposition using A^T = Q*R. Section 4.1 of ISBN 3-540-12760-7
	
		//Assigns Tangent Vector
		//Last column of Q is the tangent vector
		//Tangent1_.Resize(count);
		//Tangent2_.Resize(count);
		/*
		for(int i=0;i< count;i++)
		{
			Tangent1_[i] = Q[i][count_minus_one];
			Tangent2_[i] = Q[i][count_minus_one];
		}
		*/
		/*
		cout << "Stiffness = " << Stiffness_ << endl << endl;
		cout << "(Q*R)^T = " << (Q*R).Transpose() << endl << endl;
		cout << "Tangent1 = " << Tangent1_ << endl << endl;
		cout << "Tangent2 = " << Tangent2_ << endl << endl;
		*/   
			   break;
      }
      case 2:
      {
	 // ConsistencyCheck

	 // do nothing for now
         break;
      }
   }
}



int NewtonPCSolution::AllSolutionsFound()
{
   if (CurrentSolution_ < NumSolutions_)
   { 
      return 0;
   }
   else
   {
      return 1;
   }
   
}


double NewtonPCSolution::FindNextSolution(int &good)
{
   //Finds the next solution
	int omega=1;
	int i, j, Converge_Test ;
	int count = FirstSolution_.Dim();
	//int count_plus_one = count + 1;
	int count_minus_one = count -1;
	double Kappa, Alpha, Delta, Magnitude1, Magnitude2, temp, f, CurrentDS;
	Vector v(count,0);
	Vector u(count,0);
	Vector w(count,0);
	Vector Force1(count_minus_one,0);
	Vector Force2(count_minus_one,0);
	Vector Corrector(count, 0);
	Matrix Q(count, count);
	Matrix R(count, count_minus_one);
	
	
	//size N+1
	u = Mode_->ModeDOF();
	//NEED TO STORE TANGENT AT U
	//Bifurcation Test
	CurrentDS = CurrentDS_;
	cout << "Start: u = " << endl << setw(15) << u << endl << endl;
	cout << "CurrentDS = " << CurrentDS << endl << endl;

	Force1 = Mode_->ModeForce();
	
	cout << "Start: Tangent1 = " <<endl << setw(15) << Tangent1_ << endl << endl;
	cout << "Start: Tangent2 = " <<endl << setw(15) << Tangent2_ << endl << endl;	
	cout << " Tangent1_ * Tangent2_ = " << Tangent1_ * Tangent2_ << endl << endl;
	
	
	//Tangent of size N+1
	if(Tangent1_ * Tangent2_ < 0)
	{
		cout << "Bifurcation Point detected" << endl << endl;
		omega = -omega; 
		Tangent1_ = Tangent2_ * omega;
	}
	else
	{
	Tangent1_ = Tangent2_;
	}
	
	cout << " Tangent1_ = "<< endl << setw(15) << Tangent1_ << endl << endl;
	cout << "START OF FIND NEXT SOLUTION " << endl << endl << endl;
	
	
	//Starts solver
	do
	{
		//Obtains predicted solution
		v = u + CurrentDS * Tangent1_;
		
		cout << "CurrentDS = " << CurrentDS << endl << endl << endl;
		//cout << "Tangent1 = " <<endl <<setw(15) << Tangent2_ << endl << endl;
		
		cout << "v = " <<endl << setw(15) << v << endl << endl;		
		
		//Sets state to predicted point
		Mode_->SetModeDOF(v);
		Force1 = Mode_->ModeForce();
		
		//cout << "Force1 = " <<endl << setw(15) <<  Force1 << endl << endl;
		/*
		cout << "Rows of Modestiffness () = " << Mode_->ModeStiffness().Rows() << endl << endl;
		cout << "Rows of stiffness_ = " << Stiffness_.Rows() << endl << endl;
		cout << "Rows of Modestiffness () = " << Mode_->ModeStiffness().Cols() << endl << endl;
		cout << "Rows of stiffness_ = " << Stiffness_.Cols() << endl << endl;	
		*/
		Stiffness_ = Mode_->ModeStiffness();
		
		//QR Decomposition of Stiffness Matrix
	
		QR(Stiffness_, Q, R, 1);
		
		
		//cout << "Stiffness = " << endl << setw(15) << Stiffness_ << endl << endl;
		//cout << "(Q*R)^T = " << endl << setw(15) <<(Q*R).Transpose() << endl << endl;
		
		
		for(int i=0;i< count;i++)
		{
			//Tangent2_[i] = Q[i][count_minus_one]*omega;
			//Why omega above?
			Tangent2_[i] = Q[i][count_minus_one]*omega;
		}
		
		/*
		cout << "Q = " <<endl<<setw(15) << Q << endl << endl;
		cout << "R = " <<endl<<setw(15) << R << endl << endl;
		
		cout << "Tangent2 = " <<endl <<setw(15) << Tangent2_ << endl << endl;
		cout << "Tangent1 = " <<endl <<setw(15) << Tangent2_ << endl << endl;
		*/

		Corrector = MoorePenrose(Q,R, Force1);
		
		//cout << "Corrector = " << endl << setw(15) <<  Corrector << endl << endl;
		
		Magnitude1 = Corrector.Norm();
		//cout << "Magnitude1 = " << setprecision(20) <<Magnitude1 << endl << endl;
		//cout << setprecision(10);
		
		//CORRECTOR LOOP STARTS HERE
		
		do
		{			
			Converge_Test = 0;
			
			w =  v - Corrector;
			
			cout << "START OF CORRECTOR LOOP "  << endl << endl << endl;
			
			cout << " w = " <<setw(15)<< w << endl << endl;
			
			Mode_->SetModeDOF(w);
	
			Force2 = Mode_->ModeForce();
			
			//cout << "Force2 = " <<endl << setw(15) <<  Force2 << endl << endl;
			//cout << "Force2.Norm() = "<< Force2.Norm() << endl << endl;
		
			Magnitude2 = MoorePenrose(Q,R, Force2).Norm();
			
			//cout << "Magnitude1 = " << setprecision(20) << Magnitude1 << endl << endl;
			//cout << "Magnitude2 = " << Magnitude2 << endl << endl;
	
			//checks parameters for steplength adaptation
			
			Kappa = sqrt((Magnitude2/ Magnitude1)/cont_rate_nom_);
			Alpha = sqrt(acos(Tangent1_ * Tangent2_)/alpha_nom_);
			Delta = sqrt(Magnitude1/delta_nom_);
			
			cout << setprecision(20);
			cout << "kappa = " << Kappa << endl << endl;
			cout << "Alpha = " << Alpha << endl << endl;
			cout << "Delta = " << Delta << endl << endl;
			cout << setprecision(10) ;
			
																			
			temp = max(Kappa, Alpha);
			f = max(temp, Delta);
			
			//cout << "f = " <<f<<  endl << endl;
		
			temp = min(f,2.0);
			f = max(temp, 0.5);
			
			cout << "f = " << f << endl << endl;
			
			CurrentDS = CurrentDS/f;
						
			//cout << "CurrentDS_ /f = " << CurrentDS_ << endl << endl;			
			
			if(f >= 2)
			{
				//cout << "STEPLENGTH TOO LARGE " << endl << endl << endl;
				break;
			}
			else
			{
				v = w;
				
				cout << "STEPLENGTH OKAY " << endl << endl << endl;
				/*
				cout << " STEPLENGTH OKAY: f<2 , v = " << endl<< setw(15) <<  v << endl << endl;
				
				cout << "steplength okay: Force2.Norm() = " <<setprecision(20) <<  Force2.Norm() <<endl << endl;
				cout << "Steplength okay: magnitude1 = " << Magnitude1 << endl << endl;
				cout << "Convergence criteria = " << Converge_ << endl << endl;
				cout << setprecision(10) << endl;
				*/
				
				//Inefficient? When already converged, state is at desired mode so no need to calculate Tangent2_, etc.
				if (Force2.Norm() <= Converge_ && Corrector.Norm() <= Converge_)
				{
					Converge_Test = 1;
				}
				
				if (Converge_Test == 0)
				{
					Mode_->SetModeDOF(v);
	
					Force1 = Mode_->ModeForce();
					cout << "HAS NOT CONVERGED " << endl << endl << endl;
					
					//cout << "has not converged, Force1 = " << endl << setw(15) << Force1 << endl << endl;
					
					Stiffness_ = Mode_->ModeStiffness();
					
					//QR Decomposition of Stiffness Matrix
	
					QR(Stiffness_, Q, R, 1);
		
					//cout << "has not converged, Stiffness = " << endl << setw(15) <<Stiffness_ << endl << endl;
					//cout << "(Q*R)^T = " << endl << setw(15) <<(Q*R).Transpose() << endl << endl;
					
					/*
					for(int i=0;i< count;i++)
					{
						Tangent2_[i] = Q[i][count_minus_one];
					}
					*/
					//cout << "Q = " <<endl<<setw(15) << Q << endl << endl;
					///cout << "R = " <<endl<<setw(15) << R << endl << endl;
					
					//cout << "has not converged, Tangent2 = " << endl << setw(15) << Tangent2_ << endl << endl;
								
					
					Corrector = MoorePenrose(Q,R, Force1);
					
					//cout << " has not converged "  << endl << endl;
					//Magnitude1 = Corrector.Norm();
		
					//cout << "has not converged, Corrector = " << endl << setw(15) << Corrector << endl << endl;
					//cout << " has not converged, Magnitude1 = " << Magnitude1 << endl << endl;

				}
				
			}
			
		}
		while (Converge_Test != 1);
	}
	while (f >= 2);
	
	Stiffness_ = Mode_->ModeStiffness();

	//QR Decomposition of Stiffness Matrix
	
	QR(Stiffness_, Q, R, 1);
	
	//cout << "HAS CONVERGED, Stiffness = " << Stiffness_ << endl << endl;
	//cout << "(Q*R)^T = " << (Q*R).Transpose() << endl << endl;
					
	for(int i=0;i< count;i++)
	{
		Tangent2_[i] = Q[i][count_minus_one];
	}
	
	//cout << "HAS CONVERGED, Tangent2 = " << endl << setw(15) << Tangent2_ << endl << endl<< endl << endl<<endl << endl<< endl << endl<< endl << endl<<endl << endl;
	
	cout << "HAS CONVERGED " << endl << endl << endl;	
	
	if ((ClosedLoopStart_ >= 0) && (CurrentSolution_ > ClosedLoopStart_) &&
       ((Mode_->ModeDOF() - FirstSolution_).Norm() < CurrentDS_))
   {
      // We are done -- set currentsolution to numsolutions
      cerr << "Closed Loop detected at Solution # " << CurrentSolution_
	   << " --- Terminating!" << endl;

      CurrentSolution_ = NumSolutions_;
   }
   else
   {
      CurrentSolution_++;
   }

}


int NewtonPCSolution::BisectAlert(Lattice *Lat,char *datafile,const char *prefix,int Width,
				  fstream &out)
{
   // for the moment, do nothing
   return 1;
}

//Returns the vector H'+ * H 
Vector NewtonPCSolution::MoorePenrose(const Matrix& Q, const Matrix& R, const Vector& Force)
{
	double sum;
	int i,j;
	int k=0;
	int Size = Force.Dim()+1;
	int count_minus_one = Size -1;
	Vector Corrector(Size, 0);
	Vector y(Size, 0);
	
	for ( i = 0; i < count_minus_one; i++)
			{
				sum = 0;
				for ( j = 0; j < k;j++)
				{
					sum +=R[j][i] * y[j];
				}
		
				y[i] = (Force[i] - sum)/R[i][i];
	
				k++;
			}
	
	Corrector = Q * y;
	
	return Corrector;
	
}

