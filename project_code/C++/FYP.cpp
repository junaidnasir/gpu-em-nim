# include<iostream>
# include<cmath>
using namespace std;
const int pi = 3.14;
void main()
{
	int SIZE=1001;
    int maxTime = 1001;
	int SourceSelect = 1; // 0=Sinosoidal, 1=Gauassian
	if (SourceSelect == 0)
		maxTime = 5001;
	// Constants
	int c = 3e8;
 /* Courant Number (Accuracy) Sc
	% ideal Condition --> Sc= c*delt/delx = 1
	% f=3Ghz, lambda=c/f=0.1m, for 4 wavelengths, dx=0.4/(maxtime=1000)*/
	int PulseWidth = 800;
	int f = 3e9;
	float w = 2 * pi * f;    // omega
	float k0 = w/c     ; // free space wave number constant
	float lambda = c / f;
	float delx = (4*lambda) / SIZE;
	// so dt=dx/c=1.333e-12
	float delt = delx / c;
	float Sc = c * delt / delx;
	int epsilonr = 1;
	int mur = 1;
    // Incident and Reflected Waves Variables
	int **Eincident;
	Eincident = new int* [1] ; //1 x maxtime
	Eincident[0] = new int [maxTime];
	for(int j=0; j<maxTime; j++)
		Eincident[0][j] = 0;
	
	int **Etransmitted;
	Etransmitted = new int* [1] ; //1 x maxtime
	Etransmitted[0] = new int [maxTime];
	for(int j=0; j<maxTime; j++)
		Etransmitted[0][j] = 0;

	int **Etemp;
	Etemp = new int* [1] ; //1 x maxtime
		Etemp[0] = new int [maxTime];
		for(int j=0; j<maxTime; j++)
			Etemp[0][j] = 0;

	// refractive index variables
	int Z1 = 750;
	float z1 = Z1*delx;
	int Z2 = 760;
	float z2 = Z2*delx;
	int fspan = 100; // Points to plot in frequency domain

	// record Electric field at 750 (maxTime x 1)
	int *Exz1 =new int [maxTime];
	for (int i=0; i<maxTime ; i++)
		Exz1[i] = 0;

	// Exz2 = zeros(maxTime,1); // record electric field at 760
	int *Exz2 =new int [maxTime];
	for (int i=0; i<maxTime ; i++)
		Exz2[i] = 0;
	
	for (int medium=1; medium<=2; medium++)
	{
    // Temp Variable
    //ez=zeros(1,SIZE);
		float **ez;
		ez = new float* [1] ; //1 x maxtime
		for (int i=0; i<1 ; i++)
			ez[i] = new float [SIZE];
		for (int i=0; i<1 ; i++)
		{
			for(int j=0; j<SIZE; j++)
				ez[i][j] = 0;
		}
		//hy=zeros(1,SIZE-1);
		float **hy;
		hy = new float* [1] ; //1 x maxtime
		for (int i=0; i<1 ; i++)
			hy[i] = new float [SIZE-1];
		for (int i=0; i<1 ; i++)
		{
			for(int j=0; j<SIZE-1; j++)
				hy[i][j] = 0;
		}
		int mm = 0;
		int ez1q = 0;
		int ez2q = 0;
		int ezmq = 0;
		int ezm1q = 0;
		// Medium Specifications
		//mu=1.2566e-006*ones(1,SIZE);   %permeability of free sapce
		float **mu;
		mu = new float* [1] ; //1 x maxtime
		for (int i=0; i<1 ; i++)
			mu[i] = new float [SIZE];
		for(int j=0; j<SIZE; j++)
				mu[0][j] = 1.2566e-006;

		float **epsilon;
		epsilon = new float* [1] ; //1 x maxtime
		epsilon[0] = new float [SIZE];

		if (medium==1)
		{
			//epsilon=8.8542e-012*ones(1,SIZE); // free space
			for(int j=0; j<SIZE; j++)
				epsilon[0][j] = 8.8542e-012;
		}
		else
		{
			//epsilon=[8.8542e-012*ones(1,SIZE-500) 1.7708e-011*ones(1,500)]; // half medium
			int j;
			for(j=0; j<SIZE-500; j++)
				epsilon[0][j] = 8.8542e-012;
			for(int k=j ;k<SIZE ;k++)
				epsilon[0][j] = 1.7708e-011;

		}
		for (int qTime=0; qTime<(maxTime-1); qTime++)
		{
			//        Update Magnetic field
	 
			//hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) * (delt/(delx*mu(mm)));
			for(int nn=0; nn<(SIZE-1); nn++)
			{
				hy[0][nn] = hy[0][nn] + (ez[0][nn+1] - ez[0][nn]) * (delt/(delt*mu[0][nn]));
			}
			//         Update Electrical filed
	        
			//ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * (delt/(delx*epsilon(mm))) ;
			for(int nn=1; nn<(SIZE-1); nn++)
			{
				ez[0][nn] = ez[0][nn] + (hy[0][nn] - hy[0][nn-1]) * (delt/(delt*epsilon[0][nn]));
			}
			
	        Etemp[0][qTime]= ez[0][SIZE-498]; //just after boundary of medium
	        if (SourceSelect==0)
			//        Source node (hard coded)
			    ez[0][1] = ez[0][1] + (sin(2*pi*(qTime+1)*f*delt)*Sc);
			else
			    ez[0][1] = ez[0][1] + exp(-(qTime+1 - 30) * (qTime+1 - 30) / (PulseWidth/4));
			
		    //         Absorbing Boundary Conditions
			 //ez(1)=ez2q+(ez(2)-ez1q)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
	        ez[0][0] = ez2q+(ez[0][1]-ez1q)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );
			//ez(SIZE)=ezm1q+(ez(SIZE-1)-ezmq)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
	        ez[0][SIZE-1] = ezm1q+(ez[0][SIZE-1-1]-ezmq)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );
		    //         Saving q-1 (pervious step time values) for boundary Conditions
			ez2q = ez[0][1]; //ez2q=ez(2);
			ez1q = ez[0][0]; //ez1q=ez(1);
			ezmq = ez[0][SIZE-1]; //ezmq=ez(SIZE);
			ezm1q= ez[0][SIZE-2]; //ezm1q=ez(SIZE-1);
			Exz1[qTime] = ez[0][Z1-1]; //   Exz1(qTime)=ez(Z1);
	        Exz2[qTime] = ez[0][Z2-1];
	        //snapshot code here
		} //end qTime

		if (medium==1)
			Eincident = Etemp;
		else
			Etransmitted = Etemp;
	}
	
	
}