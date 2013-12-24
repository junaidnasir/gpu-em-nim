# include<iostream>
# include<cmath>
using namespace std;
const float pi = 3.14;
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
    // Incident and Refelected Waves Variables
	int *Eincident = new int[maxTime] ; //1 x maxtime
	for(int j=0; j<maxTime; j++)
		Eincident[j] = 0;
	
	int *Etransmitted = new int [maxTime] ; //1 x maxtime
	for(int j=0; j<maxTime; j++)
		Etransmitted[j] = 0;

	int *Etemp = new int [maxTime];
	for(int j=0; j<maxTime; j++)
		Etemp[j] = 0;

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
	
	float *mu = new float [SIZE];
	for(int j=0; j<SIZE; j++)
			mu[j] = 1.2566e-006;

	float *epsilon = new float [SIZE];

	float *ez = new float [SIZE];
	for(int j=0; j<SIZE; j++)
		ez[j] = 0;

	//hy=zeros(1,SIZE-1);
	float *hy = new float [SIZE-1];
	for(int j=0; j<SIZE-1; j++)
		hy[j] = 0;
	
	for (int medium=1; medium<=2; medium++)
	{
    // Temp Variable
    //ez=zeros(1,SIZE);
		int mm = 0;
		int ez1q = 0;
		int ez2q = 0;
		int ezmq = 0;
		int ezm1q = 0;
		// Medium Specifications
		//mu=1.2566e-006*ones(1,SIZE);   %permeability of free sapce
		
		if (medium==1)
		{
			//epsilon=8.8542e-012*ones(1,SIZE); // free space
			for(int j=0; j<SIZE; j++)
				epsilon[j] = 8.8542e-012;
		}
		else
		{
			//epsilon=[8.8542e-012*ones(1,SIZE-500) 1.7708e-011*ones(1,500)]; // half medium
			int j;
			for(j=0; j<SIZE-500; j++)
				epsilon[j] = 8.8542e-012;
			for(int k=j ;k<SIZE ;k++)
				epsilon[j] = 1.7708e-011;

		}
		for (int qTime=0; qTime<(maxTime-1); qTime++)
		{
			//        Update Magnetic field
	 
			//hy(mm) = hy(mm) + (ez(mm + 1) - ez(mm)) * (delt/(delx*mu(mm)));
			for(int nn=0; nn<(SIZE-1); nn++)
			{
				hy[nn] = hy[nn] + (ez[nn+1] - ez[nn]) * (delt/(delt*mu[nn]));
			}
			//         Update Electrical filed
	        
			//ez(mm) = ez(mm) + (hy(mm) - hy(mm - 1)) * (delt/(delx*epsilon(mm))) ;
			for(int nn=1; nn<(SIZE-1); nn++)
			{
				ez[nn] = ez[nn] + (hy[nn] - hy[nn-1]) * (delt/(delt*epsilon[nn]));
			}
			
	        Etemp[qTime]= ez[SIZE-498]; //just after boundary of medium
	        if (SourceSelect==0)
			//        Source node (hard coded)
			    ez[1] = ez[1] + (sin(2*pi*(qTime+1)*f*delt)*Sc);
			else
			    ez[1] = ez[1] + exp(-(qTime+1 - 30) * (qTime+1 - 30) / (PulseWidth/4));
			
		    //         Absorbing Boundary Conditions
			 //ez(1)=ez2q+(ez(2)-ez1q)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
	        ez[0] = ez2q+(ez[1]-ez1q)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );
			//ez(SIZE)=ezm1q+(ez(SIZE-1)-ezmq)*(((Sc/(mur*(epsilonr))^0.5)-1)/((Sc/(mur*(epsilonr))^0.5)+1));
	        ez[SIZE-1] = ezm1q+(ez[SIZE-1-1]-ezmq)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );
		    //         Saving q-1 (pervious step time values) for boundary Conditions
			ez2q = ez[1]; //ez2q=ez(2);
			ez1q = ez[0]; //ez1q=ez(1);
			ezmq = ez[SIZE-1]; //ezmq=ez(SIZE);
			ezm1q= ez[SIZE-2]; //ezm1q=ez(SIZE-1);
			Exz1[qTime] = ez[Z1-1]; //   Exz1(qTime)=ez(Z1);
	        Exz2[qTime] = ez[Z2-1];
	        //snapshot code here
		} //end qTime

		if (medium==1)
		{
			for(int i=0;i<maxTime;i++)
			Eincident[i] = Etemp[i];
		}
		else
		{
			for(int i=0;i<maxTime;i++)
			Etransmitted[i] = Etemp[i];
		}
		
	} //end medium loop

	// -------- Memory Deallocaion -----------
	
	delete[] Eincident;

	delete[] Etransmitted;

	delete[] Etemp;
	
	delete[] mu;

	delete[] epsilon;

	delete[] ez;

	delete[] hy;
	
	delete[] Exz1;

	delete[] Exz2;
}