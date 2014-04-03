#include <iostream>
#include <fstream>
#include <sstream>
#include <windows.h>
#include "Timer.h"

using namespace std;
const double pi = 3.14;
__int64 tStart=0;
__int64 tEnd=0;
__int64 tDelta=0;
double time;
bool tPaused;

void StartTimer();
void StopTimer();
void ResetTimer();
float GetElapsedTime();
void main()
{
	
	for (int a=22;a<25;a++)
	{
	tStart = GetTimeus64();
	// -------- Save to file Variables -------- 
	fstream snapshot;
	std::string filename ;
	std::stringstream stream;
	stream<<"results";
	CreateDirectory(stream.str().c_str(), NULL) ;		//create directory of results
	double temp=pow(2.,a);
	int SIZE=temp;
	cout<<SIZE<<endl;
    int maxTime = 1024;
	int SourceSelect = 1; 			// 0=Sinosoidal, 1=Gauassian
	/*cout<<"----Select Source----"<<endl;
	cout<<"0)Sinosoidal 1) Gauassian"<<endl;
	do{
		cin>>SourceSelect;
	}while(SourceSelect != 1 && SourceSelect != 0);
	if (SourceSelect == 0)
		maxTime = 5001;
*/
	// -------- Constants -------- 
	double c = 3e8;

 	/* Courant Number (Accuracy) Sc
	% ideal Condition --> Sc= c*delt/delx = 1
	% f=3Ghz, lambda=c/f=0.1m, for 4 wavelengths, dx=0.4/(maxtime=1000)*/
	double PulseWidth = 800;
	double f = 3e9;
	double w = 2 * pi * f;    		// omega
	double k0 = w/c     ; 			// free space wave number constant
	double lambda = c / f;
	double delx = (4*lambda) / SIZE;
	double delt = delx / c;
	double Sc = c * delt / delx;
	int epsilonr = 1;
	int mur = 1;

    // -------- Incident and Refelected Waves Variables -------- 
	double *Eincident = new double[maxTime] ; 
	for(int j=0; j<maxTime; j++)
		Eincident[j] = 0;
	
	double *Etransmitted = new double [maxTime] ;
	for(int j=0; j<maxTime; j++)
		Etransmitted[j] = 0;

	double *Etemp = new double [maxTime];
	for(int j=0; j<maxTime; j++)
		Etemp[j] = 0;

	// -------- refractive index variables -------- 
	int Z1 = 750;
	double z1 = Z1*delx;
	int Z2 = 760;
	double z2 = Z2*delx;

	double *Exz1 =new double [maxTime];
	for (int i=0; i<maxTime ; i++)
		Exz1[i] = 0;

	double *Exz2 =new double [maxTime];
	for (int i=0; i<maxTime ; i++)
		Exz2[i] = 0;
	
	double *mu = new double [SIZE];
	for(int j=0; j<SIZE; j++)
			mu[j] = 1.2566e-006;

	double *epsilon = new double [SIZE];

	double *ez = new double [SIZE];
	for(int j=0; j<SIZE; j++)
		ez[j] = 0;

	double *hy = new double [SIZE-1];
	for(int j=0; j<SIZE-1; j++)
		hy[j] = 0;
	
	for (int medium=1; medium<=2; medium++)
	{
    	// Temp Variable
		int mm = 0;
		double ez1q = 0;
		double ez2q = 0;
		double ezmq = 0;
		double ezm1q = 0;
		
		// -------- Medium Specifications -------- 
		if (medium==1)
		{
//			cout<<"Calculating Wave propagation in Free Space"<<endl;
			for(int j=0; j<SIZE; j++)
				epsilon[j] = 8.8542e-012;   						//epsilon=8.8542e-012*ones(1,SIZE); // free space
		}
		else
		{
//			cout<<"Calculating Wave propagation in denser Medium"<<endl;
			int j;
			for(j=0; j<SIZE-500; j++)								//epsilon=[8.8542e-012*ones(1,SIZE-500) 1.7708e-011*ones(1,500)]; // half medium
				epsilon[j] = 8.8542e-012;
			for(int k=j ;k<SIZE ;k++)
				epsilon[j] = 1.7708e-011;
		}

		// -------- Main Loop -------- 
		for (int qTime=0; qTime<(maxTime-1); qTime++)
		{
			for(int nn=0; nn<(SIZE-1); nn++)
			{
				hy[nn] = hy[nn] + (ez[nn+1] - ez[nn]) * (delt/(delx*mu[nn]));				//Update Magnetic field
			}

			for(int nn=1; nn<(SIZE-1); nn++)
			{
				ez[nn] = ez[nn] + (hy[nn] - hy[nn-1]) * (delt/(delx*epsilon[nn]));			//Update Electrical field
			}

			if (SourceSelect==0)															//Source node
			    ez[1] = ez[1] + (sin(2*pi*(qTime)*f*delt)*Sc);
			else
			    ez[1] = ez[1] + exp((-(qTime+1 - 30) * (qTime - 30)) / (PulseWidth/4));

	        Etemp[qTime]= ez[SIZE-498]; 													//Save ez after boundary
		    
		    // -------- Absorbing Boundary Conditions -------- 
	        ez[0] = ez2q+(ez[1]-ez1q)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );	
	        ez[SIZE-1] = ezm1q+(ez[SIZE-1-1]-ezmq)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );
		    
		    // -------- Saving pervious step time values -------- 
			ez2q = ez[1]; 
			ez1q = ez[0]; 
			ezmq = ez[SIZE-1];
			ezm1q= ez[SIZE-2];
			Exz1[qTime] = ez[Z1-1];
	        Exz2[qTime] = ez[Z2-1];

	        // -------- Saving to file -------- 
			stream.str(std::string());   						// clear stringstream
			stream<<"./results/"<<"Efield"<<medium<<"_"<<qTime<<".jd";   		// concatenate
			filename = stream.str();		 					// copy string
			snapshot.open(filename, ios::out|ios::binary);
			for (mm = 0; mm < SIZE; mm++)
				snapshot.write((char *)&ez[mm],sizeof(double));
			snapshot.close();
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
	
//	cout<<"Writing Values to files"<<endl;
	stream.str(std::string());
	stream<<"./results/"<<"Eincident"<<".jd";
	filename = stream.str();
	snapshot.open(filename, ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Eincident[mm],(sizeof(double)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"Etransmitted"<<".jd";
	filename = stream.str();
	snapshot.open(filename, ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Etransmitted[mm],(sizeof(double)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"Exz1"<<".jd";
	filename = stream.str();
	snapshot.open(filename, ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Exz1[mm],(sizeof(double)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"Exz2"<<".jd";
	filename = stream.str();
	snapshot.open(filename, ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Exz2[mm],(sizeof(double)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"maxTime"<<".jd";
	filename = stream.str();
	snapshot.open(filename, ios::out|ios::binary);
	snapshot.write((char *)&maxTime,(sizeof(int)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"data"<<".jd";
	filename = stream.str();
	snapshot.open(filename, ios::out|ios::binary);
	snapshot.write((char *)&delt,(sizeof(double)));
	snapshot.write((char *)&k0,(sizeof(double)));
	snapshot.write((char *)&z1,(sizeof(double)));
	snapshot.write((char *)&z2,(sizeof(double)));
	snapshot.close();

	cout<<endl;
	//cout<<"Enter anything to exit"<<endl;
	//cin>>SourceSelect;

	tEnd = GetTimeus64(); // Get ending time.
	time=((double)(tEnd-tStart))/(1000000.) ;
	cout << "Time taken = " << time << " seconds." << endl;
	ResetTimer();
	//cin>>SourceSelect;
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
	tStart=0;
	tEnd=0;
	}
	cin>>time;
}

void StartTimer()
{
	if (tPaused == true)
	{
		tStart = GetTimeus64();
		tPaused = false;
	}
}
void StopTimer()
{
	if (tPaused == false)
	{
		tEnd = GetTimeus64();
		tDelta += tEnd - tStart;
		tStart = tEnd;
		tPaused = true;
	}
}
void ResetTimer()
{
	if (tPaused == true)
		tStart = tEnd;
	else
		tStart = GetTimeus64();

	tDelta = 0UL;
}
float GetElapsedTime()
{
	if (tPaused == false)
		tEnd = GetTimeus64();

	return ((double)(tEnd-tStart+tDelta))/(1000000.);
}