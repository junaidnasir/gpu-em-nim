
#define PRECISION float

/*
#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
*/

__kernel void hy_kernel(__global PRECISION *hy, __global PRECISION *ez, __global PRECISION *mu, const PRECISION delt, const PRECISION delx, const int SIZE) 
{
	unsigned int i = get_global_id(0);
  if(i < SIZE)
    hy[i] = hy[i] + (ez[i+1] - ez[i]) * (delt/(delx*mu[i]));
}
__kernel void ez_kernel(__global PRECISION *hy, __global PRECISION *ez, __global PRECISION *epsilon, const PRECISION delt, const PRECISION delx, const int SIZE,__global PRECISION *Etemp,__global PRECISION *Exz1,__global PRECISION *Exz2,const int qTime,const PRECISION mur,const PRECISION epsilonr,const PRECISION Sc,PRECISION ez1q,PRECISION ez2q,PRECISION ezmq,PRECISION ezm1q,const PRECISION pi,const PRECISION f,const int SourceSelect,const PRECISION PulseWidth,const int Z1,const int Z2)
{
	unsigned int i = get_global_id(0);
  if(i > 0 && i < SIZE)
    ez[i] = ez[i] + (hy[i] - hy[i-1]) * (delt/(delx*epsilon[i]));

  if (SourceSelect==0)	
  {														//Source node
	ez[1] = ez[1] + (sin(2*pi*(qTime)*f*delt)*Sc);
	}
  else
  {
	ez[1] = ez[1] + exp((-(qTime+1 - 30) * (qTime - 30)) / (PulseWidth/4));
	}
  Etemp[qTime]= ez[SIZE-(SIZE/2)+2]; 													//Save ez after boundary
  // -------- Absorbing Boundary Conditions -------- 
  ez[0] = ez2q+(ez[1]-ez1q)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );
  ez[SIZE-1] = ezm1q+(ez[SIZE-1-1]-ezmq)*( ((Sc/pow(mur*epsilonr,0.5))-1 ) / ((Sc/pow(mur*epsilonr,0.5))+1) );	  
  
  barrier(CLK_GLOBAL_MEM_FENCE|CLK_LOCAL_MEM_FENCE);    
  // -------- Saving pervious step time values -------- 
  ez2q = ez[1]; 
  ez1q = ez[0]; 
  ezmq = ez[SIZE-1];
  ezm1q= ez[SIZE-2];
  Exz1[qTime] = ez[Z1-1];
  Exz2[qTime] = ez[Z2-1];
}
