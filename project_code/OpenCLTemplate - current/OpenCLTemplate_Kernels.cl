
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
__kernel void ez_kernel(__global PRECISION *hy, __global PRECISION *ez, __global PRECISION *epsilon, const PRECISION delt, const PRECISION delx, const int SIZE)
{
	unsigned int i = get_global_id(0);
  if(i > 0 && i < SIZE)
    ez[i] = ez[i] + (hy[i] - hy[i-1]) * (delt/(delx*epsilon[i]));
}
