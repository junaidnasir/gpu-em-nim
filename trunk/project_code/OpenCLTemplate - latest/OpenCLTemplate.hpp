#ifndef OpenCLTemplate_H_
#define OpenCLTemplate_H_

#define PRECISION float
//#define PRECISION double

#include <Timer.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CL/cl.h>

using namespace std;

class COpenCLTemplate
{
private:
	int maxTime;
	int SIZE;
	
    fstream snapshot;
	std::string filename ;
	std::stringstream stream;

	PRECISION pi;
    PRECISION c;
	PRECISION PulseWidth;
	PRECISION f;
	PRECISION w;
	PRECISION k0;
	PRECISION lambda;
	PRECISION delx;
	PRECISION delt;
	PRECISION Sc;
	PRECISION epsilonr;
	PRECISION mur;

	int SourceSelect;
	PRECISION ez1q; 
    PRECISION ez2q;
	PRECISION ezmq;
	PRECISION ezm1q;

	int qTime;
	int Z1;
	int Z2;
	PRECISION* Eincident;
	PRECISION* Etransmitted;
	PRECISION* Etemp;
	PRECISION* Exz1;
	PRECISION* Exz2;

	PRECISION* mu;
	PRECISION* epsilon;
	PRECISION* ez;
	PRECISION* hy;
	// Scalar.
	

	// Device data arrays.
    cl_mem hy_gpu;
    cl_mem ez_gpu;
    cl_mem mu_gpu;
    cl_mem epsilon_gpu;
	cl_mem Etemp_gpu;
	cl_mem Exz1_gpu;
	cl_mem Exz2_gpu;
	// Scalar.
	// ============ OPENCL related parameters ===========
	// OPENCL context/device/program
	cl_context context;
	cl_device_id *devices;
	cl_command_queue commandQueue;
	cl_program program;
	cl_kernel hykernel;
	cl_kernel ezkernel;
	// ==================================================

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;
	__int64 tDelta;
	bool tPaused;

public:
	COpenCLTemplate(int=1000, int=1024);	// Default width of problem is 256 and multiplier is 2.

	// Memory allocation and initialisation.
	int AllocateMemoryCPU();
	int InitialiseCPU();
	int InitialiseCL();			// Search and allocate a device.
	int AllocateMemoryGPU();
	int InitialiseCLKernelsGPU(); // Build/attach kernels to respective kernel functions and set arguments.
	int RunCLKernels();

	// Complete run encapsulating all the sub-functions.
	int CompleteRun();

	// Timing.
	void StartTimer();
	void StopTimer();
	void ResetTimer();
	PRECISION GetElapsedTime();

	std::string convertToString(const char * filename);

	int SafeCall(cl_int, const char []=NULL);

	int CleanupCPU();
	int CleanupCL();
	int CleanupGPU();
	~COpenCLTemplate ();
};
#endif // #ifndef OpenCLTemplate_H_