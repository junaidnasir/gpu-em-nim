#ifndef OpenCLTemplate_H_
#define OpenCLTemplate_H_

#define PRECISION float

#include <Timer.h>
#include <string>
#include <CL/cl.h>

class COpenCLTemplate
{
private:
	// Size of problem domain
//	const unsigned int Width;
	const unsigned int maxTime;
	// Host data arrays.
	PRECISION* input;
	PRECISION* output;
	PRECISION* arr;

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
	const PRECISION Multiplier; //variables here----------
	

	// Device data arrays.
	cl_mem d_input;
	cl_mem d_output;
	cl_mem d_arr; //buffer in gpu

	cl_mem d_Eincident;
	cl_mem d_Etransmitted;
	cl_mem d_Etemp;
	cl_mem d_Exz1;
	cl_mem d_Exz2;
	cl_mem d_mu;
	cl_mem d_epsilon;
	cl_mem d_ez;
	cl_mem d_hy;
	// Scalar.
	// ============ OPENCL related parameters ===========
	// OPENCL context/device/program
	cl_context context;
	cl_device_id *devices;
	cl_command_queue commandQueue;
	cl_program program;
	cl_kernel kernel;
	cl_kernel kernel2;
	// ==================================================

	// Timer variables.
	__int64 tStart;
	__int64 tEnd;
	__int64 tDelta;
	bool tPaused;

public:
	COpenCLTemplate(unsigned int=256U, const PRECISION=2.0);	// Default width of problem is 256 and multiplier is 2.

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
