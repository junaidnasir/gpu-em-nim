#include <OpenCLTemplate.hpp>
#include <Timer.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>
//#include <windows.h>
using namespace std;

COpenCLTemplate::COpenCLTemplate(int maxTime, int SIZE): maxTime(maxTime), SIZE(SIZE), tStart(0LL), tEnd(0LL), tDelta(0LL), tPaused(true)
{
	// Printing simulation size.
	cout << "Startin Up " << endl;
}
// Allocate memory for data arrays.
int COpenCLTemplate::AllocateMemoryCPU()
{
	Eincident =     new PRECISION[maxTime];
	Etransmitted =  new PRECISION[maxTime];
	Etemp =         new PRECISION[maxTime]; 
    Exz1 =          new PRECISION[maxTime];
	Exz2 =          new PRECISION[maxTime];

	mu =            new PRECISION[SIZE];
	epsilon =       new PRECISION[SIZE];
	ez =            new PRECISION[SIZE];
	hy =            new PRECISION[SIZE];

	return 0;
}
// Initialise CPU data.
int COpenCLTemplate::InitialiseCPU()
{

	stream<<"results";
	//CreateDirectory(stream.str().c_str(), NULL) ;		//create directory of results
    SourceSelect=1;
	pi = 3.14;
    c = 3e8;
	PulseWidth = 800;
	f = 3e9;
	w = 2 * pi * f;    		// omega
	k0 = w/c     ; 			// free space wave number constant
	lambda = c / f;
	delx = (4*lambda) / SIZE;
	delt = delx / c;
	Sc = c * delt / delx;
	epsilonr = 1;
	mur = 1;
	ez1q = 0; 
    ez2q = 0;
	ezmq = 0;
	ezm1q = 0;
	Z1=750;
	Z2=760;
	qTime=0;
	for (unsigned int i=0; i<maxTime; i++) {
		Eincident[i] = 0.; 
		Etransmitted[i] = 0.; 
		Etemp[i] = 0.;
		Exz1[i] = 0.;
		Exz2[i] = 0.;
  }

	for (unsigned int i=0; i<SIZE; i++) {
		ez[i] = 0.;
		hy[i] = 0.;
		mu[i] = 1.2566e-006;
		epsilon[i] = 8.8542e-012;
  }

	return 0;
}
int COpenCLTemplate::InitialiseCL()
{
	cl_int status = 0;
	size_t deviceListSize;

	/*
	* Have a look at the available platforms and pick either
	* the AMD one if available or a reasonable default.
	*/

	cl_uint numPlatforms;
	cl_platform_id platform = NULL;
	SafeCall(clGetPlatformIDs(0, NULL, &numPlatforms), "Error: Getting Platforms. (clGetPlatformsIDs)");

	char AMDPlatform[] = "Advanced Micro Devices, Inc.";
	char nVidiaPlatform[] = "NVIDIA Corporation";
	char *SelectedPlatform = NULL;

	char choice = '1';
/*	cout << "Choose a platform: " << endl;
	cout << "[1] Advanced Micro Devices, Inc. (default)" << endl;
	cout << "[2] NVIDIA Corporation" << endl;
	cout << ">>";
	StopTimer();
	cin >> choice;
	StartTimer();
*/
	if (choice == '1')
		SelectedPlatform = AMDPlatform;
	else if (choice == '2')
		SelectedPlatform = nVidiaPlatform;
	else
	{
		cout << "Reverting to default platform..." << endl;
		SelectedPlatform = AMDPlatform;
	}

	cout << "Detecting platforms..." << endl;
	cout << "Available platforms are: " << endl;
	if(numPlatforms > 0)
	{
		cl_platform_id* platforms = new cl_platform_id[numPlatforms];
		SafeCall(clGetPlatformIDs(numPlatforms, platforms, NULL), "Error: Getting Platform Ids. (clGetPlatformsIDs)");

		for(unsigned int i=0; i < numPlatforms; ++i)
		{
			char pbuff[100];
			SafeCall(clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(pbuff), pbuff, NULL), "Error: Getting Platform Info.(clGetPlatformInfo)");

			cout << "Platform " << i << " : " << pbuff << endl;
			if(!strcmp(pbuff, SelectedPlatform))
				platform = platforms[i];
		}
		delete platforms;
	}

	if(NULL == platform)
	{
		cout << "Selected platform not found so Exiting Application." << endl;
		return 1;
	}

	/*
	* If we could find our platform, use it. Otherwise use just available platform.
	*/
	cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };

	/////////////////////////////////////////////////////////////////
	// Create an OpenCL context
	/////////////////////////////////////////////////////////////////
	cl_device_type type;

/*	cout << "Emulate GPU run on CPU?" << endl;
	cout << "[1] Yes" << endl;
	cout << "[2] No (default)" << endl;
	cout << ">>";
	StopTimer();
	cin >> choice;
	StartTimer();
*/
	choice='2';
	if (choice == '1')
	{
		if(!strcmp(AMDPlatform, SelectedPlatform))
			cout << "Running on CPU with GPU emulation..." << endl;
		else
			cout << "Warning: Selected platform does not support GPU emulation on CPU." << endl;

		type = CL_DEVICE_TYPE_CPU;
	}
	else
	{
		cout << "Running on GPU..." << endl;
		type = CL_DEVICE_TYPE_GPU;
	}

	context = clCreateContextFromType(cps, type, NULL, NULL, &status);
	SafeCall(status, "Error: Creating Context. (clCreateContextFromType)");

	/* First, get the size of device list data */
	SafeCall(clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize), "Error: Getting Context Info (device list size, clGetContextInfo)");

	/////////////////////////////////////////////////////////////////
	// Detect OpenCL devices
	/////////////////////////////////////////////////////////////////
	devices = new cl_device_id[deviceListSize/sizeof(cl_device_id)];
	SafeCall(!devices, "Error: No devices found.");

	/* Now, get the device list data */
	SafeCall(clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceListSize, devices, NULL), "Error: Getting Context Info (device list, clGetContextInfo)");

	char platformVendor[1024];
	SafeCall(clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, sizeof(platformVendor), platformVendor, NULL), "clGetPlatformInfo failed");
	cout << "Selected Platform Vendor : " << platformVendor << endl;

	// Get number of devices available 
	cl_uint deviceCount = 0;
	SafeCall(clGetDeviceIDs(platform, type, 0, NULL, &deviceCount), "clGetDeviceIDs failed");

	cl_device_id* deviceIds = (cl_device_id*)malloc(sizeof(cl_device_id) * deviceCount);
	SafeCall(!deviceIds, "Failed to allocate memory(deviceIds)");

	// Get device ids
	SafeCall(clGetDeviceIDs(platform, type, deviceCount, deviceIds, NULL), "clGetDeviceIDs failed");

	cout << "Available devices are: " << endl;
	// Print device index and device names
	for(cl_uint i = 0; i < deviceCount; ++i)
	{
		char deviceName[1024];
		SafeCall(clGetDeviceInfo(deviceIds[i], CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL), "clGetDeviceInfo failed");
		cout << "Device " << i << " : " << deviceName <<" Device ID is "<<deviceIds[i]<< endl;
	}
	free(deviceIds);
	/////////////////////////////////////////////////////////////////
	// Create an OpenCL command queue
	/////////////////////////////////////////////////////////////////
	cout << "Running on Device 0..." << endl;
	commandQueue = clCreateCommandQueue(context, devices[0], CL_QUEUE_PROFILING_ENABLE, &status);
	SafeCall(status, "Creating Command Queue. (clCreateCommandQueue)");

	return 0;
}
int COpenCLTemplate::AllocateMemoryGPU()
{
	cl_int status;

	hy_gpu = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*SIZE, hy, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create input buffer");

	ez_gpu = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*SIZE, ez, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");
	
	mu_gpu = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*SIZE, mu, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");

	epsilon_gpu = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*SIZE, epsilon, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");

	Etemp_gpu = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*maxTime, Etemp, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");

	Exz1_gpu = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*maxTime, Exz1, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");

	Exz2_gpu = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*maxTime, Exz2, &status);
	SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");

	return 0;
}
int COpenCLTemplate::InitialiseCLKernelsGPU()
{
	int status;
	/////////////////////////////////////////////////////////////////
	// Load CL file, build CL program object, create CL kernel object
	/////////////////////////////////////////////////////////////////
	const char *filename = "OpenCLTemplate_Kernels.cl";
	string sourceStr = convertToString(filename);
	const char *source = sourceStr.c_str();
	size_t sourceSize[] = {strlen(source)};

	program = clCreateProgramWithSource( context, 1, &source, sourceSize, &status);
	SafeCall(status, "Error: Loading Binary into cl_program (clCreateProgramWithBinary)\n");

	// Create a cl program executable for all the devices specified.
	status = clBuildProgram(program, 1, devices, NULL, NULL, NULL);
	if(status == CL_BUILD_PROGRAM_FAILURE)
	{
		cl_int logStatus;
		char *buildLog = NULL;
		size_t buildLogSize = 0;
		logStatus = clGetProgramBuildInfo (program, devices[0], CL_PROGRAM_BUILD_LOG, buildLogSize, buildLog, &buildLogSize);
		SafeCall(logStatus, "clGetProgramBuildInfo failed.");
		buildLog = new char[buildLogSize];
		SafeCall(!buildLog, "Failed to allocate host memory. (buildLog)");
		memset(buildLog, 0, buildLogSize);
		logStatus = clGetProgramBuildInfo (program, devices[0], CL_PROGRAM_BUILD_LOG, buildLogSize, buildLog, NULL);
		if (logStatus != CL_SUCCESS)
		{
			cout << "clGetProgramBuildInfo failed." << endl;
			free(buildLog);
			return -1;
		}
		// Displaying build log in case of errors.
		cout << " \n\t\t\tBUILD LOG\n";
		cout << " ************************************************\n";
		cout << buildLog << endl;
		cout << " ************************************************\n";
		delete []buildLog;
	}

	// Attach kernel objects to respective kernel functions.
	hykernel = clCreateKernel(program, "hy_kernel", &status);
	SafeCall(status, "Error: Creating Kernel from program. (hy_kernel)");
	
	ezkernel = clCreateKernel(program, "ez_kernel", &status);
	SafeCall(status, "Error: Creating Kernel from program. (ez_kernel)");

	// ====== Set appropriate arguments to the kernel ======
	SafeCall(clSetKernelArg(hykernel, 0, sizeof(cl_mem), (void*)&hy_gpu), "Error: Setting kernel argument 'hy'");
	SafeCall(clSetKernelArg(hykernel, 1, sizeof(cl_mem), (void*)&ez_gpu), "Error: Setting kernel argument 'ez'");
	SafeCall(clSetKernelArg(hykernel, 2, sizeof(cl_mem), (void*)&mu_gpu), "Error: Setting kernel argument 'mu'");
	SafeCall(clSetKernelArg(hykernel, 3, sizeof(PRECISION), (void*)&delt), "Error: Setting kernel argument 'delt'");
	SafeCall(clSetKernelArg(hykernel, 4, sizeof(PRECISION), (void*)&delx), "Error: Setting kernel argument 'delx'");
	SafeCall(clSetKernelArg(hykernel, 5, sizeof(int), (void*)&SIZE), "Error: Setting kernel argument 'SIZE'");

	SafeCall(clSetKernelArg(ezkernel, 0, sizeof(cl_mem), (void*)&hy_gpu), "Error: Setting kernel argument 'hy'");
	SafeCall(clSetKernelArg(ezkernel, 1, sizeof(cl_mem), (void*)&ez_gpu), "Error: Setting kernel argument 'ez'");
	SafeCall(clSetKernelArg(ezkernel, 2, sizeof(cl_mem), (void*)&epsilon_gpu), "Error: Setting kernel argument 'mu'");
	SafeCall(clSetKernelArg(ezkernel, 3, sizeof(PRECISION), (void*)&delt), "Error: Setting kernel argument 'delt'");
	SafeCall(clSetKernelArg(ezkernel, 4, sizeof(PRECISION), (void*)&delx), "Error: Setting kernel argument 'delx'");
	SafeCall(clSetKernelArg(ezkernel, 5, sizeof(int), (void*)&SIZE), "Error: Setting kernel argument 'SIZE'");
	
	SafeCall(clSetKernelArg(ezkernel, 6, sizeof(cl_mem), (void*)&Etemp_gpu), "Error: Setting kernel argument 'Etemp'");
	SafeCall(clSetKernelArg(ezkernel, 7, sizeof(cl_mem), (void*)&Exz1_gpu), "Error: Setting kernel argument 'Exz1'");
	SafeCall(clSetKernelArg(ezkernel, 8, sizeof(cl_mem), (void*)&Exz2_gpu), "Error: Setting kernel argument 'Exz2'");
	SafeCall(clSetKernelArg(ezkernel, 9, sizeof(int), (void*)&qTime), "Error: Setting kernel argument 'qTime'");
	SafeCall(clSetKernelArg(ezkernel, 10, sizeof(PRECISION), (void*)&mur), "Error: Setting kernel argument 'mur'");
	SafeCall(clSetKernelArg(ezkernel, 11, sizeof(PRECISION), (void*)&epsilonr), "Error: Setting kernel argument 'epsilonr'");
	SafeCall(clSetKernelArg(ezkernel, 12, sizeof(PRECISION), (void*)&Sc), "Error: Setting kernel argument 'Sc'");
	SafeCall(clSetKernelArg(ezkernel, 13, sizeof(PRECISION), (void*)&ez1q), "Error: Setting kernel argument 'ez1q'");
	SafeCall(clSetKernelArg(ezkernel, 14, sizeof(PRECISION), (void*)&ez2q), "Error: Setting kernel argument 'ez2q'");
	SafeCall(clSetKernelArg(ezkernel, 15, sizeof(PRECISION), (void*)&ezmq), "Error: Setting kernel argument 'ezmq'");
	SafeCall(clSetKernelArg(ezkernel, 16, sizeof(PRECISION), (void*)&ezm1q), "Error: Setting kernel argument 'ezm1q'");
	SafeCall(clSetKernelArg(ezkernel, 17, sizeof(PRECISION), (void*)&pi), "Error: Setting kernel argument 'pi'");
	SafeCall(clSetKernelArg(ezkernel, 18, sizeof(PRECISION), (void*)&f), "Error: Setting kernel argument 'f'");
	SafeCall(clSetKernelArg(ezkernel, 19, sizeof(int), (void*)&SourceSelect), "Error: Setting kernel argument 'SourceSelect'");
	SafeCall(clSetKernelArg(ezkernel, 20, sizeof(PRECISION), (void*)&PulseWidth), "Error: Setting kernel argument 'PulseWidth'");
	SafeCall(clSetKernelArg(ezkernel, 21, sizeof(int), (void*)&Z1), "Error: Setting kernel argument 'Z1'");
	SafeCall(clSetKernelArg(ezkernel, 22, sizeof(int), (void*)&Z2), "Error: Setting kernel argument 'Z2'");
	return 0;
}
int COpenCLTemplate::RunCLKernels()
{
	//////////////////////////////	
  int mm;
  int SourceSelect = 1; 			// 0=Sinosoidal, 1=Gauassian
/*	cout<<"----Select Source----"<<endl;
	cout<<"0)Sinosoidal 1) Gauassian"<<endl;
	do{
		cin>>SourceSelect;
	}while(SourceSelect != 1 && SourceSelect != 0);
	if (SourceSelect == 0)
		maxTime = 1000;
*/
//////// my implementation
	
	// -------- refractive index variables -------- 
	PRECISION z1 = Z1*delx;
	PRECISION z2 = Z2*delx;
	////////////////////////////////////
	cl_int status;
	cl_uint maxDims;
	cl_event events[2];
	size_t globalThreads[2];
	size_t localThreads[2];
	size_t maxWorkGroupSize;
	size_t maxWorkItemSizes[3];

	// Query device capabilities. Maximum work item dimensions and the maximmum work item sizes
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), (void*)&maxWorkGroupSize, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), (void*)&maxDims, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");
	SafeCall(clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*maxDims, (void*)maxWorkItemSizes, NULL), "Error: Getting Device Info. (clGetDeviceInfo)");

  const int TILE_SIZE = 256; 
  globalThreads[0] = SIZE;//TILE_SIZE*ceil(((float)SIZE - 1)/TILE_SIZE);
	globalThreads[1] = 1;
	localThreads[0]  = 256;
	localThreads[1]  = 1;

	//cout << "Max dimensions: " << maxDims << endl;
	//cout << "Device maxWorkGroupSize = " << maxWorkGroupSize << endl;
	//cout << "Device maxWorkItemSizes = " << maxWorkItemSizes[0] << endl;
	if(localThreads[0] > maxWorkGroupSize || localThreads[0] > maxWorkItemSizes[0])
	{
		cout<<"Unsupported: Device does not support requested number of work items." << endl;
		return 1;
	}

	// Kernel timing variables.
	cl_ulong startTime, endTime;
	cl_ulong kernelExecTimeNs;
	cl_ulong kernelExecTimeNsT = 0;

	//cout << "Launching CL Kernel..." << endl;
	//cout << "Global threads: " << globalThreads[0] << "x" << globalThreads[1] << endl;
	//cout << "Local threads: " << localThreads[0] << "x" << localThreads[1] << endl;

	for (int medium=1; medium<=2; medium++)
	{  
	
		// -------- Medium Specifications -------- 
		if (medium==1)
		{
			cout<<"Calculating Wave propagation in Free Space"<<endl;
		}
		else
		{
			cout<<"Calculating Wave propagation in denser Medium"<<endl;
			int j;
			for(j=0; j<SIZE-(SIZE/2); j++)		//epsilon=[8.8542e-012*ones(1,SIZE-500) 1.7708e-011*ones(1,500)]; // half medium
				epsilon[j] = 8.8542e-012;
			for(int k=j ;k<SIZE ;k++)
				epsilon[j] = 1.7708e-011;

	    epsilon_gpu = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(PRECISION)*SIZE, epsilon, &status);
	    SafeCall(status, "Error: clCreateBuffer() cannot create output buffer");
		}

    for (int qTime=0; qTime<(maxTime-1); qTime++) {
		SafeCall(clSetKernelArg(ezkernel, 9, sizeof(int), (void*)&qTime), "Error: Setting kernel argument 'qTime'");

	    status = clEnqueueNDRangeKernel(commandQueue, hykernel, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
	    if(status != CL_SUCCESS) 
	    { 
	    	cout << "Error: Enqueueing kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
	    	if ( status == CL_INVALID_COMMAND_QUEUE ) cout << "CL_INVALID_COMMAND_QUEUE." << endl;
	    	if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
	    	if ( status == CL_INVALID_KERNEL ) cout << "CL_INVALID_KERNEL." << endl;
	    	if ( status == CL_INVALID_WORK_DIMENSION ) cout << "CL_INVALID_WORK_DIMENSION." << endl;
	    	if ( status == CL_INVALID_CONTEXT ) cout << "CL_INVALID_CONTEXT." << endl;
	    	if ( status == CL_INVALID_KERNEL_ARGS ) cout << "CL_INVALID_KERNEL_ARGS." << endl;
	    	if ( status == CL_INVALID_WORK_GROUP_SIZE ) cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
	    	if ( status == CL_INVALID_WORK_ITEM_SIZE ) cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
	    	if ( status == CL_INVALID_GLOBAL_OFFSET ) cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
	    	return 1;
	    }

	    // Wait for the kernel call to finish execution.
	    SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");
	    SafeCall(clReleaseEvent(events[0]), "Error: Release event object. (clReleaseEvent)\n");

      /////////////////////////////////////////////// 2nd kernel ///////////////////////////////////////////
	    status = clEnqueueNDRangeKernel(commandQueue, ezkernel, 2, NULL, globalThreads, localThreads, 0, NULL, &events[0]);
	    if(status != CL_SUCCESS) 
	    { 
	    	cout << "Error: Enqueueing kernel onto command queue (clEnqueueNDRangeKernel)" << endl;
	    	if ( status == CL_INVALID_COMMAND_QUEUE ) cout << "CL_INVALID_COMMAND_QUEUE." << endl;
	    	if ( status == CL_INVALID_PROGRAM_EXECUTABLE ) cout << "CL_INVALID_PROGRAM_EXECUTABLE." << endl;
	    	if ( status == CL_INVALID_KERNEL ) cout << "CL_INVALID_KERNEL." << endl;
	    	if ( status == CL_INVALID_WORK_DIMENSION ) cout << "CL_INVALID_WORK_DIMENSION." << endl;
	    	if ( status == CL_INVALID_CONTEXT ) cout << "CL_INVALID_CONTEXT." << endl;
	    	if ( status == CL_INVALID_KERNEL_ARGS ) cout << "CL_INVALID_KERNEL_ARGS." << endl;
	    	if ( status == CL_INVALID_WORK_GROUP_SIZE ) cout << "CL_INVALID_WORK_GROUP_SIZE." << endl;
	    	if ( status == CL_INVALID_WORK_ITEM_SIZE ) cout << "CL_INVALID_WORK_ITEM_SIZE." << endl;
	    	if ( status == CL_INVALID_GLOBAL_OFFSET ) cout << "CL_INVALID_GLOBAL_OFFSET." << endl;
	    	return 1;
	    }

	    // Wait for the kernel call to finish execution.
	    SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");
	    SafeCall(clReleaseEvent(events[0]), "Error: Release event object. (clReleaseEvent)\n");

      //// Copy data back to host ////
      SafeCall(clEnqueueReadBuffer(commandQueue, ez_gpu, CL_TRUE, 0,  sizeof(PRECISION)*SIZE, ez, 0, NULL, NULL), "Error reading ez back to host memory");    
      SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");
	  /////////////////////////////////////////////////////////////////////////
  

         // -------- Saving to file -------- 
			stream.str(std::string());   						// clear stringstream
			stream<<"./results/"<<"Efield"<<medium<<"_"<<qTime+1<<".jd";   		// concatenate
			filename = stream.str();		 					// copy string
			snapshot.open(filename.c_str(), ios::out|ios::binary);
			for (mm = 0; mm < SIZE; mm++)
				snapshot.write((char *)&ez[mm],sizeof(float));
			snapshot.close();
	 
      // Copy ez data to gpu
     // SafeCall(clEnqueueWriteBuffer(commandQueue, ez_gpu, CL_TRUE, 0,  sizeof(PRECISION)*SIZE, ez, 0, NULL, NULL), "Error writing ez back to GPU");    

    }
	SafeCall(clEnqueueReadBuffer(commandQueue, Etemp_gpu, CL_TRUE, 0,  sizeof(PRECISION)*maxTime, Etemp, 0, NULL, NULL), "Error reading ez back to host memory");    
    SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");
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
  }

	cout<<"Writing Values to files"<<endl;
	stream.str(std::string());
	stream<<"./results/"<<"Eincident"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Eincident[mm],(sizeof(float)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"Etransmitted"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Etransmitted[mm],(sizeof(float)));
	snapshot.close();

	SafeCall(clEnqueueReadBuffer(commandQueue, Exz1_gpu, CL_TRUE, 0,  sizeof(PRECISION)*maxTime, Exz1, 0, NULL, NULL), "Error reading Exz1 back to host memory");    
    SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");
	stream.str(std::string());
	stream<<"./results/"<<"Exz1"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Exz1[mm],(sizeof(float)));
	snapshot.close();

	SafeCall(clEnqueueReadBuffer(commandQueue, Exz2_gpu, CL_TRUE, 0,  sizeof(PRECISION)*maxTime, Exz2, 0, NULL, NULL), "Error reading Exz2 back to host memory");    
    SafeCall(clWaitForEvents(1, &events[0]), "Error: Waiting for kernel run to finish. (clWaitForEvents)");
	stream.str(std::string());
	stream<<"./results/"<<"Exz2"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	for (int mm = 0; mm < maxTime; mm++)
	snapshot.write((char *)&Exz2[mm],(sizeof(float)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"maxTime"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	snapshot.write((char *)&maxTime,(sizeof(int)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"SIZE"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	snapshot.write((char *)&SIZE,(sizeof(int)));
	snapshot.close();

	stream.str(std::string());
	stream<<"./results/"<<"data"<<".jd";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::binary);
	snapshot.write((char *)&delt,(sizeof(float)));
	snapshot.write((char *)&k0,(sizeof(float)));
	snapshot.write((char *)&z1,(sizeof(PRECISION)));
	snapshot.write((char *)&z2,(sizeof(PRECISION)));
	snapshot.close();
	
	return 0;
}
int COpenCLTemplate::CompleteRun()
{
	SafeCall(AllocateMemoryCPU(), "Error: Allocating memory on CPU.");
	SafeCall(InitialiseCPU(), "Error: Initialising data on CPU.");
	SafeCall(InitialiseCL(), "Error: Initialiasing CL.");
	SafeCall(AllocateMemoryGPU(), "Error: Allocating memory on GPU.");
	SafeCall(InitialiseCLKernelsGPU(), "Error: Copying data from CPU to GPU.");
	SafeCall(RunCLKernels(), "Error: Running kernels (GPU).");
	SafeCall(CleanupCPU(), "Error: Cleaning up CPU.");
	SafeCall(CleanupCL(), "Error: Cleaning up CL.");
	SafeCall(CleanupGPU(), "Error: Cleaning up GPU.");
  
	return 0;
}
// Converts contents of a file into a string. From OPENCL examples.
string COpenCLTemplate::convertToString(const char *filename)
{
	size_t size;
	char*  str;
	string s;
	fstream f(filename, (fstream::in | fstream::binary));

	if(f.is_open())
	{
		size_t fileSize;
		f.seekg(0, fstream::end);
		size = fileSize = (size_t)f.tellg();
		f.seekg(0, fstream::beg);

		str = new char[size+1];
		if(!str)
		{
			f.close();
			return NULL;
		}

		f.read(str, fileSize);
		f.close();
		str[size] = '\0';

		s = str;
		delete[] str;
		return s;
	}
	else
	{
		cout << "\nFile containg the kernel code(\".cl\") not found. Please copy the required file in the folder containg the executable." << endl;
		exit(1);
	}
	return NULL;
}
// Timing.
void COpenCLTemplate::StartTimer()
{
	if (tPaused == true)
	{
		tStart = GetTimeus64();
		tPaused = false;
	}
}
void COpenCLTemplate::StopTimer()
{
	if (tPaused == false)
	{
		tEnd = GetTimeus64();
		tDelta += tEnd - tStart;
		tStart = tEnd;
		tPaused = true;
	}
}
void COpenCLTemplate::ResetTimer()
{
	if (tPaused == true)
		tStart = tEnd;
	else
		tStart = GetTimeus64();

	tDelta = 0UL;
}
PRECISION COpenCLTemplate::GetElapsedTime()
{
	if (tPaused == false)
		tEnd = GetTimeus64();

	return ((double)(tEnd-tStart+tDelta))/(1000000.);
}
int COpenCLTemplate::SafeCall(cl_int Status, const char *Error)
{
	if (Status != 0)
	{
		if (Error!=NULL) cout << Error << endl;
		exit(Status);
	}
	return Status;
}
template<typename T> void DeleteArray(T *&ptr)
{
	if (ptr != NULL)
	{
		delete[] ptr;
		ptr = NULL;
	}
}
int COpenCLTemplate::CleanupCPU()
{
	DeleteArray(Eincident);
	DeleteArray(Etransmitted);
	DeleteArray(Etemp);
	DeleteArray(Exz1);
	DeleteArray(Exz2);
	DeleteArray(mu);
	DeleteArray(epsilon);
	DeleteArray(ez);
	DeleteArray(hy);

	return 0;
}
int COpenCLTemplate::CleanupCL()
{
	SafeCall(clReleaseKernel(hykernel), "Error: In clReleaseKernel");
	SafeCall(clReleaseKernel(ezkernel), "Error: In clReleaseKernel");
	SafeCall(clReleaseProgram(program), "Error: In clReleaseProgram");
	SafeCall(clReleaseCommandQueue(commandQueue), "Error: In clReleaseCommandQueue");
	SafeCall(clReleaseContext(context), "Error: In clReleaseContext");

	return 0;
}
int COpenCLTemplate::CleanupGPU()
{
	SafeCall(clReleaseMemObject(hy_gpu), "Error: clReleaseMemObject() cannot release input memory buffer");
	SafeCall(clReleaseMemObject(ez_gpu), "Error: clReleaseMemObject() cannot release output memory buffer");
	SafeCall(clReleaseMemObject(mu_gpu), "Error: clReleaseMemObject() cannot release input memory buffer");
	SafeCall(clReleaseMemObject(epsilon_gpu), "Error: clReleaseMemObject() cannot release output memory buffer");

	return 0;
}
COpenCLTemplate::~COpenCLTemplate ()
{
	// Field arrays.
	DeleteArray(Eincident);
	DeleteArray(Etransmitted);
	DeleteArray(Etemp);
	DeleteArray(Exz1);
	DeleteArray(Exz2);
	DeleteArray(mu);
	DeleteArray(epsilon);
	DeleteArray(ez);
	DeleteArray(hy);
}
