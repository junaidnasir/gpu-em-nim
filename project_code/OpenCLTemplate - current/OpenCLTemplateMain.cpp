#include <OpenCLTemplate.hpp>
#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{

  int SIZE = 1000;
  int maxTime = 1024;

	COpenCLTemplate OpenCLTemplateSim(maxTime, SIZE);
	// ================== Simulation ================
	
	OpenCLTemplateSim.CompleteRun(); // Complete GPU run.
	
	cout << "Total time taken = " << OpenCLTemplateSim.GetElapsedTime() << " seconds." << endl;

	return 0;
}
