#include <OpenCLTemplate.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

int main(int argc, char * argv[])
{

	fstream snapshot;
	std::string filename ;
	std::stringstream stream;

  float SIZE = 0;
  float maxTime = 1024;
  float time=0;
  int i=0;
  while(i<=10)
  {
	  SIZE += 1024;
	COpenCLTemplate OpenCLTemplateSim(maxTime, SIZE);
	// ================== Simulation ================
	OpenCLTemplateSim.StartTimer();
	OpenCLTemplateSim.CompleteRun(); // Complete GPU run.
	OpenCLTemplateSim.StopTimer();
	cout << "Total time taken = " << OpenCLTemplateSim.GetElapsedTime() << " seconds." << endl;
	time=OpenCLTemplateSim.GetElapsedTime();
	stream.str(std::string());
	stream<<"./results/"<<"time"<<".txt";
	filename = stream.str();
	snapshot.open(filename.c_str(), ios::out|ios::app);
	snapshot.write((char *)&SIZE,(sizeof(int)));
	snapshot.write((char *)&maxTime,(sizeof(int)));
	snapshot.write((char *)&time,(sizeof(float)));
	snapshot.close();
	i++;
  }
	return 0;
}
