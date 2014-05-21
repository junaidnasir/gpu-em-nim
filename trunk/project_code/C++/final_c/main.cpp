#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
//#include <windows.h>
using namespace std;
struct student{
char name;
char course;
int age;
double year;
};

int main()
{
	int s;
	student *s1 ;//= {'b', 'a', 22, 2002};
	s1= new student;
	s1->name='a';
cout << "sizeof s1 "<<s1->name<<endl;
//cout << "sizeof sptr "<<sizeof(sptr)<<endl;

//cout << "sizeof *sptr "<<sizeof(*sptr)<<endl;
	cin>>s;
return 0;
}

