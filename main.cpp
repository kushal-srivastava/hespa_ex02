#include<iostream>
#include<vector>
#include<fstream>

int main()
{
	double timestep_length = 0.01;
	double time_end = 1.0;
	double epsilon = 5.0;
	double sigma = 1.0;
	unsigned int N = 1000;
	std::vector<double> x;
	x.resize(3*N);
	std::vector<double> v;
	v.resize(3*N);
	std::vector<double> a;
	a.resize(3*N);
	std::vector<double> F;
	F.resize(3*N);
	
	
	//kernel for calculation of force
	cudaMalloc()


}
