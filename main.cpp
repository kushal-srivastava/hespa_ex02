#include<iostream>
#include<vector>
#include<fstream>

int main()
{
	//parameters to be read from the input file
	double timestep_length = 0.01;
	double time_end = 1.0;
	double epsilon = 5.0;
	double sigma = 1.0;
	unsigned int N = 1000;

	//initialization of the variables
	std::vector<double> x;
	x.resize(3*N);
	std::vector<double> v;
	v.resize(3*N);
	std::vector<double> a;
	a.resize(3*N);
	std::vector<double> F;
	F.resize(3*N);
	std::vector<double> F_old;

	//reading the data from the input file for test

	//initializing pointers for cudaMalloc
	double *d_x;
	double *d_v;
	double *d_a;
	double *d_F;
	double *d_F_old;
	
	
	//kernel for calculation of force
	cudaMalloc((void**)&d_x, (3*N*sizeof(double)));
	cudaMalloc((void**)&d_v, 3*N*sizeof(double));
	cudaMalloc((void**)&d_a, 3*N*sizeof(double));
	cudaMalloc((void**)&d_F, 3*N*sizeof(double));
	cudaMalloc((void**)&d_F_old, 3*N*sizeof(double));

	//memcopy from host to device
	cudaMemcpy(d_x,x,3*N*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(d_a,a,3*N*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(d_v,v,3*N*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(d_F,F,3*N*sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(d_F_old,F_old,3*N*sizeof(double),cudaMemcpyHostToDevice);

	//call kern_Force for calculating the force between particles


	//call kern_Vel_Verlet for updating the velocities
	
	//memcopy from device back to host
	cudaMemcpy(x,d_x,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(a,d_a,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(v,d_v,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(F,d_F,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(F_old,d_F_old,3*N*sizeof(double),cudaMemcpyHostToDevice);

	//free device memory after simulation
	cudaFree(x_d);
	cudaFree(v_d);
	cudaFree(F_d);
	cudaFree(F_old_d);
	cudaFree(a_d);
	cudaFree(

	return 0;
}
