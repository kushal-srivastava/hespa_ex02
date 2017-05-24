#include<iostream>
#include<vector>
#include<fstream>
#include<math.h>

__global__ void kern_Force(     double *epsilon,
				double *sigma,
				double *d_F_old,
				double *d_F,
				double *d_x){
	////+++++++++++++++++++++++++++++++++++++++++
	int x_index = threadIdx.x + blockIdx.x*blockDim.x;
	int y_index = threadIdx.y + blockIdx.y*blockDim.y;
	int z_index = threadIdx.z + blockIdx.z*blockDim.z;
	
    for(k=0;k<N;++k){
        for(i=0;i<N;++i)
        {
            
            for(j=0;j<N;j++{
                        
                
            }
        
        }
    }

}
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
	std::vector<double> d_x;
	std::vector<double> d_v;
	std::vector<double> d_a;
	std::vector<double> d_F;
	std::vector<double> d_F_old;
	
	
	
	//kernel for calculation of force
	cudaMalloc((void**)&d_x, sizeof(std::vector<double>));
	d_x.resize(N*3);
	cudaMalloc((void**)&d_v, sizeof(std::vector<double>));
	d_v.resize(N*3);
	cudaMalloc((void**)&d_a, sizeof(std::vector<double>));
	d_a.resize(N*3);
	cudaMalloc((void**)&d_F, sizeof(std::vector<double>));
	d_F.resize(3*N);
	cudaMalloc((void**)&d_F_old, sizeof(std::vector<double>));
	d_F_old.resize(3*N);

	//memcopy from host to device
	cudaMemcpy(d_x,x,3*N*sizeof(std::vector<double>),cudaMemcpyHostToDevice);
	cudaMemcpy(d_a,a,3*N*sizeof(std::vector<double>),cudaMemcpyHostToDevice);
	cudaMemcpy(d_v,v,3*N*sizeof(std::vector<double>),cudaMemcpyHostToDevice);
	cudaMemcpy(d_F,F,3*N*sizeof(std::vector<double>),cudaMemcpyHostToDevice);
	cudaMemcpy(d_F_old,F_old,3*N*sizeof(std::vector<double>),cudaMemcpyHostToDevice);

	//call kern_Force for calculating the force between particles
	

	//call kern_Vel_Verlet for updating the velocities
	
	//memcopy from device back to host
	cudaMemcpy(x,d_x,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(a,d_a,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(v,d_v,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(F,d_F,3*N*sizeof(double),cudaMemcpyHostToDevice);
        cudaMemcpy(F_old,d_F_old,3*N*sizeof(double),cudaMemcpyHostToDevice);

	//free device memory after simulation
	cudaFree(d_x);
	cudaFree(d_v);
	cudaFree(d_a);
	cudaFree(d_F_old);
	cudaFree(d_F);
	

	return 0;
}
