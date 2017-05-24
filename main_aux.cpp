#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>

/*__global__ void kern_Force(     double *epsilon,
				double *sigma,
				double *d_F_old,
				double *d_F,
				double *d_x){
	////+++++++++++++++++++++++++++++++++++++++++
	int x_index = threadIdx.x + blockIdx.x*blockDim.x;
	int y_index = threadIdx.y + blockIdx.y*blockDim.y;
	int z_index = threadIdx.z + blockIdx.z*blockDim.z;
    double  d_x = 0,
            abs_d_x = 0;
    
        

}*/
int main()
{
	//parameters to be read from the input file
	double timestep_length = 0.01;
	double time_end = 1.0;
	double epsilon = 5.0;
	double sigma = 1.0;
	unsigned int N = 2;
    
    
    double *mass  = new double[N];
    double *pos_x = new double[N];
    double *pos_y = new double[N];
    double *pos_z = new double[N];
    double *vel_x = new double[N];
    double *vel_y = new double[N];
    double *vel_z = new double[N];
        
    double *F_x = new double[N];
    double *F_y = new double[N];
    double *F_z = new double[N];
    double *F_old_x = new double[N];
    double *F_old_y = new double[N];
    double *F_old_z = new double[N];
    
    double d_x = 0,
            abs_d_x = 0,
            temp = 0,
            t = 0;
    
    double *pos[3] = {pos_x, pos_y, pos_z} ;
    double *vel[3] = {vel_x, vel_y, vel_z};
    double *Fdim[3] = {F_x, F_y,F_z};
    double *Fdim_old[3] = {F_old_x, F_old_y,F_old_z};
    //calculating force
    for(auto dim = 0;dim<3;++dim)
    {
        for(auto index = 0;index<N;++index)
        {
            for(auto z=0;z<N;++z)    {
                    for(auto y=0;y<N;++y)
                    {
                        for(auto x=0;x<N;++x){
                            
                            if(index != z*N*N + y*N + x)
                            {
                                d_x = pos[dim][index] - pos[dim][z*N*N + y*N + x];
                                abs_d_x = fabs(d_x);
                                temp = temp + (((24*epsilon/(abs_d_x * abs_d_x)) * pow(sigma/abs_d_x,6) * (2*pow(sigma/abs_d_x,6) - 1)) * d_x); 
                            }
                            
                                
                    }}}
            //final force for one particle
            //std::cout<<"temp_: "<<temp_;
            
            Fdim[dim][index] = temp;
            temp = 0;
        }
    }
    //force update complete
    
    do
    {
        
//position update and parallely copy force to force_old
            for(auto dim = 0;dim<3;++dim)
        {
            for(auto index = 0;index<N;++index)
            {
                for(auto z=0;z<N;++z)    {
                        for(auto y=0;y<N;++y)
                        {
                            for(auto x=0;x<N;++x){
                                
                                if(index != z*N*N + y*N + x)
                                {
                                    pos[dim][index] = pos[dim][index] + timestep_length * (vel[dim][index]) + ((timestep_length*timestep_length/(2*mass[index])) * (Fdim[dim][index]));
                                    
                                    
                                }
                                Fdim_old[dim][index] = Fdim[dim][index];
                                    
                        }}}

            }
        }
       
///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
    for(auto dim = 0;dim<3;++dim)
    {
        for(auto index = 0;index<N;++index)
        {
            for(auto z=0;z<N;++z)    {
                    for(auto y=0;y<N;++y)
                    {
                        for(auto x=0;x<N;++x){
                            
                            if(index != z*N*N + y*N + x)
                            {
                                d_x = pos[dim][index] - pos[dim][z*N*N + y*N + x];
                                abs_d_x = fabs(d_x);
                                temp = temp + (((24*epsilon/(abs_d_x * abs_d_x)) * pow(sigma/abs_d_x,6) * (2*pow(sigma/abs_d_x,6) - 1)) * d_x); 
                            }
                            
                                
                    }}}
            //final force for one particle
            //std::cout<<"temp_: "<<temp_;
            
            Fdim[dim][index] = temp;
            temp = 0;
        }
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(auto dim = 0;dim<3;++dim)
        {
            for(auto index = 0;index<N;++index)
            {
                for(auto z=0;z<N;++z)    {
                        for(auto y=0;y<N;++y)
                        {
                            for(auto x=0;x<N;++x){
                                
                                
                                vel[dim][index] = vel[dim][index] + timestep_length * (Fdim[dim][index] + Fdim_old[dim][index])/(2*mass[index]);
                                                                       
                                
                                
                                    
                        }}}

            }
        }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        t = t + 0.01;
    }while(t<time_end);
				
    
    delete(pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,F_x,F_y,F_z,F_old_x,F_old_y,F_old_z);
    
	return 0;
}

