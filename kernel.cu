#include<cuda_runtime.h>
typedef double real;

__global__ void Force_1D(type_T h,
		real *d_pos_x,
		real *d_pos_y,
		real *d_pos_z,
		real *d_vel_x,
		real *d_vel_y,
		real *d_vel_z,
		real *d_F_x,
		real *d_F_y,
		real *d_F_z,
		real *d_F_old_x,
		real *d_F_old_y,
		real *d_F_old_z,
		real sigma,
		real epsilon,
		unsigned int N)
{


	
	
    	real d = 0.0,
            d_2 = 0.0,
            abs_d = 0.0,
            x_i = 0.0,
            y_i = 0.0,
            dx,dy,dz,
            z_i = 0.0,
            c1 = 0.0,
            t_pow = 0.0,
            sig_abs = 0.0,
            tempx = 0.0,
            tempy = 0.0,
            tempz = 0.0;
	int id = threadIdx.x + blockIdx.x*blockDim.x;
    	c1 = 24 * epsilon;
    
            x_i = d_pos_x[id];
            y_i = d_pos_y[id];
            z_i = d_pos_z[id];
            
            
              for(auto j=0;j<N;++j){
                            
                            if(i != j)
                            {
                                
                                d_2 = (x_i - d_pos_x[j]) * (x_i - d_pos_x[j]) + (y_i - d_pos_y[j])*(y_i - d_pos_y[j]) + (z_i - d_pos_z[j]) * (z_i - d_pos_z[j]);
                                //d = (x_i - pos_x[j])  + (y_i - pos_y[j]) + (z_i - pos_z[j]);
                                //std::cout<< i<<"\t" <<j<<"\t"<<"d: "<< d <<"\n"; 
                                d = sqrt(d_2);
                                //std::cout<< i<<"\t" <<j<<"\t"<<"d_2: "<< d_2 <<"\n"; 
                                //abs_d = fabs(d);
                                dx = x_i - d_pos_x[j];
                                dy = y_i - d_pos_y[j];
                                dz = z_i - d_pos_z[j];
                                //std::cout<< i<<"\t" <<j<<"\t"<<"abs_d: "<< abs_d <<"\n"; 
                                assert(d != 0);
                                sig_abs = sigma/d;
                                t_pow = pow(sig_abs,6);
                                //std::cout<< i<<"\t" <<j<<"\t"<<"weird calc: "<< ((c1/(d_2) * t_pow * (2*t_pow - 1)) * d) <<"\t" << "c1: " << c1<<"\n"; 
                                tempx = tempx + ((c1/(d_2) * t_pow * (2*t_pow - 1)) * dx); 
                                tempy = tempy + ((c1/(d_2) * t_pow * (2*t_pow - 1)) * dy); 
                                tempz = tempz + ((c1/(d_2) * t_pow * (2*t_pow - 1)) * dz); 
                                //std::cout<< i<<"\t" <<j<<"\t"<<"temp: "<< temp <<"\n";
                            }
                                
                                
                   }
               d_F_x[i] = tempx; 
               d_F_y[i] = tempy; 
               d_F_z[i] = tempz; 
               //std::cout <<"updated F: " << i << "\t" << F[i]<<"\n";
               /*tempx = 0;
               tempy = 0;
               tempz = 0;*/
                            
}


__global__ void pos_update_1D(
		real *d_pos_x,
		real *d_pos_y,
		real *d_pos_z,
		real *d_vel_x,
		real *d_vel_y,
		real *d_vel_z,
		real *d_F_x,
		real *d_F_y,
		real *d_F_z,
		real *d_F_old_x,
		real *d_F_old_y,
		real *d_F_old_z,
		real timestep_length){

			int i = threadIdx.x + blockDim.x * blockIdx.x;
		        real del_T = timestep_length*timestep_length;
                            d_pos_x[i] = d_pos_x[i] + timestep_length * (d_vel_x[i]) + ((del_T/(2*mass[i])) * (d_F_x[i]));
                            d_pos_y[i] = d_pos_y[i] + timestep_length * (d_vel_y[i]) + ((del_T/(2*mass[i])) * (d_F_y[i]));
                            d_pos_z[i] = d_pos_z[i] + timestep_length * (d_vel_z[i]) + ((del_T/(2*mass[i])) * (d_F_z[i]));
                            //std::cout << i <<"\t" << pos_x[i] << "\t" << pos_y[i] << "\t" << pos_z[i] <<"\n";
                            d_F_old_x[i] = d_F_x[i];
                            d_F_old_y[i] = d_F_y[i];
                            d_F_old_z[i] = d_F_z[i];

}

__global__ void vel_update_1D(
			real *d_vel_x,
		real *d_vel_y,
		real *d_vel_z,
		real *d_F_x,
		real *d_F_y,
		real *d_F_z,
		real *d_F_old_x,
		real *d_F_old_y,
		real *d_F_old_z,
		real timestep_length
		real *mass){

		int i = threadIdx.x + blockDim.x * blockIdx.x;
		vel_x[i] = vel_x[i] + timestep_length * 0.5 * (F_x[i] + F_old_x[i])/mass[i];
                vel_y[i] = vel_y[i] + timestep_length * 0.5 * (F_y[i] + F_old_y[i])/mass[i];
                vel_z[i] = vel_z[i] + timestep_length * 0.5 * (F_z[i] + F_old_z[i])/mass[i];

}















