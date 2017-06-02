#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<assert.h>

typedef double real;

std::string part_input_file,part_out_name_base,vtk_out_name_base;

real timestep_length,time_end,epsilon,sigma;

int part_out_freq,vtk_out_freq,cl_workgroup_1dsize;

void fileread(std::string file);
void force_update(real* pos_x, real* pos_y, real* pos_z, real* F_x, real* F_y, real* F_z, real sigma, real epsilon, unsigned int N);
int main(int argc,char *argv[])
{
        std::cout.precision(4);
        std::string para_file;
        if(argc==1)
        {
            std::cout<< "Para file name not given" << std::endl;
            return 0;
        }
        para_file = argv[1];
        //std::cout<< para_file<< " " << argc << std::endl; //parameter file name

        fileread(para_file); //rad parameter file

        unsigned int N=0; //Number of particles 
        //*********************************************************************
        ///Read from parameter file
        //*********************************************************************
        std::ifstream f;
        f.open(part_input_file); //read input file
        f>>N;
        //std::cout<< "Number of particles: " << N << std::endl;

        real *mass  = new real[N];
        real *pos_x = new real[N];
        real *pos_y = new real[N];
        real *pos_z = new real[N];
        real *vel_x = new real[N];
        real *vel_y = new real[N];
        real *vel_z = new real[N];
                
        real *F_x = new real[N];
        real *F_y = new real[N];
        real *F_z = new real[N];
        real *F_old_x = new real[N];
        real *F_old_y = new real[N];
        real *F_old_z = new real[N];
        int count = 0;

        //d: particle distance between ith and jth particle
        //d_2: sqare of d.
        //abs_d_x: absolute of (d_x)
        //t: time step increment
        //temp: for reduction of force
        
        real d = 0.0,
            d_2 = 0.0,
            abs_d = 0,
            temp = 0,
            t = 0;
        int iter =0;
        while (true) {
            f >> mass[iter] >> pos_x[iter] >> pos_y[iter]>> pos_z[iter] >> vel_x[iter] >> vel_y[iter] >> vel_z[iter];
            if( f.eof() ) break;
            //std::cout<< mass[i] << " " <<  pos_x[i] <<" "<< pos_y[i] << " " << pos_z[i] << " " << vel_x[i] << " " << vel_y[i] <<" " << vel_z[i]<< std::endl;
            ++iter;
        }
        f.close();
        //*********************************************************************
        //*********************************************************************
        //*********************************************************************
        std::string vtk_file =" ";
        std::ofstream vtk;

        real del_T = timestep_length*timestep_length;
        
        force_update(pos_x,pos_y,pos_z,F_x, F_y, F_z,sigma,epsilon,N);
        
            do
            {
                
                    //position update and parallely copy force to force_old
                            for(auto i = 0;i<N;++i)
                            {
                            pos_x[i] = pos_x[i] + timestep_length * (vel_x[i]) + ((del_T/(2*mass[i])) * (F_x[i]));
                            pos_y[i] = pos_y[i] + timestep_length * (vel_y[i]) + ((del_T/(2*mass[i])) * (F_y[i]));
                            pos_z[i] = pos_z[i] + timestep_length * (vel_z[i]) + ((del_T/(2*mass[i])) * (F_z[i]));
                            //std::cout << i <<"\t" << pos_x[i] << "\t" << pos_y[i] << "\t" << pos_z[i] <<"\n";
                            F_old_x[i] = F_x[i];
                            F_old_y[i] = F_y[i];
                            F_old_z[i] = F_z[i];
                            }
                            
                    //Force update
                           
                           force_update(pos_x,pos_y,pos_z,F_x,F_y,F_z,sigma,epsilon,N); 
                                
                    //calculate velocity 
                            
                            
                            for(auto i = 0;i<N;++i)
                            {
                            assert(mass[i] != 0);
                            vel_x[i] = vel_x[i] + timestep_length * 0.5 * (F_x[i] + F_old_x[i])/mass[i];
                            vel_y[i] = vel_y[i] + timestep_length * 0.5 * (F_y[i] + F_old_y[i])/mass[i];
                            vel_z[i] = vel_z[i] + timestep_length * 0.5 * (F_z[i] + F_old_z[i])/mass[i];
                            
                            }
                        
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        vtk_file  = "tmp/" + vtk_out_name_base + std::to_string(count) +".vtk";
                //        std::cout << vtk_file << std::endl;
                        vtk.open(vtk_file);
                        vtk << "# vtk DataFile Version 4.0" << "\n" << "hesp visualization file" << "\n" << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n" << "POINTS "<<N<<" double" << "\n";
                        vtk<< std::fixed;
                        for(int j =0; j<N; ++j)
                            vtk<<pos_x[j] << " "<< pos_y[j]<< " " << pos_z[j] << "\n";
                        vtk << "CELLS 0 0" << "\n";
                        vtk << "CELL_TYPES 0" << "\n";
                        vtk << "POINT_DATA "<<N << "\n";
                        vtk << "SCALARS m double" << "\n";
                        vtk << "LOOKUP_TABLE default" << "\n";
                        vtk<< std::fixed;
                        for(int j =0; j<N; ++j)
                            vtk<< mass[j] << "\n";
                        vtk <<"VECTORS v double" << "\n";
                        vtk<<  std::fixed;
                        for(int j =0; j<N; ++j)
                            vtk<<vel_x[j] << " "<< vel_y[j]<< " " << vel_y[j] << "\n";
                        vtk.close();
                        vtk_file =" ";
                            
                            
                            
                            count++;
                            t = t + 0.01;
            }while(t<time_end);
                        

        delete(mass);
        delete(pos_x);
        delete(pos_y);
        delete(pos_z);
        delete(vel_x);
        delete(vel_y);
        delete(vel_z);
        delete(F_x);
        delete(F_y);
        delete(F_z);
        delete(F_old_x);
        delete(F_old_y);
        delete(F_old_z);


        return 0;
}

       

void fileread(std::string file){
    std::ifstream ff;
    std::string para_name;
    std::string value;
    ff.open(file);
    std::string::size_type sz; 
    for(int i =0; i<10; ++i)
    {
        ff>>para_name >> value;
        if(para_name=="part_input_file"){
            part_input_file = value;
            //std::cout<< "part_input_file " << part_input_file<< std::endl;
        }
        else if(para_name=="timestep_length"){
            timestep_length= std::stod  (value,&sz);
            //std::cout<< "timestep_length " << timestep_length<< std::endl;
        }
        else if(para_name=="time_end"){
           time_end= std::stod  (value,&sz);
            //std::cout<< "time_end " << time_end<< std::endl; 
        }
        else if (para_name=="sigma"){
            sigma=std::stod  (value,&sz);
            //std::cout<< "sigma " << sigma<< std::endl;
        }
        else if (para_name=="epsilon"){
            epsilon=std::stod  (value,&sz);;
            //std::cout<< "epsilon " << epsilon<< std::endl;
        }
        else if (para_name=="part_out_freq"){
            part_out_freq=std::stoi  (value,&sz);
            //std::cout<< "part_out_freq " << part_out_freq<< std::endl;
            
        }
           else if (para_name=="part_out_name_base"){
            part_out_name_base=value;
            //std::cout<< "part_out_name_base " << part_out_name_base<< std::endl;
               
        }
        else if (para_name=="vtk_out_freq")
        {
            vtk_out_freq=std::stoi  (value,&sz);
            //std::cout<< "vtk_out_freq " << vtk_out_freq<< std::endl;
        }
        else if (para_name=="vtk_out_name_base"){
            vtk_out_name_base=value;
            //std::cout<< "vtk_out_name_base " << vtk_out_name_base<< std::endl;
        }
        else if (para_name=="cl_workgroup_1dsize"){
            cl_workgroup_1dsize=std::stoi  (value,&sz);
            //std::cout<< "cl_workgroup_1dsize " << cl_workgroup_1dsize<< std::endl;
        }
    }
    ff.close();
}

void force_update(real* pos_x, real* pos_y, real* pos_z, real* F_x, real* F_y, real* F_z, real sigma, real epsilon, unsigned int N)
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
    c1 = 24 * epsilon;
    //calculating force
    // index: particle index
    
        for(auto i = 0;i<N;++i)
        {
            x_i = pos_x[i];
            y_i = pos_y[i];
            z_i = pos_z[i];
            
            
              for(auto j=0;j<N;++j){
                            
                            if(i != j)
                            {
                                
                                d_2 = (x_i - pos_x[j]) * (x_i - pos_x[j]) + (y_i - pos_y[j])*(y_i - pos_y[j]) + (z_i - pos_z[j]) * (z_i - pos_z[j]);
                                //d = (x_i - pos_x[j])  + (y_i - pos_y[j]) + (z_i - pos_z[j]);
                                //std::cout<< i<<"\t" <<j<<"\t"<<"d: "<< d <<"\n"; 
                                d = sqrt(d_2);
                                //std::cout<< i<<"\t" <<j<<"\t"<<"d_2: "<< d_2 <<"\n"; 
                                //abs_d = fabs(d);
                                dx = x_i - pos_x[j];
                                dy = y_i - pos_y[j];
                                dz = z_i - pos_z[j];
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
               F_x[i] = tempx; 
               F_y[i] = tempy; 
               F_z[i] = tempz; 
               //std::cout <<tempx<<"\t"<<tempy<<"\t"<<tempz<<std::endl;// << i << "\t" << F[i]<<"\n";
               tempx = 0;
               tempy = 0;
               tempz = 0;
        }
        
    
}
