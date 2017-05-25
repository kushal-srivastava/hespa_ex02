#include<iostream>
#include<fstream>
#include<string>
#include<cmath>

typedef double real;

std::string part_input_file,part_out_name_base,vtk_out_name_base;

real timestep_length,time_end,epsilon,sigma;

int part_out_freq,vtk_out_freq,cl_workgroup_1dsize;

void fileread(std::string file);

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

real d_x = 0,
            abs_d_x = 0,
            temp = 0,
            t = 0;
int i =0;
while (true) {
    f >> mass[i] >> pos_x[i] >> pos_y[i]>> pos_z[i] >> vel_x[i] >> vel_y[i] >> vel_z[i];
    if( f.eof() ) break;
    //std::cout<< mass[i] << " " <<  pos_x[i] <<" "<< pos_y[i] << " " << pos_z[i] << " " << vel_x[i] << " " << vel_y[i] <<" " << vel_z[i]<< std::endl;
    ++i;
}
f.close();
//*********************************************************************
//*********************************************************************
//*********************************************************************
std::string vtk_file =" ";
std::ofstream vtk;

real *pos[3] = {pos_x, pos_y, pos_z} ;
real *vel[3] = {vel_x, vel_y, vel_z};
real *Fdim[3] = {F_x, F_y,F_z};
real *Fdim_old[3] = {F_old_x, F_old_y,F_old_z};
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
            vtk_file  = "tmp/" + vtk_out_name_base + std::to_string(count) +".vtk";
    std::cout << vtk_file << std::endl;
    vtk.open(vtk_file);
    vtk << "# vtk DataFile Version 4.0" << "\n" << "hesp visualization file" << "\n" << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n" << "POINTS 2 double" << "\n";
    vtk<< std::fixed;
    for(int j =0; j<N; ++j)
        vtk<<pos_x[j] << " "<< pos_y[j]<< " " << pos_z[j] << "\n";
    vtk << "CELLS 0 0" << "\n";
    vtk << "CELL_TYPES 0" << "\n";
    vtk << "POINT_DATA 2" << "\n";
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
				







//*********************************************************************
///vtk file generation
//*********************************************************************
    /*
std::string vtk_file =" ";
std::ofstream vtk;
for (int i =0; i < 1 ; ++i){
    vtk_file  = vtk_out_name_base + std::to_string(i) +".vtk";
    std::cout << vtk_file << std::endl;
    vtk.open(vtk_file);
    vtk << "# vtk DataFile Version 4.0" << "\n" << "hesp visualization file" << "\n" << "ASCII" << "\n" << "DATASET UNSTRUCTURED_GRID" << "\n" << "POINTS 2 double" << "\n";
    vtk<< std::fixed;
    for(int j =0; j<N; ++j)
        vtk<<pos_x[j] << " "<< pos_y[j]<< " " << pos_z[j] << "\n";
    vtk << "CELLS 0 0" << "\n";
    vtk << "CELL_TYPES 0" << "\n";
    vtk << "POINT_DATA 2" << "\n";
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
}

*/

//*********************************************************************
//*********************************************************************
//*********************************************************************



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
            sigma=std::stoi  (value,&sz);;
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

