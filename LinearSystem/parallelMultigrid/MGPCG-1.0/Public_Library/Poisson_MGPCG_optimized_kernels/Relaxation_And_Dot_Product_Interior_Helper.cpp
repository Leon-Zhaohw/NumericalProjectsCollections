//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Relaxation_And_Dot_Product_Interior_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Relaxation_And_Dot_Product_Interior,(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,const T* const b,T* const delta,const unsigned char* const bit_writemask),(x_size_input,u,b,delta,bit_writemask))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Relaxation_And_Dot_Product_Interior_Compute_Delta_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,partition_number;
    Relaxation_And_Dot_Product_Interior_Compute_Delta_Size_Specific_Thread_Helper(
        Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int partition_number_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),partition_number(partition_number_input) {}
    void Run(){obj->Compute_Delta_X_Range(xmin,xmax,partition_number);}
};
template<class T,int y_size,int z_size>
struct Relaxation_And_Dot_Product_Interior_Apply_Delta_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,partition_number;
    Relaxation_And_Dot_Product_Interior_Apply_Delta_Size_Specific_Thread_Helper(
        Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int partition_number_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),partition_number(partition_number_input) {}
    void Run(){obj->Apply_Delta_X_Range(xmin,xmax,partition_number);}
};
}
template<class T,int y_size,int z_size> double Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) return Run(); // Run serially for small problems

    u_dot_b_partial_results=new double[number_of_partitions];
 
    int number_of_x_blocks=x_size/x_block_size;
    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_x_blocks%number_of_partitions,partition-1)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*partition+std::min(number_of_x_blocks%number_of_partitions,partition);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Relaxation_And_Dot_Product_Interior_Compute_Delta_Size_Specific_Thread_Helper<T,y_size,z_size>* task=
            new Relaxation_And_Dot_Product_Interior_Compute_Delta_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,partition-1);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();

    int number_of_coarse_x_blocks=coarse_x_size/x_block_size;
    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_coarse_x_blocks/number_of_partitions)*(partition-1)+std::min(number_of_coarse_x_blocks%number_of_partitions,partition-1)+1;
        int last_block_of_partition=(number_of_coarse_x_blocks/number_of_partitions)*partition+std::min(number_of_coarse_x_blocks%number_of_partitions,partition);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Relaxation_And_Dot_Product_Interior_Apply_Delta_Size_Specific_Thread_Helper<T,y_size,z_size>* task=
            new Relaxation_And_Dot_Product_Interior_Apply_Delta_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,partition-1);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();

    double u_dot_b=0;
    for(int i=0;i<number_of_partitions;i++) u_dot_b+=u_dot_b_partial_results[i];
    delete[] u_dot_b_partial_results;u_dot_b_partial_results=0;
    return u_dot_b;
}
//#####################################################################
// Function Compute_Delta_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>::
Compute_Delta_X_Range(const int xmin,const int xmax,const int partition_number)
{    
    const T one_sixth=1./6.;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {
            int index=i*x_shift+j*y_shift+k*z_shift;

            delta[index]=one_sixth*(
                u[index+x_plus_one_shift]
                +u[index+x_minus_one_shift]
                +u[index+y_plus_one_shift]
                +u[index+y_minus_one_shift]
                +u[index+z_plus_one_shift]
                +u[index+z_minus_one_shift])
                -u[index];
        }
}
//#####################################################################
// Function Apply_Delta_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_And_Dot_Product_Interior_Size_Specific_Helper<T,y_size,z_size>::
Apply_Delta_X_Range(const int xmin,const int xmax,const int partition_number)
{    
    const T two_thirds=2./3.;
    const T minus_one_ninth=-1./9.;
    double local_u_dot_b=0;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=coarse_y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=coarse_z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {    
            int index=2*(i*x_shift+j*y_shift+k*z_shift)-x_shift-y_shift-z_shift;
            int coarse_index=i*coarse_x_shift+j*coarse_y_shift+k*coarse_z_shift;
            
            if(bit_writemask[coarse_index] & 0x01){
                u[index]+=two_thirds*delta[index];
                u[index]+=minus_one_ninth*b[index];}
	    local_u_dot_b+=u[index]*b[index];
            if(bit_writemask[coarse_index] & 0x02){
                u[index+z_shift]+=two_thirds*delta[index+z_shift];
                u[index+z_shift]+=minus_one_ninth*b[index+z_shift];}
	    local_u_dot_b+=u[index+z_shift]*b[index+z_shift];
            if(bit_writemask[coarse_index] & 0x04){
                u[index+y_shift]+=two_thirds*delta[index+y_shift];
                u[index+y_shift]+=minus_one_ninth*b[index+y_shift];}
	    local_u_dot_b+=u[index+y_shift]*b[index+y_shift];
            if(bit_writemask[coarse_index] & 0x08){
                u[index+y_shift+z_shift]+=two_thirds*delta[index+y_shift+z_shift];
                u[index+y_shift+z_shift]+=minus_one_ninth*b[index+y_shift+z_shift];}
	    local_u_dot_b+=u[index+y_shift+z_shift]*b[index+y_shift+z_shift];
            if(bit_writemask[coarse_index] & 0x10){
                u[index+x_shift]+=two_thirds*delta[index+x_shift];
                u[index+x_shift]+=minus_one_ninth*b[index+x_shift];}
	    local_u_dot_b+=u[index+x_shift]*b[index+x_shift];
            if(bit_writemask[coarse_index] & 0x20){
                u[index+x_shift+z_shift]+=two_thirds*delta[index+x_shift+z_shift];
                u[index+x_shift+z_shift]+=minus_one_ninth*b[index+x_shift+z_shift];}
	    local_u_dot_b+=u[index+x_shift+z_shift]*b[index+x_shift+z_shift];
            if(bit_writemask[coarse_index] & 0x40){
                u[index+x_shift+y_shift]+=two_thirds*delta[index+x_shift+y_shift];
                u[index+x_shift+y_shift]+=minus_one_ninth*b[index+x_shift+y_shift];}
	    local_u_dot_b+=u[index+x_shift+y_shift]*b[index+x_shift+y_shift];
            if(bit_writemask[coarse_index] & 0x80){
                u[index+x_shift+y_shift+z_shift]+=two_thirds*delta[index+x_shift+y_shift+z_shift];
                u[index+x_shift+y_shift+z_shift]+=minus_one_ninth*b[index+x_shift+y_shift+z_shift];}
	    local_u_dot_b+=u[index+x_shift+y_shift+z_shift]*b[index+x_shift+y_shift+z_shift];
        }

    u_dot_b_partial_results[partition_number]=local_u_dot_b;
}
//#####################################################################
template class Relaxation_And_Dot_Product_Interior_Helper<float>;
template class Relaxation_And_Dot_Product_Interior_Helper<double>;
