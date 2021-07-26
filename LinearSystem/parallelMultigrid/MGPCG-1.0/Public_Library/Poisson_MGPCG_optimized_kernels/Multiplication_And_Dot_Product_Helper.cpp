//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Multiplication_And_Dot_Product_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Multiplication_And_Dot_Product,(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u,T* const v,const T* const diagonal_part,const T scale_factor),(x_size_input,u,v,diagonal_part,scale_factor))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Multiplication_And_Dot_Product_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Multiplication_And_Dot_Product_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,partition_number;
    Multiplication_And_Dot_Product_Size_Specific_Thread_Helper(
        Multiplication_And_Dot_Product_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int partition_number_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),partition_number(partition_number_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax,partition_number);}
};
}
template<class T,int y_size,int z_size> double Multiplication_And_Dot_Product_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    u_dot_v_partial_results=new double[number_of_partitions];

    int number_of_x_blocks=x_size/x_block_size;
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
	Multiplication_And_Dot_Product_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Multiplication_And_Dot_Product_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,partition);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();

    double u_dot_v=0;
    for(int i=0;i<number_of_partitions;i++) u_dot_v+=u_dot_v_partial_results[i];
    delete[] u_dot_v_partial_results;u_dot_v_partial_results=0;

    return u_dot_v;    
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T,int y_size,int z_size> void Multiplication_And_Dot_Product_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax,const int partition_number)
{   
    double local_u_dot_v=0;

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {
            int index=i*x_shift+j*y_shift+k*z_shift;
	    if(diagonal_part[index]==0) continue;
	    
	    v[index]=diagonal_part[index]*u[index];
	    if(diagonal_part[index+x_plus_one_shift])
		v[index]+=u[index+x_plus_one_shift];
	    if(diagonal_part[index+x_minus_one_shift])
		v[index]+=u[index+x_minus_one_shift];
	    if(diagonal_part[index+y_plus_one_shift])
		v[index]+=u[index+y_plus_one_shift];
	    if(diagonal_part[index+y_minus_one_shift])
		v[index]+=u[index+y_minus_one_shift];
	    if(diagonal_part[index+z_plus_one_shift])
		v[index]+=u[index+z_plus_one_shift];
	    if(diagonal_part[index+z_minus_one_shift])
		v[index]+=u[index+z_minus_one_shift];
	    v[index]*=scale_factor;
	    local_u_dot_v+=v[index]*u[index];
	}
    
    u_dot_v_partial_results[partition_number]=local_u_dot_v;
}
//#####################################################################
template class Multiplication_And_Dot_Product_Helper<float>;
template class Multiplication_And_Dot_Product_Helper<double>;

