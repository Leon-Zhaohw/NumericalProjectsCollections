//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
#include <limits>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema,(const int x_size_input,const int y_size_input,const int z_size_input,const T& c1,const T* const u1,T* const result,double &sum,T& minimum,T& maximum),(x_size_input,c1,u1,result,sum,minimum,maximum))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,partition_number;
    Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Thread_Helper(
        Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int partition_number_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),partition_number(partition_number_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax,partition_number);}
};
}
template<class T,int y_size,int z_size> void Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    sum_partial_results=new double[number_of_partitions];minimum_partial_results=new T[number_of_partitions];maximum_partial_results=new T[number_of_partitions];

    int number_of_x_blocks=x_size/x_block_size;
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
	Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,partition);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();
    
    sum=0;minimum=std::numeric_limits<T>::max();maximum=-std::numeric_limits<T>::max();
    
    for(int i=0;i<number_of_partitions;i++){
	sum+=sum_partial_results[i];
        if(minimum_partial_results[i]<minimum) minimum=minimum_partial_results[i];
        if(maximum_partial_results[i]>maximum) maximum=maximum_partial_results[i];}

    delete[] sum_partial_results;sum_partial_results=0;
    delete[] minimum_partial_results;minimum_partial_results=0;
    delete[] maximum_partial_results;maximum_partial_results=0;
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T,int y_size,int z_size> void Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax,const int partition_number)
{   
    double local_sum=0;
    T local_minimum=std::numeric_limits<T>::max(),local_maximum=-std::numeric_limits<T>::max();

    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {    
            int index=i*x_shift+j*y_shift+k*z_shift;
	    result[index]+=c1*u1[index];

            local_sum+=result[index];
            if(result[index]<local_minimum) local_minimum=result[index];
            if(result[index]>local_maximum) local_maximum=result[index];
	}

    sum_partial_results[partition_number]=local_sum;
    minimum_partial_results[partition_number]=local_minimum;
    maximum_partial_results[partition_number]=local_maximum;
}
//#####################################################################
template class Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper<float>;
template class Scalar_Multiply_And_Accumulate_And_Compute_Sum_And_Extrema_Helper<double>;
