//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Reduce_Maxabs_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Reduce_Maxabs,(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u),(x_size_input,u))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Reduce_Maxabs_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Reduce_Maxabs_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,partition_number;
    Reduce_Maxabs_Size_Specific_Thread_Helper(
        Reduce_Maxabs_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int partition_number_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),partition_number(partition_number_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax,partition_number);}
};
}
template<class T,int y_size,int z_size> T Reduce_Maxabs_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    maxabs_partial_results=new T[number_of_partitions];
    int number_of_x_blocks=x_size/x_block_size;
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
	Reduce_Maxabs_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Reduce_Maxabs_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,partition);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();

    T maxabs=T();
    for(int i=0;i<number_of_partitions;i++) maxabs=std::max(maxabs,maxabs_partial_results[i]);
    delete[] maxabs_partial_results;maxabs_partial_results=0;
    return maxabs;
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T,int y_size,int z_size> void Reduce_Maxabs_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax,const int partition_number)
{   
    T local_maxabs=T();
    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {    
            int index=i*x_shift+j*y_shift+k*z_shift;
	    local_maxabs=std::max(local_maxabs,u[index]<0?-u[index]:u[index]);
	}

    maxabs_partial_results[partition_number]=local_maxabs;
}
//#####################################################################
template class Reduce_Maxabs_Helper<float>;
template class Reduce_Maxabs_Helper<double>;
