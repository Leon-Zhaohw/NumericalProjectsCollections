//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Combined_Saxpy_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Combined_Saxpy,(const int x_size_input,const int y_size_input,const int z_size_input,T* const x_input,T* const p_input,const T* const z_input,const T alpha_input,const T beta_input),(x_size_input,x_input,p_input,z_input,alpha_input,beta_input))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Combined_Saxpy_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Combined_Saxpy_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Combined_Saxpy_Size_Specific_Thread_Helper(
        Combined_Saxpy_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax);}
};
}
template<class T,int y_size,int z_size> void Combined_Saxpy_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    int number_of_x_blocks=x_size/x_block_size;
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
	Combined_Saxpy_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Combined_Saxpy_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();    
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Combined_Saxpy_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax)
{   
    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size)
 
        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++)
        {    
            int index=i*x_shift+j*y_shift+k*z_shift;
            x[index]+=alpha*p[index];
            p[index]=z[index]+beta*p[index];
	}
}
//#####################################################################
template class Combined_Saxpy_Helper<float>;
template class Combined_Saxpy_Helper<double>;

