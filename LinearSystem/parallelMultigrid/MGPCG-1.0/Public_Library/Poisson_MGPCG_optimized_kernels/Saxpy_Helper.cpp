//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis,Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Saxpy_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
#include <stdlib.h>
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Saxpy,(const int x_size_input,const int y_size_input,const int z_size_input,const T& c1,const T* const u1,T* const result),(x_size_input,c1,u1,result))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Saxpy_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Saxpy_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax;
    Saxpy_Size_Specific_Thread_Helper(
        Saxpy_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax);}
};
}
template<class T,int y_size,int z_size> void Saxpy_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    int number_of_x_blocks=x_size/x_block_size;
    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
	Saxpy_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Saxpy_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Saxpy_Size_Specific_Helper<T,y_size,z_size>::
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
	    result[index]+=c1*u1[index];
        }
}
//#####################################################################
template class Saxpy_Helper<float>;
template class Saxpy_Helper<double>;

