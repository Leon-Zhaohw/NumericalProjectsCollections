//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Residual_Boundary_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Residual_Boundary,(const int x_size_input,const int y_size_input,const int z_size_input,const T* const u,const T* const b,T* const r,const T* const diagonal_part,const int* const boundary_index,const int number_of_indices),(x_size_input,u,b,r,diagonal_part,boundary_index,number_of_indices))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Residual_Boundary_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Residual_Boundary_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int index_start;
    const int index_end;
    Residual_Boundary_Size_Specific_Thread_Helper(
        Residual_Boundary_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int index_start_input, const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};
}
template<class T,int y_size,int z_size> void Residual_Boundary_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems

    for(int partition=1;partition<=number_of_partitions;partition++){
	int first_index_of_partition=(number_of_indices/number_of_partitions)*(partition-1)+std::min(number_of_indices%number_of_partitions,partition-1);
	int last_index_of_partition=(number_of_indices/number_of_partitions)*partition+std::min(number_of_indices%number_of_partitions,partition);
	Residual_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Residual_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>(this,first_index_of_partition,last_index_of_partition);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();
}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
template<class T,int y_size,int z_size> void Residual_Boundary_Size_Specific_Helper<T,y_size,z_size>::
Run_Index_Range(const int first_index, const int last_index)
{   
    for(int i=first_index;i<last_index;i++){
	int index=boundary_index[i];
	if(diagonal_part[index]!=0)
	    r[index]=b[index]-diagonal_part[index]*u[index]
                -u[index+x_plus_one_shift]
                -u[index+x_minus_one_shift]
                -u[index+y_plus_one_shift]
                -u[index+y_minus_one_shift]
                -u[index+z_plus_one_shift]
                -u[index+z_minus_one_shift];
	else
	    r[index]=T();}
}
//#####################################################################
template class Residual_Boundary_Helper<float>;
template class Residual_Boundary_Helper<double>;

