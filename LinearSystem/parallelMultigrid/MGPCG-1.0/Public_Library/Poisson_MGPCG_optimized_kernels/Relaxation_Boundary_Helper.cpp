//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Relaxation_Boundary_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Relaxation_Boundary,(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,const T* const b,const T* const one_over_diagonal_part,const int* const boundary_index, const int* const block_start, const int* const block_end, const int number_of_red_blocks, const int number_of_black_blocks,const int loops,const bool reverse_order),(x_size_input,u,b,one_over_diagonal_part,boundary_index,block_start,block_end,number_of_red_blocks,number_of_black_blocks,loops,reverse_order))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Relaxation_Boundary_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Relaxation_Boundary_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int first_block;
    const int last_block;
    Relaxation_Boundary_Size_Specific_Thread_Helper(
        Relaxation_Boundary_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int first_block_input, const int last_block_input)
        :obj(obj_input),first_block(first_block_input),last_block(last_block_input) {}
    void Run(){obj->Run_Block_Range(first_block,last_block);}
};
}
template<class T,int y_size,int z_size> void Relaxation_Boundary_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems

    if(!reverse_order){
        for(int partition=1;partition<=number_of_partitions;partition++){
            int first_block_of_partition=(number_of_red_blocks/number_of_partitions)*(partition-1)+std::min(number_of_red_blocks%number_of_partitions,partition-1);
            int last_block_of_partition=(number_of_red_blocks/number_of_partitions)*partition+std::min(number_of_red_blocks%number_of_partitions,partition)-1;
            Relaxation_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>* task=
                new Relaxation_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>(this,first_block_of_partition,last_block_of_partition);
            pthread_queue->Queue(task);}
        pthread_queue->Wait();}
    
    for(int partition=1;partition<=number_of_partitions;partition++){
        int first_block_of_partition=(number_of_black_blocks/number_of_partitions)*(partition-1)+std::min(number_of_black_blocks%number_of_partitions,partition-1) +number_of_red_blocks;
        int last_block_of_partition=(number_of_black_blocks/number_of_partitions)*partition+std::min(number_of_black_blocks%number_of_partitions,partition)-1 +number_of_red_blocks;
	Relaxation_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>* task=
            new Relaxation_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>(this,first_block_of_partition,last_block_of_partition);
	pthread_queue->Queue(task);}
    pthread_queue->Wait();
    
    if(reverse_order){
	for(int partition=1;partition<=number_of_partitions;partition++){
            int first_block_of_partition=(number_of_red_blocks/number_of_partitions)*(partition-1)+std::min(number_of_red_blocks%number_of_partitions,partition-1);
            int last_block_of_partition=(number_of_red_blocks/number_of_partitions)*partition+std::min(number_of_red_blocks%number_of_partitions,partition)-1;
            Relaxation_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>* task=
                new Relaxation_Boundary_Size_Specific_Thread_Helper<T,y_size,z_size>(this,first_block_of_partition,last_block_of_partition);
            pthread_queue->Queue(task);}
        pthread_queue->Wait();}
}
//#####################################################################
// Function Run_Block_Range
//#####################################################################
template<class T,int y_size,int z_size> void Relaxation_Boundary_Size_Specific_Helper<T,y_size,z_size>::
Run_Block_Range(const int first_block,const int last_block)
{
    if(reverse_order)
        for(int block=last_block;block>=first_block;block--){
            int first_index=block_start[block],last_index=block_end[block];
            for(int loop=0;loop<loops;loop++)
                for(int i=last_index-1;i>=first_index;i--){
                    int index=boundary_index[i];
                    u[index]=(b[index]
                        -u[index+x_plus_one_shift]
                        -u[index+x_minus_one_shift]
                        -u[index+y_minus_one_shift]
                        -u[index+y_plus_one_shift]
                        -u[index+z_plus_one_shift]
                        -u[index+z_minus_one_shift])
                        *one_over_diagonal_part[i];}}
    else
        for(int block=first_block;block<=last_block;block++){
            int first_index=block_start[block],last_index=block_end[block];
            for(int loop=0;loop<loops;loop++)
                for(int i=first_index;i<last_index;i++){
                    int index=boundary_index[i];
                    u[index]=(b[index]
                        -u[index+x_plus_one_shift]
                        -u[index+x_minus_one_shift]
                        -u[index+y_minus_one_shift]
                        -u[index+y_plus_one_shift]
                        -u[index+z_plus_one_shift]
                        -u[index+z_minus_one_shift])
                        *one_over_diagonal_part[i];}}
}
//#####################################################################
template class Relaxation_Boundary_Helper<float>;
template class Relaxation_Boundary_Helper<double>;
