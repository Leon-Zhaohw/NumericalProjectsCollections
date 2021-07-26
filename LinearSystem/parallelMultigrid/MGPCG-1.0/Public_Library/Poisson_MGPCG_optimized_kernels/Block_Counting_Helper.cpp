//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis, Aleka McAdams
// This file is governed by the license contained in the accompanying file COPYRIGHT.txt.
//#####################################################################
#include "Block_Counting_Helper.h"
#include "PTHREAD_QUEUE.h"
#include "Instantiation_Helpers.h"
using namespace PhysBAM;
extern PTHREAD_QUEUE* pthread_queue;
//#####################################################################
// Constructor/Desctructor of non-size-specific class
//#####################################################################
CONSTRUCTOR_INSTANTIATION_HELPERX(Block_Counting,(const int x_size_input,const int y_size_input,const int z_size_input,const unsigned char* const is_boundary_bytemask,const unsigned char* const is_extended_boundary_bytemask,int* const boundary_nodes_in_block,int* const extended_boundary_nodes_in_block,
	int& total_red_boundary_blocks_input,int& total_black_boundary_blocks_input,int& total_extended_boundary_blocks_input,
	int& total_red_boundary_indices_input,int& total_black_boundary_indices_input,int& total_extended_boundary_indices_input),(x_size_input,is_boundary_bytemask,is_extended_boundary_bytemask,boundary_nodes_in_block,extended_boundary_nodes_in_block,
	    total_red_boundary_blocks_input,total_black_boundary_blocks_input,total_extended_boundary_blocks_input,
	    total_red_boundary_indices_input,total_black_boundary_indices_input,total_extended_boundary_indices_input))
//#####################################################################
// Function Run_Parallel
//#####################################################################
namespace{
template<class T,int y_size,int z_size>
struct Block_Counting_Size_Specific_Thread_Helper:public PhysBAM::PTHREAD_QUEUE::TASK
{
    Block_Counting_Size_Specific_Helper<T,y_size,z_size>* const obj;
    const int xmin,xmax,partition;
    Block_Counting_Size_Specific_Thread_Helper(Block_Counting_Size_Specific_Helper<T,y_size,z_size>* const obj_input,const int xmin_input,const int xmax_input,const int partition_input)
        :obj(obj_input),xmin(xmin_input),xmax(xmax_input),partition(partition_input) {}
    void Run(){obj->Run_X_Range(xmin,xmax,partition);}
};
}
template<class T,int y_size,int z_size> void Block_Counting_Size_Specific_Helper<T,y_size,z_size>::
Run_Parallel(const int number_of_partitions)
{
    if(x_size*y_size*z_size<=32*32*32) {Run();return;} // Run serially for small problems


    red_boundary_blocks_partial_results = new int[number_of_partitions];
    black_boundary_blocks_partial_results = new int[number_of_partitions];
    extended_boundary_blocks_partial_results = new int[number_of_partitions];
    red_boundary_indices_partial_results = new int[number_of_partitions];
    black_boundary_indices_partial_results = new int[number_of_partitions];
    extended_boundary_indices_partial_results = new int[number_of_partitions];


    for(int partition=0;partition<number_of_partitions;partition++){
        int first_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition)+std::min(number_of_x_blocks%number_of_partitions,partition)+1;
        int last_block_of_partition=(number_of_x_blocks/number_of_partitions)*(partition+1)+std::min(number_of_x_blocks%number_of_partitions,partition+1);
        int xmin=(first_block_of_partition-1)*x_block_size+1;
        int xmax=last_block_of_partition*x_block_size;
        Block_Counting_Size_Specific_Thread_Helper<T,y_size,z_size>* task=new Block_Counting_Size_Specific_Thread_Helper<T,y_size,z_size>(this,xmin,xmax,partition);
        pthread_queue->Queue(task);}
    pthread_queue->Wait();

    
    total_red_boundary_blocks = 0;
    total_black_boundary_blocks = 0;
    total_extended_boundary_blocks=0;
    total_red_boundary_indices=0;
    total_black_boundary_indices=0;
    total_extended_boundary_indices=0;
    for(int partition=0;partition<number_of_partitions;partition++){
	total_red_boundary_blocks += red_boundary_blocks_partial_results[partition];
	total_black_boundary_blocks += black_boundary_blocks_partial_results[partition];
	total_extended_boundary_blocks+=extended_boundary_blocks_partial_results[partition];
	total_red_boundary_indices+=red_boundary_indices_partial_results[partition];
	total_black_boundary_indices+=black_boundary_indices_partial_results[partition];
	total_extended_boundary_indices+=extended_boundary_indices_partial_results[partition];
    }
	
    delete[] red_boundary_blocks_partial_results; red_boundary_blocks_partial_results=0;
    delete[] black_boundary_blocks_partial_results; black_boundary_blocks_partial_results=0;
    delete[] extended_boundary_blocks_partial_results; extended_boundary_blocks_partial_results=0;
    delete[] red_boundary_indices_partial_results; red_boundary_indices_partial_results=0;
    delete[] black_boundary_indices_partial_results; black_boundary_indices_partial_results=0;
    delete[] extended_boundary_indices_partial_results; extended_boundary_indices_partial_results=0;
}
//#####################################################################
// Function Run_X_Range
//#####################################################################
template<class T,int y_size,int z_size> void Block_Counting_Size_Specific_Helper<T,y_size,z_size>::
Run_X_Range(const int xmin,const int xmax,const int partition)
{    
    
    int local_rbb=0,local_bbb=0,local_ebb=0,local_rbi=0,local_bbi=0,local_ebi=0;
    for(int block_i=xmin;block_i<=xmax;block_i+=x_block_size)
    for(int block_j=1;block_j<=y_size;block_j+=y_block_size)
    for(int block_k=1;block_k<=z_size;block_k+=z_block_size){
 
        int bi=(block_i-1)/x_block_size;
        int bj=(block_j-1)/y_block_size;
        int bk=(block_k-1)/z_block_size;
        
        int block_index=(bi*number_of_y_blocks+bj)*number_of_z_blocks+bk;
	int block_color=(bi+bj+bk)%2;
            
        boundary_nodes_in_block[block_index]=0;
	extended_boundary_nodes_in_block[block_index]=0;
	

        for(int i=block_i;i<block_i+x_block_size;i++)
        for(int j=block_j;j<block_j+y_block_size;j++)
        for(int k=block_k;k<block_k+z_block_size;k++){

            int index=i*x_shift+j*y_shift+k*z_shift;
            if(is_boundary_bytemask[index]) boundary_nodes_in_block[block_index]++;
            if(is_extended_boundary_bytemask[index]) extended_boundary_nodes_in_block[block_index]++;
        }

	local_ebb+=extended_boundary_nodes_in_block[block_index]?1:0;
	local_ebi+=extended_boundary_nodes_in_block[block_index];

	if(block_color){
	    local_bbb+=boundary_nodes_in_block[block_index]?1:0;
	    local_bbi+=boundary_nodes_in_block[block_index];
	}else{	    
	    local_rbb+=boundary_nodes_in_block[block_index]?1:0;
	    local_rbi+=boundary_nodes_in_block[block_index];
	}
	
    }
    red_boundary_blocks_partial_results[partition]=local_rbb;
    black_boundary_blocks_partial_results[partition]=local_bbb;
    extended_boundary_blocks_partial_results[partition]=local_ebb;
    red_boundary_indices_partial_results[partition]=local_rbi;
    black_boundary_indices_partial_results[partition]=local_bbi;
    extended_boundary_indices_partial_results[partition]=local_ebi;
}
//#####################################################################
template class Block_Counting_Helper<void>;

